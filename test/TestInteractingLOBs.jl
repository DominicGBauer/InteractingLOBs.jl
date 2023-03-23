# -*- coding: utf-8 -*-
using Revise
using Plots
using RelevanceStacktrace
using StatsBase 
using Distributions
using ProgressMeter
using Statistics
using Random
import Random:rand

using InteractingLOBs

# +
struct Spl <: Sampleable{Univariate, Continuous}
    p::Float64
end  

function rand(s::Spl)
    if rand()<s.p
        return rand(Normal(0.0, 1.0))
    else
        return 3
    end
end

# +
# Configuration Arguments
num_paths = 3

L = 200     # real system width (e.g. 200 meters)
M = 400     # divided into M pieces , 400

T = 20      # simulation runs until real time T (e.g. 80 seconds)
p₀ = 238.75 #this is the mid_price at t=0 

# Free-Parameters for gaussian version
D = 1.0 # real diffusion constant e.g. D=1 (meters^2 / second), 1
σ = 1.0 

ν = 0.5
α_slob = 40.0
α_lob = 0.0
α = α_lob

#dist = Normal(0.0,1.0)
#dist = TDist(1)
dist = Spl(1);
# -

Δx = L / M  # real gap between simulation points 
Δt = (Δx^2) / (2.0 * D)  # real time seperation between simulation points

# +
λ = 1.0
μ = 0.1 

mySourceTerm = SourceTerm(λ, μ);

# +
# coupling:
a = 13.0  #gap between stocks before at full strength: strong is 0.3
b = 1.0   #weighting of interaction term: strong is 2
c = 1.2   #skew factor: strong is 2

myCouplingTerm = CouplingTerm(μ, a, b, c);

# +
RealStartTime = 10 # when, in real time, to kick the system
SimStartTime = to_simulation_time(RealStartTime,Δt) # convert to simulation time
SimEndTime = SimStartTime + 1 # when to stop kicking, in simulation time
Position = 0
Volume = 8;

# If position == -x where x>=0, then put it x above the mid price each time

# +
myRLPusher1 = RLPushTerm(SimStartTime,SimEndTime,Position,Volume,true)

lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, σ, ν, α, dist, 
    mySourceTerm, myCouplingTerm, myRLPusher1);

# +
myRLPusher2 = RLPushTerm(SimStartTime,SimEndTime,Position,Volume,false)

lob_model² = SLOB(num_paths, T, p₀, M, L, D, σ, ν, α, dist,
    mySourceTerm, myCouplingTerm, myRLPusher2);

# +
lob_densities¹, sources¹, couplings¹, rl_pushes¹, raw_price_paths¹, sample_price_paths¹, P⁺s¹, P⁻s¹, Ps¹, 
lob_densities², sources², couplings², rl_pushes², raw_price_paths², sample_price_paths², P⁺s², P⁻s², Ps² =

InteractOrderBooks(lob_model¹,lob_model², -1, true) ;
# -

to_sim(x) = to_simulation_time(x,lob_model¹.Δt)
to_real(x) = to_real_time(x,lob_model¹.Δt)

# +
total_steps = to_sim(T)
total_length = to_sim(T)
step = floor(Int,total_length/total_steps)

# the below just ensure we see the full graph (100-myp)% of the time
myp = 10
max_y¹ = percentile( [maximum(lob_densities¹[:,i,1]) for i in 1:length(lob_densities¹[1,:,1])] , 100-myp)
min_y¹ = percentile( [minimum(lob_densities¹[:,i,1]) for i in 1:length(lob_densities¹[1,:,1])] , myp)
max_y² = percentile( [maximum(lob_densities²[:,i,1]) for i in 1:length(lob_densities²[1,:,1])] , 100-myp)
min_y² = percentile( [minimum(lob_densities²[:,i,1]) for i in 1:length(lob_densities²[1,:,1])] , myp)


function plot_price_path(s, r, lob_model, raw_price_paths, sample_price_paths)
    plt = plot((0:s-1).*lob_model.Δt, raw_price_paths[1:s,1],color=1,w=0.6) ;
    for path in 2:lob_model.num_paths
        plot!((0:s-1).*lob_model.Δt, raw_price_paths[1:s,path],color=5,w=0.6) ;
    end
    plot!(0:r-1, sample_price_paths[1:r,1],color=1,w=2.7) ;
    plot!(legend=false, ylab="Price", xlab="Time") ;   
    return plt
end

function plot_density_visual(s, r, lob_model, lob_densities, couplings, sources, rl_pushes)
    plt = plot(lob_model.x, lob_densities[:,s,1], color=1,label="Density"); 
    plot!(lob_model.x, couplings[:,s,1], color=2, label="Coupling") ;
    plot!(lob_model.x, sources[:,s,1], color=3, label="Source") ;
    plot!(lob_model.x, rl_pushes[:,s,1], color=4, label="RL") ;
    #plot!(lob_model.x, x -> 0,color="black", primary=false) ; #draw horizontal line
    plot!( legend=:bottomleft, title="time=$r", xlab="Price", ylab="LOB&Source",ylim=[min_y¹,max_y¹],xlim=[150,340])#xlim=[150,340]) ;
    return plt
end


range = 1:step:(total_length-step)
#range = [SimStartTime:SimStartTime+to_sim(1);]

anim = @animate for s = range           #s is the time in simulation time
    r = to_real_time(s, lob_model¹.Δt)  #r is the time in real time
    
    #global prev_r, plt1, plt2, plt3, plt4, plt5, plt6, sim_times, real_times
    
    l = @layout [a d; b e; c f]
    
    plt1 = plot_price_path(s,r,lob_model¹,raw_price_paths¹,sample_price_paths¹)
    plt3 = plot_price_path(s,r,lob_model²,raw_price_paths²,sample_price_paths²)
    plt5 = plot_price_path(s,r,lob_model²,raw_price_paths¹-raw_price_paths²,sample_price_paths¹-sample_price_paths²)
    
    plt2 = plot_density_visual(s, r, lob_model¹, lob_densities¹, couplings¹, sources¹, rl_pushes¹)
    plt4 = plot_density_visual(s, r, lob_model², lob_densities², couplings², sources², rl_pushes²)
    plt6 = plot_density_visual(s, r, lob_model², lob_densities¹+lob_densities², 
                                                 couplings¹+couplings², 
                                                 sources¹+sources², 
                                                 rl_pushes¹+rl_pushes²)
    
    plot(plt1, plt2, plt3, plt4, plt5, plt6 ,layout=l,size=(1000,1000))
end

gif(anim, "/tmp/LOB.gif", fps=20*length(range)/200)
#gif(anim, "~/Desktop/Masters/GoodPics/x.gif", fps=8)
#gif(anim, "~/Desktop/Masters/GoodPics/LotsOfStuff.gif", fps=16)
# -

sum(abs.(lob_densities¹[:,10,1]))

real_times

# +
observed_price_path = sample_price_paths¹[1:end-1,1];
observed_log_returns = diff(log.(observed_price_path[:,1]));

data_stylized_facts = StylizedFacts.StylizedFactsPlot(observed_price_path);
# -

StylizedFacts.plot_all_stylized_facts(data_stylized_facts,(1000,1200))

raw_price_paths¹

Revise.revise()

# +
#middle = 3:length(sample_price_paths²[:,1])
#price_changes = (sample_price_paths²[middle.-1,:] - sample_price_paths²[middle.-2,:])[:]

# +
#Plot Log Returns Hist
#dist = fit(Normal, observed_log_returns)
#plt = histogram(observed_log_returns, normalize=true, legend=false)
#plot!(plt, x -> pdf(dist, x), xlim=xlims(), legend=false)
#
#StylizedFacts.plot_log_returns(data_stylized_facts, "")
#StylizedFacts.plot_qq_log_returns(data_stylized_facts, "")
#StylizedFacts.plot_acf_order_flow(market_data_stylized_facts, "")
#StylizedFacts.plot_acf_log_returns(data_stylized_facts, "")
#StylizedFacts.plot_acf_abs_log_returns(data_stylized_facts, "")

# +
# save_object("OneRun.jld2",(lob_densities¹, sources¹, couplings¹, rl_pushes¹, raw_price_paths¹, sample_price_paths¹, P⁺s¹, P⁻s¹, Ps¹, 
# lob_densities², sources², couplings², rl_pushes², raw_price_paths², sample_price_paths², P⁺s², P⁻s², Ps²))
# temp = load_object("OneRun.jld2")

# +
# combine all different paths into one large vector
#diff(log.(sample_price_paths¹[1:end-1,:]),dims=1)[:]

# +
#observed_summary_stats = get_summary_stats(observed_log_returns) #is in AdaptiveABC
