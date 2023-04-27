# -*- coding: utf-8 -*-
using Revise
using Plots
using RelevanceStacktrace
using StatsBase 
using Distributions
using ProgressMeter
using Statistics
using Random
using CurveFit
import Random:rand
using SpecialFunctions

using InteractingLOBs

Revise.revise()

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
#num_paths = 2
num_paths = 3#30

L = 200     # real system width (e.g. 200 meters)
M = 400     # divided into M pieces , 400

T = 5000      # simulation runs until real time T (e.g. 80 seconds)
#T = 30   # simulation runs until real time T (e.g. 80 seconds)
p₀ = 0.0#238.75 #this is the mid_price at t=0 

# Free-Parameters for gaussian version
D = 0.5 # real diffusion constant e.g. D=1 (meters^2 / second), 1
σ = 1.0 

ν = 0.8#0.8
α_slob = 40.0
α_lob = 0.0
α = α_lob
γ = 0.7
r = 0.5 #proportion of time in which it jumps left or right

dist = Normal(0.0,1.0)
#dist = TDist(1)
#dist = Spl(1);
# -

Δx = L / M  # real gap between simulation points 
#Δt = (Δx^2) / (2.0 * D) / 10 #* (2.0/3.0) # real time seperation between simulation points
Δt = (r * (Δx^2) / (2.0 * D))^(1/γ)
(Δt,to_simulation_time(T,Δt),num_paths*to_simulation_time(T,Δt)) #about 2GB RAM per 100K, i.e. can only do about 1.8 million

# +
λ = 1.0
μ = 0.1 

mySourceTerm = SourceTerm(λ, μ, true);

# +
# coupling:
a = 13.0  #gap between stocks before at full strength: strong is 0.3
b = 1.0   #weighting of interaction term: strong is 2
c = 1.2   #skew factor: strong is 2

myCouplingTerm = CouplingTerm(μ, a, b, c, true);

# +
RealStartTime = 0 # when, in real time, to kick the system
SimStartTime = to_simulation_time(RealStartTime,Δt)-2 # convert to simulation time
SimEndTime = SimStartTime + 10 # when to stop kicking, in simulation time
Position = 200
Volume = -8;

# If position == -x where x>=0, then put it x above the mid price each time

# +
myRLPusher1 = RLPushTerm(SimStartTime,SimEndTime,Position,Volume,false)

lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, σ, ν, α, r, γ, dist, 
    mySourceTerm, myCouplingTerm, myRLPusher1,true);
# -

#to_simulation_time(T,Δt)
lob_model¹.SK_DP

# +
myRLPusher2 = RLPushTerm(SimStartTime,SimEndTime,Position,Volume,false)

lob_model² = SLOB(num_paths, T, p₀, M, L, D, σ, ν, α, r, γ, dist,
    mySourceTerm, myCouplingTerm, myRLPusher2,true);

# +
lob_densities¹, sources¹, couplings¹, rl_pushes¹, raw_price_paths¹, sample_price_paths¹, P⁺s¹, P⁻s¹, Ps¹, V¹,
lob_densities², sources², couplings², rl_pushes², raw_price_paths², sample_price_paths², P⁺s², P⁻s², Ps², V² = 
    [nothing for _ in 1:20]

GC.gc()

lob_densities¹, sources¹, couplings¹, rl_pushes¹, raw_price_paths¹, sample_price_paths¹, P⁺s¹, P⁻s¹, Ps¹, V¹,
lob_densities², sources², couplings², rl_pushes², raw_price_paths², sample_price_paths², P⁺s², P⁻s², Ps², V², broke_points =
    InteractOrderBooks(lob_model¹,lob_model², -1, true) ;

# +
#function NameAll(a)
#    return  (lob_densities¹ = a[1], sources¹ = a[2], couplings¹ = a[3], rl_pushes¹ = a[4], 
#            raw_price_paths¹ = a[5], sample_price_paths¹ = a[6], P⁺s¹ = a[7], P⁻s¹ = a[8], Ps¹ = a[9], V¹ = a[10],
#            lob_densities² = a[11], sources² = a[12], couplings² = a[13], rl_pushes² = a[14], 
#            raw_price_paths² = a[15], sample_price_paths² = a[16], P⁺s² = a[17], P⁻s² = a[18], Ps² = a[19], V² = a[20]) 
#end;
# -

to_sim(x) = to_simulation_time(x,lob_model¹.Δt)
to_real(x) = to_real_time(x,lob_model¹.Δt)

# +
total_steps = min(to_sim(T),100)
total_length = to_sim(T)
step = floor(Int,total_length/total_steps)

# the below just ensure we see the full graph (100-myp)% of the time
myp = 10
max_y¹ = percentile( [maximum(lob_densities¹[:,i,1]) for i in 1:length(lob_densities¹[1,:,1])] , 100-myp)
min_y¹ = percentile( [minimum(lob_densities¹[:,i,1]) for i in 1:length(lob_densities¹[1,:,1])] , myp)
max_y² = percentile( [maximum(lob_densities²[:,i,1]) for i in 1:length(lob_densities²[1,:,1])] , 100-myp)
min_y² = percentile( [minimum(lob_densities²[:,i,1]) for i in 1:length(lob_densities²[1,:,1])] , myp)

#x_axis_width = #4

l = @layout [a d; b e; c f];


# +
function plot_price_path(s, r, lob_model, raw_price_paths, sample_price_paths)
    plt = plot((0:s-1).*lob_model.Δt, raw_price_paths[1:s,1],color=1,w=0.6) ;
    for path in 2:lob_model.num_paths
        plot!((0:s-1).*lob_model.Δt, raw_price_paths[1:s,path],color=5,w=0.6) ;
    end
    plot!(0:r-1, sample_price_paths[1:r,1],color=1,w=2.7) ;
    plot!(legend=false, ylab="Price", xlab="Time") ;   
    return plt
end

function plot_density_visual(s, r, lob_model, lob_densities, couplings, sources, rl_pushes, raw_price_paths, 
                                                    plot_raw_price=true,x_axis_width = L/2)
    x_axis  = [p₀-x_axis_width,p₀+x_axis_width]
    path_to_plot = 1
    @assert path_to_plot <= lob_model.num_paths
    plt = plot(lob_model.x, lob_densities[:,s,path_to_plot], color=1,label="Density"); 
    plot!(lob_model.x, couplings[:,s,path_to_plot], color=2, label="Coupling") ;
    plot!(lob_model.x, sources[:,s,path_to_plot], color=3, label="Source") ;
    plot!(lob_model.x, rl_pushes[:,s,path_to_plot], color=4, label="RL") ;
    if (plot_raw_price)
        if (isinteger(s*lob_model.Δt))
            mycol = :black
        else 
            mycol = 1
        end
        scatter!([raw_price_paths[s,path_to_plot]],[0],label="Midprice",markercolor=mycol, markersize=3,markerstrokewidth=0.5)
    end
    #plot!(lob_model.x, x -> 0,color="black", primary=false) ; #draw horizontal line
    #plot!( legend=:bottomleft, title="time=$r", xlab="Price", ylab="LOB&Source",ylim=[min_y¹,max_y¹],xlim=x_axis)
    #plot!( legend=:bottomleft, title="time=$r", xlab="Price", ylab="LOB&Source",ylim=[-20,1],xlim=x_axis)
    plot!( legend=:bottomleft, title="time=$r", xlab="Price", ylab="LOB&Source",xlim=x_axis)
    return plt
end;

# +
range = 1:step:(total_length-step)
#range = [SimStartTime:SimStartTime+to_sim(1)*2;]
#range = [SimStartTime]

p_outer = Progress(length(range),dt=0.1)

anim = @animate for s = range           #s is the time in simulation time
    r = to_real_time(s, lob_model¹.Δt)  #r is the time in real time
    #s = to_simulation_time(r, lob_model¹.Δt)  #s is the time in real time
    
    #global prev_r, plt1, plt2, plt3, plt4, plt5, plt6, sim_times, real_times
    
    
    plt1 = plot_price_path(s,r,lob_model¹,raw_price_paths¹,sample_price_paths¹)
    plt3 = plot_price_path(s,r,lob_model²,raw_price_paths²,sample_price_paths²)
    plt5 = plot_price_path(s,r,lob_model²,raw_price_paths¹-raw_price_paths²,sample_price_paths¹-sample_price_paths²)
    
    plt2 = plot_density_visual(s, r, lob_model¹, lob_densities¹, couplings¹, sources¹, rl_pushes¹,raw_price_paths¹)
    plt4 = plot_density_visual(s, r, lob_model², lob_densities², couplings², sources², rl_pushes²,raw_price_paths²)
    plt6 = plot_density_visual(s, r, lob_model², lob_densities¹+lob_densities², 
                                                 couplings¹+couplings², 
                                                 sources¹+sources², 
                                                 rl_pushes¹+rl_pushes²,raw_price_paths¹,false)
    
    plot(plt1, plt2, plt3, plt4, plt5, plt6 ,layout=l,size=(1000,1000))
    
    next!(p_outer)
end

gif(anim, "/tmp/LOB.gif", fps=20*length(range)/200)
#gif(anim, "~/Desktop/Masters/StillWorking/Random_Walk_And_Coupling_For_Alpha_Is_0.7.gif", fps=20*length(range)/200)
# -

how_many = 12
#l2 = @layout [a d; b e; c f; ];
p_arr1 = Array{Plots.Plot{Plots.GRBackend},1}(undef,how_many)
for i in [1:how_many;]
p_arr1[i] = plot_density_visual(SimStartTime-2+i, to_real(SimStartTime-2+i), lob_model¹, lob_densities¹, couplings¹, sources¹, rl_pushes¹,raw_price_paths¹,true,10)
end
plot(p_arr1...,size=(1000,1000))
#gif(anim, "~/Desktop/Masters/GoodPics/ForM2.gif", fps=16)
#png("/home/derickdiana/Desktop/Masters/GoodPics/ForM2.png")

SimStartTime

# +
observed_price_path = sample_price_paths¹[1:end-1,1];
observed_log_returns = diff(log.(observed_price_path[:,1]));

data_stylized_facts = StylizedFacts.StylizedFactsPlot(observed_price_path);

# +
StylizedFacts.plot_all_stylized_facts(data_stylized_facts,(1000,1200))

#png("/home/derickdiana/Desktop/Masters/GoodPics/DefinitelyNoTails.png")

# +
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

# +
function myvar(p,x,Δx)
    mu = mymean(p,x,Δx)
    Ex2 = 0
    for i in 1:length(p)
        Ex2 += p[i]*x[i]^2*Δx
    end
    return Ex2 - mu^2
end

function mymean(p,x,Δx)
    sum = 0
    for i in 1:length(p)
        sum += p[i]*x[i]*Δx
    end
    return sum
end

test = lob_densities¹[:,end-5,2]
test = test/sum(test)
t_mean = mymean(test ,lob_model¹.x,lob_model¹.Δx)
t_var = myvar(test ,lob_model¹.x,lob_model¹.Δx)
print(t_mean," ",t_var^0.5)
plot(lob_model¹.x,test)
plot!([t_mean],seriestype="vline")
plot!([t_mean - t_var^0.5],seriestype="vline")
plot!([t_mean + t_var^0.5],seriestype="vline")

# +
l = length(lob_densities¹[1,:,1]) # number of time steps

variances = zeros(Float64,l-1)
sums = zeros(Float64,l-1)
for t in 2:l 
    temp = lob_densities¹[:,t,1]/sum(lob_densities¹[:,t,1]) #across all x values, at time "t",  on path 1
    sums[t-1] = sum(temp)
    variances[t-1] = myvar(temp,lob_model¹.x,lob_model¹.Δx)
end

l2 =  @layout [a;b];
sub = 7:length(sums)-1
plt1 = plot(variances[sub],label="Variance")
plt2 = plot(sums[sub],label=string("Area underneath has range ",maximum(sums[sub]) - minimum(sums[sub])),color="red")
plot(plt1, plt2, layout=l2, size=(700,700),xlab="Time")
#png("/home/derickdiana/Desktop/Masters/StillWorking/VarianceAndSumsForAlphaIs0.5.png")

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

# +
# function fractional_choose(n,k)
#     if k==1
#         temp = 0
#     else 
#         temp = 0
#     end
    
#     if (n==0)
#         return 0.0 + temp
#     end
    
#     return gamma(n+1)/(gamma(n-k+1)*gamma(k+1)) - 1/(n+1)*1/beta(n-k+1,k+1)#*(-1)^k + temp
# end
# n = 200
# for m in 1:(n-1)
#     print(fractional_choose(1-0.8,n-m)," ")
# end
