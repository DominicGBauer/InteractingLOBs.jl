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
using JLD2

using InteractingLOBs

Revise.revise()

# +
#to_sim(x) = to_simulation_time(x,lob_model¹.Δt)
#to_real(x) = to_real_time(x,lob_model¹.Δt)
# -

# # General working

# +
# Configuration Arguments
num_paths = 1#30

L = 200     # real system width (e.g. 200 meters)
M = 400     # divided into M pieces , 400

T = 10000  # simulation runs until real time T (e.g. 80 seconds)
p₀ = 230.0  #this is the mid_price at t=0  238.75 

# Free-Parameters for gaussian version
D = 0.5 # real diffusion constant e.g. D=1 (meters^2 / second), 1
α = 0.0 # legacy, no longer used

ν = 10.0 #removal rate
γ = 1.0 #fraction of derivative (1 is normal diffusion, less than 1 is D^{1-γ} derivative on the RHS)

# Source term:
λ = 1.0 #
μ = 0.1 #

mySourceTerm = SourceTerm(λ, μ, true);

# Coupling term:
a = 13.0  #gap between stocks before at full strength: strong is 0.3
b = 1.0   #weighting of interaction term: strong is 2
c = 1.2   #skew factor: strong is 2

myCouplingTerm = CouplingTerm(μ, a, b, c, true);

# My randomness term
σ = 1.0 #variance in randomness
r = 0.5 #proportion of time in which it jumps left or right

myRandomnessTerm = RandomnessTerm(σ,r,true)

Δx = L / M  # real gap between simulation points 
Δt = (r * (Δx^2) / (2.0 * D))^(1/γ)

# RL Stuff:
RealStartTime = 2 # when, in real time, to kick the system
SimStartTime = to_simulation_time(RealStartTime,Δt)-2 # convert to simulation time
SimEndTime = SimStartTime + 10 # when to stop kicking, in simulation time
Position = 200
Volume = -8; # If position == -x where x>=0, then put it x above the mid price each time

myRLPusher1 = RLPushTerm(SimStartTime,SimEndTime,Position,Volume,false)

myRLPusher2 = RLPushTerm(SimStartTime,SimEndTime,Position,Volume,false)

lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ, 
    mySourceTerm, myCouplingTerm, myRLPusher1, myRandomnessTerm);

lob_model² = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
    mySourceTerm, myCouplingTerm, myRLPusher2, myRandomnessTerm);

print((Δt,to_simulation_time(T,Δt),num_paths*to_simulation_time(T,Δt))) #about 2GB RAM per 100K, i.e. can only do about 1.8 million
lob_model¹.SK_DP

# +
# clear everything pointed to by the dictionary then garbage collect. If it wasn't assigned yet it will left you know
try clear_double_dict(Dat) catch e print("Not initialized") end
GC.gc()

Dat = InteractOrderBooks([lob_model¹,lob_model²], -1, true) ;
#Dat = InteractOrderBooks([lob_model¹], -1, true) ;

# +
# slob = lob_model¹
# num_time_steps = to_simulation_time(T,Δt)

# lob_densities¹ =    zeros(Float64, slob.M+1, num_time_steps + 1, 2, num_paths)
# sources¹ =          zeros(Float64, slob.M+1, num_time_steps + 1, 2, num_paths)
# couplings¹ =        zeros(Float64, slob.M+1, num_time_steps + 1, 2, num_paths)
# rl_pushes¹ =        zeros(Float64, slob.M+1, num_time_steps + 1, 2, num_paths)

# raw_price_paths¹ =  ones(Float64,           num_time_steps + 1, 2, num_paths)
# sample_price_paths¹ =  ones(Float64,           slob.T         + 1, 2, num_paths)

# lob_densities² =    zeros(Float64, slob.M+1, num_time_steps + 1, 2, num_paths)
# sources² =          zeros(Float64, slob.M+1, num_time_steps + 1, 2, num_paths)
# couplings² =        zeros(Float64, slob.M+1, num_time_steps + 1, 2, num_paths)
# rl_pushes² =        zeros(Float64, slob.M+1, num_time_steps + 1, 2, num_paths)

# raw_price_paths² =  ones(Float64,           num_time_steps + 1, 2, num_paths)
# sample_price_paths² =  ones(Float64,           slob.T         + 1, 2, num_paths)


# for path_num in 1:num_paths
#     lob_densities¹[:,:,path_num] = Dat[path_num][1].lob_densities[:,:]
#     sources¹[:,:,path_num] = Dat[path_num][1].sources[:,:,:]
#     couplings¹[:,:,path_num] = Dat[path_num][1].couplings[:,:,:]
#     rl_pushes¹[:,:,path_num] = Dat[path_num][1].rl_pushes[:,:,:]
#     raw_price_paths¹[:,path_num] = Dat[path_num][1].raw_price_paths[:,:]
#     sample_price_paths¹[:,path_num] = Dat[path_num][1].obs_price_paths[:,:]

#     lob_densities²[:,:,path_num] = Dat[path_num][2].lob_densities[:,:,:]
#     sources²[:,:,path_num] = Dat[path_num][2].sources[:,:,:]
#     couplings²[:,:,path_num] = Dat[path_num][2].couplings[:,:,:]
#     rl_pushes²[:,:,path_num] = Dat[path_num][2].rl_pushes[:,:,:]
#     raw_price_paths²[:,path_num] = Dat[path_num][2].raw_price_paths[:,:]
#     sample_price_paths²[:,path_num] = Dat[path_num][2].obs_price_paths[:,:];
# end

# +
# the below just ensure we see the full graph (100-myp)% of the time
myp = 10
path_to_plot = 1
num_time_steps = to_simulation_time(T,Δt)

max_y¹ = percentile( [maximum(Dat[path_to_plot][1].lob_densities[:,i]) for i in 1:num_time_steps] , 100-myp)
min_y¹ = percentile( [minimum(Dat[path_to_plot][1].lob_densities[:,i]) for i in 1:num_time_steps] , myp)
if length(Dat[1])>1
    max_y² = percentile( [maximum(Dat[path_to_plot][2].lob_densities[:,i]) for i in 1:num_time_steps] , 100-myp)
    min_y² = percentile( [minimum(Dat[path_to_plot][2].lob_densities[:,i]) for i in 1:num_time_steps] , myp);
end

#x_axis_width = #4

l = @layout [a d; b e; c f];

function plot_price_path(s, r, lob_num, Dat, diff=false)
    lob_model = Dat[1][lob_num].slob
    
    plt = plot()
    for path in 1:lob_model.num_paths
        raw_price_paths = Dat[path][lob_num].raw_price_paths[1:s]
        if diff
            raw_price_paths .-=  Dat[path][3-lob_num].raw_price_paths[1:s]
        end
        plot!((0:s-1).*lob_model.Δt,raw_price_paths ,color=5,w=0.6) ;
        
    end
    
    obs_price_paths = Dat[path_to_plot][lob_num].obs_price_paths[1:r]
    if diff
        obs_price_paths .-=  Dat[path_to_plot][3-lob_num].obs_price_paths[1:r]
    end
    plot!(0:r-1, obs_price_paths,color=1,w=2.7) ;
    
    plot!(legend=false, ylab="Price", xlab="Time") ;   
    return plt
end

function plot_density_visual(s, r, lob_num, Dat, dosum=false, plot_raw_price=true, x_axis_width = L/2)
    lob_model = Dat[1][lob_num].slob
    x_axis  = [lob_model.p₀-x_axis_width,lob_model.p₀+x_axis_width]
    
    lob_densities = Dat[path_to_plot][lob_num].lob_densities[:,s]
    couplings = Dat[path_to_plot][lob_num].couplings[:,s]
    sources = Dat[path_to_plot][lob_num].sources[:,s]
    rl_pushes = Dat[path_to_plot][lob_num].rl_pushes[:,s]
    
    if dosum
        lob_densities .+= Dat[path_to_plot][3-lob_num].lob_densities[:,s]
        couplings .+= Dat[path_to_plot][3-lob_num].couplings[:,s]
        sources .+= Dat[path_to_plot][3-lob_num].sources[:,s]
        rl_pushes .+= Dat[path_to_plot][3-lob_num].rl_pushes[:,s]
    end
    
    plt = plot(lob_model.x, lob_densities , color=1,label="Density"); 
    plot!(lob_model.x, couplings, color=2, label="Coupling") ;
    plot!(lob_model.x, sources, color=3, label="Source") ;
    plot!(lob_model.x, rl_pushes, color=4, label="RL") ;
    
    if (plot_raw_price)
        if (isinteger(s*lob_model.Δt))
            mycol = :black
        else 
            mycol = 1
        end
        scatter!([Dat[path_to_plot][lob_num].raw_price_paths[s]],[0],label="Midprice",markercolor=mycol, markersize=3,markerstrokewidth=0.5)
    end
    #plot!(lob_model.x, x -> 0,color="black", primary=false) ; #draw horizontal line
    #plot!( legend=:bottomleft, title="time=$r", xlab="Price", ylab="LOB&Source",ylim=[min_y¹,max_y¹],xlim=x_axis)
    #plot!( legend=:bottomleft, title="time=$r", xlab="Price", ylab="LOB&Source",ylim=[-20,1],xlim=x_axis)
    plot!( legend=:bottomleft, title="time=$r", xlab="Price", ylab="LOB&Source",xlim=x_axis)
    return plt
end;

# +
#total_steps = min(to_simulation_time(T,Δt),100)
total_steps = 10
total_length = to_simulation_time(T,Δt)
step = floor(Int,total_length/total_steps)

range = 1:step:(total_length-step)
#range = [SimStartTime:SimStartTime+to_sim(1)*2;]
#range = [SimStartTime]

p_outer = Progress(length(range),dt=0.1)

anim = @animate for s = range           #s is the time in simulation time
    r = to_real_time(s, lob_model¹.Δt)  #r is the time in real time
    #s = to_simulation_time(r, lob_model¹.Δt)  #s is the time in real time
    
    plt1 = plot_price_path(s,r,1,Dat)
    if length(Dat[1])>1
        plt3 = plot_price_path(s,r,2,Dat)
        plt5 = plot_price_path(s,r,1,Dat,true)
    end
    
    plt2 = plot_density_visual(s, r, 1, Dat)
    if length(Dat[1])>1
        plt4 = plot_density_visual(s, r, 2, Dat)
        plt6 = plot_density_visual(s, r, 1, Dat, true, false)
    end
    
    if length(Dat[1])>1
        plot(plt1, plt2, plt3, plt4, plt5, plt6 ,layout=l,size=(1000,1000))
    else
        plot(plt1,plt2,size=(1000,400))
    end
    
    next!(p_outer)
end

gif(anim, "/tmp/LOB.gif", fps=20*length(range)/200)
#gif(anim, "~/Desktop/Masters/StillWorking/Random_Walk_And_Coupling_For_Alpha_Is_0.7.gif", fps=20*length(range)/200)
# -

# # RANDOM KICKS PRICE IMPACT STYLIZED FACTS NO FRACTIONAL

# ## Plot price impacts

# +
# Configuration Arguments
num_paths = 50#30

L = 200     # real system width (e.g. 200 meters)
M = 400     # divided into M pieces , 400

p₀ = 230.0  #this is the mid_price at t=0  238.75 

# Free-Parameters for gaussian version
D = 0.5 # real diffusion constant e.g. D=1 (meters^2 / second), 1
α = 0.0 # legacy, no longer used

ν = 12.0 #removal rate
γ = 1.0 #fraction of derivative (1 is normal diffusion, less than 1 is D^{1-γ} derivative on the RHS)

# Source term:
λ = 1.0 #
μ = 0.1 #
mySourceTerm = SourceTerm(λ, μ, true);

# Coupling term:
a = 13.0  #gap between stocks before at full strength: strong is 0.3
b = 1.0   #weighting of interaction term: strong is 2
c = 1.2   #skew factor: strong is 2

myCouplingTerm = CouplingTerm(μ, a, b, c, false);

# My randomness term
σ = 1.0 #variance in randomness
r = 0.5 #proportion of time in which it jumps left or right
myRandomnessTerm = RandomnessTerm(σ,r,true)

Δx = L / M  # real gap between simulation points 
Δt = (r * (Δx^2) / (2.0 * D))^(1/γ)

# RL Stuff:
T = 10
RealStartTime = 8 # when, in real time, to kick the system
SimStartTime = to_simulation_time(RealStartTime,Δt)-2 # convert to simulation time
Position = 0
Volume = 20

myRLPusher = RLPushTerm(SimStartTime,SimStartTime+1,Position,Volume,true)

lob_model = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ, 
    mySourceTerm, myCouplingTerm, myRLPusher, myRandomnessTerm);

print((Δt,to_simulation_time(T,Δt),num_paths*to_simulation_time(T,Δt))) #about 2GB RAM per 100K, i.e. can only do about 1.8 million
lob_model.SK_DP
# -

# check something actually happens for one example
if true
    try clear_double_dict(Dat) catch e print("Not initialized") end
    GC.gc()

    Dat = InteractOrderBooks([lob_model], -1, true) ;

    how_many = 12
    p_arr1 = Array{Plots.Plot{Plots.GRBackend},1}(undef,how_many)
    for i in [1:how_many;]
        p_arr1[i] = plot_density_visual(SimStartTime-2+i, to_real_time(SimStartTime-2+i,Δt), 1, Dat,false, true, 10)
    end
    plot(p_arr1...,size=(1000,1000))
    #png("/home/derickdiana/Desktop/Masters/OfficialOne/DraftOnePics/SinglePointDiffusionStepByStep.png")
    plot!()
end

# +
volumes = [1:50;] 
volumes = volumes.*2

vi_len = length(volumes)
mean_price_impacts = ones(Float64,vi_len)
var_price_impacts = ones(Float64,vi_len)

p_outer = Progress(vi_len,dt=0.1)

l = 1#to_simulation_time(1,Δt)

price_impact = zeros(Float64,num_paths)

for Volume in 1:vi_len
    myRLPusher = RLPushTerm(SimStartTime,SimStartTime+l,Position,volumes[Volume],true)
    
    lob_model = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
        mySourceTerm, myCouplingTerm, myRLPusher,myRandomnessTerm);

    try clear_double_dict(Dat) catch e print("Not initialized") end
    GC.gc()

    Dat = InteractOrderBooks([lob_model,lob_model], -1, false) ;
    
    price_impact = zeros(Float64,num_paths)
    for path_num in 1:num_paths
        price_impact[path_num] = Dat[path_num][1].raw_price_paths[SimStartTime+l] - Dat[path_num][1].raw_price_paths[SimStartTime]
    end
    
    mean_price_impacts[Volume] = mean(price_impact)
    var_price_impacts[Volume] = var(price_impact)
    
    next!(p_outer)
end
# -

mean_price_impacts = .- mean_price_impacts;

# +
x = volumes
a,b = log_fit(x,mean_price_impacts)
#c,d = power_fit(x,mean_price_impacts)

#scatter(x,mean_price_impacts,label="data: p(t+1)-p(t)",ms=3, ma=1,yerr=var_price_impacts)
scatter(x,mean_price_impacts,label="data: p(t+1)-p(t)",ms=3, ma=1)
plot!(x,mean_price_impacts,ribbon=var_price_impacts,alpha=0,fillalpha=0.4,fillcolor="blue",label="")
plot!(x,a.+b.*log.(x),label="Log fit",w=1.5,color="red")
#plot!(x,c.*((x).^d),label="Power fit",w=1.5,color="green")
plot!(xlabel="Volume",ylabel="Price impact i.e. p(t+1)-p(t)")
#png("/home/derickdiana/Desktop/Masters/OfficialOne/DraftOnePics/PriceImpactRandomKicksNoFractional.png")
plot!()
# -
# ## Plot stylised facts

# +
# Configuration Arguments
num_paths = 1#30

L = 200     # real system width (e.g. 200 meters)
M = 400     # divided into M pieces , 400

T = 5000  # simulation runs until real time T (e.g. 80 seconds)
p₀ = 230.0  #this is the mid_price at t=0  238.75 

# Free-Parameters for gaussian version
D = 0.5 # real diffusion constant e.g. D=1 (meters^2 / second), 1
α = 0.0 # legacy, no longer used

ν = 10.0 #removal rate
γ = 1.0 #fraction of derivative (1 is normal diffusion, less than 1 is D^{1-γ} derivative on the RHS)

# Source term:
λ = 1.0 #
μ = 0.1 #

mySourceTerm = SourceTerm(λ, μ, true);

# Coupling term:
a = 13.0  #gap between stocks before at full strength: strong is 0.3
b = 1.0   #weighting of interaction term: strong is 2
c = 1.2   #skew factor: strong is 2

myCouplingTerm = CouplingTerm(μ, a, b, c, false);

# My randomness term
σ = 1.5 #variance in randomness
r = 0.5 #proportion of time in which it jumps left or right
β = 0.5 #probability of being the value of the previous lag or mean reversion strength
lag = 10 #lag
do_random_walk = true #behave like a random walk
myRandomnessTerm = RandomnessTerm(σ,r,β,lag,do_random_walk,true)

Δx = L / M  # real gap between simulation points 
Δt = (r * (Δx^2) / (2.0 * D))^(1/γ)

# RL Stuff:
RealStartTime = 2 # when, in real time, to kick the system
SimStartTime = to_simulation_time(RealStartTime,Δt)-2 # convert to simulation time
SimEndTime = SimStartTime + 10 # when to stop kicking, in simulation time
Position = 200
Volume = -8; # If position == -x where x>=0, then put it x above the mid price each time

myRLPusher = RLPushTerm(SimStartTime,SimEndTime,Position,Volume,false)

lob_model = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ, 
    mySourceTerm, myCouplingTerm, myRLPusher, myRandomnessTerm);

print((Δt,to_simulation_time(T,Δt),num_paths*to_simulation_time(T,Δt))) #about 2GB RAM per 100K, i.e. can only do about 1.8 million
lob_model.SK_DP

# +
# clear everything pointed to by the dictionary then garbage collect. If it wasn't assigned yet it will left you know
try clear_double_dict(Dat) catch e print("Not initialized") end
GC.gc()

Dat = InteractOrderBooks([lob_model], -1, true) ;
#Dat = InteractOrderBooks([lob_model¹], -1, true) ;
# -

histogram(Dat[1][1].V)

# +
observed_price_path = Dat[1][1].raw_price_paths[1:end-1];
observed_log_returns = diff(log.(observed_price_path[:,1]));

data_stylized_facts = StylizedFacts.StylizedFactsPlot(observed_price_path);

StylizedFacts.plot_all_stylized_facts(data_stylized_facts,(1000,1200))

#png("/home/derickdiana/Desktop/Masters/OfficialOne/DraftOnePics/StylizedFactsRandomKicksNoFractional.png")
plot!()
# -

# # RANDOM KICKS PRICE IMPACT STYLIZED FACTS FRACTIONAL

# ## Plot stylized facts for gamma is 0.8

# +
# Configuration Arguments
T = 10000 
num_paths = 1 

L = 200     # real system width (e.g. 200 meters)
M = 400     # divided into M pieces , 400

p₀ = 230.0  #this is the mid_price at t=0  238.75 

# Free-Parameters for gaussian version
D = 0.5 # real diffusion constant e.g. D=1 (meters^2 / second), 1
σ = 1.0 #variance in randomness
α = 0.0 # legacy, no longer used

ν = 12.0 #removal rate
γ = 0.8 #fraction of derivative (1 is normal diffusion, less than 1 is D^{1-γ} derivative on the RHS)
r = 0.5 #proportion of time in which it jumps left or right

dist = Normal(0.0,σ) #dist = TDist(1) #dist = Spl(1);

Δx = L / M  # real gap between simulation points 
#Δt = (Δx^2) / (2.0 * D) / 10 #* (2.0/3.0) # real time seperation between simulation points
Δt = (r * (Δx^2) / (2.0 * D))^(1/γ)
print((Δt,to_simulation_time(T,Δt),num_paths*to_simulation_time(T,Δt))) #about 2GB RAM per 100K, i.e. can only do about 1.8 million
lob_model¹.SK_DP

# Source term:
λ = 1.0
μ = 0.1 

mySourceTerm = SourceTerm(λ, μ, true);

# Coupling term:
myCouplingTerm = CouplingTerm(0.0, 0.0, 0.0, 0.0, false);

# RL Stuff:
RealStartTime = 6 # when, in real time, to kick the system
SimStartTime = to_simulation_time(RealStartTime,Δt)-2 # convert to simulation time
Position = 0
# Volume set below

myRLPusher1 = RLPushTerm(SimStartTime,SimEndTime,Position,Volume,false)

lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, σ, ν, α, r, γ, dist, 
    mySourceTerm, myCouplingTerm, myRLPusher1,true);

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
observed_price_path = sample_price_paths¹[1:end-1,1];
observed_log_returns = diff(log.(observed_price_path[:,1]));

data_stylized_facts = StylizedFacts.StylizedFactsPlot(observed_price_path);

StylizedFacts.plot_all_stylized_facts(data_stylized_facts,(1000,1200))

#png("/home/derickdiana/Desktop/Masters/OfficialOne/DraftOnePics/StylizedFactsRandomKicksFractionalAlphaIsPoint8.png")
plot!()
# -

# ## Price impact fixed delay and duration

# +
# Configuration Arguments
num_paths = 1

L = 200     # real system width (e.g. 200 meters)
M = 400     # divided into M pieces , 400

p₀ = 230.0  #this is the mid_price at t=0  238.75 

# Free-Parameters for gaussian version
D = 0.5 # real diffusion constant e.g. D=1 (meters^2 / second), 1
σ = 1.0 #variance in randomness
α = 0.0 # legacy, no longer used

ν = 14.0 #removal rate
γ = 0.6 #fraction of derivative (1 is normal diffusion, less than 1 is D^{1-γ} derivative on the RHS)
r = 0.5 #proportion of time in which it jumps left or right

dist = Normal(0.0,σ) #dist = TDist(1) #dist = Spl(1);

# Source term:
λ = 1.0
μ = 0.1 

mySourceTerm = SourceTerm(λ, μ, true);

# Coupling term:
myCouplingTerm = CouplingTerm(0.0, 0.0, 0.0, 0.0, false);

myRLPusher¹ = RLPushTerm(SimStartTime,SimStartTime+1,Position,-8,false)
lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, σ, ν, α, r, γ, dist,
    mySourceTerm, myCouplingTerm, myRLPusher¹,true);

myRLPusher² = RLPushTerm(SimStartTime,SimStartTime+1,Position,-8,false)
lob_model² = SLOB(num_paths, T, p₀, M, L, D, σ, ν, α, r, γ, dist,
    mySourceTerm, myCouplingTerm, myRLPusher²,true);

# -

Δx = L / M  # real gap between simulation points 
Δt = (r * (Δx^2) / (2.0 * D))^(1/γ)
print((Δt,to_simulation_time(T,Δt),num_paths*to_simulation_time(T,Δt))) #about 2GB RAM per 100K, i.e. can only do about 1.8 million
lob_model¹.SK_DP

# +
lob_densities¹, sources¹, couplings¹, rl_pushes¹, raw_price_paths¹, sample_price_paths¹, P⁺s¹, P⁻s¹, Ps¹, V¹,
lob_densities², sources², couplings², rl_pushes², raw_price_paths², sample_price_paths², P⁺s², P⁻s², Ps², V² = 
    [nothing for _ in 1:20]

GC.gc()

lob_densities¹, sources¹, couplings¹, rl_pushes¹, raw_price_paths¹, sample_price_paths¹, P⁺s¹, P⁻s¹, Ps¹, V¹,
lob_densities², sources², couplings², rl_pushes², raw_price_paths², sample_price_paths², P⁺s², P⁻s², Ps², V², broke_points =
    InteractOrderBooks(lob_model¹,lob_model², -1, true) ;

l = @layout [a d; b e; c f];
print(size(raw_price_paths²))



total_steps = 10#min(to_simulation_time(T,Δt),100)
total_length = to_simulation_time(T,Δt)
step = floor(Int,total_length/total_steps)

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

# +
volume_indices = [1:50;] 
volume_indices = volume_indices./2
gamma_indices = [1.0,0.9,0.8,0.7,0.6]

vi_len = length(volume_indices)
gi_len = length(gamma_indices)
#mean_price_impacts = [0.0 for x in volume_indices]
#var_price_impacts = [0.0 for x in volume_indices]
#mean_log_price_impacts = [0.0 for x in volume_indices]
#var_log_price_impacts = [0.0 for x in volume_indices]
mean_price_impacts_frac = ones(Float64,gi_len,vi_len)
var_price_impacts_frac = ones(Float64,gi_len,vi_len)

Δt_largest = (r * (Δx^2) / (2.0 * D))^(1/1)

for Gamma in 1:gi_len 
    p_outer = Progress(vi_len,dt=0.1)
    
    γ = gamma_indices[Gamma]
    Δt = (r * (Δx^2) / (2.0 * D))^(1/γ)

    T = 10
    RealStartTime = 8 # when, in real time, to kick the system. Starts late to give system time to relax. Must be latest of all the gammas to try
    SimStartTime = to_simulation_time(RealStartTime,Δt)-2 # convert to simulation time
    
    print(2)
    
    for Volume in 1:vi_len
        volume = volume_indices[Volume]
        
        myRLPusher¹ = RLPushTerm(SimStartTime,SimStartTime+2,Position,volume,true)
        lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, σ, ν, α, r, γ, dist,
            mySourceTerm, myCouplingTerm, myRLPusher¹,true);

        myRLPusher² = RLPushTerm(SimStartTime,SimStartTime+2,Position,volume,true)
        lob_model² = SLOB(num_paths, T, p₀, M, L, D, σ, ν, α, r, γ, dist,
            mySourceTerm, myCouplingTerm, myRLPusher²,true);
        

        lob_densities¹, sources¹, couplings¹, rl_pushes¹, raw_price_paths¹, sample_price_paths¹, P⁺s¹, P⁻s¹, Ps¹, V¹,
        lob_densities², sources², couplings², rl_pushes², raw_price_paths², sample_price_paths², P⁺s², P⁻s², Ps², V² = 
            [nothing for _ in 1:20]

        GC.gc()

        lob_densities¹, sources¹, couplings¹, rl_pushes¹, raw_price_paths¹, sample_price_paths¹, P⁺s¹, P⁻s¹, Ps¹, V¹,
        lob_densities², sources², couplings², rl_pushes², raw_price_paths², sample_price_paths², P⁺s², P⁻s², Ps², V² =
        InteractOrderBooks(lob_model¹,lob_model², -1, false) ;

        #price_impact¹ = sample_price_paths¹[RealStartTime+1,:] - sample_price_paths¹[RealStartTime,:]
        #price_impact² = sample_price_paths²[RealStartTime+1,:] - sample_price_paths²[RealStartTime,:]
        #price_impact = cat(price_impact¹,price_impact²,dims=1)
        #print(size(raw_price_paths¹))
        price_impact¹ = raw_price_paths¹[SimStartTime+2,:] - raw_price_paths¹[SimStartTime,:]
        price_impact² = raw_price_paths²[SimStartTime+2,:] - raw_price_paths²[SimStartTime,:]
        price_impact = cat(price_impact¹,price_impact²,dims=1)
        mean_price_impacts_frac[Gamma,Volume] = mean(price_impact)
        var_price_impacts_frac[Gamma,Volume] = var(price_impact)


        #log_price_impact¹ = log.(sample_price_paths¹[RealStartTime,:]) - log.(sample_price_paths¹[RealStartTime-1,:])
        #log_price_impact² = log.(sample_price_paths²[RealStartTime,:]) - log.(sample_price_paths²[RealStartTime-1,:])
        #log_price_impact  = cat(log_price_impact¹,log_price_impact²,dims=1)

        #mean_log_price_impacts[Volume] = mean(log_price_impact)
        #var_log_price_impacts[Volume] = var(log_price_impact)

        next!(p_outer)
    end
end
# -

mean_price_impacts_frac = .- mean_price_impacts_frac;

# +
x = volume_indices
colors = ["blue","red","green","orange","purple"]
plot()
for Gamma in 1:gi_len
    a,b = log_fit(x,mean_price_impacts_frac[Gamma,:])
    c,d = power_fit(x,mean_price_impacts_frac[Gamma,:])

    #scatter(x,mean_price_impacts,label="data: p(t+1)-p(t)",ms=3, ma=1,yerr=var_price_impacts)
    scatter!(x,mean_price_impacts_frac[Gamma,:],label=string("data: gamma = ",gamma_indices[Gamma]),ms=3, ma=1,color=colors[Gamma])
    plot!(x,mean_price_impacts_frac[Gamma,:],ribbon=var_price_impacts_frac[Gamma,:],alpha=0,fillalpha=0.4,fillcolor=colors[Gamma],label="")
    plot!(x,a.+b.*log.(x),label="Log fit",w=1.5,color=colors[Gamma])
    plot!(x,c.*((x).^d),label="Power fit",w=1.5,color=colors[Gamma],linestyle=:dash)
    plot!(xlabel="Volume",ylabel="Price impact i.e. p(t+1)-p(t)")
    #png("/home/derickdiana/Desktop/Masters/OfficialOne/DraftOnePics/PriceImpactRandomKicksFractional0.7.png")
end

plot!(size=(1000,1000))
png("/home/derickdiana/Desktop/Masters/OfficialOne/DraftOnePics/PriceImpactFractionalFixedDelayAndDuration.png")
plot!(size=(1000,1000))
# -
# ## Price impact adaptive delay and duration

# +
volume_indices = [1:50;] 
volume_indices = volume_indices./2
gamma_indices = [1.0,0.9,0.8,0.7,0.6]

vi_len = length(volume_indices)
gi_len = length(gamma_indices)
#mean_price_impacts = [0.0 for x in volume_indices]
#var_price_impacts = [0.0 for x in volume_indices]
#mean_log_price_impacts = [0.0 for x in volume_indices]
#var_log_price_impacts = [0.0 for x in volume_indices]
mean_price_impacts_frac = ones(Float64,gi_len,vi_len)
var_price_impacts_frac = ones(Float64,gi_len,vi_len)

Δt_largest = (r * (Δx^2) / (2.0 * D))^(1/1)

for Gamma in 1:gi_len 
    p_outer = Progress(vi_len,dt=0.1)
    
    γ = gamma_indices[Gamma]
    Δt = (r * (Δx^2) / (2.0 * D))^(1/γ)

    T = 10
    RealStartTime = 8 # when, in real time, to kick the system. Starts late to give system time to relax. Must be latest of all the gammas to try
    SimStartTime = to_simulation_time(RealStartTime,Δt)-2 # convert to simulation time
    
    l = to_simulation_time(Δt_largest,Δt)
    print(l)
    
    for Volume in 1:vi_len
        volume = volume_indices[Volume]
        
        myRLPusher¹ = RLPushTerm(SimStartTime,SimStartTime+l,Position,volume,true)
        lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, σ, ν, α, r, γ, dist,
            mySourceTerm, myCouplingTerm, myRLPusher¹,true);

        myRLPusher² = RLPushTerm(SimStartTime,SimStartTime+l,Position,volume,true)
        lob_model² = SLOB(num_paths, T, p₀, M, L, D, σ, ν, α, r, γ, dist,
            mySourceTerm, myCouplingTerm, myRLPusher²,true);
        

        lob_densities¹, sources¹, couplings¹, rl_pushes¹, raw_price_paths¹, sample_price_paths¹, P⁺s¹, P⁻s¹, Ps¹, V¹,
        lob_densities², sources², couplings², rl_pushes², raw_price_paths², sample_price_paths², P⁺s², P⁻s², Ps², V² = 
            [nothing for _ in 1:20]

        GC.gc()

        lob_densities¹, sources¹, couplings¹, rl_pushes¹, raw_price_paths¹, sample_price_paths¹, P⁺s¹, P⁻s¹, Ps¹, V¹,
        lob_densities², sources², couplings², rl_pushes², raw_price_paths², sample_price_paths², P⁺s², P⁻s², Ps², V² =
        InteractOrderBooks(lob_model¹,lob_model², -1, false) ;

        #price_impact¹ = sample_price_paths¹[RealStartTime+1,:] - sample_price_paths¹[RealStartTime,:]
        #price_impact² = sample_price_paths²[RealStartTime+1,:] - sample_price_paths²[RealStartTime,:]
        #price_impact = cat(price_impact¹,price_impact²,dims=1)
        #print(size(raw_price_paths¹))
        price_impact¹ = raw_price_paths¹[SimStartTime+l,:] - raw_price_paths¹[SimStartTime,:]
        price_impact² = raw_price_paths²[SimStartTime+l,:] - raw_price_paths²[SimStartTime,:]
        price_impact = cat(price_impact¹,price_impact²,dims=1)
        mean_price_impacts_frac[Gamma,Volume] = mean(price_impact)
        var_price_impacts_frac[Gamma,Volume] = var(price_impact)


        #log_price_impact¹ = log.(sample_price_paths¹[RealStartTime,:]) - log.(sample_price_paths¹[RealStartTime-1,:])
        #log_price_impact² = log.(sample_price_paths²[RealStartTime,:]) - log.(sample_price_paths²[RealStartTime-1,:])
        #log_price_impact  = cat(log_price_impact¹,log_price_impact²,dims=1)

        #mean_log_price_impacts[Volume] = mean(log_price_impact)
        #var_log_price_impacts[Volume] = var(log_price_impact)

        next!(p_outer)
    end
end
# -

mean_price_impacts_frac = .- mean_price_impacts_frac;

# +
x = volume_indices
colors = ["blue","red","green","orange","purple"]
plot()
for Gamma in 1:gi_len
    a,b = log_fit(x,mean_price_impacts_frac[Gamma,:])
    c,d = power_fit(x,mean_price_impacts_frac[Gamma,:])

    #scatter(x,mean_price_impacts,label="data: p(t+1)-p(t)",ms=3, ma=1,yerr=var_price_impacts)
    scatter!(x,mean_price_impacts_frac[Gamma,:],label=string("data: gamma = ",gamma_indices[Gamma]),ms=3, ma=1,color=colors[Gamma])
    plot!(x,mean_price_impacts_frac[Gamma,:],ribbon=var_price_impacts_frac[Gamma,:],alpha=0,fillalpha=0.4,fillcolor=colors[Gamma],label="")
    plot!(x,a.+b.*log.(x),label="Log fit",w=1.5,color=colors[Gamma])
    plot!(x,c.*((x).^d),label="Power fit",w=1.5,color=colors[Gamma],linestyle=:dash)
    plot!(xlabel="Volume",ylabel="Price impact i.e. p(t+1)-p(t)")
    #png("/home/derickdiana/Desktop/Masters/OfficialOne/DraftOnePics/PriceImpactRandomKicksFractional0.7.png")
end

plot!(size=(1000,1000))
png("/home/derickdiana/Desktop/Masters/OfficialOne/DraftOnePics/PriceImpactFractionalAdaptiveDelayAndDuration.png")
plot!(size=(1000,1000))
# -
# ## Price impact adaptive delay only no random component

# +
volume_indices = [1:200;] 
volume_indices = volume_indices./2
gamma_indices = [1.0,0.9,0.8,0.7,0.6]

vi_len = length(volume_indices)
gi_len = length(gamma_indices)
#mean_price_impacts = [0.0 for x in volume_indices]
#var_price_impacts = [0.0 for x in volume_indices]
#mean_log_price_impacts = [0.0 for x in volume_indices]
#var_log_price_impacts = [0.0 for x in volume_indices]
mean_price_impacts_frac = ones(Float64,gi_len,vi_len)
var_price_impacts_frac = ones(Float64,gi_len,vi_len)

Δt_largest = (r * (Δx^2) / (2.0 * D))^(1/1)

for Gamma in 1:gi_len 
    p_outer = Progress(vi_len,dt=0.1)
    
    γ = gamma_indices[Gamma]
    Δt = (r * (Δx^2) / (2.0 * D))^(1/γ)

    T = 10
    RealStartTime = 8 # when, in real time, to kick the system. Starts late to give system time to relax. Must be latest of all the gammas to try
    SimStartTime = to_simulation_time(RealStartTime,Δt)-2 # convert to simulation time
    
    l = to_simulation_time(Δt_largest,Δt)
    print(l)
    
    for Volume in 1:vi_len
        volume = volume_indices[Volume]
        
        myRLPusher¹ = RLPushTerm(SimStartTime,SimStartTime+1,Position,volume,true)
        lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, σ, ν, α, r, γ, dist,
            mySourceTerm, myCouplingTerm, myRLPusher¹,false);

        myRLPusher² = RLPushTerm(SimStartTime,SimStartTime+1,Position,volume,true)
        lob_model² = SLOB(num_paths, T, p₀, M, L, D, σ, ν, α, r, γ, dist,
            mySourceTerm, myCouplingTerm, myRLPusher²,false);
        

        lob_densities¹, sources¹, couplings¹, rl_pushes¹, raw_price_paths¹, sample_price_paths¹, P⁺s¹, P⁻s¹, Ps¹, V¹,
        lob_densities², sources², couplings², rl_pushes², raw_price_paths², sample_price_paths², P⁺s², P⁻s², Ps², V² = 
            [nothing for _ in 1:20]

        GC.gc()

        lob_densities¹, sources¹, couplings¹, rl_pushes¹, raw_price_paths¹, sample_price_paths¹, P⁺s¹, P⁻s¹, Ps¹, V¹,
        lob_densities², sources², couplings², rl_pushes², raw_price_paths², sample_price_paths², P⁺s², P⁻s², Ps², V² =
        InteractOrderBooks(lob_model¹,lob_model², -1, false) ;

        #price_impact¹ = sample_price_paths¹[RealStartTime+1,:] - sample_price_paths¹[RealStartTime,:]
        #price_impact² = sample_price_paths²[RealStartTime+1,:] - sample_price_paths²[RealStartTime,:]
        #price_impact = cat(price_impact¹,price_impact²,dims=1)
        #print(size(raw_price_paths¹))
        price_impact¹ = raw_price_paths¹[SimStartTime+l,:] - raw_price_paths¹[SimStartTime,:]
        price_impact² = raw_price_paths²[SimStartTime+l,:] - raw_price_paths²[SimStartTime,:]
        price_impact = cat(price_impact¹,price_impact²,dims=1)
        mean_price_impacts_frac[Gamma,Volume] = mean(price_impact)
        var_price_impacts_frac[Gamma,Volume] = var(price_impact)


        #log_price_impact¹ = log.(sample_price_paths¹[RealStartTime,:]) - log.(sample_price_paths¹[RealStartTime-1,:])
        #log_price_impact² = log.(sample_price_paths²[RealStartTime,:]) - log.(sample_price_paths²[RealStartTime-1,:])
        #log_price_impact  = cat(log_price_impact¹,log_price_impact²,dims=1)

        #mean_log_price_impacts[Volume] = mean(log_price_impact)
        #var_log_price_impacts[Volume] = var(log_price_impact)

        next!(p_outer)
    end
end
# -

mean_price_impacts_frac = .- mean_price_impacts_frac;

# +
x = volume_indices
colors = ["blue","red","green","orange","purple"]
plot()
for Gamma in 1:gi_len
    a,b = log_fit(x,mean_price_impacts_frac[Gamma,:])
    #c,d = power_fit(x,mean_price_impacts_frac[Gamma,:])

    scatter!(x,mean_price_impacts_frac[Gamma,:],label=string("data: gamma = ",gamma_indices[Gamma]),ms=2,markerstrokewidth=0.3,  ma=1,color=colors[Gamma])
    plot!(x,mean_price_impacts_frac[Gamma,:],ribbon=var_price_impacts_frac[Gamma,:],alpha=0,fillalpha=0.4,fillcolor=colors[Gamma],label="")
    plot!(x,a.+b.*log.(x),label=string("Log fit: ",round(a,digits=2)," + ",round(b,digits=2),"log(x)"),w=1.5,color=colors[Gamma])
    #plot!(x,c.*((x).^d),label="Power fit",w=1.5,color=colors[Gamma],linestyle=:dash)
    plot!(xlabel="Volume",ylabel="Price impact i.e. p(t+1)-p(t)")
end

plot!(size=(1000,1000))
#png("/home/derickdiana/Desktop/Masters/OfficialOne/DraftOnePics/PriceImpactFractionalAdaptiveDelayOnlyNoRandom.png")
plot!(size=(1000,1000))
# -
# ## Price impact adaptive delay only

# +
volume_indices = [1:200;] 
volume_indices = volume_indices./2
gamma_indices = [1.0,0.9,0.8,0.7,0.6]

vi_len = length(volume_indices)
gi_len = length(gamma_indices)
#mean_price_impacts = [0.0 for x in volume_indices]
#var_price_impacts = [0.0 for x in volume_indices]
#mean_log_price_impacts = [0.0 for x in volume_indices]
#var_log_price_impacts = [0.0 for x in volume_indices]
mean_price_impacts_frac = ones(Float64,gi_len,vi_len)
var_price_impacts_frac = ones(Float64,gi_len,vi_len)

Δt_largest = (r * (Δx^2) / (2.0 * D))^(1/1)

for Gamma in 1:gi_len 
    p_outer = Progress(vi_len,dt=0.1)
    
    γ = gamma_indices[Gamma]
    Δt = (r * (Δx^2) / (2.0 * D))^(1/γ)

    T = 10
    RealStartTime = 8 # when, in real time, to kick the system. Starts late to give system time to relax. Must be latest of all the gammas to try
    SimStartTime = to_simulation_time(RealStartTime,Δt)-2 # convert to simulation time
    
    l = to_simulation_time(Δt_largest,Δt)
    print(l)
    
    for Volume in 1:vi_len
        volume = volume_indices[Volume]
        
        myRLPusher¹ = RLPushTerm(SimStartTime,SimStartTime+l,Position,volume,true)
        lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, σ, ν, α, r, γ, dist,
            mySourceTerm, myCouplingTerm, myRLPusher¹,false);

        myRLPusher² = RLPushTerm(SimStartTime,SimStartTime+l,Position,volume,true)
        lob_model² = SLOB(num_paths, T, p₀, M, L, D, σ, ν, α, r, γ, dist,
            mySourceTerm, myCouplingTerm, myRLPusher²,false);
        

        lob_densities¹, sources¹, couplings¹, rl_pushes¹, raw_price_paths¹, sample_price_paths¹, P⁺s¹, P⁻s¹, Ps¹, V¹,
        lob_densities², sources², couplings², rl_pushes², raw_price_paths², sample_price_paths², P⁺s², P⁻s², Ps², V² = 
            [nothing for _ in 1:20]

        GC.gc()

        lob_densities¹, sources¹, couplings¹, rl_pushes¹, raw_price_paths¹, sample_price_paths¹, P⁺s¹, P⁻s¹, Ps¹, V¹,
        lob_densities², sources², couplings², rl_pushes², raw_price_paths², sample_price_paths², P⁺s², P⁻s², Ps², V² =
        InteractOrderBooks(lob_model¹,lob_model², -1, false) ;

        #price_impact¹ = sample_price_paths¹[RealStartTime+1,:] - sample_price_paths¹[RealStartTime,:]
        #price_impact² = sample_price_paths²[RealStartTime+1,:] - sample_price_paths²[RealStartTime,:]
        #price_impact = cat(price_impact¹,price_impact²,dims=1)
        #print(size(raw_price_paths¹))
        price_impact¹ = raw_price_paths¹[SimStartTime+1,:] - raw_price_paths¹[SimStartTime,:]
        price_impact² = raw_price_paths²[SimStartTime+1,:] - raw_price_paths²[SimStartTime,:]
        price_impact = cat(price_impact¹,price_impact²,dims=1)
        mean_price_impacts_frac[Gamma,Volume] = mean(price_impact)
        var_price_impacts_frac[Gamma,Volume] = var(price_impact)


        #log_price_impact¹ = log.(sample_price_paths¹[RealStartTime,:]) - log.(sample_price_paths¹[RealStartTime-1,:])
        #log_price_impact² = log.(sample_price_paths²[RealStartTime,:]) - log.(sample_price_paths²[RealStartTime-1,:])
        #log_price_impact  = cat(log_price_impact¹,log_price_impact²,dims=1)

        #mean_log_price_impacts[Volume] = mean(log_price_impact)
        #var_log_price_impacts[Volume] = var(log_price_impact)

        next!(p_outer)
    end
end
# -

mean_price_impacts_frac = .- mean_price_impacts_frac;

# +
x = volume_indices
colors = ["blue","red","green","orange","purple"]
plot()
for Gamma in 1:gi_len
    a,b = log_fit(x,mean_price_impacts_frac[Gamma,:])
    c,d = power_fit(x,mean_price_impacts_frac[Gamma,:])

    #scatter(x,mean_price_impacts,label="data: p(t+1)-p(t)",ms=3, ma=1,yerr=var_price_impacts)
    scatter!(x,mean_price_impacts_frac[Gamma,:],label=string("data: gamma = ",gamma_indices[Gamma]),ms=3, ma=1,color=colors[Gamma])
    plot!(x,mean_price_impacts_frac[Gamma,:],ribbon=var_price_impacts_frac[Gamma,:],alpha=0,fillalpha=0.4,fillcolor=colors[Gamma],label="")
    plot!(x,a.+b.*log.(x),label="Log fit",w=1.5,color=colors[Gamma])
    plot!(x,c.*((x).^d),label="Power fit",w=1.5,color=colors[Gamma],linestyle=:dash)
    plot!(xlabel="Volume",ylabel="Price impact i.e. p(t+1)-p(t)")
    #png("/home/derickdiana/Desktop/Masters/OfficialOne/DraftOnePics/PriceImpactRandomKicksFractional0.7.png")
end

plot!(size=(1000,1000))
#png("/home/derickdiana/Desktop/Masters/OfficialOne/DraftOnePics/PriceImpactFractionalAdaptiveDelayOnly.png")
plot!(size=(1000,1000))
# -
mean_price_impacts

# ## Stylised facts 

# +
# Configuration Arguments
T = 10000 
num_paths = 1 

L = 200     # real system width (e.g. 200 meters)
M = 400     # divided into M pieces , 400

p₀ = 230.0  #this is the mid_price at t=0  238.75 

# Free-Parameters for gaussian version
D = 0.5 # real diffusion constant e.g. D=1 (meters^2 / second), 1
σ = 1.0 #variance in randomness
α = 0.0 # legacy, no longer used

ν = 6.0 #removal rate
γ = 1.0 #fraction of derivative (1 is normal diffusion, less than 1 is D^{1-γ} derivative on the RHS)
r = 0.5 #proportion of time in which it jumps left or right

dist = Normal(0.0,σ) #dist = TDist(1) #dist = Spl(1);

# Source term:
λ = 1.0
μ = 0.1 

mySourceTerm = SourceTerm(λ, μ, true);

# Coupling term:
myCouplingTerm = CouplingTerm(0.0, 0.0, 0.0, 0.0, false);

# RL Stuff:
RealStartTime = 6 # when, in real time, to kick the system
SimStartTime = to_simulation_time(RealStartTime,Δt)-2 # convert to simulation time
Position = 0
# Volume set below

myRLPusher1 = RLPushTerm(SimStartTime,SimEndTime,Position,Volume,false)

lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, σ, ν, α, r, γ, dist, 
    mySourceTerm, myCouplingTerm, myRLPusher1,true);

myRLPusher2 = RLPushTerm(SimStartTime,SimEndTime,Position,Volume,false)

lob_model² = SLOB(num_paths, T, p₀, M, L, D, σ, ν, α, r, γ, dist,
    mySourceTerm, myCouplingTerm, myRLPusher2,true);
# -

Δx = L / M  # real gap between simulation points 
#Δt = (Δx^2) / (2.0 * D) / 10 #* (2.0/3.0) # real time seperation between simulation points
Δt = (r * (Δx^2) / (2.0 * D))^(1/γ)
print((Δt,to_simulation_time(T,Δt),num_paths*to_simulation_time(T,Δt))) #about 2GB RAM per 100K, i.e. can only do about 1.8 million
lob_model¹.SK_DP

# +
lob_densities¹, sources¹, couplings¹, rl_pushes¹, raw_price_paths¹, sample_price_paths¹, P⁺s¹, P⁻s¹, Ps¹, V¹,
lob_densities², sources², couplings², rl_pushes², raw_price_paths², sample_price_paths², P⁺s², P⁻s², Ps², V² = 
    [nothing for _ in 1:20]

GC.gc()

lob_densities¹, sources¹, couplings¹, rl_pushes¹, raw_price_paths¹, sample_price_paths¹, P⁺s¹, P⁻s¹, Ps¹, V¹,
lob_densities², sources², couplings², rl_pushes², raw_price_paths², sample_price_paths², P⁺s², P⁻s², Ps², V², broke_points =
    InteractOrderBooks(lob_model¹,lob_model², -1, true) ;

# +
observed_price_path = sample_price_paths¹[1:end-1,1];
observed_log_returns = diff(log.(observed_price_path[:,1]));

data_stylized_facts = StylizedFacts.StylizedFactsPlot(observed_price_path);

StylizedFacts.plot_all_stylized_facts(data_stylized_facts,(1000,1200))

#png("/home/derickdiana/Desktop/Masters/OfficialOne/DraftOnePics/StylizedFactsRandomKicksNoFractional.png")
plot!()
# -

# # SINGLE POINT DIFFUSION

# ## Plot spike variance for normal diffusion

# +
# Configuration Arguments
num_paths = 1#30

L = 200     # real system width (e.g. 200 meters)
M = 400     # divided into M pieces , 400

p₀ = 230.0  #this is the mid_price at t=0  238.75 

# Free-Parameters for gaussian version
D = 0.5 # real diffusion constant e.g. D=1 (meters^2 / second), 1
α = 0.0 # legacy, no longer used

ν = 6.0 #removal rate
γ = 1.0 #fraction of derivative (1 is normal diffusion, less than 1 is D^{1-γ} derivative on the RHS)

# Source term:

mySourceTerm = SourceTerm(0.0, 0.0, false);

myCouplingTerm = CouplingTerm(0.0, 0.0, 0.0, 0.0, false);

# My randomness term
σ = 6.0 #variance in randomness
r = 0.5 #proportion of time in which it jumps left or right

myRandomnessTerm = RandomnessTerm(σ,r,false)

Δx = L / M  # real gap between simulation points 
Δt = (r * (Δx^2) / (2.0 * D))^(1/γ)

# RL Stuff:
T = 10
RealStartTime = 8 # when, in real time, to kick the system
SimStartTime = to_simulation_time(RealStartTime,Δt)-2 # convert to simulation time
Position = Int(M/2)
Volume = -80; # If position == -x where x>=0, then put it x above the mid price each time

myRLPusher = RLPushTerm(SimStartTime,SimStartTime+1,Position,Volume,true)

lob_model = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ, 
    mySourceTerm, myCouplingTerm, myRLPusher, myRandomnessTerm);

print((Δt,to_simulation_time(T,Δt),num_paths*to_simulation_time(T,Δt))) #about 2GB RAM per 100K, i.e. can only do about 1.8 million
lob_model¹.SK_DP

# +
# clear everything pointed to by the dictionary then garbage collect. If it wasn't assigned yet it will left you know
try clear_double_dict(Dat) catch e print("Not initialized") end
GC.gc()

Dat = InteractOrderBooks([lob_model], -1, true) ;
#Dat = InteractOrderBooks([lob_model¹], -1, true) ;
# -

how_many = 12
p_arr1 = Array{Plots.Plot{Plots.GRBackend},1}(undef,how_many)
for i in [1:how_many;]
    p_arr1[i] = plot_density_visual(SimStartTime-2+i, to_real_time(SimStartTime-2+i,Δt), 1, Dat,false, true, 10)
end
plot(p_arr1...,size=(1000,1000))
#png("/home/derickdiana/Desktop/Masters/OfficialOne/DraftOnePics/SinglePointDiffusionStepByStep.png")
plot!()

# +
observed_price_path = sample_price_paths¹[1:end-1,1];
#observed_log_returns = diff(log.(observed_price_path[:,1]));

#data_stylized_facts = StylizedFacts.StylizedFactsPlot(observed_price_path);

# +
function myvar(p,x,Δx)
    mu = mymean(p,x,Δx)
    Ex2 = 0
    for i in 1:length(p)
        Ex2 += p[i]*x[i]^2
    end
    return Ex2 - mu^2
end

function mymean(p,x,Δx)
    sum = 0
    for i in 1:length(p)
        sum += p[i]*x[i]
    end
    return sum
end

test = lob_densities¹[:,end-5,1]
test = test/sum(test)
t_mean = mymean(test ,lob_model¹.x,lob_model¹.Δx)
t_var = myvar(test ,lob_model¹.x,lob_model¹.Δx)
print(t_mean," ",t_var^0.5)
plot(lob_model¹.x,test)
plot!([t_mean],seriestype="vline",color="red")
plot!([t_mean - t_var^0.5],seriestype="vline",color="green")
plot!([t_mean + t_var^0.5],seriestype="vline",color="green")

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
sub = SimStartTime:length(sums)-1
plt1 = plot(variances[sub],label="Variance")
plt2 = plot(sums[sub],label=string("Area underneath has range ",maximum(sums[sub]) - minimum(sums[sub])),color="red",ylim=[0.99999,1.00001])
plot(plt1, plt2, layout=l2, size=(700,700),xlab="Time")
png("/home/derickdiana/Desktop/Masters/OfficialOne/DraftOnePics/SinglePointVarianceAlphaIs.png")
plot!()
# -

# ## Plot spike variance for various gamma

# +
gamma_indices = [1.0,0.9,0.8,0.7]
gi_len = length(gamma_indices)
T = 100
l = to_simulation_time(T,Δt)#length(lob_densities¹[1,:,1]) # number of time steps
sums = ones(Float64,gi_len,l-1)
variances = ones(Float64,gi_len,l-1)

for Gamma in 1:gi_len
    # Configuration Arguments
    num_paths = 1

    L = 200     # real system width (e.g. 200 meters)
    M = 400     # divided into M pieces , 400

    p₀ = 230.0  #this is the mid_price at t=0  238.75 

    # Free-Parameters for gaussian version
    D = 0.5 # real diffusion constant e.g. D=1 (meters^2 / second), 1
    σ = 1.0 #variance in randomness
    α = 0.0 # legacy, no longer used

    ν = 12.0#14.0 #removal rate
    γ = gamma_indices[Gamma] #fraction of derivative (1 is normal diffusion, less than 1 is D^{1-γ} derivative on the RHS)
    r = 0.5 #proportion of time in which it jumps left or right

    dist = Normal(0.0,σ) #dist = TDist(1) #dist = Spl(1);

    Δx = L / M  # real gap between simulation points 
    #Δt = (Δx^2) / (2.0 * D) / 10 #* (2.0/3.0) # real time seperation between simulation points
    Δt = (r * (Δx^2) / (2.0 * D))^(1/γ)

    RealStartTime = 1 # when, in real time, to kick the system. Starts late to give system time to relax. Must be latest of all the gammas to try
    SimStartTime = to_simulation_time(RealStartTime,Δt)-2 # convert to simulation time

    mySourceTerm = SourceTerm(0.0, 0.0, false);

    # Coupling term:
    myCouplingTerm = CouplingTerm(0.0, 0.0, 0.0, 0.0, false);
    Position = Int(M/2)

    myRLPusher¹ = RLPushTerm(SimStartTime,SimStartTime+1,Position,-8,true)
    lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, σ, ν, α, r, γ, dist,
        mySourceTerm, myCouplingTerm, myRLPusher¹,false);

    myRLPusher² = RLPushTerm(SimStartTime,SimStartTime+1,Position,-8,false)
    lob_model² = SLOB(num_paths, T, p₀, M, L, D, σ, ν, α, r, γ, dist,
        mySourceTerm, myCouplingTerm, myRLPusher²,false);
    
    lob_densities¹, sources¹, couplings¹, rl_pushes¹, raw_price_paths¹, sample_price_paths¹, P⁺s¹, P⁻s¹, Ps¹, V¹,
    lob_densities², sources², couplings², rl_pushes², raw_price_paths², sample_price_paths², P⁺s², P⁻s², Ps², V² = 
        [nothing for _ in 1:20]

    GC.gc()

    lob_densities¹, sources¹, couplings¹, rl_pushes¹, raw_price_paths¹, sample_price_paths¹, P⁺s¹, P⁻s¹, Ps¹, V¹,
    lob_densities², sources², couplings², rl_pushes², raw_price_paths², sample_price_paths², P⁺s², P⁻s², Ps², V², broke_points =
        InteractOrderBooks(lob_model¹,lob_model², -1, true) ;
    

    for t in 2:l 
        temp = lob_densities¹[:,t,1]/sum(lob_densities¹[:,t,1]) #across all x values, at time "t",  on path 1
        sums[Gamma,t-1] = sum(temp)
        variances[Gamma,t-1] = myvar(temp,lob_model¹.x,lob_model¹.Δx)
    end

    
end

l2 =  @layout [a;b];
sub = SimStartTime:l-1
plot()
for Gamma in 1:gi_len
    slope = round((variances[Gamma,sub[end]]-variances[Gamma,sub[1]])/(l),digits=4)
    plot!(variances[Gamma,sub],label=string("Variance with alpha=",gamma_indices[Gamma]," has slope ",slope))
    #plt2 = plot(sums[Gamma,sub],label=string("Area underneath has range ",maximum(sums[sub]) - minimum(sums[sub])),color="red",ylim=[0.99999,1.00001])
    #plot!(plt1, plt2, layout=l2, size=(700,700),xlab="Time")
end
png("/home/derickdiana/Desktop/Masters/OfficialOne/DraftOnePics/SinglePointVarianceForChangingAlpha.png")
plot!()
# -

# # Junk

# + active=""
# Revise.revise()
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
# -



# +
lob_densities¹, sources¹, couplings¹, rl_pushes¹, raw_price_paths¹, sample_price_paths¹, P⁺s¹, P⁻s¹, Ps¹, V¹,
lob_densities², sources², couplings², rl_pushes², raw_price_paths², sample_price_paths², P⁺s², P⁻s², Ps², V² = 
    [nothing for _ in 1:20]

GC.gc()

lob_densities¹, sources¹, couplings¹, rl_pushes¹, raw_price_paths¹, sample_price_paths¹, P⁺s¹, P⁻s¹, Ps¹, V¹,
lob_densities², sources², couplings², rl_pushes², raw_price_paths², sample_price_paths², P⁺s², P⁻s², Ps², V², broke_points =
    InteractOrderBooks(lob_model¹,lob_model², -1, true) ;

l = @layout [a d; b e; c f];

total_steps = 10#min(to_sim(T),100)
total_length = to_sim(T)
step = floor(Int,total_length/total_steps)

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
