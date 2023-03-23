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
using CurveFit

using InteractingLOBs

# +
# Configuration Arguments
# T = ... set below
# num_paths = ... set below

L = 200     # real system width (e.g. 200 meters)
M = 400     # divided into M pieces  
p₀ = 238.75 #this is the mid_price at t=0 

# Free-Parameters for gaussian version
D = 1.0 # real diffusion constant e.g. D=1 (meters^2 / second)
σ = 1.0 

ν = 0.5
α_slob = 40.0
α_lob = 0.0
α = α_lob

dist = Normal(0.0,1.0);
#dist = TDist(1)
# -

Δx = L / M  # real gap between simulation points 
Δt = (Δx^2) / (2.0 * D); # real time seperation between simulation points

# +
λ = 1.0
μ = 0.1 

mySourceTerm = SourceTerm(λ, μ);
# -

myCouplingTerm = CouplingTerm(0.0,0.0,0.0,0.0,false);

# +
T = 4     # simulation runs until real time T (e.g. 80 seconds)
num_paths = 400 # number of paths per call to the simulation function

RealStartTime = 2
SimStartTime = to_simulation_time(RealStartTime,Δt) #at real time 10, kick the system
Position = 0
#Amount =  ... set below

# If position == -x where x>=0, then put it x above the mid price each time

# -

(RealStartTime,SimStartTime,to_simulation_time(T,Δt))

# +
volume_indices = [1:60;] #try volumes from 1 to volume_count

mean_price_impacts = [0.0 for x in volume_indices]
var_price_impacts = [0.0 for x in volume_indices]
mean_log_price_impacts = [0.0 for x in volume_indices]
var_log_price_impacts = [0.0 for x in volume_indices]

p_outer = Progress(volume_count,dt=0.1)
scale_volume = 0.5

for Volume in volume_indices
    myRLPusher¹ = RLPushTerm(SimStartTime,SimStartTime+1,Position,scale_volume*Volume,true)
    lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, σ, ν, α, dist,
        mySourceTerm, myCouplingTerm, myRLPusher¹);
    
    myRLPusher² = RLPushTerm(SimStartTime,SimStartTime+1,Position,scale_volume*Volume,true)
    lob_model² = SLOB(num_paths, T, p₀, M, L, D, σ, ν, α, dist,
        mySourceTerm, myCouplingTerm, myRLPusher²);
    
    lob_densities¹, sources¹, couplings¹, rl_pushes¹, raw_price_paths¹, sample_price_paths¹, P⁺s¹, P⁻s¹, Ps¹, 
    lob_densities², sources², couplings², rl_pushes², raw_price_paths², sample_price_paths², P⁺s², P⁻s², Ps² =
    InteractOrderBooks(lob_model¹,lob_model², -1, false) ;
    
    #price_impact¹ = sample_price_paths¹[RealStartTime+1,:] - sample_price_paths¹[RealStartTime,:]
    #price_impact² = sample_price_paths²[RealStartTime+1,:] - sample_price_paths²[RealStartTime,:]
    #price_impact = cat(price_impact¹,price_impact²,dims=1)
    price_impact¹ = raw_price_paths¹[SimStartTime+1,:] - raw_price_paths¹[SimStartTime,:]
    price_impact² = raw_price_paths²[SimStartTime+1,:] - raw_price_paths²[SimStartTime,:]
    price_impact = cat(price_impact¹,price_impact²,dims=1)
    mean_price_impacts[Volume-volume_indices[1]+1] = mean(price_impact)
    var_price_impacts[Volume-volume_indices[1]+1] = var(price_impact)
    
    
    #log_price_impact¹ = log.(sample_price_paths¹[RealStartTime,:]) - log.(sample_price_paths¹[RealStartTime-1,:])
    #log_price_impact² = log.(sample_price_paths²[RealStartTime,:]) - log.(sample_price_paths²[RealStartTime-1,:])
    #log_price_impact  = cat(log_price_impact¹,log_price_impact²,dims=1)
    
    #mean_log_price_impacts[Volume] = mean(log_price_impact)
    #var_log_price_impacts[Volume] = var(log_price_impact)
    
    next!(p_outer)
end

# +
x = scale_volume.*volume_indices
a,b = log_fit(x,mean_price_impacts)
#c,d = power_fit(x,mean_price_impacts)
var_price_impacts[9] = var_price_impacts[2]

scatter(x,mean_price_impacts,label="data: p(t+1)-p(t)",ms=3, ma=1, yerr=var_price_impacts)
plot!(x,a.+b.*log.(x),label="Log fit",w=2)
#plot!(x,c.*((x).^d),label="Power fit",w=2)
plot!(xlabel="Volume",ylabel="Price impact i.e. p(t+1)-p(t)")
#png("/home/derickdiana/Desktop/Masters/GoodPics/PriceDiffImpactGivenVolume.png")

# +
x = [1:volume_count;]
a,b = log_fit(x,mean_log_price_impacts)
c,d = power_fit(x,mean_log_price_impacts)

scatter(x,mean_log_price_impacts,label="data: log(p(t+1))-log(p(t))",ms=3, ma=1, yerr = var_log_price_impacts)
plot!(x,a.+b.*log.(x),label="Log fit")
plot!(x,c.*((x).^d),label="Power fit")
plot!(xlabel="Volume",ylabel="Price impact i.e. log(p(t+1))-log(p(t))")
#png("/home/derickdiana/Desktop/Masters/GoodPics/PriceLogImpactGivenVolume.png")
# -

Revise.revise()


