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
γ = 0.5
r = 0.5

dist = Normal(0.0,1.0);
#dist = TDist(1)
# -

Δx = L / M  # real gap between simulation points 
Δt = (r * (Δx^2) / (2.0 * D))^(1/γ)

# +
λ = 1.0
μ = 0.1 

mySourceTerm = SourceTerm(λ, μ, true);
# -

myCouplingTerm = CouplingTerm(0.0,0.0,0.0,0.0,false);

# +
T = 4     # simulation runs until real time T (e.g. 80 seconds)
num_paths = 200 # number of paths per call to the simulation function

RealStartTime = 2 # when in real time to kick the system
SimStartTime = to_simulation_time(RealStartTime,Δt) #when in simulation time to kick the system
Position = 0
#Amount =  ... set below

# If position == -x where x>=0, then put it x above the mid price each time

# -

(Δt,to_simulation_time(T,Δt),num_paths*to_simulation_time(T,Δt)) #about 2GB RAM per 100K, i.e. can only do about 1.8 million

# +
volume_indices = [1:20;] #try volumes from 1 to volume_count

mean_price_impacts = [0.0 for x in volume_indices]
var_price_impacts = [0.0 for x in volume_indices]
mean_log_price_impacts = [0.0 for x in volume_indices]
var_log_price_impacts = [0.0 for x in volume_indices]

p_outer = Progress(volume_indices[end],dt=0.1)
scale_volume = 1.5  #volume_indices must be counting numbers i.e. 1,2,3..., 
                    #so stretch it out with this if you need to count in, say, 10s by settings scale_volume=10

l = 1

for Volume in volume_indices
    myRLPusher¹ = RLPushTerm(SimStartTime,SimStartTime+1,Position,scale_volume*Volume,true)
    lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, σ, ν, α, r, γ, dist,
        mySourceTerm, myCouplingTerm, myRLPusher¹,false);
    
    myRLPusher² = RLPushTerm(SimStartTime,SimStartTime+1,Position,scale_volume*Volume,true)
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
    price_impact¹ = raw_price_paths¹[SimStartTime+l,:] - raw_price_paths¹[SimStartTime,:]
    price_impact² = raw_price_paths²[SimStartTime+l,:] - raw_price_paths²[SimStartTime,:]
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
# -

mean_price_impacts = .- mean_price_impacts

# +
x = scale_volume.*volume_indices
a,b = log_fit(x,mean_price_impacts)
#c,d = power_fit(x,mean_price_impacts)

scatter(x,mean_price_impacts,label="data: p(t+1)-p(t)",ms=3, ma=1,yerr=var_price_impacts)
plot!(x,a.+b.*log.(x),label="Log fit",w=2)
#plot!(x,c.*((x).^d),label="Power fit",w=2)
plot!(xlabel="Volume",ylabel="Price impact i.e. p(t+1)-p(t)")
#png("/home/derickdiana/Desktop/Masters/StillWorking/Price_Impact_For_Alpha_Is_0.7_No_Random_Component.png")
# -
png("/home/derickdiana/Desktop/Masters/StillWorking/Price_Impact_For_Alpha_Is_0.5_No_Random_Component.png")


Revise.revise()
