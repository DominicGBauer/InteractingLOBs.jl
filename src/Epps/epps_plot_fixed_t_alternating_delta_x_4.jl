# L = 200
# M = 400
# Δx = 0.5

using InteractingLOBs

include("../setup.jl")
include("./epps.jl")

num_paths = 30#30

L = 200     # real system width (e.g. 200 meters)
M = 400     # divided into M pieces , 400

T = 2000  # simulation runs until real time T (e.g. 80 seconds)
p₀ = 230.0  #this is the mid_price at t=0  238.75

# Free-Parameters for gaussian version
D = 0.5 # real diffusion constant e.g. D=1 (meters^2 / second), 1
α = 0.0 # legacy, no longer used

ν = 14.0 #removal rate
γ = 1.0 #fraction of derivative (1 is normal diffusion, less than 1 is D^{1-γ} derivative on the RHS)

# Source term:
λ = 1.0 #
μ = 0.1 #

mySourceTerm = SourceTerm(λ, μ, true);

# Coupling term:
a = 6.0  #gap between stocks before at full strength: strong is 0.3
b = 1.5   #weighting of interaction term: strong is 2
c = 1.2   #skew factor: strong is 2

myCouplingTerm = CouplingTerm(μ, a, b, c, true);

# My randomness term
σ = 1.0 #variance in randomness
r = 0.5 #proportion of time in which it jumps left or right
β = 0.0 #probability of being the value of the previous lag or mean reversion strength
lag = 10 #lag
do_random_walk = false #behave like a random walk
myRandomnessTerm = RandomnessTerm(σ,r,β,lag,do_random_walk,true)


Δx = L / M  # real gap between simulation points
Δt = 0.03

# RL Stuff:
RealStartTime = 50 # when, in real time, to kick the system
SimStartTime = to_simulation_time(RealStartTime,Δt)-2 # convert to simulation time
SimEndTime = SimStartTime + 3 # when to stop kicking, in simulation time
Position = 200
Volume = -8; # If position == -x where x>=0, then put it x above the mid price each time

myRLPusher1 = RLPushTerm(SimStartTime,SimEndTime,Position,Volume,true)

myRLPusher2 = RLPushTerm(SimStartTime,SimEndTime,Position,Volume,false)

lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
    mySourceTerm, myCouplingTerm, myRLPusher1, myRandomnessTerm);

lob_model² = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
    mySourceTerm, myCouplingTerm, myRLPusher2, myRandomnessTerm);

# total_steps = 10
# total_length = to_simulation_time(T,Δt)
# step = floor(Int,total_length/total_steps)

# range = 1:step:(total_length-step)

r = to_real_time(7113, lob_model¹.Δt)  #r is the time in real time
s = to_simulation_time(r, lob_model¹.Δt)  #s is the time in real time

Data = InteractOrderBooks([lob_model¹,lob_model²], -1, true);
# path1 = Data[1][1].raw_price_paths[1:s]
# path2 = Data[1][2].raw_price_paths[1:s]
# index_vector = 0:1.0:(size(path2)[1]-1)
# epps_data = hcat(index_vector, path1, path2)

# # uncomment to rerun
# epps_4 = Empirical(epps_data)
# save("Computed Data/EppsCorrection/Empirical_delta_x_1_over_2.png.jld", "epps", epps_4)

# ComputedResults = load("Computed Data/EppsCorrection/Empirical_delta_x_1_over_2.png.jld")
# epps_4 = ComputedResults["epps"]

# Plots
dt = collect(1:1:400)
q = quantile.(TDist(m-1), [0.975])

(average_epps_mean, average_epps_value, m) = generate_epps_plots_values(Data)

p4_x = plot(dt, average_epps_mean, legend = false,dpi=300, ribbon=(q .* std(average_epps_value, dims = 2) * 0.001), fillalpha=.15, title="Δt=0.03 Δx=1/2")

xlabel!(p4_x, L"\Delta t\textrm{[sec]}")
ylabel!(p4_x, L"\rho_{\Delta t}^{ij}")

savefig(p4_x, "Plots/Epps/Epps_delta_x_1_over_2.png")
