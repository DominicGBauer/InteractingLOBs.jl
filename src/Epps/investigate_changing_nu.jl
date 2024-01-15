# L = 100
# M = 400
# Δx = 0.375

using InteractingLOBs

include("../setup.jl")
include("./epps.jl")
include("./generate_plot.jl")

num_paths = 10#30

L = 200     # real system width (e.g. 200 meters)
M = 400     # divided into M pieces , 400

T = 2000  # simulation runs until real time T (e.g. 80 seconds)
p₀ = 230.0  #this is the mid_price at t=0  238.75

# Free-Parameters for gaussian version
D = 0.5 # real diffusion constant e.g. D=1 (meters^2 / second), 1
α = 0.0 # legacy, no longer used

ν = 3.0 #removal rate
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
myRandomnessTerm = RandomnessTerm(σ, r, β, lag, do_random_walk, true)


Δx = L / M  # real gap between simulation points
Δt = (r * (Δx^2) / (2.0 * D))^(1 / γ) #0.0625

# RL Stuff:
RealStartTime = 50 # when, in real time, to kick the system
SimStartTime = to_simulation_time(RealStartTime, Δt) - 2 # convert to simulation time
SimEndTime = SimStartTime + 3 # when to stop kicking, in simulation time
Position = 200
Volume = -8; # If position == -x where x>=0, then put it x above the mid price each time

myRLPusher1 = RLPushTerm(SimStartTime, SimEndTime, Position, Volume, true)

myRLPusher2 = RLPushTerm(SimStartTime, SimEndTime, Position, Volume, false)

lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
  mySourceTerm, myCouplingTerm, myRLPusher1, myRandomnessTerm);

lob_model² = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
  mySourceTerm, myCouplingTerm, myRLPusher2, myRandomnessTerm);

r = to_real_time(14401, lob_model¹.Δt)  #r is the time in real time
s = to_simulation_time(r, lob_model¹.Δt)  #s is the time in real time

Data = InteractOrderBooks([lob_model¹, lob_model²], -1, true);


ν = 5.0
lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
  mySourceTerm, myCouplingTerm, myRLPusher1, myRandomnessTerm);

lob_model² = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
  mySourceTerm, myCouplingTerm, myRLPusher2, myRandomnessTerm);


Data1 = InteractOrderBooks([lob_model¹, lob_model²], -1, true);

ν = 8.0
lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
  mySourceTerm, myCouplingTerm, myRLPusher1, myRandomnessTerm);

lob_model² = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
  mySourceTerm, myCouplingTerm, myRLPusher2, myRandomnessTerm);
Data2 = InteractOrderBooks([lob_model¹, lob_model²], -1, true);

ν = 12.0
lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
  mySourceTerm, myCouplingTerm, myRLPusher1, myRandomnessTerm);

lob_model² = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
  mySourceTerm, myCouplingTerm, myRLPusher2, myRandomnessTerm);
Data3 = InteractOrderBooks([lob_model¹, lob_model²], -1, true);

ν = 15.0
lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
  mySourceTerm, myCouplingTerm, myRLPusher1, myRandomnessTerm);

lob_model² = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
  mySourceTerm, myCouplingTerm, myRLPusher2, myRandomnessTerm);
Data4 = InteractOrderBooks([lob_model¹, lob_model²], -1, true);

ν = 18.0
lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
  mySourceTerm, myCouplingTerm, myRLPusher1, myRandomnessTerm);

lob_model² = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
  mySourceTerm, myCouplingTerm, myRLPusher2, myRandomnessTerm);
Data5 = InteractOrderBooks([lob_model¹, lob_model²], -1, true);


# Plots
dt = collect(1:1:400)

(average_epps_mean, average_epps_value, m) = generate_epps_plots_values(Data)
(average_epps_mean1, average_epps_value1, m) = generate_epps_plots_values(Data1)
(average_epps_mean2, average_epps_value2, m) = generate_epps_plots_values(Data2)
(average_epps_mean3, average_epps_value3, m) = generate_epps_plots_values(Data3)
(average_epps_mean4, average_epps_value4, m) = generate_epps_plots_values(Data4)
(average_epps_mean5, average_epps_value5, m) = generate_epps_plots_values(Data5)
q = quantile.(TDist(m - 1), [0.975])

p1_x = plot(dt, average_epps_mean, dpi=300, fillalpha=0.15, label="ν=3")
plot!(dt, average_epps_mean1, dpi=300, fillalpha=0.15, label="ν=6")
plot!(dt, average_epps_mean2, dpi=300, fillalpha=0.15, label="ν=9")
plot!(dt, average_epps_mean3, dpi=300, fillalpha=0.15, label="ν=12")
plot!(dt, average_epps_mean4, dpi=300, fillalpha=0.15, label="ν=15")
plot!(dt, average_epps_mean5, dpi=300, fillalpha=0.15, label="ν=18")

xlabel!(p1_x, L"\Delta t\textrm{[sec]}")
ylabel!(p1_x, L"\rho_{\Delta t}^{ij}")

savefig(p1_x, "Plots/Epps/Chaning_Nu.png")
