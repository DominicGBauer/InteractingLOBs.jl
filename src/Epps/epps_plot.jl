using InteractingLOBs, LaTeXStrings

include("./epps.jl")

num_paths = 2#30

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
Δt = (r * (Δx^2) / (2.0 * D))^(1/γ)

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

r = to_real_time(14401, lob_model¹.Δt)  #r is the time in real time
s = to_simulation_time(r, lob_model¹.Δt)  #s is the time in real time

Data = InteractOrderBooks([lob_model¹,lob_model²], -1, true);
index_vector = 0:1.0:(size(Data[1][1].raw_price_paths[1:s])[1]-1)


function generate_epps_plots_values(simulation_data)
    epps_value = [zeros(M, 1) for i = 1:num_paths]
    epps_mean = [zeros(M, 1) for i = 1:num_paths]
    total_epps_mean = zeros(size(epps_mean[1])[1])
    total_epps_value = zeros(size(epps_value[1])[1])
    m = 0

    for i in 1:num_paths
        path1 = simulation_data[i][1].raw_price_paths[1:s]
        path2 = simulation_data[i][2].raw_price_paths[1:s]
        index_vector = 0:1.0:(size(simulation_data[i][2].raw_price_paths[1:s])[1]-1)
        epps_data = hcat(index_vector, path1, path2)
        epps = Empirical(epps_data)
        m = size(epps_data)[1]
        # Save and Load
        # save("Computed Data/EppsCorrection/Empirical$i.jld", "epps$i", epps[i])
        # ComputedResults = load("Computed Data/EppsCorrection/Empirical$i.jld")
        # epps_result = ComputedResults["epps$i"]
        epps_value[i] = epps[1]
        epps_mean[i] = mean(epps[1], dims=2)
        total_epps_mean = total_epps_mean .+ epps_mean[i]
        total_epps_value = total_epps_value .+ epps_value[i]
    end

    average_epps_mean = total_epps_mean ./ num_paths
    average_epps_value = total_epps_value ./ num_paths

    return (average_epps_mean, average_epps_value, m)
end

# (average_epps_mean, average_epps_value, m) = generate_epps_plots_values(Data)

# # Plots
# dt = collect(1:1:400)
# q = quantile.(TDist(m-1), [0.975])

# p1 = plot(dt, average_epps_mean, legend = false,dpi=300, ribbon=(q .* std(average_epps_value, dims = 2) * 0.001), fillalpha=.15)
# xlabel!(p1, L"\Delta t\textrm{[sec]}")
# ylabel!(p1, L"\rho_{\Delta t}^{ij}")

# savefig(p1, "Plots/Epps/Epps.png")

# p1 = plot(dt, epps_mean[1], legend = false,dpi=300, fillalpha=.15)
# plot!(dt, epps_mean[2], legend = false,dpi=300, fillalpha=.15)
# plot!(dt, epps_mean[3], legend = false,dpi=300, fillalpha=.15)
# plot!(dt, epps_mean[4], legend = false,dpi=300, fillalpha=.15)
# plot!(dt, epps_mean[5], legend = false,dpi=300, fillalpha=.15)

# p1 = plot(dt, epps_mean[1], legend = false,dpi=300, fillalpha=.15)
# plot!(dt, epps_mean[2], legend = false,dpi=300, fillalpha=.15)
# plot!(dt, epps_mean[3], legend = false,dpi=300, fillalpha=.15)
# plot!(dt, epps_mean[4], legend = false,dpi=300, fillalpha=.15)
# plot!(dt, epps_mean[5], legend = false,dpi=300, fillalpha=.15)
# plot(dt, total_epps_mean, legend = false,dpi=300, ribbon=(q .* std(total_epps_value, dims = 2) .* 0.001), line=(2, [:solid]), fillalpha=.15, color=:red)

# function get_indicative_path(data)
#     total_prices_path_1 = zeros(size(data[1][1].raw_price_paths[1:s]))
#     total_prices_path_2 = zeros(size(data[1][1].raw_price_paths[1:s]))
#     for i in 1:num_paths
#         total_prices_path_1 = total_prices_path_1 + Data[i][1].raw_price_paths[1:s]
#         total_prices_path_2 = total_prices_path_1 + Data[i][2].raw_price_paths[1:s]
#     end
#     average_prices_path_1 = total_prices_path_1 ./ num_paths
#     average_prices_path_2 = total_prices_path_2 ./ num_paths

#     return (average_prices_path_1, average_prices_path_2)
# end
# (path1, path2) = get_indicative_path(Data)
# epps_data = hcat(index_vector, path1, path2)

# epps = Empirical(epps_data)

# Save and Load
# save("Computed Data/EppsCorrection/Empirical.jld", "epps", epps)

# ComputedResults = load("Computed Data/EppsCorrection/Empirical.jld")
# epps = ComputedResults["epps"]

# Plots
# dt = collect(1:1:400)
# m = size(epps_data)[1]
# q = quantile.(TDist(m-1), [0.975])

# epps_mean = mean(epps[1], dims=2)
# p1 = plot(dt, epps_mean, legend = false,dpi=300, ribbon=(q .* std(epps[1], dims = 2) * 0.00001), fillalpha=.15)

# p1 = plot(dt, epps_mean, legend = :bottomright,dpi=300, label = L"\textrm{Measured}", ribbon=(movingaverage * 0.05))
# plot!(movingaverage_upper, color=:red, label = L"\textrm{LOL}")
# plot!(movingaverage_lower, color=:red, label = L"\textrm{SMA -2% lower bound}")
# p1 = plot(dt, mean(SBKFSR[1], dims=2), ribbon=(q .* std(SBKFSR[1], dims = 2)), fillalpha=.15, legend = :topright, color = :red, line=(1, [:solid]), label = L"\textrm{Measured}", marker=([:+ :d],1,0,stroke(2,:red)), dpi = 300, ylims = (-0.1, 1.55))
# plot!(p1, dt, mean(epps[2], dims=2), ribbon=(q .* std(epps[2], dims = 2)), fillalpha=.15, color = :blue, line=(1, [:solid]), label = L"\textrm{Flat trade correction}", marker=([:x :d],1,0,stroke(2,:blue)))
# plot!(p1, dt, mean(epps[3], dims=2), ribbon=(q .* std(epps[3], dims = 2)), fillalpha=.15, color = :green, line=(1, [:solid]), label = L"\textrm{Overlap correction}", marker=([:circle :d],1,0,stroke(2,:green)))
# hline!(p1, [mean(epps[4])], ribbon=(q .* std(epps[4], dims = 1)), fillalpha=.15, color = :brown, line=(1, [:dash]), label = L"\textrm{HY}")
# xlabel!(p1, L"\Delta t\textrm{[sec]}")
# ylabel!(p1, L"\rho_{\Delta t}^{ij}")

# savefig(p1, "Plots/Epps/Epps.png")
