using DifferentialEquations, Plots

include("./epps.jl")
include("../setup_variables.jl")

function coupling_inner1(x, p¹, p²)
    f(y) = -(μ * (y)) * exp(-(μ * (y))^2) #y is a temporary variable
    g = 1 + tanh(abs(p² - p¹) / a) * c

    if (p² > p¹)
        if x > p¹
            return b * f(1 / g * (x - p¹))
        else
            return b * g * f(x - p¹)
        end
    else
        if x < p¹
            return b * f(1 / g * (x - p¹))
        else
            return b * g * f(x - p¹)
        end
    end
end

W1 = WienerProcess(0.0, 1.0, 1.0)
dt = 0.1
W1.dt = dt
u = nothing;
p = nothing; # for state-dependent distributions
calculate_step!(W1, dt, u, p)
for i in 1:14401
    accept_step!(W1, dt, u, p)
end

W2 = WienerProcess(0.0, 1.0, 1.0)
dt = 0.1
W2.dt = dt
calculate_step!(W2, dt, u, p)
for i in 1:14401
    accept_step!(W2, dt, u, p)
end

coupling_inner_W1 = [coupling_inner1(x, W1.u[400], W2.u[400]) for x in W1.u]
coupling_inner_W1 = coupling_inner_W1 .+ abs(minimum(coupling_inner_W1)) .+ 0.001
coupling_inner_W2 = [coupling_inner1(x, W1.u[400], W2.u[400]) for x in W2.u]
coupling_inner_W2 = coupling_inner_W2 .+ abs(minimum(coupling_inner_W2)) .+ 0.001

index_vector = 0:1.0:(size(W1)[1]-1)

epps_data_random_walk = hcat(index_vector, coupling_inner_W1, coupling_inner_W2)
epps_random_walk = Empirical(epps_data_random_walk)

# # Save and Load
save("Computed Data/EppsCorrection/Empirical_Random_Walk.jld", "epps", epps_random_walk)

ComputedResults = load("Computed Data/EppsCorrection/Empirical_Random_Walk.jld")
epps_random_walk = ComputedResults["epps"]

epps_mean_random_walk = mean(epps_random_walk[1], dims=2)
p1 = plot(1:size(epps_mean_random_walk)[1], epps_mean_random_walk, legend=false, dpi=300, fillalpha=0.15)

xlabel!(p1, L"\Delta t\textrm{[sec]}")
ylabel!(p1, L"\rho_{\Delta t}^{ij}")
ylims!(-0.01, 0.01)


savefig(p1, "Plots/Epps/Epps_Random_Walk.png")
