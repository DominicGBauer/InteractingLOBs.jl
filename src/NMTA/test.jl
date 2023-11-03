using InteractingLOBs, Optimization, OptimizationOptimJL, CSV, DataFrames
include("../setup.jl")

r = to_real_time(14401, lob_model¹.Δt)  #r is the time in real time
s = to_simulation_time(r, lob_model¹.Δt)  #s is the time in real time
real_data = vec(Matrix(CSV.read("Data/Original_Price_Bars_2300.csv", DataFrame)))



lob_model¹ = SLOB(2, T, p₀, M, L, D, ν, α, γ, mySourceTerm, myCouplingTerm, myRLPusher1, myRandomnessTerm);
lob_model² = SLOB(2, T, p₀, M, L, D, ν, α, γ, mySourceTerm, myCouplingTerm, myRLPusher2, myRandomnessTerm);
Data = InteractOrderBooks([lob_model¹,lob_model²], -1, true);
path1 = Data[1][1].raw_price_paths[1:s]

variables = [M, L, D, ν, γ]

r = to_real_time(14401, lob_model¹.Δt)  #r is the time in real time
s = to_simulation_time(r, lob_model¹.Δt)  #s is the time in real time

function sum_of_squared_errors(variables, p)
    lob_model¹ = SLOB(2, 2299, p₀, variables[1], variables[2], variables[3], variables[4], α, variables[5], mySourceTerm, myCouplingTerm, myRLPusher1, myRandomnessTerm);
    lob_model² = SLOB(2, 2299, p₀, variables[1], variables[2], variables[3], variables[4], α, variables[5], mySourceTerm, myCouplingTerm, myRLPusher2, myRandomnessTerm);
    Data = InteractOrderBooks([lob_model¹,lob_model²], -1, true);
    generated_data = Data[1][1].raw_price_paths[1:s]
    return sum((real_data .- generated_data).^2)
end


# rosenbrock(x, p) = (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2

x0 = zeros(2299)

# p = [1.0, 100.0]

prob = Optimization.OptimizationProblem(sum_of_squared_errors, x0, variables)

# sol = solve(prob, Optim.NelderMead())
