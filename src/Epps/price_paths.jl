using Plots

include("../setup.jl")

myRandomnessTerm = RandomnessTerm(σ, r, β, lag, do_random_walk, true)
myCouplingTerm = CouplingTerm(μ, a, b, c, true);

myRLPusher1 = RLPushTerm(SimStartTime, SimEndTime, Position, Volume, true)
myRLPusher2 = RLPushTerm(SimStartTime, SimEndTime, Position, Volume, false)

lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ, mySourceTerm, myCouplingTerm, myRLPusher1, myRandomnessTerm);

lob_model² = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
    mySourceTerm, myCouplingTerm, myRLPusher2, myRandomnessTerm);

Data = InteractOrderBooks([lob_model¹, lob_model²], -1, true);
path1 = Data[1][1].obs_price_paths
path2 = Data[2][1].obs_price_paths
plot1 = plot(path1, label="Price path A", dpi=300)
plot!(path2, label="Price path B")
xlabel!("t")
ylabel!("p")
savefig(plot1, "Plots/PricePath/PricePaths.png",)

# myRandomnessTerm = RandomnessTerm(σ,r,β,lag,do_random_walk,true)
# myCouplingTerm = CouplingTerm(μ, a, b, c, true);

# myRLPusher1 = RLPushTerm(SimStartTime,SimEndTime,Position,Volume,true)
# myRLPusher2 = RLPushTerm(SimStartTime,SimEndTime,Position,Volume,false)

# lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ, mySourceTerm, myCouplingTerm, myRLPusher1, myRandomnessTerm);

# lob_model² = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
#     mySourceTerm, myCouplingTerm, myRLPusher2, myRandomnessTerm);

# Data = InteractOrderBooks([lob_model¹, lob_model²], -1, true);
# path_with_coupling = Data[1][1].obs_price_paths
# plot2 = plot(path_with_coupling,label="Price path", dpi=300)
# xlabel!("t")
# ylabel!("p")
# savefig(plot2, "Plots/PricePath/PathWithCoupling.png")
