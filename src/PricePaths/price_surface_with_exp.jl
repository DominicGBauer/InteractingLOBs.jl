include("../setup.jl")

RealStartTime = 50 # when, in real time, to kick the system
SimStartTime = to_simulation_time(RealStartTime, Δt) - 10 # convert to simulation time
SimEndTime = SimStartTime + 3 # when to stop kicking, in simulation time
Position = 200
Volume = -8; # If position == -x where x>=0, then put it x above the mid price each time

myRLPusher1 = RLPushTerm(SimStartTime, SimEndTime, Position, Volume, true)

myRLPusher2 = RLPushTerm(SimStartTime, SimEndTime, Position, Volume, false)

lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, ν, γ,
  mySourceTerm, myCouplingTerm, myRLPusher1, myRandomnessTerm, do_exp_dist_times=true);

lob_model² = SLOB(num_paths, T, p₀, M, L, D, ν, γ,
  mySourceTerm, myCouplingTerm, myRLPusher2, myRandomnessTerm, do_exp_dist_times=true);

Data = InteractOrderBooks([lob_model¹, lob_model²], -1, true);


function quick_plot(Dat, sim_start_time, how_many=12; center_at_midprice=true, dosum=false, plot_raw_price=true, x_axis_width=-1, center=-2, scale_down=true, marker_size=1.6, overall_alpha=[1.0, 1.0, 1.0, 1.0, 1.0])
  if x_axis_width == -1
    x_axis_width = 10#lob_models[1].L/8
  end

  p_arr1 = Array{Plots.Plot{Plots.GRBackend},1}(undef, how_many)
  for i in [1:how_many;]
    p_arr1[i] = plot_density_visual(Dat, sim_start_time - 2 + i, 1; do_left=false, do_right=false)
    savefig(p_arr1[i], "Plots/PriceSurfaceWithExp/PriceSurfaceOneShockWithExp$i.png")
  end
  p1 = plot(p_arr1..., size=(1200, 1000), dpi=300)

  return p1
end

quick_plot(Data, SimStartTime, 9; dosum=false, plot_raw_price=true, x_axis_width=10)

path1 = Data[1][1].obs_price_paths
path2 = Data[2][1].obs_price_paths
plot1 = plot(path1, label="Price path A", dpi=300)
plot!(path2, label="Price path B")
xlabel!("t")
ylabel!("p")
savefig(plot1, "Plots/PricePath/PricePathsWithExp.png",)
