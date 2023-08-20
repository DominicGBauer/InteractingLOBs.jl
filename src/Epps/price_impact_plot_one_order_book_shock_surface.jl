include("../setup.jl")

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

function quick_plot(lob_models,sim_start_time,how_many=12; center_at_midprice=true, dosum=false, plot_raw_price=true, x_axis_width = -1, center = -2, scale_down=true, marker_size=1.6,overall_alpha=[1.0,1.0,1.0,1.0,1.0])
    try clear_double_dict(Dat) catch e print("Not initialized") end
    GC.gc()

    Dat = InteractOrderBooks(lob_models, -1, true) ;

    if x_axis_width==-1
        x_axis_width=10#lob_models[1].L/8
    end

    p_arr1 = Array{Plots.Plot{Plots.GRBackend},1}(undef,how_many)
    for i in [1:how_many;]
        p_arr1[i] = plot_density_visual(sim_start_time-2+i, to_real_time(sim_start_time-2+i,Δt), 1, Dat;
            dosum=dosum, plot_raw_price=plot_raw_price, x_axis_width = x_axis_width, center = center, scale_down=scale_down, marker_size=marker_size, overall_alpha=overall_alpha )#,overall_alpha=[0.4,1.0,0.4,0.4,0.4,0.4])
    end
    p1 = plot(p_arr1...,size=(1200,1000), dpi=100)
    savefig(p1, "Plots/Epps/PriceSurfaceOneShock.png")

    return Dat
end

quick_plot([lob_model¹, lob_model²],SimStartTime,9; dosum=false,plot_raw_price=true,x_axis_width=10)
