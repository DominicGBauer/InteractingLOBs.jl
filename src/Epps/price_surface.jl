using InteractingLOBs, Plots

path_to_plot = 1

function plot_density_visual(
    s,
    lob_num,
    Dat;
    dosum=false,
    plot_raw_price=true,
    x_axis_width = -1,
    center = -2,
    marker_size=1.6,
    overall_alpha=[1.0,1.0,1.0,1.0,1.0]
)

    lob_model = Dat[1][lob_num].slob

    # the below just ensure we see the full graph (100-myp)% of the time
    #myp = 10
    #num_time_steps = to_simulation_time(T,Δt)

    #max_y¹ = percentile( [maximum(Dat[path_to_plot][1].lob_densities[:,i]) for i in 1:num_time_steps] , 100-myp)
    #min_y¹ = percentile( [minimum(Dat[path_to_plot][1].lob_densities[:,i]) for i in 1:num_time_steps] , myp)
    #if length(Dat[1])>1
    #    max_y² = percentile( [maximum(Dat[path_to_plot][2].lob_densities[:,i]) for i in 1:num_time_steps] , 100-myp)
    #    min_y² = percentile( [minimum(Dat[path_to_plot][2].lob_densities[:,i]) for i in 1:num_time_steps] , myp);
    #end

    mult = Dat[1][lob_num].slob.Δt

    if center == -1
        center = Dat[1][lob_num].raw_price_paths[s]
    end

    if center == -2
        center = lob_model.p₀
    end

    if x_axis_width==-1
        x_axis_width = Dat[1][lob_num].slob.L/8
    end

    if Dat[1][lob_num].slob.old_way
        removal_fraction = - Dat[1][lob_num].slob.nu * Dat[1][lob_num].slob.Δt
    else
        removal_fraction = exp(-Dat[1][lob_num].slob.nu * Dat[1][lob_num].slob.Δt) - 1
    end

    x_axis  = [center-x_axis_width,center+x_axis_width]

    lob_densities = Dat[path_to_plot][lob_num].lob_densities[:,s]
    couplings = Dat[path_to_plot][lob_num].couplings[:,s]
    sources = Dat[path_to_plot][lob_num].sources[:,s]
    rl_pushes = Dat[path_to_plot][lob_num].rl_pushes[:,s]

    if dosum
        lob_densities .+= Dat[path_to_plot][3-lob_num].lob_densities[:,s]
        couplings .+= Dat[path_to_plot][3-lob_num].couplings[:,s]
        sources .+= Dat[path_to_plot][3-lob_num].sources[:,s]
        rl_pushes .+= Dat[path_to_plot][3-lob_num].rl_pushes[:,s]
    end

    plt = plot(lob_model.x,                lob_densities, color=1,label="Density",alpha = overall_alpha[1]);
    plot!(lob_model.x,               mult.*couplings, color=2, label="Coupling",alpha = overall_alpha[2]) ;
    plot!(lob_model.x,               mult.*sources, color=3, label="Source",alpha = overall_alpha[3]) ;
    plot!(lob_model.x,   removal_fraction.*lob_densities, color=5, label="Removal",alpha = overall_alpha[4]) ;
    plot!(lob_model.x,               mult.*rl_pushes, color=4, label="RL",alpha = overall_alpha[5]) ;

    scatter!(lob_model.x,                lob_densities, color=1,label="",markersize=marker_size,markerstrokewidth=0.0,alpha = overall_alpha[1]);
    scatter!(lob_model.x,               mult.*couplings, color=2, label="",markersize=marker_size,markerstrokewidth=0.0,alpha = overall_alpha[2]) ;
    scatter!(lob_model.x,               mult.*sources, color=3, label="",markersize=marker_size,markerstrokewidth=0.0,alpha = overall_alpha[3]) ;
    scatter!(lob_model.x,   removal_fraction.*lob_densities, color=5, label="",markersize=marker_size,markerstrokewidth=0.0,alpha = overall_alpha[4]) ;
    scatter!(lob_model.x,               mult.*rl_pushes, color=4, label="",markersize=marker_size,markerstrokewidth=0.0,alpha = overall_alpha[5]) ;

    if (plot_raw_price)
        if (isinteger(s*lob_model.Δt))
            mycol = :black
        else
            mycol = 1
    end

        scatter!([Dat[path_to_plot][lob_num].raw_price_paths[s]],[0],label="Midprice",markercolor=mycol, markersize=3,markerstrokewidth=0.5)
    end

    real_time = round(lob_model.Δt * s,digits=3)
    #plot!(lob_model.x, x -> 0,color="black", primary=false) ; #draw horizontal line
    #plot!( legend=:bottomleft, title="time=$r", xlab="Price", ylab="LOB&Source",ylim=[min_y¹,max_y¹],xlim=x_axis)
    #plot!( legend=:bottomleft, title="time=$r", xlab="Price", ylab="LOB&Source",ylim=[-20,1],xlim=x_axis)
    plot!( legend=:bottomleft, title="time=$real_time", xlab="Price", ylab="LOB&Source",xlim=x_axis)
    return plt
end;


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

lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ, mySourceTerm, myCouplingTerm, myRLPusher1, myRandomnessTerm);

lob_model² = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
    mySourceTerm, myCouplingTerm, myRLPusher2, myRandomnessTerm);

Data = InteractOrderBooks([lob_model¹, lob_model²], -1, true);

# s = 14401

plot_density_visual(s, 1, Data)
