# -*- coding: utf-8 -*-
using Revise
using Plots
using RelevanceStacktrace
using StatsBase 
using Distributions
using ProgressMeter
using Statistics
using Random  
using CurveFit
import Random:rand
using SpecialFunctions
using JLD2
using Combinatorics

using InteractingLOBs

Revise.revise()

# # USEFUL FUNCTIONS

# +
function get_sums(lob_densities;absolute=true)
    l = size(lob_densities)[2]
    sums = zeros(Float64, l)
    
    if absolute 
        my_sum = (x) -> sum(abs.(x))
    else
        my_sum = (x) -> sum(x)
    end
    
    for t in 1:l
        sums[t] = my_sum(lob_densities[:,t])
    end
    
    return sums
end

function save_png(folder_name,name)
    png(string("/home/derickdiana/Desktop/Masters/",folder_name,"/",name,".png"))
end

function plot_sums(Dat;path_num=1,slob_num=1,new_plot=true,absolute=true)
    if new_plot
        plot()
    end
    sums = get_sums(Dat[path_num][slob_num].lob_densities;absolute=absolute)
    slob = Dat[path_num][slob_num].slob

    scatter!(sums,markersize=1.4,label="Area",size=(1400,500))
    vline!([slob.rl_push_term.StartTime],label="Kick time")
    if absolute
        hline!([sums[slob.rl_push_term.StartTime-1]+slob.Δt*slob.rl_push_term.Amount],label="Theoretical kick volume")
    else
        hline!([-slob.Δt*slob.rl_push_term.Amount],label="Theoretical kick volume")
    end
    
    source_sums = get_sums(Dat[path_num][slob_num].sources;absolute=absolute)
    hline!([source_sums[3]/slob.nu],label="Theoretical equilibrium volume")
    
    plot!(xlabel="Simulation steps",ylabel="Signed area under system")
    plot!(size=(1200,500))
end

path_to_plot = 1
l = @layout [a d; b e; c f];

function plot_price_path(s, r, lob_num, Dat, diff=false)
    lob_model = Dat[1][lob_num].slob
    
    plt = plot()
    for path in 1:lob_model.num_paths
        raw_price_paths = Dat[path][lob_num].raw_price_paths[1:s]
        if diff
            raw_price_paths .-=  Dat[path][3-lob_num].raw_price_paths[1:s]
        end
        plot!((0:s-1).*lob_model.Δt,raw_price_paths ,color=5,w=0.6) ;
        
    end
    
    obs_price_paths = Dat[path_to_plot][lob_num].obs_price_paths[1:r]
    if diff
        obs_price_paths .-=  Dat[path_to_plot][3-lob_num].obs_price_paths[1:r]
    end
    plot!(0:r-1, obs_price_paths,color=1,w=2.7) ;
    
    plot!(legend=false, ylab="Price", xlab="Time") ;   
    return plt
end

function plot_density_visual(s, r, lob_num, Dat; 
                dosum=false, plot_raw_price=true, x_axis_width = -1, center = -2, scale_down=true, marker_size=1.6,overall_alpha=[1.0,1.0,1.0,1.0,1.0])
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



function get_second_derivative(x,temp)
    middle = 2:length(temp)-1
    return ((temp[middle.+1]-temp[middle])./(x[middle.+1]-x[middle])-(temp[middle]-temp[middle.-1])./(x[middle]-x[middle.-1]))./(x[middle.+1]-x[middle.-1])
end;

function fit_and_plot_price_impact(volumes,mean_price_impacts,var_price_impacts,labels;
                                                sub=-1,do_kinks=true,colors=-1,do_power_fit=false,xticks=-1,new_plot=false,forplot=(),
                                                do_ribbon=true,do_horiz_log_shift=false,do_log_fit=true,do_log_plot=false,do_vert_log_shift=true)
    if new_plot
        plot()
    end
    
    volumes_original = volumes[:]
    
    if do_log_plot
        volumes = log.(volumes.+1)[:]
    end
    
    #labels create indices which are usually the different gammas
    
    li_len = length(labels)
    vi_len = length(volumes)
    
    if sub==-1
        sub = 1:vi_len
    end
    
    if colors==-1 #change later!
        if li_len<=5
            colors = ["blue","red","green","orange","purple"]
        else
            colors = 1:li_len
        end
    end
    
    shift = do_horiz_log_shift ? 1 : 0
    
    for ind in 1:li_len
        if do_log_fit
            a,b = log_fit( volumes_original[sub].+shift ,  mean_price_impacts[sub,ind] )
            a = do_vert_log_shift ? a : 0
            y_line_log_fit = a.+b.*log.(volumes_original[sub].+shift)

            plot!(         volumes[sub] , y_line_log_fit ,
                    label=string("Log fit: ",round(a,digits=2)," + ",round(b,digits=2),"log(x+1)"),w=1.5,color=colors[ind])
        end
        
        scatter!(      volumes[sub],   mean_price_impacts[sub,ind],
                label=labels[ind],ms=1.5,markerstrokewidth=0.1,  ma=1,color=colors[ind])
        
        if do_ribbon
            plot!(     volumes[sub],   mean_price_impacts[sub,ind],  ribbon=var_price_impacts[sub,ind].^0.5,alpha=0,
                fillalpha=0.4,fillcolor=colors[ind],label="")
        end
        
        if do_power_fit
            c,d = power_fit(volumes[sub],mean_price_impacts[sub,ind])
            
            plot!(volumes[sub],        c.*((volumes[sub]).^d),label="Power fit",
                        w=1.5,color=colors[ind],linestyle=:dash)
        end

            
        if xticks!=-1
            plot!(xticks=xticks)
        end


        if do_kinks
            vol_scale = (volumes[end]-volumes[1])/20
            impact_scale = (maximum(mean_price_impacts[sub,ind],dims=1)[1]-minimum(mean_price_impacts[sub,ind],dims=1)[1])/30
            
            second_deriv = get_second_derivative(volumes_original[sub],mean_price_impacts[sub,ind])
            kink_position = findnext(x->x>0,second_deriv,1)
            
            while !(kink_position===nothing)
                target_x, target_y = (volumes[sub][kink_position+1],mean_price_impacts[sub[kink_position+1],ind])
                quiver!([target_x+vol_scale],[target_y-impact_scale],quiver=([-vol_scale],[impact_scale]),color=colors[ind])
                scatter!([target_x],[target_y],markershape=:star5,color=colors[ind],
                    label=string("Kink at position ",round(target_x,digits=2)),markerstrokewidth=0.1,  ma=1)
                
                kink_position = findnext(x->x>0,second_deriv,kink_position+2)
            end
        end

    end
    
    my_xlabel = do_log_plot ? "log(Volume)" : Volume
    plot!(xlabel=my_xlabel,ylabel="Price impact i.e. p(t+1)-p(t)";forplot...)
    plot!(size=(1000,1000))
end

function quick_plot(lob_models,sim_start_time,how_many=12; center_at_midprice=true,
        dosum=false, plot_raw_price=true, x_axis_width = -1, center = -2, scale_down=true, marker_size=1.6,overall_alpha=[1.0,1.0,1.0,1.0,1.0])
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
    plot(p_arr1...,size=(1200,1000))
    
    return Dat
end

function calculate_price_impacts(volumes,inputs,get_set; measured_slob=1) #lob_sets must have the form ((lob_a1,lob_b1),(lob_a2,lob_b2)) where the inner brackets are meant
                                                                                #to be passed to InteractOrderBooks together
    
    vol_len = length(volumes)
    input_len = length(inputs)

    mean_price_impacts = ones(Float64,vol_len,input_len)
    var_price_impacts  = ones(Float64,vol_len,input_len)

    for Ind in 1:input_len 
        p_outer = Progress(vol_len,dt=0.1)
        
        #for Volume in 1:vol_len
        Threads.@threads for Volume in 1:vol_len  
            lob_models, sim_start_time, l  = get_set(volumes[Volume],inputs[Ind])

            #try clear_double_dict(Dat) catch e print("Not initialized") end
            #GC.gc()

            Dat = InteractOrderBooks(lob_models, -1, false) ;

            num_paths = lob_models[1].num_paths
            price_impact = zeros(Float64,num_paths)
            for path_num in 1:num_paths
                price_impact[path_num] = Dat[path_num][measured_slob].raw_price_paths[sim_start_time+l] - Dat[path_num][measured_slob].raw_price_paths[sim_start_time]
            end

            mean_price_impacts[Volume,Ind] = mean(price_impact)
            var_price_impacts[Volume,Ind] = var(price_impact)

            next!(p_outer)

        end
    end

    mean_price_impacts = .- mean_price_impacts;
    
    return (mean_price_impacts,var_price_impacts)
end

function myvar(p,x,Δx)
    mu = mymean(p,x,Δx)
    Ex2 = 0
    for i in 1:length(p)
        Ex2 += p[i]*x[i]^2
    end
    return Ex2 - mu^2
end

function mymean(p,x,Δx)
    sum = 0
    for i in 1:length(p)
        sum += p[i]*x[i]
    end
    return sum
end

# -

function myfunc(a,args...;b=1,argz...)
    return args[1]
end
function myfunc2(a,args...;b=1,argz...)
    return myfunc(9,args...)
end
myfunc2(8,2,3)

# # GENERAL WORKING

# +
# Configuration Arguments
num_paths = 10#30

L = 200     # real system width (e.g. 200 meters)
M = 400     # divided into M pieces , 400

T = 2000  # simulation runs until real time T (e.g. 80 seconds)
p₀ = 230.0  #this is the mid_price at t=0  238.75 

# Free-Parameters for gaussian version
D = 0.5 # real diffusion constant e.g. D=1 (meters^2 / second), 1
α = 0.0 # legacy, no longer used

ν = 2.0 #removal rate
γ = 1.0 #fraction of derivative (1 is normal diffusion, less than 1 is D^{1-γ} derivative on the RHS)

# Source term:
λ = 3.0 #
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

print((Δt,to_simulation_time(T,Δt),num_paths*to_simulation_time(T,Δt))) #about 2GB RAM per 100K, i.e. can only do about 1.8 million
lob_model¹.SK_DP

# +
# clear everything pointed to by the dictionary then garbage collect. If it wasn't assigned yet it will left you know
try clear_double_dict(Dat) catch e print("Not initialized") end
GC.gc()

Dat = InteractOrderBooks([lob_model¹,lob_model²], -1, true) ;
#Dat = InteractOrderBooks([lob_model¹], -1, true) ;

# +
#total_steps = min(to_simulation_time(T,Δt),100)
total_steps = 10
total_length = to_simulation_time(T,Δt)
step = floor(Int,total_length/total_steps)

range = 1:step:(total_length-step)
#range = [SimStartTime:SimStartTime+to_sim(1)*2;]
#range = [SimStartTime]

p_outer = Progress(length(range),dt=0.1)

l = @layout [a d; b e; c f];

anim = @animate for s = myrange           #s is the time in simulation time
    r = to_real_time(s, lob_model¹.Δt)  #r is the time in real time
    #s = to_simulation_time(r, lob_model¹.Δt)  #s is the time in real time
    
    plt1 = plot_price_path(s,r,1,Dat)
    if length(Dat[1])>1
        plt3 = plot_price_path(s,r,2,Dat)
        plt5 = plot_price_path(s,r,1,Dat,true)
    end
    
    plt2 = plot_density_visual(s, r, 1, Dat)
    if length(Dat[1])>1
        plt4 = plot_density_visual(s, r, 2, Dat)
        plt6 = plot_density_visual(s, r, 1, Dat; dosum=true, plot_raw_price=false)
    end
    
    if length(Dat[1])>1
        plot(plt1, plt2, plt3, plt4, plt5, plt6 ,layout=l,size=(1000,1000))
    else
        plot(plt1,plt2,size=(1000,1000))
    end
    
    next!(p_outer)
end

gif(anim, "/tmp/LOB.gif", fps=20*length(myrange)/200)
#gif(anim, "~/Desktop/Masters/StillWorking/Random_Walk_And_Coupling_For_Alpha_Is_0.7.gif", fps=20*length(myrange)/200)
# -

# # DIFFERENT CORRELATIONS OF STYLIZED FACTS

# +
# Configuration Arguments
num_paths = 1#30

L = 200     # real system width (e.g. 200 meters)
M = 400     # divided into M pieces , 400

T = 10000  # simulation runs until real time T (e.g. 80 seconds)
p₀ = 230.0  #this is the mid_price at t=0  238.75 

# Free-Parameters for gaussian version
D = 0.5 # real diffusion constant e.g. D=1 (meters^2 / second), 1
α = 0.0 # legacy, no longer used

ν = 10.0 #removal rate
γ = 1.0 #fraction of derivative (1 is normal diffusion, less than 1 is D^{1-γ} derivative on the RHS)

# Source term:
λ = 1.0 #
μ = 0.1 #

mySourceTerm = SourceTerm(λ, μ, true);

# Coupling term:
a = 13.0  #gap between stocks before at full strength: strong is 0.3
b = 1.0   #weighting of interaction term: strong is 2
c = 1.2   #skew factor: strong is 2

myCouplingTerm = CouplingTerm(μ, a, b, c, true);



Δx = L / M  # real gap between simulation points 
Δt = (r * (Δx^2) / (2.0 * D))^(1/γ)

# RL Stuff:
RealStartTime = 2 # when, in real time, to kick the system
SimStartTime = to_simulation_time(RealStartTime,Δt)-2 # convert to simulation time
SimEndTime = SimStartTime + 10 # when to stop kicking, in simulation time
Position = 200
Volume = -8; # If position == -x where x>=0, then put it x above the mid price each time

myRLPusher1 = RLPushTerm(SimStartTime,SimEndTime,Position,Volume,false)

myRLPusher2 = RLPushTerm(SimStartTime,SimEndTime,Position,Volume,false)

print((Δt,to_simulation_time(T,Δt),num_paths*to_simulation_time(T,Δt))) #about 2GB RAM per 100K, i.e. can only do about 1.8 million
lob_model¹.SK_DP
# -

if false
    #total_steps = min(to_simulation_time(T,Δt),100)
    total_steps = 10
    total_length = to_simulation_time(T,Δt)
    step = floor(Int,total_length/total_steps)

    myrange = 1:step:(total_length-step)
    #myrange = [SimStartTime:SimStartTime+to_sim(1)*2;]
    #myrange = [SimStartTime]

    p_outer = Progress(length(myrange),dt=0.1)

    l = @layout [a d; b e; c f];

    anim = @animate for s = myrange           #s is the time in simulation time
        r = to_real_time(s, lob_model¹.Δt)  #r is the time in real time
        #s = to_simulation_time(r, lob_model¹.Δt)  #s is the time in real time

        plt1 = plot_price_path(s,r,1,Dat)
        if length(Dat[1])>1
            plt3 = plot_price_path(s,r,2,Dat)
            plt5 = plot_price_path(s,r,1,Dat,true)
        end

        plt2 = plot_density_visual(s, r, 1, Dat)
        if length(Dat[1])>1
            plt4 = plot_density_visual(s, r, 2, Dat)
            plt6 = plot_density_visual(s, r, 1, Dat, true, false)
        end

        if length(Dat[1])>1
            plot(plt1, plt2, plt3, plt4, plt5, plt6 ,layout=l,size=(1000,1000))
        else
            plot(plt1,plt2,size=(1000,1000))
        end

        next!(p_outer)
    end

    gif(anim, "/tmp/LOB.gif", fps=20*length(myrange)/200)
end

# +
# My randomness term
σ = 1.0 #variance in randomness
r = 0.5 #proportion of time in which it jumps left or right
β = 0.0 #probability of being the value of the previous lag or mean reversion strength
lag = 10 #lag
do_random_walk = false #behave like a random walk
myRandomnessTerm = RandomnessTerm(σ,r,β,lag,do_random_walk,true)
T = 20000

lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ, 
    mySourceTerm, myCouplingTerm, myRLPusher1, myRandomnessTerm);

lob_model² = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
    mySourceTerm, myCouplingTerm, myRLPusher2, myRandomnessTerm);
# clear everything pointed to by the dictionary then garbage collect. If it wasn't assigned yet it will left you know
try clear_double_dict(Dat) catch e print("Not initialized") end
GC.gc()

Dat = InteractOrderBooks([lob_model¹,lob_model²], -1, true) ;
#Dat = InteractOrderBooks([lob_model¹], -1, true) ;

# +
observed_price_path = Dat[1][1].raw_price_paths[1:end-1,1];
observed_log_returns = diff(log.(observed_price_path[:,1]));

data_stylized_facts = StylizedFacts.StylizedFactsPlot(observed_price_path);

StylizedFacts.plot_all_stylized_facts(data_stylized_facts,(1000,1200))

#png("/home/derickdiana/Desktop/Masters/Reworked/StylizedFactsRandomKicksNoFractional.png")
#plot!()

# +
# My randomness term
σ = 1.0 #variance in randomness
r = 0.5 #proportion of time in which it jumps left or right
β = 0.6 #probability of being the value of the previous lag or mean reversion strength
lag = 10 #lag
do_random_walk = true #behave like a random walk
myRandomnessTerm = RandomnessTerm(σ,r,β,lag,do_random_walk,true)
T = 20000

lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ, 
    mySourceTerm, myCouplingTerm, myRLPusher1, myRandomnessTerm);

lob_model² = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
    mySourceTerm, myCouplingTerm, myRLPusher2, myRandomnessTerm);

# clear everything pointed to by the dictionary then garbage collect. If it wasn't assigned yet it will left you know
try clear_double_dict(Dat) catch e print("Not initialized") end
GC.gc()

Dat = InteractOrderBooks([lob_model¹,lob_model²], -1, true) ;
#Dat = InteractOrderBooks([lob_model¹], -1, true) ;

# +
observed_price_path = Dat[1][1].raw_price_paths[1:end-1,1];
observed_log_returns = diff(log.(observed_price_path[:,1]));

data_stylized_facts = StylizedFacts.StylizedFactsPlot(observed_price_path);

StylizedFacts.plot_all_stylized_facts(data_stylized_facts,(1000,1200))

#png("/home/derickdiana/Desktop/Masters/Reworked/StylizedFactsRandomWalkNoFractional.png")
#plot!()

# +
# My randomness term
σ = 1.0 #variance in randomness
r = 0.5 #proportion of time in which it jumps left or right
β = 0.9 #probability of being the value of the previous lag or mean reversion strength
lag = 10 #lag
do_random_walk = true #behave like a random walk
myRandomnessTerm = RandomnessTerm(σ,r,β,lag,do_random_walk,true)
T = 20000

lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ, 
    mySourceTerm, myCouplingTerm, myRLPusher1, myRandomnessTerm);

lob_model² = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
    mySourceTerm, myCouplingTerm, myRLPusher2, myRandomnessTerm);

# clear everything pointed to by the dictionary then garbage collect. If it wasn't assigned yet it will left you know
try clear_double_dict(Dat) catch e print("Not initialized") end
GC.gc()

Dat = InteractOrderBooks([lob_model¹,lob_model²], -1, true) ;
#Dat = InteractOrderBooks([lob_model¹], -1, true) ;

# +
observed_price_path = Dat[1][1].raw_price_paths[1:end-1,1];
observed_log_returns = diff(log.(observed_price_path[:,1]));

data_stylized_facts = StylizedFacts.StylizedFactsPlot(observed_price_path);

StylizedFacts.plot_all_stylized_facts(data_stylized_facts,(1000,1200))

#png("/home/derickdiana/Desktop/Masters/Reworked/StylizedFactsStrongerRandomWalkNoFractional.png")
plot!()
# -

# # STYLIZED FACTS

# ## Plot stylised facts for no fractional

# ## Stylised facts for fractional

# +
# Configuration Arguments
num_paths = 1#30

L = 200     # real system width (e.g. 200 meters)
M = 400     # divided into M pieces , 400

T = 10000  # simulation runs until real time T (e.g. 80 seconds)
p₀ = 230.0  #this is the mid_price at t=0  238.75 

# Free-Parameters for gaussian version
D = 0.5 # real diffusion constant e.g. D=1 (meters^2 / second), 1
α = 0.0 # legacy, no longer used

ν = 10.0 #removal rate
γ = 1.0 #fraction of derivative (1 is normal diffusion, less than 1 is D^{1-γ} derivative on the RHS)

# Source term:
λ = 1.0 #
μ = 0.1 #

mySourceTerm = SourceTerm(λ, μ, true);

# Coupling term:
a = 13.0  #gap between stocks before at full strength: strong is 0.3
b = 1.0   #weighting of interaction term: strong is 2
c = 1.2   #skew factor: strong is 2

myCouplingTerm = CouplingTerm(μ, a, b, c, false);

# My randomness term
σ = 0.8 #variance in randomness
r = 0.3 #proportion of time in which it jumps left or right
β = 0.0 #probability of being the value of the previous lag or mean reversion strength
lag = 10 #lag
do_random_walk = true #behave like a random walk
myRandomnessTerm = RandomnessTerm(σ,r,β,lag,do_random_walk,true)

Δx = L / M  # real gap between simulation points 
Δt = (r * (Δx^2) / (2.0 * D))^(1/γ)

# RL Stuff:
RealStartTime = 2 # when, in real time, to kick the system
SimStartTime = to_simulation_time(RealStartTime,Δt)-2 # convert to simulation time
SimEndTime = SimStartTime + 10 # when to stop kicking, in simulation time
Position = 200
Volume = -8; # If position == -x where x>=0, then put it x above the mid price each time

myRLPusher = RLPushTerm(SimStartTime,SimEndTime,Position,Volume,false)

lob_model = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ, 
    mySourceTerm, myCouplingTerm, myRLPusher, myRandomnessTerm);

print((Δt,to_simulation_time(T,Δt),num_paths*to_simulation_time(T,Δt))) #about 2GB RAM per 100K, i.e. can only do about 1.8 million
lob_model.SK_DP

# +
# clear everything pointed to by the dictionary then garbage collect. If it wasn't assigned yet it will left you know
try clear_double_dict(Dat) catch e print("Not initialized") end
GC.gc()

Dat = InteractOrderBooks([lob_model], -1, true) ;
#Dat = InteractOrderBooks([lob_model¹], -1, true) ;

# +
observed_price_path = Dat[1][1].raw_price_paths[1:end-1];
observed_log_returns = diff(log.(observed_price_path[:,1]));

data_stylized_facts = StylizedFacts.StylizedFactsPlot(observed_price_path);

StylizedFacts.plot_all_stylized_facts(data_stylized_facts,(1000,1200))

#png("/home/derickdiana/Desktop/Masters/OfficialOne/DraftOnePics/StylizedFactsRandomKicksNoFractional.png")
plot!()
# -

# # MICHAELS VS NEW

# +
# Configuration Arguments
num_paths = 1#50#30

L = 200     # real system width (e.g. 200 meters)
M = 400     # divided into M pieces , 400

p₀ = 230.0  #this is the mid_price at t=0  238.75 

# Free-Parameters for gaussian version
D = 0.5#0.5/8 # real diffusion constant e.g. D=1 (meters^2 / second), 1
α = 0.0 # legacy, no longer used

#ν = 3.0 #removal rate
#γ = 1.0 #fraction of derivative (1 is normal diffusion, less than 1 is D^{1-γ} derivative on the RHS)

ν = 0.1#1.0#3.0 #removal rate
ν = 1.0#3.0 #removal rate
γ = 1.0 #fraction of derivative (1 is normal diffusion, less than 1 is D^{1-γ} derivative on the RHS)

# Source term:
λ = 1.0 #
μ = 0.1 #
mySourceTerm = SourceTerm(λ, μ, true);

# Coupling term:
a = 0.01  #gap between stocks before at full strength: strong is 0.3
b = 2.0   #weighting of interaction term: strong is 2
c = 2.0   #skew factor: strong is 2

myCouplingTerm = CouplingTerm(μ, a, b, c, true);

# My randomness term
σ = 1.0 #variance in randomness
r = 0.5 #proportion of time in which it jumps left or right
β = 0.0 #probability of being the value of the previous lag or mean reversion strength
lag = 10 #lag
do_random_walk = false #behave like a random walk
myRandomnessTerm = RandomnessTerm(σ,r,β,lag,do_random_walk,false)

Δx = L / M  # real gap between simulation points 
Δt = (r * (Δx^2) / (2.0 * D))^(1/γ)

# RL Stuff:
T = 20
RealKickStartTime = 8 # when, in real time, to kick the system
SimKickStartTime = to_simulation_time(RealKickStartTime,Δt)-2 # convert to simulation time
Position = -2 
Volume = 10

myRLPusherPush = RLPushTerm(SimKickStartTime,SimKickStartTime+1,Position,Volume,true)
myRLPusherNoPush = RLPushTerm(SimKickStartTime,SimKickStartTime+1,Position,Volume,false)

lob_model_push = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
        mySourceTerm, myCouplingTerm, myRLPusherPush,myRandomnessTerm,shift=-1,old_way=true);

lob_model_no_push = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
        mySourceTerm, myCouplingTerm, myRLPusherNoPush,myRandomnessTerm,shift=-1,old_way=true);

print((Δt,to_simulation_time(T,Δt),num_paths*to_simulation_time(T,Δt))) #about 2GB RAM per 100K, i.e. can only do about 1.8 million
lob_model_push.SK_DP
# -

# check something actually happens for one example
if true
    myCouplingTerm = CouplingTerm(μ, a, b, c, false);
    
    Volume = 10

    myRLPusherPush = RLPushTerm(SimKickStartTime,SimKickStartTime+1,Position,Volume,true)
    myRLPusherNoPush = RLPushTerm(SimKickStartTime,SimKickStartTime+1,Position,Volume,false)
    
    lob_model_push = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
        mySourceTerm, myCouplingTerm, myRLPusherPush,myRandomnessTerm,shift=-1,michaels_way=true);
    
    lob_model_no_push = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
        mySourceTerm, myCouplingTerm, myRLPusherNoPush,myRandomnessTerm,shift=-1,michaels_way=true);
    
    quick_plot([lob_model_push,lob_model_no_push],SimKickStartTime,12)
    
    #png("/home/derickdiana/Desktop/Masters/Reworked/KickTheSytemWithRemoval.png")
    #png("/home/derickdiana/Desktop/Masters/Reworked/NumericalInstabilityForLargeRemovalRate.png")
    plot!()
end


plot_sums(Dat)

# +
gammas = [1.0]#[1.0,0.9,0.8,0.7,0.6]#[1.0,0.9,0.8,0.7,0.6]
volumes = range(1,100,length=100)

function get_set_inner(volume,γ,do_michaels_way,ν)
    Δt = (r * (Δx^2) / (2.0 * D * myscale))^(1/γ)
        
    l = Int(round(to_simulation_time(1,Δt)/3,digits=0))

    RealKickStartTime = 8 # when, in real time, to kick the system
    SimKickStartTime = to_simulation_time(RealKickStartTime,Δt)-2 # convert to simulation time

    #myCouplingTerm = CouplingTerm(μ, a, b, c, false);

    myRLPusherPush = RLPushTerm(SimKickStartTime,SimKickStartTime+1,Position,volume,true)
    myRLPusherNoPush = RLPushTerm(SimKickStartTime,SimKickStartTime+1,Position,volume,false)

    lob_model_push = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
        mySourceTerm, myCouplingTerm, myRLPusherPush,myRandomnessTerm, michaels_way = do_michaels_way);

    lob_model_no_push = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
        mySourceTerm, myCouplingTerm, myRLPusherNoPush,myRandomnessTerm, michaels_way = do_michaels_way);
    
    return ([lob_model_push, lob_model_no_push], SimKickStartTime, l)
    
end

function get_set(volume,γ)
    return get_set_inner(volume,γ,true,2.0)
end
(mean_price_impacts_frac_no_random_michaels_way,var_price_impacts_frac_no_random_michaels_way) = calculate_price_impacts(volumes,gammas,get_set)

function get_set(volume,γ)
    return get_set_inner(volume,γ,false,2.0)
end
(mean_price_impacts_frac_no_random_new_way,var_price_impacts_frac_no_random_new_way) = calculate_price_impacts(volumes,gammas,get_set)

function get_set(volume,γ)
    return get_set_inner(volume,γ,true,2.0*6/10)
end
(mean_price_impacts_frac_no_random_nu_adjust_michaels_way,var_price_impacts_frac_no_random_nu_adjust_michaels_way) = calculate_price_impacts(volumes,gammas,get_set);


# +
for Gamma in 1:(gi_len)
    fit_and_plot_price_impact(volumes,mean_price_impacts_frac_no_random_michaels_way,var_price_impacts_frac_no_random_michaels_way,["Michaels with nu=2"];colors=["blue"],new_plot=true)
    fit_and_plot_price_impact(volumes,mean_price_impacts_frac_no_random_new_way,var_price_impacts_frac_no_random_new_way,["New way with nu=2"];colors=["red"])
    fit_and_plot_price_impact(volumes,mean_price_impacts_frac_no_random_nu_adjust_michaels_way,var_price_impacts_frac_no_random_nu_adjust_michaels_way,["Michaels moved closer to new way by setting nu = 1.3"];colors=["green"])
end
plot!()

#png("/home/derickdiana/Desktop/Masters/Reworked/PriceImpactFractionalAdaptiveDelayOnlyNoRandom.png")
#png("/home/derickdiana/Desktop/Masters/Reworked/PriceImpactMichaelsCode.png")
#png("/home/derickdiana/Desktop/Masters/Reworked/75ScaleWithRandom.png")
#png("/home/derickdiana/Desktop/Masters/Reworked/ThreeWays.png")
# -
# # NO RANDOM KICKS PRICE IMPACTS

# +
# Configuration Arguments
num_paths = 1#50#30

L = 200     # real system width (e.g. 200 meters)
M = 400     # divided into M pieces , 400

p₀ = 230.0  #this is the mid_price at t=0  238.75 

# Free-Parameters for gaussian version
D = 0.5#0.5/8 # real diffusion constant e.g. D=1 (meters^2 / second), 1
α = 0.0 # legacy, no longer used

#ν = 3.0 #removal rate
#γ = 1.0 #fraction of derivative (1 is normal diffusion, less than 1 is D^{1-γ} derivative on the RHS)

ν = 1.0#1.0#3.0 #removal rate
γ = 1.0 #fraction of derivative (1 is normal diffusion, less than 1 is D^{1-γ} derivative on the RHS)

# Source term:
λ = 1.0 #
μ = 0.1 #
mySourceTerm = SourceTerm(λ, μ, true);

# Coupling term:
a = 0.01  #gap between stocks before at full strength: strong is 0.3
b = 2.0   #weighting of interaction term: strong is 2
c = 2.0   #skew factor: strong is 2

myCouplingTerm = CouplingTerm(μ, a, b, c, false);

# My randomness term
σ = 1.0 #variance in randomness
r = 0.5 #proportion of time in which it jumps left or right
β = 0.0 #probability of being the value of the previous lag or mean reversion strength
lag = 10 #lag
do_random_walk = false #behave like a random walk
myRandomnessTerm = RandomnessTerm(σ,r,β,lag,do_random_walk,false)

Δx = L / M  # real gap between simulation points 
Δt = (r * (Δx^2) / (2.0 * D))^(1/γ)

# RL Stuff:
T = 20
RealKickStartTime = 8 # when, in real time, to kick the system
SimKickStartTime = to_simulation_time(RealKickStartTime,Δt)-2 # convert to simulation time
Position = -1 
Volume = 10#10

myRLPusherPush = RLPushTerm(SimKickStartTime,SimKickStartTime+1,Position,Volume,true)
myRLPusherNoPush = RLPushTerm(SimKickStartTime,SimKickStartTime+1,Position,Volume,false)

lob_model_push = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
        mySourceTerm, myCouplingTerm, myRLPusherPush,myRandomnessTerm,shift=-1,old_way=true);

lob_model_no_push = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
        mySourceTerm, myCouplingTerm, myRLPusherNoPush,myRandomnessTerm,shift=-1,old_way=true);

print((Δt,to_simulation_time(T,Δt),num_paths*to_simulation_time(T,Δt))) #about 2GB RAM per 100K, i.e. can only do about 1.8 million
lob_model_push.SK_DP
# -

# check something actually happens for one example
if true
    myCouplingTerm = CouplingTerm(μ, a, b, c, false);
    
    Volume = 100

    myRLPusherPush = RLPushTerm(SimKickStartTime,SimKickStartTime+1,Position,Volume,true)
    myRLPusherNoPush = RLPushTerm(SimKickStartTime,SimKickStartTime+1,Position,Volume,false)
    
    lob_model_push = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
        mySourceTerm, myCouplingTerm, myRLPusherPush,myRandomnessTerm,shift=-1,michaels_way=true);
    
    lob_model_no_push = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
        mySourceTerm, myCouplingTerm, myRLPusherNoPush,myRandomnessTerm,shift=-1,michaels_way=true);
    
    Dat = quick_plot([lob_model_push,lob_model_no_push],SimKickStartTime,12)
    
    #png("/home/derickdiana/Desktop/Masters/Reworked/KickTheSytemWithRemoval.png")
    #png("/home/derickdiana/Desktop/Masters/Reworked/NumericalInstabilityForLargeRemovalRate.png")
    plot!()
end


plot_sums(Dat)

# ## Price impact for different D, gamma and nu

# +
function get_set(volume,combined_slice)
    D = combined_slice[1]
    ν = combined_slice[2]
    γ = combined_slice[3]
    return get_set_inner(volume,γ,D,ν)
end

function get_set_inner(volume,γ,D,ν)
    Δt = (r * (Δx^2) / (2.0 * D))^(1/γ)
    #l = Int(round(to_simulation_time(1,Δt)/3,digits=0))
    l = 50
    SimKickStartTime = to_simulation_time(RealKickStartTime,Δt)-2 # convert to simulation time

    myRLPusherPush   = RLPushTerm(SimKickStartTime, SimKickStartTime+1, Position, volume,  true)
    myRLPusherNoPush = RLPushTerm(SimKickStartTime, SimKickStartTime+1, Position, volume, false)

    lob_model_push    = SLOB(num_paths, RealKickStartTime+to_real_time(l,Δt)+1, p₀, M, L, D, ν, α, γ,
        mySourceTerm, myCouplingTerm, myRLPusherPush,   myRandomnessTerm);

    lob_model_no_push = SLOB(num_paths, RealKickStartTime+to_real_time(l,Δt)+1, p₀, M, L, D, ν, α, γ,
        mySourceTerm, myCouplingTerm, myRLPusherNoPush, myRandomnessTerm);
    
    return ([lob_model_push, lob_model_no_push], SimKickStartTime, l)
    
end;
# -


# ### Change Ds

# +
volumes = range(1,100,length=100)

#what to try
nus = [2.0]
ds = [0.2,0.5,1.0,1.5]
gammas = [1.0]
#all combinations of the above
combined = collect(Iterators.product(nus,ds,gammas))[:]

(mean_price_impacts,var_price_impacts) = calculate_price_impacts(volumes,  combined,   get_set)
# -
my_labels = map(c -> string("Data: nu=",c[1]," and D=",c[2]),combined)
for Gamma in 1:length(gammas)
    fit_and_plot_price_impact(volumes,mean_price_impacts,var_price_impacts,my_labels;colors=["red","green","blue","purple"],new_plot=true,
        forplot=(legend=:bottomright,yticks=range(0,1,step=L/M)))
end
plot!()

# ### Change nu's

# +
volumes = range(1,100,length=100)

#what to try
ds = [1.0]
nus = [0.8,1.5,2.0,2.5]
gammas = [1.0]
#all combinations of the above
combined = collect(Iterators.product(nus,ds,gammas))[:]

(mean_price_impacts,var_price_impacts) = calculate_price_impacts(volumes,  combined,   get_set)
# -

my_labels = map(c -> string("Data: nu=",c[1]," and D=",c[2]),combined)
for Gamma in 1:length(gammas)
    fit_and_plot_price_impact(volumes,mean_price_impacts,var_price_impacts,my_labels;colors=["red","green","blue","purple"],new_plot=true,
        forplot=(legend=:bottomright,yticks=range(-2,10,step=L/M)))
end
plot!()
# ### Change gamma's

# ### Change gammas

# +
volumes = range(1,5000,length=100)

#what to try
ds = [1.0]
nus = [1.0]
gammas = [1.0,0.9,0.8,0.7]
#all combinations of the above
combined = collect(Iterators.product(nus,ds,gammas))[:]

(mean_price_impacts,var_price_impacts) = calculate_price_impacts(volumes,  combined,   get_set)
# -
my_labels = map(c -> string("Data: gamma=",c[3]),combined)
for Gamma in 1:length(gammas)
    fit_and_plot_price_impact(volumes,mean_price_impacts,var_price_impacts,my_labels;colors=["red","green","blue","purple"],new_plot=true,
        forplot=(legend=:bottomright,))
end
plot!()

# ## Price impact for different delays

# +
function get_set(volume,combined_slice)
    l = combined_slice[1]
    return get_set_inner(volume,l)
end

function get_set_inner(volume,l)
    L = 200
    M = 400
    ν = 0.5
    Δt = (r * (Δx^2) / (2.0 * D))^(1/γ)
    #l = Int(round(to_simulation_time(1,Δt)/3,digits=0))
    SimKickStartTime = to_simulation_time(RealKickStartTime,Δt)-2 # convert to simulation time

    myRLPusherPush   = RLPushTerm(SimKickStartTime, SimKickStartTime+1, -1, volume,  true)
    myRLPusherNoPush = RLPushTerm(SimKickStartTime, SimKickStartTime+1, -1, volume, false)

    lob_model_push    = SLOB(num_paths, RealKickStartTime+to_real_time(l,Δt)+1, p₀, M, L, D, ν, α, γ,
        mySourceTerm, myCouplingTerm, myRLPusherPush,   myRandomnessTerm);

    lob_model_no_push = SLOB(num_paths, RealKickStartTime+to_real_time(l,Δt)+1, p₀, M, L, D, ν, α, γ,
        mySourceTerm, myCouplingTerm, myRLPusherNoPush, myRandomnessTerm);
    
    return ([lob_model_push, lob_model_no_push], SimKickStartTime, l)
    
end
# -


# #### Structure of kinks

# +
#volumes = range(0.1,10000,length=100)
#volumes = exp.(range(log(0.1),log(1000),length=100))
volumes = exp.(range(log(1),log(20000),length=10000))

#what to try
ls = [4,5]#1:5
#all combinations of the above
combined = collect(Iterators.product(ls))[:]

(mean_price_impacts,var_price_impacts) = calculate_price_impacts(volumes,  combined,   get_set)
# +
my_labels = map(l -> string("Data: l=",l[1]),combined)

fit_and_plot_price_impact(volumes,mean_price_impacts,var_price_impacts,my_labels;
    new_plot=true,do_shift=true,forplot=(legend=:bottomright,yticks=range(0,5,step=1/2)),do_log_fit=true,do_log_plot=true,do_log_shift=true)
hline!([1/2*mi for mi in 1:6],color="black",linewidth=0.3)


#save_png("Reworked","SubstructureToKinksLog")
plot!()
# -

# #### Many different delays

# +
#volumes = range(0.1,10000,length=100)
#volumes = exp.(range(log(0.1),log(100),length=1000))
volumes = exp.(range(-1,13,length=1000))
#volumes = exp.(range(log(0.1),log(1000),length=1000))
#volumes = exp.(range(log(1),log(10000),length=10000))

#what to try
ls = 1:7
#all combinations of the above
combined = collect(Iterators.product(ls))[:]

(mean_price_impacts,var_price_impacts) = calculate_price_impacts(volumes,  combined,   get_set)
# +
my_labels = map(l -> string("Data: l=",l[1]),combined)

fit_and_plot_price_impact(volumes,mean_price_impacts,var_price_impacts,my_labels;
    new_plot=true,forplot=(legend=:topleft,yticks=range(0,5,step=1/2),
    xticks=range(0,100,step=1)),do_log_fit=true,do_log_plot=true,do_vert_log_shift=true,do_horiz_log_shift=true)
hline!([1/2*mi for mi in 0:6],color="black",linewidth=0.3)
#save_png("Reworked","KinksAsFunctionOfDelayLog")

plot!()
# -

# ## Price impact for different divisions of space

# +
function get_set(volume,combined_slice)
    M = combined_slice[1]
    return get_set_inner(volume,M)
end

function get_set_inner(volume,M)
    D = 0.7
    Δt = (r * (Δx^2) / (2.0 * D))^(1/γ)
    #l = Int(round(to_simulation_time(1,Δt)/3,digits=0))
    #l = 3
    l = 100
    SimKickStartTime = to_simulation_time(RealKickStartTime,Δt)-2 # convert to simulation time

    myRLPusherPush   = RLPushTerm(SimKickStartTime, SimKickStartTime+1, Position, volume,  true)
    myRLPusherNoPush = RLPushTerm(SimKickStartTime, SimKickStartTime+1, Position, volume, false)

    lob_model_push    = SLOB(num_paths, RealKickStartTime + to_real_time(l,Δt) + 1, p₀, M, L, D, ν, α, γ,
        mySourceTerm, myCouplingTerm, myRLPusherPush,   myRandomnessTerm);

    lob_model_no_push = SLOB(num_paths, RealKickStartTime + to_real_time(l,Δt) + 1, p₀, M, L, D, ν, α, γ,
        mySourceTerm, myCouplingTerm, myRLPusherNoPush, myRandomnessTerm);
    
    return ([lob_model_push, lob_model_no_push], SimKickStartTime, l)
    
end


# +
#volumes = range(0.1,10000,length=100)
volumes = exp.(range(log(1),log(10),length=100))

#what to try
ms = [400,800,1600]#[100,200,400,800]
#all combinations of the above
combined = collect(Iterators.product(ms))[:]

(mean_price_impacts,var_price_impacts) = calculate_price_impacts(volumes,  combined,   get_set)
# +
my_labels = map(l -> string("Data: l=",l[1]),combined)

fit_and_plot_price_impact(volumes,mean_price_impacts,var_price_impacts,my_labels;new_plot=true,do_shift=true,
    forplot=(legend=:bottomright,))

plot!()
# -

# ## Price impact for different couplings

# ## Stylized facts

# +
# Configuration Arguments
T = 10000 
num_paths = 1 

L = 200     # real system width (e.g. 200 meters)
M = 400     # divided into M pieces , 400

p₀ = 230.0  #this is the mid_price at t=0  238.75 

# Free-Parameters for gaussian version
D = 0.5 # real diffusion constant e.g. D=1 (meters^2 / second), 1
α = 0.0 # legacy, no longer used

ν = 5.0 #removal rate
γ = 1.0 #fraction of derivative (1 is normal diffusion, less than 1 is D^{1-γ} derivative on the RHS)

dist = Normal(0.0,σ) #dist = TDist(1) #dist = Spl(1);

# Source term:
λ = 1.0
μ = 0.1 

mySourceTerm = SourceTerm(λ, μ, true);

# Coupling term:
myCouplingTerm = CouplingTerm(0.0, 0.0, 0.0, 0.0, false);

Δx = L / M  # real gap between simulation points 
Δt = (r * (Δx^2) / (2.0 * D))^(1/γ)

# My randomness term
σ = 1.0 #variance in randomness
r = 0.5 #proportion of time in which it jumps left or right
β = 0.0 #probability of being the value of the previous lag or mean reversion strength
lag = 10 #lag
do_random_walk = false #behave like a random walk
myRandomnessTerm = RandomnessTerm(σ,r,β,lag,do_random_walk,true)

# RL Stuff:
RealStartTime = 6 # when, in real time, to kick the system
SimStartTime = to_simulation_time(RealStartTime,Δt)-2 # convert to simulation time
Position = 0
# Volume set below

myRLPusher = RLPushTerm(SimStartTime,SimEndTime,Position,Volume,false)

lob_model = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ, 
    mySourceTerm, myCouplingTerm, myRLPusher, myRandomnessTerm);

print((Δt,to_simulation_time(T,Δt),num_paths*to_simulation_time(T,Δt))) #about 2GB RAM per 100K, i.e. can only do about 1.8 million
lob_model.SK_DP

# +
# clear everything pointed to by the dictionary then garbage collect. If it wasn't assigned yet it will left you know
try clear_double_dict(Dat) catch e print("Not initialized") end
GC.gc()

Dat = InteractOrderBooks([lob_model¹,lob_model²], -1, true) ;
#Dat = InteractOrderBooks([lob_model¹], -1, true) ;

# +
observed_price_path = Dat[1][1].obs_price_paths[1:end-1,1];
observed_log_returns = diff(log.(observed_price_path[:,1]));

data_stylized_facts = StylizedFacts.StylizedFactsPlot(observed_price_path);

StylizedFacts.plot_all_stylized_facts(data_stylized_facts,(1000,1200))

#png("/home/derickdiana/Desktop/Masters/OfficialOne/DraftOnePics/StylizedFactsRandomKicksNoFractional.png")
plot!()
# -

# # SINGLE POINT DIFFUSION

# ## Plot spike variance for various gamma

# +
# Configuration Arguments
num_paths = 1#30

L = 200     # real system width (e.g. 200 meters)
M = 400     # divided into M pieces , 400

p₀ = 230.0  #this is the mid_price at t=0  238.75 

# Free-Parameters for gaussian version
D = 0.4 # real diffusion constant e.g. D=1 (meters^2 / second), 1
α = 0.0 # legacy, no longer used

ν = 0.0 #removal rate
γ = 0.4 #fraction of derivative (1 is normal diffusion, less than 1 is D^{1-γ} derivative on the RHS)

# Source term:

mySourceTerm = SourceTerm(0.0, 0.0, false);

myCouplingTerm = CouplingTerm(0.0, 0.0, 0.0, 0.0, false);

# My randomness term
σ = 6.0 #variance in randomness
r = 0.5 #proportion of time in which it jumps left or right
β = 0.0 #probability of being the value of the previous lag or mean reversion strength
lag = 10 #lag
do_random_walk = false #behave like a random walk
myRandomnessTerm = RandomnessTerm(σ,r,β,lag,do_random_walk,false)


Δx = L / M  # real gap between simulation points 
Δt = (r * (Δx^2) / (2.0 * D))^(1/γ)

# RL Stuff:
RealKickStartTime = 1 # when, in real time, to kick the system
SimKickStartTime = to_simulation_time(RealKickStartTime,Δt)-2 # convert to simulation time
Position = Int(M/2)
Volume = -80; # If position == -x where x>=0, then put it x above the mid price each time

myRLPusher = RLPushTerm(SimKickStartTime,SimKickStartTime+1,Position,Volume,true)

T = 1000
lob_model = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ, 
    mySourceTerm, myCouplingTerm, myRLPusher, myRandomnessTerm,kernel_cut_off=0.000001);

RealKickStartTime = to_real_time(lob_model.cut_off,Δt) # when, in real time, to kick the system
SimKickStartTime = to_simulation_time(RealKickStartTime,Δt)-2 # convert to simulation time
Position = Int(M/2)
Volume = -80; # If position == -x where x>=0, then put it x above the mid price each time

myRLPusher = RLPushTerm(SimKickStartTime,SimKickStartTime+1,Position,Volume,true)

T = RealKickStartTime+1
lob_model = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ, 
    mySourceTerm, myCouplingTerm, myRLPusher, myRandomnessTerm,kernel_cut_off=0.000001);

MinRealKickStartTime = RealKickStartTime# when, in real time, to kick the system. Starts late to give system time to relax. Must be latest of all the gammas to try
print((Δt,to_simulation_time(T,Δt),num_paths*to_simulation_time(T,Δt))) #about 2GB RAM per 100K, i.e. can only do about 1.8 million
lob_model.SK_DP

MinRealKickStartTime
# -

quick_plot([lob_model],SimKickStartTime,12; dosum=false,plot_raw_price=true,x_axis_width=10)

plot_sums(Dat)

test = Dat[1][1].lob_densities[:,SimKickStartTime+5]
t_mean = mymean(test/sum(test) ,lob_model.x,lob_model.Δx)
t_var = myvar(test/sum(test) ,lob_model.x,lob_model.Δx)
print(t_mean," ",t_var^0.5)
plot(lob_model.x,test/sum(test),label="Data")
plot!(xlim=[t_mean-10,t_mean+10])
plot!([t_mean],seriestype="vline",color="red",label="mean")
plot!([t_mean - t_var^0.5],seriestype="vline",color="green",label="1 st.dev.")
plot!([t_mean + t_var^0.5],seriestype="vline",color="green",label="")

# +
# works with kernel_cut_off =0.000001
#gamma_indices = [1.0,0.9,0.8,0.7,0.6,0.5,0.4]
#ν = 0.0 #removal rate
#γ above set to 0.6
# ... up to the fact that the original thing is not linear

gamma_indices = [1.0]#,0.9,0.8,0.7,0.6,0.5,0.4]
gi_len = length(gamma_indices)

delta_ts = map(γ->(r * (Δx^2) / (2.0 * D))^(1/γ),gamma_indices)

how_long = 20
min_Δt = minimum(delta_ts,dims=1)[1]
sim_how_long = to_simulation_time(how_long,min_Δt)

lob_model = SLOB(num_paths, T, p₀, M, L, D, ν, α, minimum(gamma_indices,dims=1)[1], 
    mySourceTerm, myCouplingTerm, myRLPusher, myRandomnessTerm,kernel_cut_off=0.000001);
MinRealKickStartTime = to_real_time(lob_model.cut_off,min_Δt) + 1# when, in real time, to kick the system. Starts late to give system time to relax. Must be latest of all the gammas to try


variances = zeros(Float64,gi_len,sim_how_long-1)
sums = zeros(Float64,gi_len,sim_how_long-1)

for Gamma in 1:gi_len
    γ = gamma_indices[Gamma] #fraction of derivative (1 is normal diffusion, less than 1 is D^{1-γ} derivative on the RHS)
    Δt = delta_ts[Gamma]
    
    RealKickStartTime = max(26,MinRealKickStartTime)
    SimKickStartTime = to_simulation_time(RealKickStartTime,Δt)-2 # convert to simulation time

    myRLPusher = RLPushTerm(SimKickStartTime,SimKickStartTime+1,Position,-8,true)
    
    lob_model = SLOB(num_paths, RealKickStartTime+how_long+1, p₀, M, L, D, ν, α, γ,
        mySourceTerm, myCouplingTerm, myRLPusher, myRandomnessTerm; shift=-1 , kernel_cut_off=0.000001 )
    
    try clear_double_dict(Dat) catch e print("Not initialized") end
    GC.gc()

    Dat = InteractOrderBooks([lob_model], -1, true) 
    
    sim_how_long = to_simulation_time(how_long,delta_ts[Gamma])
    for t in 1:sim_how_long-1
        temp = Dat[1][1].lob_densities[:,SimKickStartTime+t] #across all x values, at time "t",  on path 1
        sums[Gamma,t] = sum(temp)
        variances[Gamma,t] = myvar(temp./sum(temp),lob_model.x,lob_model.Δx)
    end
    
end

# +
l2 =  @layout [a;b]
plot()
#l = to_simulation_time(1,maximum(delta_ts[Gamma],dims=1)[1])
for Gamma in 1:gi_len
    l = to_simulation_time(how_long,delta_ts[Gamma])
    sub = 2:l-1
    real_time = sub.*delta_ts[Gamma]
    if ν!=0
        slope = round((variances[Gamma,sub[end]]-variances[Gamma,sub[1]])/(l*delta_ts[Gamma]),digits=4)
        plot!(real_time,variances[Gamma,sub],label=string("Variance with alpha=",gamma_indices[Gamma]," has slope ",slope),color=Gamma)
    else
        a,b = power_fit(real_time,variances[Gamma,sub])
        scatter!(real_time,variances[Gamma,sub],label=string("Variance with alpha=",gamma_indices[Gamma]),color=Gamma,markersize=1,markerstrokewidth=0.2)
        plot!(real_time,(2*D)/gamma(1+gamma_indices[Gamma]).*real_time.^gamma_indices[Gamma],label=string("Theoretical fit for alpha=",gamma_indices[Gamma]," is y=",round((2*D)/gamma(1+gamma_indices[Gamma]),digits=2),"t^",gamma_indices[Gamma]),ls=:solid,color=Gamma)
        plot!(real_time,a.*real_time.^b,label=string("Numeric fit for alpha=",gamma_indices[Gamma]," is y=",round(a,digits = 2),"t^",round(b,digits=2)),ls=:dash,color=Gamma)
    end
        
    #plt2 = plot(sums[Gamma,sub],label=string("Area underneath has myrange ",maximum(sums[sub]) - minimum(sums[sub])),color="red",ylim=[0.99999,1.00001])
    #plot!(plt1, plt2, layout=l2, size=(700,700),xlab="Time")
end

plot!(xlab="Real time (t)",ylab="Variance of spike",size=(800,800))
#png("/home/derickdiana/Desktop/Masters/Reworked/DiffusionForManyGammas.png")
plot!()
# -

# # COMPARE DIFFERENT VERSIONS

# +
# Configuration Arguments
num_paths = 1#30

L = 200     # real system width (e.g. 200 meters)
M = 400     # divided into M pieces , 400

p₀ = 230.0  #this is the mid_price at t=0  238.75 

# Free-Parameters for gaussian version
D = 0.5#0.5/8 # real diffusion constant e.g. D=1 (meters^2 / second), 1
α = 0.0 # legacy, no longer used

ν = 2.0 #removal rate
γ = 1.0 #fraction of derivative (1 is normal diffusion, less than 1 is D^{1-γ} derivative on the RHS)

T = 3000

# Source term:
λ = 1.0 #
μ = 0.1 #
mySourceTerm = SourceTerm(λ, μ, true);

# Coupling term:
a = 13.0  #gap between stocks before at full strength: strong is 0.3
b = 1.0   #weighting of interaction term: strong is 2
c = 1.2   #skew factor: strong is 2

myCouplingTerm = CouplingTerm(μ, a, b, c, false);

# My randomness term
σ = 1.0 #variance in randomness
r = 0.5 #proportion of time in which it jumps left or right
β = 0.0 #probability of being the value of the previous lag or mean reversion strength
lag = 10 #lag
do_random_walk = false #behave like a random walk
myRandomnessTerm = RandomnessTerm(σ,r,β,lag,do_random_walk,false)

Δx = L / M  # real gap between simulation points 
Δt = (r * (Δx^2) / (2.0 * D))^(1/γ)

# RL Stuff:
RealKickStartTime = 8 # when, in real time, to kick the system
SimKickStartTime = to_simulation_time(RealKickStartTime,Δt)-2 # convert to simulation time
Position = 0
Volume = 20

myRLPusher = RLPushTerm(SimKickStartTime,SimKickStartTime+8,Position,Volume,true)

lob_model = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ, 
    mySourceTerm, myCouplingTerm, myRLPusher, myRandomnessTerm,shift=-1,old_way=true);

print((Δt,to_simulation_time(T,Δt),num_paths*to_simulation_time(T,Δt))) #about 2GB RAM per 100K, i.e. can only do about 1.8 million
lob_model.SK_DP
# -

# check something actually happens for one example
if true
    try clear_double_dict(Dat) catch e print("Not initialized") end
    GC.gc()
    
    lob_model = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ, 
        mySourceTerm, myCouplingTerm, myRLPusher, myRandomnessTerm,michaels_way=true);
    TrulyOldWay = InteractOrderBooks([lob_model], -1, true) ;
    
    lob_model = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ, 
        mySourceTerm, myCouplingTerm, myRLPusher, myRandomnessTerm,michaels_way=false,shift=-1,old_way=true);
    DatSpecialCase = InteractOrderBooks([lob_model], -1, true) ;
    
    lob_model = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ, 
        mySourceTerm, myCouplingTerm, myRLPusher, myRandomnessTerm,michaels_way=false,shift=0,old_way=false);
    DatNewWayWithNoShift = InteractOrderBooks([lob_model], -1, true) ;
    
    lob_model = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ, 
        mySourceTerm, myCouplingTerm, myRLPusher, myRandomnessTerm,michaels_way=false,shift=-1,old_way=false);
    DatNewWayWithShift = InteractOrderBooks([lob_model], -1, true) ;

    
    sub_x = Int(3M/8):Int(5*M/8)
    how_many = 1
    p_arr1 = Array{Plots.Plot{Plots.GRBackend},1}(undef,how_many)
    start_pos =  SimKickStartTime+14
    #start_pos = to_simulation_time(T,Δt) - how_many - 1
    for i in [start_pos:(start_pos+how_many-1);]
        #p_arr1[i] = plot_density_visual(SimKickStartTime-2+i, to_real_time(SimKickStartTime-2+i,Δt), 1, Dat,false, true, 10, Dat[1][1].raw_price_paths[SimKickStartTime])
        
        plot()
        plot(TrulyOldWay[1][1].slob.x[sub_x], TrulyOldWay[1][1].lob_densities[sub_x,i], color=1,ls=:solid,lw=3,label="Old way (Michael's code exactly)");
        plot!(DatSpecialCase[1][1].slob.x[sub_x], DatSpecialCase[1][1].lob_densities[sub_x,i], color=2,ls=:dash,lw=2,label="New way analytically reduced to old (bound correction and only keep first order in Dt)")
        plot!(DatNewWayWithNoShift[1][1].slob.x[sub_x], DatNewWayWithNoShift[1][1].lob_densities[sub_x,i], color=3,ls=:dash,lw=2,label="New way with no bound correction (no bound correction and all orders in Dt)")
        plot!(DatNewWayWithShift[1][1].slob.x[sub_x], DatNewWayWithShift[1][1].lob_densities[sub_x,i], color=4,ls=:dash,lw=2,label="New way with bound correction (bound correction and all orders in Dt)")
        
        scatter!(TrulyOldWay[1][1].slob.x[sub_x], TrulyOldWay[1][1].lob_densities[sub_x,i], color=1,label="",markersize=4,markerstrokewidth=0.0);
        scatter!(DatSpecialCase[1][1].slob.x[sub_x], DatSpecialCase[1][1].lob_densities[sub_x,i], color=2,label="",markersize=4,markerstrokewidth=0.0)
        scatter!(DatNewWayWithNoShift[1][1].slob.x[sub_x], DatNewWayWithNoShift[1][1].lob_densities[sub_x,i], color=3,label="",markersize=4,markerstrokewidth=0.0)
        scatter!(DatNewWayWithShift[1][1].slob.x[sub_x], DatNewWayWithShift[1][1].lob_densities[sub_x,i], color=4,label="",markersize=4,markerstrokewidth=0.0)
        
        p_arr1[i-start_pos+1] = plot!()
    end
    plot(p_arr1...,size=(1200,1000))
    
    #png("/home/derickdiana/Desktop/Masters/Reworked/CompareDifferentSchemesWithoutRandomness.png")
    plot!()
end


# # Junk

# +
# slob = lob_model¹
# num_time_steps = to_simulation_time(T,Δt)

# lob_densities¹ =    zeros(Float64, slob.M+1, num_time_steps + 1, 2, num_paths)
# sources¹ =          zeros(Float64, slob.M+1, num_time_steps + 1, 2, num_paths)
# couplings¹ =        zeros(Float64, slob.M+1, num_time_steps + 1, 2, num_paths)
# rl_pushes¹ =        zeros(Float64, slob.M+1, num_time_steps + 1, 2, num_paths)

# raw_price_paths¹ =  ones(Float64,           num_time_steps + 1, 2, num_paths)
# sample_price_paths¹ =  ones(Float64,           slob.T         + 1, 2, num_paths)

# lob_densities² =    zeros(Float64, slob.M+1, num_time_steps + 1, 2, num_paths)
# sources² =          zeros(Float64, slob.M+1, num_time_steps + 1, 2, num_paths)
# couplings² =        zeros(Float64, slob.M+1, num_time_steps + 1, 2, num_paths)
# rl_pushes² =        zeros(Float64, slob.M+1, num_time_steps + 1, 2, num_paths)

# raw_price_paths² =  ones(Float64,           num_time_steps + 1, 2, num_paths)
# sample_price_paths² =  ones(Float64,           slob.T         + 1, 2, num_paths)


# for path_num in 1:num_paths
#     lob_densities¹[:,:,path_num] = Dat[path_num][1].lob_densities[:,:]
#     sources¹[:,:,path_num] = Dat[path_num][1].sources[:,:,:]
#     couplings¹[:,:,path_num] = Dat[path_num][1].couplings[:,:,:]
#     rl_pushes¹[:,:,path_num] = Dat[path_num][1].rl_pushes[:,:,:]
#     raw_price_paths¹[:,path_num] = Dat[path_num][1].raw_price_paths[:,:]
#     sample_price_paths¹[:,path_num] = Dat[path_num][1].obs_price_paths[:,:]

#     lob_densities²[:,:,path_num] = Dat[path_num][2].lob_densities[:,:,:]
#     sources²[:,:,path_num] = Dat[path_num][2].sources[:,:,:]
#     couplings²[:,:,path_num] = Dat[path_num][2].couplings[:,:,:]
#     rl_pushes²[:,:,path_num] = Dat[path_num][2].rl_pushes[:,:,:]
#     raw_price_paths²[:,path_num] = Dat[path_num][2].raw_price_paths[:,:]
#     sample_price_paths²[:,path_num] = Dat[path_num][2].obs_price_paths[:,:];
# end

# + active=""
# Revise.revise()
# -

how_many = 12
#l2 = @layout [a d; b e; c f; ];
p_arr1 = Array{Plots.Plot{Plots.GRBackend},1}(undef,how_many)
for i in [1:how_many;]
p_arr1[i] = plot_density_visual(SimStartTime-2+i, to_real(SimStartTime-2+i), lob_model¹, lob_densities¹, couplings¹, sources¹, rl_pushes¹,raw_price_paths¹,true,10)
end
plot(p_arr1...,size=(1000,1000))
#gif(anim, "~/Desktop/Masters/GoodPics/ForM2.gif", fps=16)
#png("/home/derickdiana/Desktop/Masters/GoodPics/ForM2.png")

SimStartTime

# +
observed_price_path = sample_price_paths¹[1:end-1,1];
observed_log_returns = diff(log.(observed_price_path[:,1]));

data_stylized_facts = StylizedFacts.StylizedFactsPlot(observed_price_path);

# +
StylizedFacts.plot_all_stylized_facts(data_stylized_facts,(1000,1200))

#png("/home/derickdiana/Desktop/Masters/GoodPics/DefinitelyNoTails.png")

# +
function myvar(p,x,Δx)
    mu = mymean(p,x,Δx)
    Ex2 = 0
    for i in 1:length(p)
        Ex2 += p[i]*x[i]^2*Δx
    end
    return Ex2 - mu^2
end

function mymean(p,x,Δx)
    sum = 0
    for i in 1:length(p)
        sum += p[i]*x[i]*Δx
    end
    return sum
end

test = lob_densities¹[:,end-5,2]
test = test/sum(test)
t_mean = mymean(test ,lob_model¹.x,lob_model¹.Δx)
t_var = myvar(test ,lob_model¹.x,lob_model¹.Δx)
print(t_mean," ",t_var^0.5)
plot(lob_model¹.x,test)
plot!([t_mean],seriestype="vline")
plot!([t_mean - t_var^0.5],seriestype="vline")
plot!([t_mean + t_var^0.5],seriestype="vline")

# +
l = length(lob_densities¹[1,:,1]) # number of time steps

variances = zeros(Float64,l-1)
sums = zeros(Float64,l-1)
for t in 2:l 
    temp = lob_densities¹[:,t,1]/sum(lob_densities¹[:,t,1]) #across all x values, at time "t",  on path 1
    sums[t-1] = sum(temp)
    variances[t-1] = myvar(temp,lob_model¹.x,lob_model¹.Δx)
end

l2 =  @layout [a;b];
sub = 7:length(sums)-1
plt1 = plot(variances[sub],label="Variance")
plt2 = plot(sums[sub],label=string("Area underneath has myrange ",maximum(sums[sub]) - minimum(sums[sub])),color="red")
plot(plt1, plt2, layout=l2, size=(700,700),xlab="Time")
#png("/home/derickdiana/Desktop/Masters/StillWorking/VarianceAndSumsForAlphaIs0.5.png")

# +
#middle = 3:length(sample_price_paths²[:,1])
#price_changes = (sample_price_paths²[middle.-1,:] - sample_price_paths²[middle.-2,:])[:]

# +
#Plot Log Returns Hist
#dist = fit(Normal, observed_log_returns)
#plt = histogram(observed_log_returns, normalize=true, legend=false)
#plot!(plt, x -> pdf(dist, x), xlim=xlims(), legend=false)
#
#StylizedFacts.plot_log_returns(data_stylized_facts, "")
#StylizedFacts.plot_qq_log_returns(data_stylized_facts, "")
#StylizedFacts.plot_acf_order_flow(market_data_stylized_facts, "")
#StylizedFacts.plot_acf_log_returns(data_stylized_facts, "")
#StylizedFacts.plot_acf_abs_log_returns(data_stylized_facts, "")

# +
# save_object("OneRun.jld2",(lob_densities¹, sources¹, couplings¹, rl_pushes¹, raw_price_paths¹, sample_price_paths¹, P⁺s¹, P⁻s¹, Ps¹, 
# lob_densities², sources², couplings², rl_pushes², raw_price_paths², sample_price_paths², P⁺s², P⁻s², Ps²))
# temp = load_object("OneRun.jld2")

# +
# combine all different paths into one large vector
#diff(log.(sample_price_paths¹[1:end-1,:]),dims=1)[:]

# +
#observed_summary_stats = get_summary_stats(observed_log_returns) #is in AdaptiveABC

# +
# function fractional_choose(n,k)
#     if k==1
#         temp = 0
#     else 
#         temp = 0
#     end
    
#     if (n==0)
#         return 0.0 + temp
#     end
    
#     return gamma(n+1)/(gamma(n-k+1)*gamma(k+1)) - 1/(n+1)*1/beta(n-k+1,k+1)#*(-1)^k + temp
# end
# n = 200
# for m in 1:(n-1)
#     print(fractional_choose(1-0.8,n-m)," ")
# end
# -



# +
lob_densities¹, sources¹, couplings¹, rl_pushes¹, raw_price_paths¹, sample_price_paths¹, P⁺s¹, P⁻s¹, Ps¹, V¹,
lob_densities², sources², couplings², rl_pushes², raw_price_paths², sample_price_paths², P⁺s², P⁻s², Ps², V² = 
    [nothing for _ in 1:20]

GC.gc()

lob_densities¹, sources¹, couplings¹, rl_pushes¹, raw_price_paths¹, sample_price_paths¹, P⁺s¹, P⁻s¹, Ps¹, V¹,
lob_densities², sources², couplings², rl_pushes², raw_price_paths², sample_price_paths², P⁺s², P⁻s², Ps², V², broke_points =
    InteractOrderBooks(lob_model¹,lob_model², -1, true) ;

l = @layout [a d; b e; c f];

total_steps = 10#min(to_sim(T),100)
total_length = to_sim(T)
step = floor(Int,total_length/total_steps)

myrange = 1:step:(total_length-step)
#myrange = [SimStartTime:SimStartTime+to_sim(1)*2;]
#myrange = [SimStartTime]

p_outer = Progress(length(myrange),dt=0.1)

anim = @animate for s = myrange           #s is the time in simulation time
    r = to_real_time(s, lob_model¹.Δt)  #r is the time in real time
    #s = to_simulation_time(r, lob_model¹.Δt)  #s is the time in real time
    
    #global prev_r, plt1, plt2, plt3, plt4, plt5, plt6, sim_times, real_times
    
    
    plt1 = plot_price_path(s,r,lob_model¹,raw_price_paths¹,sample_price_paths¹)
    plt3 = plot_price_path(s,r,lob_model²,raw_price_paths²,sample_price_paths²)
    plt5 = plot_price_path(s,r,lob_model²,raw_price_paths¹-raw_price_paths²,sample_price_paths¹-sample_price_paths²)
    
    plt2 = plot_density_visual(s, r, lob_model¹, lob_densities¹, couplings¹, sources¹, rl_pushes¹,raw_price_paths¹)
    plt4 = plot_density_visual(s, r, lob_model², lob_densities², couplings², sources², rl_pushes²,raw_price_paths²)
    plt6 = plot_density_visual(s, r, lob_model², lob_densities¹+lob_densities², 
                                                 couplings¹+couplings², 
                                                 sources¹+sources², 
                                                 rl_pushes¹+rl_pushes²,raw_price_paths¹,false)
    
    plot(plt1, plt2, plt3, plt4, plt5, plt6 ,layout=l,size=(1000,1000))
    
    next!(p_outer)
end

gif(anim, "/tmp/LOB.gif", fps=20*length(myrange)/200)
#gif(anim, "~/Desktop/Masters/StillWorking/Random_Walk_And_Coupling_For_Alpha_Is_0.7.gif", fps=20*length(myrange)/200)
# -

# ## Plot stylized facts for gamma is 0.8

# +
# Configuration Arguments
T = 10000 
num_paths = 1 

L = 200     # real system width (e.g. 200 meters)
M = 400     # divided into M pieces , 400

p₀ = 230.0  #this is the mid_price at t=0  238.75 

# Free-Parameters for gaussian version
D = 0.5 # real diffusion constant e.g. D=1 (meters^2 / second), 1
σ = 1.0 #variance in randomness
α = 0.0 # legacy, no longer used

ν = 12.0 #removal rate
γ = 0.8 #fraction of derivative (1 is normal diffusion, less than 1 is D^{1-γ} derivative on the RHS)
r = 0.5 #proportion of time in which it jumps left or right

dist = Normal(0.0,σ) #dist = TDist(1) #dist = Spl(1);

Δx = L / M  # real gap between simulation points 
#Δt = (Δx^2) / (2.0 * D) / 10 #* (2.0/3.0) # real time seperation between simulation points
Δt = (r * (Δx^2) / (2.0 * D))^(1/γ)
print((Δt,to_simulation_time(T,Δt),num_paths*to_simulation_time(T,Δt))) #about 2GB RAM per 100K, i.e. can only do about 1.8 million
lob_model¹.SK_DP

# Source term:
λ = 1.0
μ = 0.1 

mySourceTerm = SourceTerm(λ, μ, true);

# Coupling term:
myCouplingTerm = CouplingTerm(0.0, 0.0, 0.0, 0.0, false);

# RL Stuff:
RealStartTime = 6 # when, in real time, to kick the system
SimStartTime = to_simulation_time(RealStartTime,Δt)-2 # convert to simulation time
Position = -1
# Volume set below

myRLPusher1 = RLPushTerm(SimStartTime,SimEndTime,Position,Volume,false)

lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, σ, ν, α, r, γ, dist, 
    mySourceTerm, myCouplingTerm, myRLPusher1,true);

myRLPusher2 = RLPushTerm(SimStartTime,SimEndTime,Position,Volume,false)

lob_model² = SLOB(num_paths, T, p₀, M, L, D, σ, ν, α, r, γ, dist,
    mySourceTerm, myCouplingTerm, myRLPusher2,true);

# +
lob_densities¹, sources¹, couplings¹, rl_pushes¹, raw_price_paths¹, sample_price_paths¹, P⁺s¹, P⁻s¹, Ps¹, V¹,
lob_densities², sources², couplings², rl_pushes², raw_price_paths², sample_price_paths², P⁺s², P⁻s², Ps², V² = 
    [nothing for _ in 1:20]

GC.gc()

lob_densities¹, sources¹, couplings¹, rl_pushes¹, raw_price_paths¹, sample_price_paths¹, P⁺s¹, P⁻s¹, Ps¹, V¹,
lob_densities², sources², couplings², rl_pushes², raw_price_paths², sample_price_paths², P⁺s², P⁻s², Ps², V², broke_points =
    InteractOrderBooks(lob_model¹,lob_model², -1, true) ;

# +
observed_price_path = sample_price_paths¹[1:end-1,1];
observed_log_returns = diff(log.(observed_price_path[:,1]));

data_stylized_facts = StylizedFacts.StylizedFactsPlot(observed_price_path);

StylizedFacts.plot_all_stylized_facts(data_stylized_facts,(1000,1200))

#png("/home/derickdiana/Desktop/Masters/OfficialOne/DraftOnePics/StylizedFactsRandomKicksFractionalAlphaIsPoint8.png")
plot!()
# -

# ## Price impact fixed delay and duration

# +
# Configuration Arguments
num_paths = 1

L = 200     # real system width (e.g. 200 meters)
M = 400     # divided into M pieces , 400

p₀ = 230.0  #this is the mid_price at t=0  238.75 

# Free-Parameters for gaussian version
D = 0.5 # real diffusion constant e.g. D=1 (meters^2 / second), 1
σ = 1.0 #variance in randomness
α = 0.0 # legacy, no longer used

ν = 14.0 #removal rate
γ = 0.6 #fraction of derivative (1 is normal diffusion, less than 1 is D^{1-γ} derivative on the RHS)
r = 0.5 #proportion of time in which it jumps left or right

dist = Normal(0.0,σ) #dist = TDist(1) #dist = Spl(1);

# Source term:
λ = 1.0
μ = 0.1 

mySourceTerm = SourceTerm(λ, μ, true);

# Coupling term:
myCouplingTerm = CouplingTerm(0.0, 0.0, 0.0, 0.0, false);

myRLPusher¹ = RLPushTerm(SimStartTime,SimStartTime+1,Position,-8,false)
lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, σ, ν, α, r, γ, dist,
    mySourceTerm, myCouplingTerm, myRLPusher¹,true);

myRLPusher² = RLPushTerm(SimStartTime,SimStartTime+1,Position,-8,false)
lob_model² = SLOB(num_paths, T, p₀, M, L, D, σ, ν, α, r, γ, dist,
    mySourceTerm, myCouplingTerm, myRLPusher²,true);

# -

Δx = L / M  # real gap between simulation points 
Δt = (r * (Δx^2) / (2.0 * D))^(1/γ)
print((Δt,to_simulation_time(T,Δt),num_paths*to_simulation_time(T,Δt))) #about 2GB RAM per 100K, i.e. can only do about 1.8 million
lob_model¹.SK_DP

# +
lob_densities¹, sources¹, couplings¹, rl_pushes¹, raw_price_paths¹, sample_price_paths¹, P⁺s¹, P⁻s¹, Ps¹, V¹,
lob_densities², sources², couplings², rl_pushes², raw_price_paths², sample_price_paths², P⁺s², P⁻s², Ps², V² = 
    [nothing for _ in 1:20]

GC.gc()

lob_densities¹, sources¹, couplings¹, rl_pushes¹, raw_price_paths¹, sample_price_paths¹, P⁺s¹, P⁻s¹, Ps¹, V¹,
lob_densities², sources², couplings², rl_pushes², raw_price_paths², sample_price_paths², P⁺s², P⁻s², Ps², V², broke_points =
    InteractOrderBooks(lob_model¹,lob_model², -1, true) ;

l = @layout [a d; b e; c f];
print(size(raw_price_paths²))



total_steps = 10#min(to_simulation_time(T,Δt),100)
total_length = to_simulation_time(T,Δt)
step = floor(Int,total_length/total_steps)

myrange = 1:step:(total_length-step)
#myrange = [SimStartTime:SimStartTime+to_sim(1)*2;]
#myrange = [SimStartTime]

p_outer = Progress(length(myrange),dt=0.1)

anim = @animate for s = myrange           #s is the time in simulation time
    r = to_real_time(s, lob_model¹.Δt)  #r is the time in real time
    #s = to_simulation_time(r, lob_model¹.Δt)  #s is the time in real time
    
    #global prev_r, plt1, plt2, plt3, plt4, plt5, plt6, sim_times, real_times
    
    
    plt1 = plot_price_path(s,r,lob_model¹,raw_price_paths¹,sample_price_paths¹)
    plt3 = plot_price_path(s,r,lob_model²,raw_price_paths²,sample_price_paths²)
    plt5 = plot_price_path(s,r,lob_model²,raw_price_paths¹-raw_price_paths²,sample_price_paths¹-sample_price_paths²)
    
    plt2 = plot_density_visual(s, r, lob_model¹, lob_densities¹, couplings¹, sources¹, rl_pushes¹,raw_price_paths¹)
    plt4 = plot_density_visual(s, r, lob_model², lob_densities², couplings², sources², rl_pushes²,raw_price_paths²)
    plt6 = plot_density_visual(s, r, lob_model², lob_densities¹+lob_densities², 
                                                 couplings¹+couplings², 
                                                 sources¹+sources², 
                                                 rl_pushes¹+rl_pushes²,raw_price_paths¹,false)
    
    plot(plt1, plt2, plt3, plt4, plt5, plt6 ,layout=l,size=(1000,1000))
    
    next!(p_outer)
end

gif(anim, "/tmp/LOB.gif", fps=20*length(myrange)/200)

# +
volume_indices = [1:50;] 
volume_indices = volume_indices./2
gamma_indices = [1.0,0.9,0.8,0.7,0.6]

vi_len = length(volume_indices)
gi_len = length(gamma_indices)
#mean_price_impacts = [0.0 for x in volume_indices]
#var_price_impacts = [0.0 for x in volume_indices]
#mean_log_price_impacts = [0.0 for x in volume_indices]
#var_log_price_impacts = [0.0 for x in volume_indices]
mean_price_impacts_frac = ones(Float64,gi_len,vi_len)
var_price_impacts_frac = ones(Float64,gi_len,vi_len)

Δt_largest = (r * (Δx^2) / (2.0 * D))^(1/1)

for Gamma in 1:gi_len 
    p_outer = Progress(vi_len,dt=0.1)
    
    γ = gamma_indices[Gamma]
    Δt = (r * (Δx^2) / (2.0 * D))^(1/γ)

    T = 10
    RealStartTime = 8 # when, in real time, to kick the system. Starts late to give system time to relax. Must be latest of all the gammas to try
    SimStartTime = to_simulation_time(RealStartTime,Δt)-2 # convert to simulation time
    
    print(2)
    
    for Volume in 1:vi_len
        volume = volume_indices[Volume]
        
        myRLPusher¹ = RLPushTerm(SimStartTime,SimStartTime+2,Position,volume,true)
        lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, σ, ν, α, r, γ, dist,
            mySourceTerm, myCouplingTerm, myRLPusher¹,true);

        myRLPusher² = RLPushTerm(SimStartTime,SimStartTime+2,Position,volume,true)
        lob_model² = SLOB(num_paths, T, p₀, M, L, D, σ, ν, α, r, γ, dist,
            mySourceTerm, myCouplingTerm, myRLPusher²,true);
        

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
        #print(size(raw_price_paths¹))
        price_impact¹ = raw_price_paths¹[SimStartTime+2,:] - raw_price_paths¹[SimStartTime,:]
        price_impact² = raw_price_paths²[SimStartTime+2,:] - raw_price_paths²[SimStartTime,:]
        price_impact = cat(price_impact¹,price_impact²,dims=1)
        mean_price_impacts_frac[Gamma,Volume] = mean(price_impact)
        var_price_impacts_frac[Gamma,Volume] = var(price_impact)


        #log_price_impact¹ = log.(sample_price_paths¹[RealStartTime,:]) - log.(sample_price_paths¹[RealStartTime-1,:])
        #log_price_impact² = log.(sample_price_paths²[RealStartTime,:]) - log.(sample_price_paths²[RealStartTime-1,:])
        #log_price_impact  = cat(log_price_impact¹,log_price_impact²,dims=1)

        #mean_log_price_impacts[Volume] = mean(log_price_impact)
        #var_log_price_impacts[Volume] = var(log_price_impact)

        next!(p_outer)
    end
end
# -

mean_price_impacts_frac = .- mean_price_impacts_frac;

# +
x = volume_indices
colors = ["blue","red","green","orange","purple"]
plot()
for Gamma in 1:gi_len
    a,b = log_fit(x,mean_price_impacts_frac[Gamma,:])
    c,d = power_fit(x,mean_price_impacts_frac[Gamma,:])

    #scatter(x,mean_price_impacts,label="data: p(t+1)-p(t)",ms=3, ma=1,yerr=var_price_impacts)
    scatter!(x,mean_price_impacts_frac[Gamma,:],label=string("data: gamma = ",gamma_indices[Gamma]),ms=3, ma=1,color=colors[Gamma])
    plot!(x,mean_price_impacts_frac[Gamma,:],ribbon=var_price_impacts_frac[Gamma,:],alpha=0,fillalpha=0.4,fillcolor=colors[Gamma],label="")
    plot!(x,a.+b.*log.(x),label="Log fit",w=1.5,color=colors[Gamma])
    plot!(x,c.*((x).^d),label="Power fit",w=1.5,color=colors[Gamma],linestyle=:dash)
    plot!(xlabel="Volume",ylabel="Price impact i.e. p(t+1)-p(t)")
    #png("/home/derickdiana/Desktop/Masters/OfficialOne/DraftOnePics/PriceImpactRandomKicksFractional0.7.png")
end

plot!(size=(1000,1000))
png("/home/derickdiana/Desktop/Masters/OfficialOne/DraftOnePics/PriceImpactFractionalFixedDelayAndDuration.png")
plot!(size=(1000,1000))
# -
# ## Price impact adaptive delay and duration

# +
volume_indices = [1:50;] 
volume_indices = volume_indices./2
gamma_indices = [1.0,0.9,0.8,0.7,0.6]

vi_len = length(volume_indices)
gi_len = length(gamma_indices)
#mean_price_impacts = [0.0 for x in volume_indices]
#var_price_impacts = [0.0 for x in volume_indices]
#mean_log_price_impacts = [0.0 for x in volume_indices]
#var_log_price_impacts = [0.0 for x in volume_indices]
mean_price_impacts_frac = ones(Float64,gi_len,vi_len)
var_price_impacts_frac = ones(Float64,gi_len,vi_len)

Δt_largest = (r * (Δx^2) / (2.0 * D))^(1/1)

for Gamma in 1:gi_len 
    p_outer = Progress(vi_len,dt=0.1)
    
    γ = gamma_indices[Gamma]
    Δt = (r * (Δx^2) / (2.0 * D))^(1/γ)

    T = 10
    RealStartTime = 8 # when, in real time, to kick the system. Starts late to give system time to relax. Must be latest of all the gammas to try
    SimStartTime = to_simulation_time(RealStartTime,Δt)-2 # convert to simulation time
    
    l = to_simulation_time(Δt_largest,Δt)
    print(l)
    
    for Volume in 1:vi_len
        volume = volume_indices[Volume]
        
        myRLPusher¹ = RLPushTerm(SimStartTime,SimStartTime+l,Position,volume,true)
        lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, σ, ν, α, r, γ, dist,
            mySourceTerm, myCouplingTerm, myRLPusher¹,true);

        myRLPusher² = RLPushTerm(SimStartTime,SimStartTime+l,Position,volume,true)
        lob_model² = SLOB(num_paths, T, p₀, M, L, D, σ, ν, α, r, γ, dist,
            mySourceTerm, myCouplingTerm, myRLPusher²,true);
        

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
        #print(size(raw_price_paths¹))
        price_impact¹ = raw_price_paths¹[SimStartTime+l,:] - raw_price_paths¹[SimStartTime,:]
        price_impact² = raw_price_paths²[SimStartTime+l,:] - raw_price_paths²[SimStartTime,:]
        price_impact = cat(price_impact¹,price_impact²,dims=1)
        mean_price_impacts_frac[Gamma,Volume] = mean(price_impact)
        var_price_impacts_frac[Gamma,Volume] = var(price_impact)


        #log_price_impact¹ = log.(sample_price_paths¹[RealStartTime,:]) - log.(sample_price_paths¹[RealStartTime-1,:])
        #log_price_impact² = log.(sample_price_paths²[RealStartTime,:]) - log.(sample_price_paths²[RealStartTime-1,:])
        #log_price_impact  = cat(log_price_impact¹,log_price_impact²,dims=1)

        #mean_log_price_impacts[Volume] = mean(log_price_impact)
        #var_log_price_impacts[Volume] = var(log_price_impact)

        next!(p_outer)
    end
end
# -

mean_price_impacts_frac = .- mean_price_impacts_frac;

# +
x = volume_indices
colors = ["blue","red","green","orange","purple"]
plot()
for Gamma in 1:gi_len
    a,b = log_fit(x,mean_price_impacts_frac[Gamma,:])
    c,d = power_fit(x,mean_price_impacts_frac[Gamma,:])

    #scatter(x,mean_price_impacts,label="data: p(t+1)-p(t)",ms=3, ma=1,yerr=var_price_impacts)
    scatter!(x,mean_price_impacts_frac[Gamma,:],label=string("data: gamma = ",gamma_indices[Gamma]),ms=3, ma=1,color=colors[Gamma])
    plot!(x,mean_price_impacts_frac[Gamma,:],ribbon=var_price_impacts_frac[Gamma,:],alpha=0,fillalpha=0.4,fillcolor=colors[Gamma],label="")
    plot!(x,a.+b.*log.(x),label="Log fit",w=1.5,color=colors[Gamma])
    plot!(x,c.*((x).^d),label="Power fit",w=1.5,color=colors[Gamma],linestyle=:dash)
    plot!(xlabel="Volume",ylabel="Price impact i.e. p(t+1)-p(t)")
    #png("/home/derickdiana/Desktop/Masters/OfficialOne/DraftOnePics/PriceImpactRandomKicksFractional0.7.png")
end

plot!(size=(1000,1000))
png("/home/derickdiana/Desktop/Masters/OfficialOne/DraftOnePics/PriceImpactFractionalAdaptiveDelayAndDuration.png")
plot!(size=(1000,1000))

# +
# Configuration Arguments
num_paths = 50#30

L = 200     # real system width (e.g. 200 meters)
M = 400     # divided into M pieces , 400

p₀ = 230.0  #this is the mid_price at t=0  238.75 

# Free-Parameters for gaussian version
D = 0.5#0.5/8 # real diffusion constant e.g. D=1 (meters^2 / second), 1
α = 0.0 # legacy, no longer used

ν = 3.0 #removal rate
γ = 1.0 #fraction of derivative (1 is normal diffusion, less than 1 is D^{1-γ} derivative on the RHS)

# Source term:
λ = 1.0 #
μ = 0.1 #
mySourceTerm = SourceTerm(λ, μ, true);

# Coupling term:
a = 0.01  #gap between stocks before at full strength: strong is 0.3
b = 2.0   #weighting of interaction term: strong is 2
c = 2.0   #skew factor: strong is 2

myCouplingTerm = CouplingTerm(μ, a, b, c, false);

# My randomness term
σ = 1.0 #variance in randomness
r = 0.5 #proportion of time in which it jumps left or right
β = 0.0 #probability of being the value of the previous lag or mean reversion strength
lag = 10 #lag
do_random_walk = false #behave like a random walk
myRandomnessTerm = RandomnessTerm(σ,r,β,lag,do_random_walk,false)

Δx = L / M  # real gap between simulation points 
Δt = (r * (Δx^2) / (2.0 * D))^(1/γ)

# RL Stuff:
T = 20
RealKickStartTime = 8 # when, in real time, to kick the system
SimKickStartTime = to_simulation_time(RealKickStartTime,Δt)-2 # convert to simulation time
Position = -1
Volume = 20

myRLPusherPush = RLPushTerm(SimKickStartTime,SimKickStartTime+1,Position,Volume,true)
myRLPusherNoPush = RLPushTerm(SimKickStartTime,SimKickStartTime+1,Position,Volume,false)

lob_model_push = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
        mySourceTerm, myCouplingTerm, myRLPusherPush,myRandomnessTerm,shift=-1,old_way=true);

lob_model_no_push = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
        mySourceTerm, myCouplingTerm, myRLPusherNoPush,myRandomnessTerm,shift=-1,old_way=true);

print((Δt,to_simulation_time(T,Δt),num_paths*to_simulation_time(T,Δt))) #about 2GB RAM per 100K, i.e. can only do about 1.8 million
lob_model_push.SK_DP
# -

print(lob_model_push.shift)
fieldnames(SLOB)

# check something actually happens for one example
if true
    myCouplingTerm = CouplingTerm(μ, a, b, c, false);
    
    Volume = 10

    myRLPusherPush = RLPushTerm(SimKickStartTime,SimKickStartTime+1,Position,Volume,true)
    myRLPusherNoPush = RLPushTerm(SimKickStartTime,SimKickStartTime+1,Position,Volume,false)
    
    lob_model_push = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
        mySourceTerm, myCouplingTerm, myRLPusherPush,myRandomnessTerm,shift=0,old_way=false);
    
    lob_model_no_push = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
        mySourceTerm, myCouplingTerm, myRLPusherNoPush,myRandomnessTerm,shift=0,old_way=false);
    
    try clear_double_dict(Dat) catch e print("Not initialized") end
    GC.gc()

    Dat = InteractOrderBooks([lob_model_push,lob_model_no_push], -1, true) ;

    how_many = 12
    p_arr1 = Array{Plots.Plot{Plots.GRBackend},1}(undef,how_many)
    for i in [1:how_many;]
        p_arr1[i] = plot_density_visual(SimKickStartTime-2+i, to_real_time(SimKickStartTime-2+i,Δt), 1, Dat,false, true, 10, Dat[1][1].raw_price_paths[SimKickStartTime])
    end
    plot(p_arr1...,size=(1200,1000))
    
    #png("/home/derickdiana/Desktop/Masters/OfficialOne/DraftOnePics/SinglePointDiffusionStepByStep.png")
    plot!()
end


# +
sums = get_sums(Dat[1][1].lob_densities)

scatter(sums,markersize=2.5,label="Area")
vline!([SimKickStartTime],label="Kick time")
hline!([-Δt*Volume],label="Theoretical kick volume")

# +
volumes = [1:70;] 
volumes = volumes.*0.03

vi_len = length(volumes)
mean_price_impacts = ones(Float64,vi_len)
var_price_impacts = ones(Float64,vi_len)

p_outer = Progress(vi_len,dt=0.1)

l = 3#Int(round(to_simulation_time(1,Δt),digits=0))

price_impact = zeros(Float64,num_paths)

for Volume in 1:vi_len
    myCouplingTerm = CouplingTerm(μ, a, b, c, false);
    
    myRLPusherPush = RLPushTerm(SimKickStartTime,SimKickStartTime+l,Position,volumes[Volume],true)
    myRLPusherNoPush = RLPushTerm(SimKickStartTime,SimKickStartTime+l,Position,volumes[Volume],false)
    
    lob_model_push = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
        mySourceTerm, myCouplingTerm, myRLPusherPush,myRandomnessTerm);
    
    lob_model_no_push = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
        mySourceTerm, myCouplingTerm, myRLPusherNoPush,myRandomnessTerm);

    try clear_double_dict(Dat) catch e print("Not initialized") end
    GC.gc()

    Dat = InteractOrderBooks([lob_model_push,lob_model_no_push], -1, false) ;
    
    price_impact = zeros(Float64,num_paths)
    for path_num in 1:num_paths
        price_impact[path_num] = Dat[path_num][1].raw_price_paths[SimKickStartTime+l] - Dat[path_num][1].raw_price_paths[SimKickStartTime]
    end
    
    mean_price_impacts[Volume] = mean(price_impact)
    var_price_impacts[Volume] = var(price_impact)
    
    next!(p_outer)
end

mean_price_impacts = .- mean_price_impacts;

# +
sub = 1:vi_len

x = volumes[sub]
y = mean_price_impacts[sub]
a,b = log_fit(x,y)  #y[i] = a + b * log(x[i])

scatter(x,y,label="data: p(t+1)-p(t)",ms=3, ma=1)
plot!(x,a.+b.*log.(x),label="Log fit",w=1.5,color="red")
plot!(x,y,ribbon=var_price_impacts[sub].^0.5,alpha=0,fillalpha=0.4,fillcolor="blue",label="")

#png("/home/derickdiana/Desktop/Masters/OfficialOne/DraftOnePics/PriceImpactRandomKicksNoFractional.png")

plot!(xlabel="Volume",ylabel="Price impact i.e. p(t+1)-p(t)")
plot!()
