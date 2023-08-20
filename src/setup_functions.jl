using Plots, CurveFit, Statistics, ProgressMeter, InteractingLOBs

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
    png(string("/Users/dominic/Desktop/personal-projects/InteractingLOBs.jl/",folder_name,"/",name,".png"))
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

    plt = plot(lob_model.x,lob_densities, color=1,label="Density",alpha = overall_alpha[1]);
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

# Shows n pictures through time
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
    p1 = plot(p_arr1...,size=(1200,1000), dpi=100)
    savefig(p1, "Plots/Epps/ThroughTime.png")

    return Dat
end

#NB
function calculate_price_impacts(volumes,inputs,get_set; measured_slob=1) #lob_sets must have the form ((lob_a1,lob_b1),(lob_a2,lob_b2)) where the inner brackets are meant
                                                                                #to be passed to InteractOrderBooks together

    vol_len = length(volumes)
    input_len = length(inputs) #encoding of all paramters we want to change

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
