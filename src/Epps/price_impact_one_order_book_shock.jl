myscale = 1
include("../setup.jl")

volumes = 1:100

function get_set_inner(volume,γ,ν)
    myRandomnessTerm = RandomnessTerm(σ,r,β,lag,do_random_walk,true)
    myCouplingTerm = CouplingTerm(μ, a, b, c, true);

    Δt = (r * (Δx^2) / (2.0 * D * myscale))^(1/γ)

    l = Int(round(to_simulation_time(1,Δt)/3,digits=0))

    RealKickStartTime = 8 # when, in real time, to kick the system
    SimKickStartTime = to_simulation_time(RealKickStartTime,Δt)-2 # convert to simulation time

    myRLPusherPush = RLPushTerm(SimKickStartTime,SimKickStartTime+1,Position,volume,true)
    myRLPusherNoPush = RLPushTerm(SimKickStartTime,SimKickStartTime+1,Position,volume,false)

    lob_model_push = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
        mySourceTerm, myCouplingTerm, myRLPusherPush,myRandomnessTerm, michaels_way = false);

    lob_model_no_push = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,
        mySourceTerm, myCouplingTerm, myRLPusherNoPush,myRandomnessTerm, michaels_way = false);

    return ([lob_model_push, lob_model_no_push], SimKickStartTime, l)
end

function get_set(volume,γ)
    return get_set_inner(volume,γ,ν)
end

function fit_and_plot_price_impact(; volumes,
    mean_price_impacts,
    var_price_impacts,
    labels,
    sub=-1,
    do_kinks=false,
    colors=-1,
    do_power_fit=false,
    xticks=-1,
    new_plot=true,
    forplot=(),
    do_ribbon=true,
    do_horiz_log_shift=false,
    do_log_fit=false,
    do_log_plot=false,
    do_vert_log_shift=true
)
    if new_plot
        p1 = plot()
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

            plot!(volumes[sub], y_line_log_fit ,label=string("Log fit: ",round(a,digits=2)," + ",round(b,digits=2),"log(x+1)"),w=1.5,color=colors[ind])
        end

        scatter!(volumes[sub], mean_price_impacts[sub,ind], label=labels[ind], ms=1.5, markerstrokewidth=0.1, ma=1, color=colors[ind])

        if do_ribbon
            plot!(volumes[sub], mean_price_impacts[sub,ind], ribbon=var_price_impacts[sub,ind].^0.5, alpha=0, fillalpha=0.4,fillcolor=colors[ind], label="")
        end

        if do_power_fit
            c,d = power_fit(volumes[sub], mean_price_impacts[sub,ind])
            plot!(volumes[sub], c.*((volumes[sub]).^d),label="Power fit", w=1.5, color=colors[ind], linestyle=:dash)
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
                quiver!([target_x+vol_scale], [target_y-impact_scale], quiver=([-vol_scale],[impact_scale]), color=colors[ind])
                scatter!([target_x], [target_y], markershape=:star5,color=colors[ind], label=string("Kink at position ",round(target_x,digits=2)), markerstrokewidth=0.1,  ma=1)

                kink_position = findnext(x->x>0,second_deriv,kink_position+2)
            end
        end

    end

    my_xlabel = do_log_plot ? "log(Volume)" : "Volume"
    plot!(xlabel=my_xlabel,ylabel="Price impact i.e. p(t+1)-p(t)";forplot...)
    plot!(size=(1000,1000), dpi=300)

    savefig(p1, "Plots/Epps/PriceImpactOneShock.png")
end


(mean_price_impacts_frac_no_random_new_way,var_price_impacts_frac_no_random_new_way) = calculate_price_impacts(volumes,[1.0],get_set)

fit_and_plot_price_impact(;
    volumes = volumes,
    mean_price_impacts = mean_price_impacts_frac_no_random_new_way,
    var_price_impacts = var_price_impacts_frac_no_random_new_way,
    labels = ["Price Impact"],
    colors=["red"]
)
