using Plots, CurveFit, Statistics, ProgressMeter, InteractingLOBs

# +
function get_sums(lob_densities; absolute=true)
    l = size(lob_densities)[2]
    sums = zeros(Float64, l)

    if absolute
        my_sum = (x) -> sum(abs.(x))
    else
        my_sum = (x) -> sum(x)
    end

    for t in 1:l
        sums[t] = my_sum(lob_densities[:, t])
    end

    return sums
end

function save_png(folder_name, name)
    png(string("/Users/dominic/Desktop/personal-projects/InteractingLOBs.jl/", folder_name, "/", name, ".png"))
end

function plot_sums(Dat; path_num=1, slob_num=1, new_plot=true, absolute=true)
    if new_plot
        plot()
    end
    sums = get_sums(Dat[path_num][slob_num].lob_densities; absolute=absolute)
    slob = Dat[path_num][slob_num].slob

    scatter!(sums, markersize=1.4, label="Area", size=(1400, 500))
    vline!([slob.rl_push_term.StartTime], label="Kick time")
    if absolute
        hline!([sums[slob.rl_push_term.StartTime-1] + slob.Δt * slob.rl_push_term.Amount], label="Theoretical kick volume")
    else
        hline!([-slob.Δt * slob.rl_push_term.Amount], label="Theoretical kick volume")
    end

    source_sums = get_sums(Dat[path_num][slob_num].sources; absolute=absolute)
    hline!([source_sums[3] / slob.nu], label="Theoretical equilibrium volume")

    plot!(xlabel="Simulation steps", ylabel="Signed area under system")
    plot!(size=(1200, 500))
end

path_to_plot = 1
l = @layout [a d; b e; c f];

function plot_price_path(s, r, lob_num, Dat, diff=false)
    lob_model = Dat[1][lob_num].slob

    plt = plot()
    for path in 1:lob_model.num_paths
        raw_price_paths = Dat[path][lob_num].raw_price_paths[1:s]
        if diff
            raw_price_paths .-= Dat[path][3-lob_num].raw_price_paths[1:s]
        end
        plot!((0:s-1) .* lob_model.Δt, raw_price_paths, color=5, w=0.6)

    end

    obs_price_paths = Dat[path_to_plot][lob_num].obs_price_paths[1:r]
    if diff
        obs_price_paths .-= Dat[path_to_plot][3-lob_num].obs_price_paths[1:r]
    end
    plot!(0:r-1, obs_price_paths, color=1, w=2.7)

    plot!(legend=false, ylab="Price", xlab="Time")
    return plt
end

function my_pad(str::String, target_length::Int; do_left=true::Bool)
    l = length(str)

    if l < target_length
        diff = target_length - l
        add_space = repeat(" ", diff)

        if do_left
            return add_space * str
        else
            return str * add_space
        end
    end

    return str
end

function plot_density_visual(Dat, s, lob_num;
    dosum=false, plot_raw_price=true, x_axis_width=-1, center=-2, shift_with_loop=false, marker_size=1.6, overall_alpha=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    path_to_plot=1, size=(1000, 1000), kw_for_plot=(), shift_back_approx=true, do_left=false, do_right=false, label_price=true, label_Δxₘ=true,
    annotate_pos=:topright)


    dat = Dat[path_to_plot][lob_num]
    lob_model = dat.slob

    if !lob_model.store_past_densities
        @warn "Should not be using this function since you did not store the densities\n"
    end

    mult = lob_model.Δt

    shift = dat.x_shifts[s]

    if center == -1
        center = dat.raw_price_paths[s]
    end

    if center == -2
        center = lob_model.p₀
    end

    if shift_with_loop
        camera_shift = shift
    else
        camera_shift = 0
    end

    if x_axis_width == -1
        x_axis_width = dat.slob.L / 15
    end


    removal_fraction = exp(-lob_model.nu * lob_model.Δt) - 1



    x_axis = [center - x_axis_width + camera_shift, center + x_axis_width + camera_shift]

    lob_densities = dat.lob_densities[:, s]
    lob_densities_L = dat.lob_densities_L[:, s]
    lob_densities_R = dat.lob_densities_R[:, s]
    couplings = dat.couplings[:, s+1]#-1 not there before
    sources = dat.sources[:, s+1]#-1 not there before
    rl_pushes = dat.rl_pushes[:, s+1]#-1 not there before
    removals = dat.lob_densities[:, s]

    if dosum
        dat2 = Dat[path_to_plot][3-lob_num]
        lob_densities .+= dat2.lob_densities[:, s]
        lob_densities_L .+= dat2.lob_densities_L[:, s]
        lob_densities_R .+= dat2.lob_densities_R[:, s]
        couplings .+= dat2.couplings[:, s+1]#-1 not there before
        sources .+= dat2.sources[:, s+1]#-1 not there before
        rl_pushes .+= dat2.rl_pushes[:, s+1]#-1 not there before
        removals .+= dat2.lob_densities[:, s]
    end

    lob_densities = lob_densities
    lob_densities_L = lob_densities_L
    lob_densities_R = lob_densities_R
    couplings = mult .* couplings
    sources = mult .* sources
    rl_pushes = mult .* rl_pushes
    removals = removal_fraction .* removals

    x_range = lob_model.x .+ shift
    common = (label="", markersize=marker_size, markerstrokewidth=0.0)


    Δx_ = dat.slob.Δxs[s-1]

    if shift_back_approx
        x_range_shifted_left = x_range .- Δx_
        x_range_shifted_right = x_range .+ Δx_
    else
        x_range_shifted_left = x_range
        x_range_shifted_right = x_range
    end

    plt =
        scatter(x_range, lob_densities, color=1, alpha=overall_alpha[1]; common...)
    if do_left
        scatter!(x_range_shifted_left, lob_densities_L, color=1, alpha=overall_alpha[2]; common...)
    end
    if do_right
        scatter!(x_range_shifted_right, lob_densities_R, color=1, alpha=overall_alpha[3]; common...)
    end
    # scatter!(x_range, couplings, color=2, alpha=overall_alpha[4]; common...)
    scatter!(x_range, sources, color=3, alpha=overall_alpha[5]; common...)
    scatter!(x_range, removals, color=5, alpha=overall_alpha[6]; common...)
    scatter!(x_range, rl_pushes, color=4, alpha=overall_alpha[7]; common...)


    x_range_dense = lob_model.x .+ shift

    if shift_back_approx
        x_range_shifted_left_dense = x_range_dense .- Δx_
        x_range_shifted_right_dense = x_range_dense .+ Δx_
    else
        x_range_shifted_left_dense = x_range_dense
        x_range_shifted_right_dense = x_range_dense
    end


    plot!(x_range_dense, lob_densities, label="φⁱ", color=1, alpha=overall_alpha[1], dpi=300)
    if do_left
        plot!(x_range_shifted_left_dense, lob_densities_L, label="φⁱ⁻¹", color=1, alpha=overall_alpha[2], style=:dash, dpi=300)
    end
    if do_right
        plot!(x_range_shifted_right_dense, lob_densities_R, label="φⁱ⁺¹", color=1, alpha=overall_alpha[3], style=:dashdotdot, dpi=300)
    end
    #plot!(x_range_dense,           couplings, color=2, label="Coupling",alpha = overall_alpha[4]) ;
    plot!(x_range_dense, sources, label="Arrivals", color=3, alpha=overall_alpha[5], dpi=300)
    plot!(x_range_dense, removals, label="Removals", color=5, alpha=overall_alpha[6], dpi=300)
    plot!(x_range_dense, rl_pushes, label="Impulse", color=4, alpha=overall_alpha[7], dpi=300)


    price_pos = dat.raw_price_paths[s]
    if (plot_raw_price)
        if (isinteger(s * lob_model.Δt))
            mycol = :black
        else
            mycol = 1
        end


        scatter!([price_pos], [0], label="p";
            markersize=3, markercolor=mycol, markerstrokewidth=0.5)
        if do_left
            scatter!([price_pos - Δx_], [0], label="p-Δxₘ"; markershape=:star4,
                markersize=5, markercolor="black", markerstrokewidth=0.5)
        end
        if do_right
            scatter!([price_pos + Δx_], [0], label="p+Δxₘ"; markershape=:star4,
                markersize=5, markercolor="black", markerstrokewidth=0.5)
        end
    end

    if lob_model.do_exp_dist_times
        real_time = round(lob_model.Δts_cum[s], digits=2)#round(lob_model.Δt * (s-1),digits=3)
    else
        real_time = round(lob_model.Δt * (s - 1), digits=3)
    end

    top_right_label = ""
    if label_price
        top_right_label *= my_pad("  p   = ", 8; do_left=false) * my_pad(string(round(price_pos, digits=2)), 8; do_left=false)
    end
    top_right_label *= "\n"
    if label_Δxₘ
        top_right_label *= my_pad("Δxₘ = ", 8; do_left=false) * my_pad(string(round(Δx_, digits=2)), 8; do_left=false)
    end
    top_right_label *= "\n"

    annotate!((annotate_pos, text(top_right_label, 8)))

    plot!(legend=:bottomleft, title="", xlab="Price at time=$real_time", ylab="Densities", xlim=x_axis, dpi=300)#, legendfontsize=17.0)
    plot!(; dpi=300, kw_for_plot...)
    return plt
end;



function get_second_derivative(x, temp)
    middle = 2:length(temp)-1
    return ((temp[middle.+1] - temp[middle]) ./ (x[middle.+1] - x[middle]) - (temp[middle] - temp[middle.-1]) ./ (x[middle] - x[middle.-1])) ./ (x[middle.+1] - x[middle.-1])
end;


function fit_and_plot_price_impact(volumes, mean_price_impacts, var_price_impacts, labels;
    sub=-1, do_kinks=true, colors=-1, do_power_fit=false, xticks=-1, new_plot=false, forplot=(),
    do_ribbon=true, do_horiz_log_shift=false, do_log_fit=true, do_log_plot=false, do_vert_log_shift=true)
    if new_plot
        plot()
    end

    volumes_original = volumes[:]

    if do_log_plot
        volumes = log.(volumes .+ 1)[:]
    end

    #labels create indices which are usually the different gammas

    li_len = length(labels)
    vi_len = length(volumes)

    if sub == -1
        sub = 1:vi_len
    end

    if colors == -1 #change later!
        if li_len <= 5
            colors = ["blue", "red", "green", "orange", "purple"]
        else
            colors = 1:li_len
        end
    end

    shift = do_horiz_log_shift ? 1 : 0

    for ind in 1:li_len
        if do_log_fit
            a, b = log_fit(volumes_original[sub] .+ shift, mean_price_impacts[sub, ind])
            a = do_vert_log_shift ? a : 0
            y_line_log_fit = a .+ b .* log.(volumes_original[sub] .+ shift)

            plot!(volumes[sub], y_line_log_fit,
                label=string("Log fit: ", round(a, digits=2), " + ", round(b, digits=2), "log(x+1)"), w=1.5, color=colors[ind])
        end

        scatter!(volumes[sub], mean_price_impacts[sub, ind],
            label=labels[ind], ms=1.5, markerstrokewidth=0.1, ma=1, color=colors[ind])

        if do_ribbon
            plot!(volumes[sub], mean_price_impacts[sub, ind], ribbon=var_price_impacts[sub, ind] .^ 0.5, alpha=0,
                fillalpha=0.4, fillcolor=colors[ind], label="")
        end

        if do_power_fit
            c, d = power_fit(volumes[sub], mean_price_impacts[sub, ind])

            plot!(volumes[sub], c .* ((volumes[sub]) .^ d), label="Power fit",
                w=1.5, color=colors[ind], linestyle=:dash)
        end


        if xticks != -1
            plot!(xticks=xticks)
        end


        if do_kinks
            vol_scale = (volumes[end] - volumes[1]) / 20
            impact_scale = (maximum(mean_price_impacts[sub, ind], dims=1)[1] - minimum(mean_price_impacts[sub, ind], dims=1)[1]) / 30

            second_deriv = get_second_derivative(volumes_original[sub], mean_price_impacts[sub, ind])
            kink_position = findnext(x -> x > 0, second_deriv, 1)

            while !(kink_position === nothing)
                target_x, target_y = (volumes[sub][kink_position+1], mean_price_impacts[sub[kink_position+1], ind])
                quiver!([target_x + vol_scale], [target_y - impact_scale], quiver=([-vol_scale], [impact_scale]), color=colors[ind])
                scatter!([target_x], [target_y], markershape=:star5, color=colors[ind],
                    label=string("Kink at position ", round(target_x, digits=2)), markerstrokewidth=0.1, ma=1)

                kink_position = findnext(x -> x > 0, second_deriv, kink_position + 2)
            end
        end

    end

    my_xlabel = do_log_plot ? "log(Volume)" : Volume
    plot!(xlabel=my_xlabel, ylabel="Price impact i.e. p(t+1)-p(t)"; forplot...)
    plot!(size=(1000, 1000))
end

# Shows n pictures through time
function quick_plot(lob_models, sim_start_time, how_many=12; center_at_midprice=true,
    dosum=false, plot_raw_price=true, x_axis_width=-1, center=-2, scale_down=true, marker_size=1.6, overall_alpha=[1.0, 1.0, 1.0, 1.0, 1.0])
    try
        clear_double_dict(Dat)
    catch e
        print("Not initialized")
    end
    GC.gc()

    Dat = InteractOrderBooks(lob_models, -1, true)

    if x_axis_width == -1
        x_axis_width = 10#lob_models[1].L/8
    end

    p_arr1 = Array{Plots.Plot{Plots.GRBackend},1}(undef, how_many)
    for i in [1:how_many;]
        p_arr1[i] = plot_density_visual(sim_start_time - 2 + i, to_real_time(sim_start_time - 2 + i, Δt), 1, Dat;
            dosum=dosum, plot_raw_price=plot_raw_price, x_axis_width=x_axis_width, center=center, scale_down=scale_down, marker_size=marker_size, overall_alpha=overall_alpha)#,overall_alpha=[0.4,1.0,0.4,0.4,0.4,0.4])
    end
    p1 = plot(p_arr1..., size=(1200, 1000), dpi=100)
    savefig(p1, "Plots/Epps/ThroughTime.png")

    return Dat
end

#NB
function calculate_price_impacts(volumes, inputs, get_set; measured_slob=1) #lob_sets must have the form ((lob_a1,lob_b1),(lob_a2,lob_b2)) where the inner brackets are meant
    #to be passed to InteractOrderBooks together

    vol_len = length(volumes)
    input_len = length(inputs) #encoding of all paramters we want to change

    mean_price_impacts = ones(Float64, vol_len, input_len)
    var_price_impacts = ones(Float64, vol_len, input_len)

    for Ind in 1:input_len
        p_outer = Progress(vol_len, dt=0.1)

        #for Volume in 1:vol_len
        Threads.@threads for Volume in 1:vol_len
            lob_models, sim_start_time, l = get_set(volumes[Volume], inputs[Ind])

            #try clear_double_dict(Dat) catch e print("Not initialized") end
            #GC.gc()

            Dat = InteractOrderBooks(lob_models, -1, false)

            num_paths = lob_models[1].num_paths
            price_impact = zeros(Float64, num_paths)
            for path_num in 1:num_paths
                price_impact[path_num] = Dat[path_num][measured_slob].raw_price_paths[sim_start_time+l] - Dat[path_num][measured_slob].raw_price_paths[sim_start_time]
            end

            mean_price_impacts[Volume, Ind] = mean(price_impact)
            var_price_impacts[Volume, Ind] = var(price_impact)

            next!(p_outer)

        end
    end

    mean_price_impacts = .-mean_price_impacts

    return (mean_price_impacts, var_price_impacts)
end

function myvar(p, x, Δx)
    mu = mymean(p, x, Δx)
    Ex2 = 0
    for i in 1:length(p)
        Ex2 += p[i] * x[i]^2
    end
    return Ex2 - mu^2
end

function mymean(p, x, Δx)
    sum = 0
    for i in 1:length(p)
        sum += p[i] * x[i]
    end
    return sum
end
