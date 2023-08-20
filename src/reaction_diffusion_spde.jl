# -*- coding: utf-8 -*-
# +
function initial_conditions_numerical(D, t, slob_num, V₀=0)
    slob¹ = D[slob_num].slob
    #slob² = D[2+(1-slob_num)].slob

    if (t == 1)
        t = t + 1
    end#

    # eq 36:
    ud = (-V₀/(2.0*slob¹.Δx) + slob¹.D/(slob¹.Δx^2)) * ones(Float64, slob¹.M)#upper diagonal
    md = ((-2.0*slob¹.D)/(slob¹.Δx^2) - slob¹.nu) * ones(Float64, slob¹.M+1) #middle diagonal
    ld = (V₀/(2.0*slob¹.Δx) + slob¹.D/(slob¹.Δx^2)) * ones(Float64, slob¹.M) #lower diagonal
    A = Tridiagonal(ld, md, ud)

    A[1,2] = 2*slob¹.D/(slob¹.Δx^2)
    A[end, end-1] = 2*slob¹.D/(slob¹.Δx^2)

    source   = slob¹.source_term(   D,
                                    slob_num, t)

    coupling = slob¹.coupling_term( D,
                                    slob_num, t)

    rl_push  = slob¹.rl_push_term(  D,
                                    slob_num, t)

    B = .-(source.+coupling.+rl_push)#B is s in (37)

    #temp
    if slob¹.source_term.do_source
        φ = A \ B
    else
        φ = source .* 0
    end

    D[slob_num].lob_densities[:,t] = φ
    D[slob_num].sources[:,t]       = source
    D[slob_num].couplings[:,t]     = coupling
    D[slob_num].rl_pushes[:,t]     = rl_push
    D[slob_num].V[t]               = V₀

    #give x such that A*x = B
    #return φ, source, coupling, rl_push, V₀

end
# -

function extract_mid_price_index(slob,lob_density)
    mid_price_ind = 2
    l = length(lob_density)
    #while (mid_price_ind<l)&&((lob_density[mid_price_ind] > 0) || (lob_density[mid_price_ind+1]>lob_density[mid_price_ind]))
    while ((mid_price_ind<l)&&(lob_density[mid_price_ind] > 0))
        mid_price_ind += 1
    end #scan array in x looking for cross over point of the mid price
    if ( mid_price_ind==l || mid_price_ind==2 )
        return -1
    end
    return mid_price_ind
end

function extract_mid_price(slob, lob_density)
    mid_price_ind = extract_mid_price_index(slob,lob_density)
    if mid_price_ind == -1
        return -1
    end

    y1 = lob_density[mid_price_ind-1]#y value (density) to the left
    y2 = lob_density[mid_price_ind]#y value (density) to the right
    x1 = slob.x[mid_price_ind-1]#x value value to the left

    mid_price = (-y1 * slob.Δx)/(y2 - y1) + x1 #solution to assuming straight approximation between
                                               #left and right point (did the math) (see page 20)
    return mid_price
end


# +
using Interpolations
using Plots

x = 1.0:1.0:10.0
#y = @. cos(x^2 / 9.0)
y = [1.0 for xi in x]
y[3] = 10.0

A = hcat(x,y)

# +
itps = Array{AbstractInterpolation}(undef,0)
xs = Array{Vector{Float64}}(undef,0)
ys = Array{Vector{Float64}}(undef,0)
labels = Array{String}(undef,0)

tfine = 1.0:0.1:10.0

push!(itps,Interpolations.scale(interpolate(A, (BSpline(Linear()))), x, 1:2))
push!(labels,"BSplit Linear Natural On Grid No Interp")

push!(itps,Interpolations.scale(interpolate(A, (BSpline(Quadratic(Natural(OnGrid()))))), x, 1:2))
push!(labels,"BSplit Quadratic Natural On Grid")

push!(itps,Interpolations.scale(interpolate(A, (BSpline(Cubic(Natural(OnGrid()))))), x, 1:2))
push!(labels,"BSplit Cubic Natural On Grid")



for i in 1:length(itps)
    push!(xs,[itps[i](t,1) for t in tfine])
    push!(ys,[itps[i](t,2) for t in tfine])
end
# -

scatter(x, y, label="knots")
for i in 1:length(itps)
    plot!(xs[i], ys[i], label=labels[i],ls=:dash)
end
plot!(size=(1100,1100))

function SibuyaKernelDP(n,slob)
    return slob.SK_DP[n]
end

function calculate_next_step_original(φ, slob, P, P⁺, P⁻,net_source) #NB, pass in only most recent time
    ###### LITERALLY Michael's code #####################################

    φ₋₁ = φ[1]
    φ₊₁ = φ[end]
    φ_next = zeros(Float64, size(φ,1))

    φ_next[1] = P⁺ * φ₋₁ + P⁻ * φ[2] + P * φ[1] -
        slob.nu * slob.Δt * φ[1] + slob.Δt * net_source[1]

    φ_next[end] = P⁻ * φ₊₁ + P⁺ * φ[end-1] + P * φ[end] -
        slob.nu * slob.Δt * φ[end] + slob.Δt * net_source[end]

    φ_next[2:end-1] = P⁺ * φ[1:end-2] + P⁻ * φ[3:end] + P * φ[2:end-1] -
        slob.nu * slob.Δt * φ[2:end-1] + slob.Δt * net_source[2:end-1]

    #####################################################################

    return φ_next
end

function calculate_next_step(φ, slob, P, P⁺, P⁻,net_source,
                                t)
    myend = size(φ,1)
    φ_next = zeros(Float64,myend)
    middle = 2:(myend-1)    #all indices excluding the first and last

    if (slob.old_way)&&(slob.shift==-1)&&(slob.cut_off==1)
        Pm = P
        removal_scale = - slob.nu * slob.Δt
    else
        Pm = P - 1
        removal_scale = exp( - slob.nu * slob.Δt)
    end

    for tₘ in max(1,t-slob.cut_off):(t-1)

        front = SibuyaKernelDP(t-tₘ,slob)

        φ_next[1] += front *
            (  P⁺ * φ[1        ,tₘ] +   #1      -1    (but not because boundary)
               P⁻ * φ[2        ,tₘ] +   #1      +1
               Pm * φ[1        ,tₘ]  )  #1      +0

        φ_next[end] += front *
            (  P⁺ * φ[end-1    ,tₘ] +   #end    -1
               P⁻ * φ[end      ,tₘ] +   #end    +1    (but not because boundary)
               Pm * φ[end      ,tₘ]  )  #end    +0

        φ_next[middle] +=  front *
            (  P⁺ * φ[middle.-1,tₘ] +   #middle -1
               P⁻ * φ[middle.+1,tₘ] +   #middle +1
               Pm * φ[middle   ,tₘ]  )  #middle +0

    end

    φ_next += removal_scale * φ[:,t-1] + slob.Δt * net_source

    return φ_next
end

# +
function intra_time_period_simulate_fractional( D,
                                                t, slob_num)
    slob = D[slob_num].slob

    P, P⁺, P⁻, V_t  = slob.randomness_term(slob, D[slob_num].V, t)

    source   = slob.source_term(    D,
                                    slob_num, t)

    coupling = slob.coupling_term(  D,
                                    slob_num, t)

    rl_push  = slob.rl_push_term(   D,
                                    slob_num, t)

    net_source = source .+ coupling .+ rl_push

    if slob.michaels_way
        ###### LITERALLY Michael's code #####################################
        φ_next = calculate_next_step_original(D[slob_num].lob_densities[:,t-1], slob, P, P⁺, P⁻, net_source)

    else
        # the newest version with fractional derivatives
        φ_next = calculate_next_step(D[slob_num].lob_densities, slob, P, P⁺, P⁻, net_source, t)

    end

    # store all those intermediate values
    D[slob_num].lob_densities[:,t] = φ_next
    D[slob_num].sources[:,t]       = source
    D[slob_num].couplings[:,t]     = coupling
    D[slob_num].rl_pushes[:,t]     = rl_push
    D[slob_num].Ps[t]              = P
    D[slob_num].P⁺s[t]             = P⁺
    D[slob_num].P⁻s[t]             = P⁻
    D[slob_num].V[t]               = V_t

end
# -


function move_back_densities(D)
    for slob_num in 1:length(D)
        D[slob_num].lob_densities[:,1] = D[slob_num].lob_densities[:,2]
        D[slob_num].sources[:,1]       = D[slob_num].sources[:,2]
        D[slob_num].couplings[:,1]     = D[slob_num].couplings[:,2]
        D[slob_num].rl_pushes[:,1]     = D[slob_num].rl_pushes[:,2]
        D[slob_num].Ps[1]              = D[slob_num].Ps[2]
        D[slob_num].P⁺s[1]             = D[slob_num].P⁺s[2]
        D[slob_num].P⁻s[1]             = D[slob_num].P⁻s[2]
        D[slob_num].V[1]               = D[slob_num].V[2]
    end
end

# +
function dtrw_solver_fractional(D) #handles change over time logic
    broke_point = -1
    slob = D[1].slob

    time_steps = to_simulation_time(slob.T, slob.Δt)


    t = 1 #initial conditions take current t to read values

    for slob_num in 1:length(D)
         D[slob_num].raw_price_paths[1] = D[slob_num].slob.p₀
         D[slob_num].obs_price_paths[1] = D[slob_num].slob.p₀
         initial_conditions_numerical( D, t, slob_num)
    end


    t = 2

    not_broke = reduce(&,[(D[l].raw_price_paths[t-1]!=-1 || !D[l].slob.source_term.do_source) for l in 1:length(D)],init=true) #only true when non of the most recent raw price values are -1


    while (t <= time_steps) && (not_broke)

        if slob.store_past_densities
            t_store = t
        else
            move_back_densities(D)
            t_store = 2
        end

        for slob_num in 1:length(D)
            intra_time_period_simulate_fractional( D, t_store, slob_num)

            D[slob_num].raw_price_paths[t] =
                                    extract_mid_price(D[slob_num].slob, D[slob_num].lob_densities[:,t])
        end


        t += 1
        not_broke = reduce(&,[(D[l].raw_price_paths[t-1]!=-1 || !D[l].slob.source_term.do_source) for l in 1:length(D)],init=true)
    end

    #check if broken
    if !(not_broke)
        print("BE")
        broke_point = t-1
    end

    for slob_num in 1:length(D)
        D[slob_num].obs_price_paths =
                    sample_mid_price_path(D[slob_num].slob, D[slob_num].raw_price_paths)
    end

end
# -
