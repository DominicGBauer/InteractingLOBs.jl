# -*- coding: utf-8 -*-
mutable struct SLOB
    num_paths::Int64 #How many paths to simulate
    T::Int64    # How many time periods it runs for
    p₀::Float64 # The initial mid price
    M::Int64    # Number of discrete price points in addition to p₀ (the central point)
    L::Float64  # The length of the price dimension in the lattice
    D::Float64  # The diffusion constant (determines Δt)
    nu::Float64 # Cancellation rate
    α::Float64  # Waiting time in terms of Exp(α)
    source_term::SourceTerm # Source term s(x,t)(μ,λ,p¹,p²) 
    coupling_term::CouplingTerm # Interaction term c(x,t)(μ,λ,p¹,p²)
    rl_push_term::RLPushTerm
    randomness_term::RandomnessTerm
    x::Array{Float64, 1}    #
    Δx::Float64 # The grid spacing which is calculated as Δx=L/M
    Δt::Float64 # The temporal grid spacing calculated as Δt=Δx*Δx/(2D)
    γ::Float64 #Degree of fractional diffusion
    SK_DP::Array{Float64,1}
    cut_off::Int64
end


function SLOB(num_paths::Int64, T::Int64, p₀::Float64,
    M::Int64, L::Real, D::Float64, nu::Float64,
    α::Float64, γ::Float64, 
    source_term::SourceTerm, coupling_term::CouplingTerm, rl_push_term::RLPushTerm, randomness_term::RandomnessTerm) 
    #translates convinient order into order expected by object itself
    
    x₀ = p₀ - 0.5*L
    x_m = p₀ + 0.5*L
    #@assert x₀ >= 0
    x = collect(Float64, range(x₀, stop=x_m, length=M+1)) #creates an array of the entries. So returns set of x points
    r = randomness_term.r
    Δx = L/M
    Δt = (r * (Δx^2) / (2.0 * D))^(1/γ) #(Δx^2) / (2.0*D)
    
    function partially_applied_Sibuya(n)
         return SibuyaKernelModified(n,γ,nu,Δt)
    end
    
    l = to_simulation_time(T,Δt)
    SK_DP = zeros(Float64,l)
    
    cut_off = 0

    # basic loop
    next_sk = γ
    SK_DP[cut_off+1] = next_sk*exp(-(cut_off+1)*nu*Δt)                                   
    cut_off += 1                                    
    
    # special case for cut_off = 2
    next_sk = (-1+γ)*(1-(2-γ)/(cut_off+1))
    SK_DP[cut_off+1] = next_sk*exp(-(cut_off+1)*nu*Δt)                                   
    
    # check if latest meets condition to be included. Only then do we increase the cutoff and calculate again
    while (cut_off<(l-1)) && (SK_DP[cut_off+1]<=-0.001)            
        cut_off += 1
        
        next_sk = next_sk * (1-(2-γ)/(cut_off+1)) 
        SK_DP[cut_off+1] = next_sk * exp(-(cut_off+1)*nu*Δt)                    
    end
   
    SK_DP = SK_DP[1:cut_off] #shorten the array to whatever cutoff was found
    
    return SLOB(num_paths, T, p₀, M, L, D, nu, α, source_term, coupling_term, rl_push_term, randomness_term, x, Δx, Δt, γ, SK_DP,cut_off)
end

function to_real_time(simulation_time, Δt::Float64) # from simulation time
    return floor(Int, (floor(Int,simulation_time - 1)+1) * Δt)
end


function to_simulation_time(real_time, Δt::Float64) #from real time 
    # if real_time = end time, this gives the total number of time indices in the simulation
    simulation_time = real_time / Δt
    return floor(Int, simulation_time + eps(simulation_time)) + 1 #how many simulation steps did it take to reach 'real_time'?
end


# +
# the above two are inverses for whatever unknown reason:
#temp = 
#        [ 
#           [to_real_time(to_simulation_time(xi,m),m) - xi for xi in [1:10000;]]
#           #[round(m,digits=2) for xi in [1:10000;]]
#        for m in rand(100)]
#
#sum(sum(temp))
# -

function get_sample_inds_by_sampling_at_integer_real_times(slob::SLOB) #sample the simulation at integer times in real time
    partially_applied_to_sim_time(real_times) = to_simulation_time(real_times, slob.Δt) 
    my_real_times = 0:slob.T
    sample_inds = partially_applied_to_sim_time.(my_real_times)
    return sample_inds
end

function sample_mid_price_path(slob::SLOB, price_path)
    sample_inds = get_sample_inds_by_sampling_at_integer_real_times(slob)
    mid_prices = price_path[sample_inds]
    return mid_prices
end

# +
#########################My own code below here########################
# -

function InteractOrderBooks(slobs::Array{SLOB,1}, seed::Int=-1, progress = false)#same but returns more stuff
    #handles change over paths logic
    
    slob = slobs[1]
    
    if (progress==true)
        p = Progress(slob.num_paths,dt=0.1)
    end
    
    # do some logging IDK
    logger = FileLogger(Dict(Logging.Info => "info.log", Logging.Error => "error.log"), append=false)
    global oldglobal = global_logger(logger)
    
    # generate a set of seeds
    if seed == -1
            seeds = Int.(rand(MersenneTwister(), UInt32, slob.num_paths))
        else
            seeds = Int.(rand(MersenneTwister(seed), UInt32, slob.num_paths))
    end
    
    
    num_time_steps = to_simulation_time(slob.T, slob.Δt) #how many num_time steps in the base simulation? i.e.Raw num_time indices
    
    # the below dimensions are ordered according to the frequency with which they appear
    
    
    outer_dict = Dict{Int64,Dict{Int64,DataPasser}}()
    for path_num in 1:slob.num_paths
        inner_dict = Dict{Int64,DataPasser}()
        for slob_num in 1:length(slobs)
            ### reserve memory
            lob_densities =    zeros(Float64, slob.M+1, num_time_steps + 1)
            sources =          zeros(Float64, slob.M+1, num_time_steps + 1)
            couplings =        zeros(Float64, slob.M+1, num_time_steps + 1)
            rl_pushes =        zeros(Float64, slob.M+1, num_time_steps + 1)

            raw_price_paths =  ones(Float64,           num_time_steps + 1)
            obs_price_paths =  ones(Float64,           slob.T         + 1)

            P⁺s =              ones(Float64,           num_time_steps    )
            P⁻s =              ones(Float64,           num_time_steps    )
            Ps =               ones(Float64,           num_time_steps    )
            V =                ones(Float64,           num_time_steps + 1)

            my_struct = DataPasser(slobs[slob_num], lob_densities, sources, couplings, rl_pushes,
                                    raw_price_paths, obs_price_paths,
                                    P⁺s, P⁻s, Ps, V)
                        #DataPasser(1)
            ### end of reserve memory
            
            inner_dict[slob_num] = my_struct
        end
        
        outer_dict[path_num] = inner_dict
        
    end
    
    #print(outer_dict)
    
    
    broke_points =     ones(Float64,             slob.num_paths)
    
    recalc = slob.α > 0.0 
    
    #print(seeds[1]," ")
    
    counter = 0
    for path_num = 1:slob.num_paths
        broke_point = -1;
        Random.seed!(seeds[path_num])
        
        @info "path $path_num with seed $(seeds[path_num])"
        
        dtrw_solver_fractional(outer_dict[path_num])
       
        counter += 1
        if (progress==true)
            next!(p)
        end
    end
    
    if (progress==true)
        finish!(p)
    end

    return   outer_dict
end

# +
outer1 = Dict(1:3 .=> 1:3)
outer2 = Dict(1:3 .=> 1:3)


outer = Dict(1:2 .=> (outer1,outer2))

length(outer)


# +
# mutable struct DataPasser
#     lob_densities::Array{Float64,2}
#     sources::Array{Float64,2}
#     couplings::Array{Float64,2}  
#     rl_pushes::Array{Float64,2}
#     raw_price_paths::Array{Float64,1}
#     obs_price_paths::Array{Float64,1}
#     P⁺s::Array{Float64,1}
#     P⁻s::Array{Float64,1}
#     Ps::Array{Float64,1} 
#     V::Array{Float64,1}
# end


# outer_dict = Dict{Int64,Dict{Int64,DataPasser}}()

# +
d =   zeros(Float64, 200, 10000, 2, 10);
function loop_over1()
    for first in 1:size(d)[1]
        for second in 1:size(d)[2]
            for third in 1:size(d)[3]
                for fourth in 1:size(d)[4]
                    d[first,second,third,fourth]
                end
            end
        end
    end
end

function loop_over2()
    for first in 1:size(d)[4]
        for second in 1:size(d)[3]
            for third in 1:size(d)[2]
                for fourth in 1:size(d)[1]
                    d[fourth,third,second,first]
                end
            end
        end
    end
end
# -

@elapsed loop_over1()

@elapsed loop_over2()

# +
a = zeros(Float64,2,3)

b = @view a[1,:]

b[3]=1

a

# +
#SLOB(dict) = SLOB(
#    dict["num_paths"], dict["T"], dict["p₀"], dict["M"],
#    dict["L"], dict["D"], dict["σ"], dict["nu"], dict["α"],
#    SourceTerm(dict["λ"], dict["μ"]),InteractionTerm(dict["κ"],dict["β"])) #parse a dictionary form of the input


# +
#SLOB(;num_paths::Int64=1, T::Int64=100, p₀::Real=100.0,
#    M::Int64=100, L::Real=50.0, D::Real=4.0, σ::Real=1.0,
#    nu::Real=0.1, α::Real=20.0, λ::Real=1.0, μ::Real=0.5, κ::Real=0.5, β::Real=0.5) =
#    SLOB(num_paths, T, p₀, M, L, D, σ, nu, α,
#    SourceTerm(λ, μ),InteractionTerm(κ, β)) #


# +
#not changed:
# function InteractOrderBooks(slob¹::SLOB,slob²::SLOB,kick,seed::Int=-1)
#     if seed == -1
#             seeds = Int.(rand(MersenneTwister(), UInt32, slob¹.num_paths))
#         else
#             seeds = Int.(rand(MersenneTwister(seed), UInt32, slob¹.num_paths))
#     end

#     time_steps = get_time_steps(slob¹.T, slob¹.Δt)

#     mid_price_paths = ones(Float64, slob.T + 1, slob.num_paths)
#     for path in 1:slob.num_paths
#         Random.seed!(seeds[path])#! means it does modify its arguments
#         mid_price_paths[:, path] = dtrw_solver(slob)#discrete time random walk solver
#     end #[rows,columns]

#     return mid_price_paths
# end


# +
###############################################Ignore for now
# function (slob::SLOB)(debug::Bool, seed::Int=-1)#same but returns more stuff
#     logger = FileLogger(Dict(Logging.Info => "info.log", Logging.Error => "error.log"), append=false)
#     global oldglobal = global_logger(logger)
#     if seed == -1
#             seeds = Int.(rand(MersenneTwister(), UInt32, slob.num_paths))
#         else
#             seeds = Int.(rand(MersenneTwister(seed), UInt32, slob.num_paths))
#     end
#     time_steps = get_time_steps(slob.T, slob.Δt)

#     raw_price_paths = ones(Float64, time_steps + 1, slob.num_paths)

#     sample_price_paths = ones(Float64, slob.T + 1, slob.num_paths)

#     lob_densities = zeros(Float64, slob.M+1, time_steps + 1, slob.num_paths)#[x_entries,time t,path i]
    
#     P⁺s = ones(Float64, time_steps, slob.num_paths)
#     P⁻s = ones(Float64, time_steps, slob.num_paths)
#     Ps = ones(Float64, time_steps, slob.num_paths)
    
#     for path in 1:slob.num_paths
#         Random.seed!(seeds[path])
#         @info "path $path with seed $(seeds[path])"
#         lob_densities[:, :, path], raw_price_paths[:, path],
#             sample_price_paths[:, path], P⁺s[:, path],
#             P⁻s[:, path], Ps[:, path] = dtrw_solver(debug, slob)
#     end

#     return lob_densities, raw_price_paths, sample_price_paths, P⁺s, P⁻s, Ps
# end

# +
###############################################Ignore for now
# function (slob::SLOB)(seed::Int=-1)
#     if seed == -1
#             seeds = Int.(rand(MersenneTwister(), UInt32, slob.num_paths))
#         else
#             seeds = Int.(rand(MersenneTwister(seed), UInt32, slob.num_paths))
#     end

#     time_steps = get_time_steps(slob.T, slob.Δt)

#     mid_price_paths = ones(Float64, slob.T + 1, slob.num_paths)
#     for path in 1:slob.num_paths
#         Random.seed!(seeds[path])#! means it does modify its arguments
#         mid_price_paths[:, path] = dtrw_solver(slob)#discrete time random walk solver
#     end #[rows,columns]

#     return mid_price_paths
# end

