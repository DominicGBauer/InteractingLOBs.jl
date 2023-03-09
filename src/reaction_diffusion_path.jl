# -*- coding: utf-8 -*-
mutable struct SLOB
    num_paths::Int64 #How many paths to simulate
    T::Int64    # How many time periods it runs for
    p₀::Float64 # The initial mid price
    M::Int64    # Number of discrete price points in addition to p₀ (the central point)
    L::Float64  # The length of the price dimension in the lattice
    D::Float64  # The diffusion constant (determines Δt)
    σ::Float64  # Standard deviation in V(0,σ)
    nu::Float64 # Cancellation rate
    α::Float64  # Waiting time in terms of Exp(α)
    source_term::SourceTerm # Source term s(x,t)(μ,λ,p¹,p²) 
    interaction_term::InteractionTerm # Interaction term c(x,t)(μ,λ,p¹,p²)
    rl_push_term::RLPushTerm
    x::Array{Float64, 1}    #
    Δx::Float64 # The grid spacing which is calculated as Δx=L/M
    Δt::Float64 # The temporal grid spacing calculated as Δt=Δx*Δx/(2D)
end


function SLOB(num_paths::Int64, T::Int64, p₀::Float64,
    M::Int64, L::Real, D::Float64, σ::Float64, nu::Float64,
    α::Float64, source_term::SourceTerm, interaction_term::InteractionTerm, rl_push_term::RLPushTerm)

    x₀ = p₀ - 0.5*L
    x_m = p₀ + 0.5*L
    @assert x₀ >= 0
    x = collect(Float64, range(x₀, stop=x_m, length=M+1)) #creates an array of the entries. So returns set of x points
    Δx = L/M
    Δt = (Δx^2) / (2.0*D)
    return SLOB(num_paths, T, p₀, M, L, D, σ, nu, α, source_term, interaction_term, rl_push_term, x, Δx, Δt)
end

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
# -


function get_time_steps(T::Int, Δt::Float64)
    time_steps =  T / Δt #returns float like 5/3=1.666666666....
    return floor(Int, time_steps + eps(time_steps)) #makes float tiny bit bigger then rounds down
end


function get_mid_price_inds(slob::SLOB) 
    get_mid_price_inds_temp(x) = get_time_steps(x, slob.Δt) #return ≈ floor(x/Δt). These are mid price indexes
    mid_price_inds = get_mid_price_inds_temp.(0:slob.T) .+ 1
end

function sample_mid_price_path(slob::SLOB, price_path)
    mid_price_inds = get_mid_price_inds(slob)
    mid_prices = price_path[mid_price_inds]
    return mid_prices
end

# +
#########################My own code below here########################
# -

function InteractOrderBooks(slob¹::SLOB,slob²::SLOB, seed::Int=-1)#same but returns more stuff
    # do some logging IDK
    logger = FileLogger(Dict(Logging.Info => "info.log", Logging.Error => "error.log"), append=false)
    global oldglobal = global_logger(logger)
    
    # generate a set of seeds
    if seed == -1
            seeds = Int.(rand(MersenneTwister(), UInt32, slob¹.num_paths))
        else
            seeds = Int.(rand(MersenneTwister(seed), UInt32, slob¹.num_paths))
    end
    
    
    time_steps = get_time_steps(slob¹.T, slob¹.Δt)

    raw_price_paths¹ = ones(Float64, time_steps + 1, slob¹.num_paths)
    raw_price_paths² = ones(Float64, time_steps + 1, slob².num_paths)

    sample_price_paths¹ = ones(Float64, slob¹.T + 1, slob¹.num_paths)
    sample_price_paths² = ones(Float64, slob².T + 1, slob².num_paths)

    lob_densities¹ = zeros(Float64, slob¹.M+1, time_steps + 1, slob¹.num_paths)#[x_entries,time t,path i]
    lob_densities² = zeros(Float64, slob².M+1, time_steps + 1, slob².num_paths)#[x_entries,time t,path i]
    
    sources¹ = zeros(Float64, slob¹.M+1, time_steps + 1, slob¹.num_paths)#[x_entries,time t,path i]
    sources² = zeros(Float64, slob².M+1, time_steps + 1, slob².num_paths)#[x_entries,time t,path i]
    
    couplings¹ = zeros(Float64, slob¹.M+1, time_steps + 1, slob¹.num_paths)#[x_entries,time t,path i]
    couplings² = zeros(Float64, slob².M+1, time_steps + 1, slob².num_paths)#[x_entries,time t,path i]
    
    rl_pushes¹ = zeros(Float64, slob¹.M+1, time_steps + 1, slob¹.num_paths)#[x_entries,time t,path i]
    rl_pushes² = zeros(Float64, slob².M+1, time_steps + 1, slob².num_paths)#[x_entries,time t,path i]
    
    P⁺s¹ = ones(Float64, time_steps, slob¹.num_paths)
    P⁻s¹ = ones(Float64, time_steps, slob¹.num_paths)
    Ps¹ = ones(Float64, time_steps, slob¹.num_paths)
    
    P⁺s² = ones(Float64, time_steps, slob².num_paths)
    P⁻s² = ones(Float64, time_steps, slob².num_paths)
    Ps² = ones(Float64, time_steps, slob².num_paths)
    
    recalc = slob¹.α > 0.0 
    
    #for path in 1:slob¹.num_paths
    counter = 0
    lk = Threads.SpinLock()
    Threads.@threads for path = 1:slob¹.num_paths
        Random.seed!(seeds[path])
        
        @info "path $path with seed $(seeds[path])"
        lob_densities¹[:, :, path], sources¹[:,:,path], couplings¹[:,:,path], rl_pushes¹[:,:,path],
            raw_price_paths¹[:, path], sample_price_paths¹[:, path], 
            P⁺s¹[:, path], P⁻s¹[:, path], Ps¹[:, path],
        lob_densities²[:, :, path], sources²[:,:,path], couplings²[:,:,path], rl_pushes²[:,:,path],
            raw_price_paths²[:, path], sample_price_paths²[:, path], 
            P⁺s²[:, path], P⁻s²[:, path], Ps²[:, path] =
        dtrw_solver(slob¹,slob², recalc)
        
        lock(lk)
        try
            counter += 1
            #Base.buffer_writes(string(counter," "))
            #print(string(counter," "))
            write(stdout,string(counter," "))
            flush(stdout)
        finally
            unlock(lk)
        end
    end

    return   lob_densities¹, sources¹, couplings¹, rl_pushes¹, raw_price_paths¹, sample_price_paths¹, P⁺s¹, P⁻s¹, Ps¹,
             lob_densities², sources², couplings², rl_pushes², raw_price_paths², sample_price_paths², P⁺s², P⁻s², Ps²
end

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

# -

Base.buffer_writes(stdout)="test"
flush(stdout)




