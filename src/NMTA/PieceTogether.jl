using DataFrames, CSV, Dates, ProgressMeter, JLD, InteractingLOBs, Distributions
using NLSolversBase: value, value!, value!!, NonDifferentiable

include("./Moments.jl")
include("./NMTA.jl")
include("../setup.jl")

mutable struct Parameters
    M::Int64
    L::Int64
    D::Float64
    ν::Float64 #removal rate
    γ::Float64 #fraction of derivative (1 is normal diffusion, less than 1 is D^{1-γ} derivative on the RHS)
	function Parameters(; M = 400, L = 200, D = 0.5, ν = 14.0, γ = 1.0) # 1000, 25000
		new(M, L, D, ν, γ)
	end
end

function Optimize(f::NonDifferentiable{Tf, Tx}, initial_x::Tx, options::Options{T} = Options(), state::NelderMeadState = InitialState(NelderMead(), f, initial_x)) where {Tx <: AbstractArray, Tf <: Real, T <: Real}
    t₀ = time() # Initial time stamp used to control early stopping by options.time_limit
    tr = OptimizationTrace{Tf}() # Store optimization trace
    thresholds = !isempty(options.ta_rounds) ? reduce(vcat, fill.(options.f_reltol, options.ta_rounds)) : zeros(options.iterations)
    tracing = options.store_trace || options.show_trace || options.extended_trace
    stopped, stopped_by_time_limit, f_increased = false, false, false
    g_converged = InitialConvergence(state, initial_x, options) # Converged if criterion is met
    # iteration = 0 # Counter
    if options.show_trace # Print header
        @printf "Iter     Function value    √(Σ(yᵢ-ȳ)²)/n \n"
        @printf "------   --------------    --------------\n"
    end
    t = time()
    Trace!(tr, state, state.iteration, options, t - t₀)
    start_time = now()
    while !g_converged && !stopped_by_time_limit && state.iteration < options.iterations
        state.iteration += 1
        println(string(thresholds[state.iteration], "      ", thresholds[state.iteration] * sum(state.f_simplex) / state.m))
        # println(run(`free -m`))
        println("Iterations = ", state.iteration)
        sleep(1)
        if rand() < options.ξ
            try
                ThresholdAccepting!(f, state, thresholds[state.iteration] * (sum(state.f_simplex) / state.m))
            catch e
                println(e)
                @error "Something went wrong" exception=(e, catch_backtrace())
                PostProcessError!(f, state)
                return OptimizationResults{Tx, Tf}(initial_x, state.x, value(f), state.iteration, state.iteration == options.iterations, options.f_reltol, g_converged, Float64(options.g_abstol), f_increased, tr, options.time_limit, t - t₀, stopped_by_time_limit, state)
            end
        else
            try
                SimplexSearch!(f, state, thresholds[state.iteration] * (sum(state.f_simplex) / state.m)) # Percentage of best solution
            catch e
                println(e)
                @error "Something went wrong" exception=(e, catch_backtrace())
                PostProcessError!(f, state)
                return OptimizationResults{Tx, Tf}(initial_x, state.x, value(f), state.iteration, state.iteration == options.iterations, options.f_reltol, g_converged, Float64(options.g_abstol), f_increased, tr, options.time_limit, t - t₀, stopped_by_time_limit, state)
            end
        end
        g_converged, f_increased = AssessConvergence(state, options)
        if tracing
            Trace!(tr, state, state.iteration, options, time() - t₀)
        end
        t = time()
        stopped_by_time_limit = t - t₀ > options.time_limit
    end
    PostProcess!(f, state)
    println("While loop time = ", now() - start_time)
    return OptimizationResults{Tx, Tf}(initial_x, state.x, value(f), state.iteration, state.iteration == options.iterations, options.f_reltol, g_converged, Float64(options.g_abstol), f_increased, tr, options.time_limit, t - t₀, stopped_by_time_limit, state)
end

#----- Generate emperical and simulated log-returns and emperical moments -----#
function GenerateEmpericalReturnsAndMoments(startTime::DateTime, endTime::DateTime)
    println("Generating returns and moments for: " * Dates.format(startTime, "yyyy-mm-ddTHH:MM:SS") * " to " * Dates.format(endTime, "yyyy-mm-ddTHH:MM:SS"))
    empericalData = CSV.File(string("Data/L1LOB.csv"), missingstring = "missing", types = Dict(:DateTime => DateTime, :Type => Symbol)) |> DataFrame
    filter!(x -> startTime <= x.DateTime, empericalData)
    filter!(x -> !ismissing(x.MidPrice), empericalData); filter!(x -> !ismissing(x.MicroPrice), empericalData)
    midPriceLogReturns = diff(log.(empericalData.MidPrice))
    microPriceLogReturns = diff(log.(empericalData.MicroPrice))
    empericalLogReturns = DataFrame(MidPriceLogReturns = midPriceLogReturns, MicroPriceLogReturns = microPriceLogReturns)
    empericalMidPriceMoments = Moments(midPriceLogReturns, midPriceLogReturns)
    empericalMicroPriceMoments = Moments(microPriceLogReturns, microPriceLogReturns)
    empericalMoments = Dict("empericalMidPriceMoments" => empericalMidPriceMoments, "empericalMicroPriceMoments" => empericalMicroPriceMoments)
    return empericalLogReturns, empericalMoments
end

# function GenerateSimulatedReturnsAndMoments()
#     simulatedPrices = Data[1][1].raw_price_paths[1:s]
#     logReturns = diff(log.(simulatedPrices))
#     moments = Moments(logReturns, logReturns)
#     return logReturns, moments
# end

date = DateTime("2019-07-08")
startTime = date + Hour(9) + Minute(1)
endTime = date + Hour(16) + Minute(50)
empiricalLogReturns, empiricalMoments = GenerateEmpericalReturnsAndMoments(startTime, endTime)
# simulatedLogReturns, simulatedMoments = GenerateSimulatedReturnsAndMoments()


#----- Create weights -----#
function MovingBlockBootstrap(logreturns::Vector{Float64}, iterations::Int64 = 1000, windowsize::Int64 = 2000)
    bootstrapmoments = fill(0.0, (iterations, 8))
    @showprogress "Computing weight matrix..." for i in 1:iterations
        indeces = rand(1:(length(logreturns) - windowsize + 1), Int(ceil(length(logreturns)/windowsize)))
        bootstrapreturns = Vector{Float64}()
        for index in indeces
            append!(bootstrapreturns, logreturns[index:(index  + windowsize - 1)])
        end
        moments = Moments(bootstrapreturns[1:length(logreturns)], logreturns)
        bootstrapmoments[i,:] = [moments.μ moments.σ moments.ks moments.hurst moments.gph moments.adf moments.garch moments.hill]
    end
    DataFrame(bootstrapmoments, Symbol.(["Mean","Std","KS","Hurst","GPH","ADF","GARCH","Hill"]))
    W = inv(cov(bootstrapmoments))
    save("Data/Calibration/W.jld", "W", W)
end
# MovingBlockBootstrap(empiricalLogReturns.MicroPriceLogReturns)

function MovingBlockBootstrapSimulated(logreturns::Vector{Float64}, empiricallogreturns::Vector{Float64}, iterations::Int64 = 1000, windowsize::Int64 = 2000)
    bootstrapmoments = fill(0.0, (iterations, 8))
    @showprogress "Computing weight matrix..." for i in 1:iterations
        indeces = rand(1:(length(logreturns) - windowsize + 1), Int(ceil(length(logreturns)/windowsize)))
        bootstrapreturns = Vector{Float64}()
        for index in indeces
            append!(bootstrapreturns, logreturns[index:(index  + windowsize - 1)])
        end
        moments = Moments(bootstrapreturns[1:length(logreturns)], empiricallogreturns)
        bootstrapmoments[i,:] = [moments.μ moments.σ moments.ks moments.hurst moments.gph moments.adf moments.garch moments.hill]
    end
    bootstrapmoments_df = DataFrame(bootstrapmoments, Symbol.(["Mean","Std","KS","Hurst","GPH","ADF","GARCH","Hill"]))
    CSV.write("Data/Calibration/SimulatedBootstrapMoments.csv", bootstrapmoments_df)
    W = inv(cov(bootstrapmoments))
    save("Data/Calibration/SimulatedW.jld", "SimulatedW", W)
end
# MovingBlockBootstrapSimulated(simulatedLogReturns, empiricalLogReturns.MicroPriceLogReturns)
# Note replications is I the number of the Monte Carlo replications to use
function WeightedSumofSquaredErrors(parameters::Parameters, replications::Int64, W::Array{Float64, 2}, empiricalmoments::Moments, empiricallogreturns::Vector{Float64})
    errormatrix = fill(0.0, (replications, 8))
    s = 2299
    for i in 1:replications
        lob_model¹ = SLOB(2, s, p₀, parameters.M, parameters.L, parameters.D, parameters.ν, α, parameters.γ, mySourceTerm, myCouplingTerm, myRLPusher1, myRandomnessTerm);
        lob_model² = SLOB(2, s, p₀, parameters.M, parameters.L, parameters.D, parameters.ν, α, parameters.γ, mySourceTerm, myCouplingTerm, myRLPusher2, myRandomnessTerm);
        Data = InteractOrderBooks([lob_model¹,lob_model²], -1, true);
        microprice = Data[1][1].raw_price_paths[1:s]
        # midprice, microprice = simulate(parameters, false, false, seed = i)
        if !isempty(microprice)
            filter!(x -> !isnan(x), microprice)
            logreturns = diff(log.(microprice))
            try
                simulatedmoments = Moments(logreturns, empiricallogreturns)
                errormatrix[i, :] = [simulatedmoments.μ-empiricalmoments.μ simulatedmoments.σ-empiricalmoments.σ simulatedmoments.ks-empiricalmoments.ks simulatedmoments.hurst-empiricalmoments.hurst simulatedmoments.gph-empiricalmoments.gph simulatedmoments.adf-empiricalmoments.adf simulatedmoments.garch-empiricalmoments.garch simulatedmoments.hill-empiricalmoments.hill]
            catch e
                println(e)
                errormatrix[i, :] = errormatrix[i - 1, :]
            end
        else
            return Inf
        end
    end
    GC.gc() # Garbage collection
    errors = mean(errormatrix, dims = 1)
    return (errors * W * transpose(errors))[1]
end

ta_rounds_arg = [5, 10, 20, 30, 35]
f_reltol_arg = [0.3, 0.2, 0.1, 0.05, 0]
initialsolution = [400, 200, 0.5, 14.0, 0.8]


function Calibrate(initialsolution::Vector{Float64}, empiricallogreturns::Vector{Float64}, empiricalmoments::Moments; f_reltol::Vector{Float64} = [0.3, 0.2, 0.1, 0], ta_rounds::Vector{Int64} = [4, 3, 2, 1], neldermeadstate = nothing)
    W = load("Data/Calibration/W.jld")["W"]
    # W = fill(1.0, 8, 8)
    objective = NonDifferentiable(x -> WeightedSumofSquaredErrors(
            Parameters(),
            5,
            W,
            empiricalmoments,
            empiricallogreturns
        ),
        initialsolution
    )
    optimizationoptions = Options(show_trace = true, store_trace = true, trace_simplex = true, extended_trace = true, iterations = sum(ta_rounds), ξ = 0.15, ta_rounds = ta_rounds, f_reltol = f_reltol)
    @time result = !isnothing(neldermeadstate) ? Optimize(objective, initialsolution, optimizationoptions, neldermeadstate) : Optimize(objective, initialsolution, optimizationoptions)
    save("Data/Calibration/OptimizationResult.jld", "result", result)
end

# @time Calibrate(initialsolution, empiricalLogReturns.MicroPriceLogReturns, empiricalMoments["empericalMicroPriceMoments"], ta_rounds = ta_rounds_arg, f_reltol = f_reltol_arg) # , neldermeadstate = neldermeadstate) [takes about 8hrs]

function PlotObjectiveConvergence(stacktrace)

    iters = iterations(stacktrace) + 1
    f = zeros(Float64, iters); g_norm = zeros(Float64, iters); f_simplex = fill(0.0, iters, 6)#; centr = fill(0.0, length(stacktrace), 5); metadata = Vector{Dict{Any, Any}}()
    i = 1
    j = 1
    for s in trace(stacktrace)
        f[i] = s.value                         # vertex with the lowest value (lowest with a tolerence in the begining)
        g_norm[i] = s.g_norm                   # √(Σ(yᵢ-ȳ)²)/n
        f_simplex[i, :] = transpose(s.metadata["simplex_values"])
        i += 1
        j += 1
    end
    # Objectives
    objectives = plot(1:iters, f, seriestype = :line, linecolor = :blue, label = "Weighted SSE objective", xlabel = "Iteration", ylabel = "Weighted SSE objective", legendfontsize = 5, fg_legend = :transparent, tickfontsize = 5, xaxis = false, xticks = false, legend = :bottomleft, guidefontsize = 7, yscale = :log10, minorticks = true, left_margin = 5Plots.mm, right_margin = 15Plots.mm, fontfamily = "Computer Modern")
    plot!(twinx(), 1:iters, g_norm, seriestype = :line, linecolor = :purple, label = "Convergence criterion", ylabel = "Convergence criterion", legend = :topright, legendfontsize = 5, fg_legend = :transparent, tickfontsize = 5, yscale = :log10, minorticks = true, guidefontsize = 7)
    savefig(objectives, "Data/Calibration/NMTAFitnessBestVertexOG.png")
    # Simplex values
    convergence = plot(1:iters, f_simplex, seriestype = :line, linecolor = [:blue :purple :green :orange :red :black :magenta], xlabel = "Iteration", ylabel = "Weighted SSE objective", legend = false, tickfontsize = 5, guidefontsize = 7, yscale = :log10, minorticks = true, fontfamily = "Computer Modern")
    savefig(convergence, "Data/Calibration/NMTAFitnessAllSimplexValuesOG.png")
end
# stacktrace = load("Data/Calibration/OptimizationResult.jld")["result"]
# PlotObjectiveConvergence(stacktrace)
#----------------------------------------

function ParameterTracePlots(stacktrace)
    iters = iterations(stacktrace) + 1
    parameters = [("M", raw"$M$"), ("L", raw"$L$"), ("D",raw"$D$"), ("ν", raw"$\nu$"), ("γ", raw"$\gamma$")]
    c = [:blue :purple :green :orange :red :black]
    for (i,param) in enumerate(parameters)
        t = fill(0.0, iters, 6)
        for (j,s) in enumerate(simplex_trace(stacktrace))
            t[j,:] = transpose(hcat(s...))[:,i]
        end
        p = plot(1:iters, t, seriestype = :line, linestyle = :dash, linecolor = c[i], xlabel = "Iteration", ylabel = last(param), legend = false, tickfontsize = 5, guidefontsize = 7, minorticks = true, fontfamily = "Computer Modern", dpi=300)
        plot!(1:iters, transpose(hcat(centroid_trace(stacktrace)...))[:,i], seriestype = :line, linestyle = :solid, linecolor = c[i], linewidth = 2, legend = false)
        savefig(p, "Data/Calibration/ParameterConvergence" * first(param) * ".png")
    end

end
# ParameterTracePlots(stacktrace)

# final_nmta_values = stacktrace.end_state.x_centroid
