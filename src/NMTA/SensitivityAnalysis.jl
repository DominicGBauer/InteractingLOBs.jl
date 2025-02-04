#=
SensitivityAnalysis:
- Julia version: 1.7.1
- Authors: Matthew Dicks, Tim Gebbie, (some code was adapted from https://github.com/IvanJericevich/IJPCTG-ABMCoinTossX)
- Function: Perform the sensitivity analysis, plot the results and save them in figures
- Structure:
    1. Sensitivity analysis
    2. Visualisations
- Examples:
    1. Sensitivity analysis
        date = DateTime("2019-07-08")
        startTime = date + Hour(9) + Minute(1)
        endTime = date + Hour(16) + Minute(50) # Hour(17)
        empericalLogReturns, empericalMoments = GenerateEmpericalReturnsAndMoments(startTime, endTime)
        NᴸₜRange = [3,6,9,12]
        NᴸᵥRange = [3,6,9,12]
        δRange = collect(range(0.01, 0.2, length = 4))
        κRange = collect(range(2, 5, length = 4))
        νRange = collect(range(2, 8, length = 4))
        σᵥRange = collect(range(0.0025, 0.025, length = 4))
        parameterCombinations = GenerateParameterCombinations(NᴸₜRange, NᴸᵥRange, δRange, κRange, νRange, σᵥRange)
        @time SensitivityAnalysis(empericalLogReturns, empericalMoments, parameterCombinations, parameterCombinationsRange) [takes about 30hrs]
    2. Visualisations
        MomentViolinPlots("MicroPrice", true); MomentViolinPlots("MidPrice", true)
        MomentInteractionSurfaces("MicroPrice", false); MomentInteractionSurfaces("MidPrice", false)
        ObjectiveInteractionSurfaces("MicroPrice", false); ObjectiveInteractionSurfaces("MidPrice", false)
        ParameterMomentCorrelationMatrix("MicroPrice", false); ParameterMomentCorrelationMatrix("MidPrice", false)
=#
ENV["JULIA_COPY_STACKS"]=1
using InteractingLOBs, ProgressMeter, CSV, Plots, DataFrames, StatsPlots, Statistics, ColorSchemes, Dates, JLD, Combinatorics, Colors, LaTeXStrings
using LinearAlgebra: diag, inv, transpose
include("./Moments.jl")

# Make sure to keep this the same as in PieceTogether.jl
mutable struct Parameters
    M::Int64
    L::Int64
    D::Float64
    ν::Float64 #removal rate
    γ::Float64 #fraction of derivative (1 is normal diffusion, less than 1 is D^{1-γ} derivative on the RHS)
	function Parameters(; M = 400.0, L = 200.0, D = 0.5, ν = 14.0, γ = 1.0) # 1000, 25000
		new(M, L, D, ν, γ)
	end
end

#----- Sensitivity analysis -----#

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

#---------------------------------------------------------------------------------------------------

#----- Generating parameter combinations -----#
function GenerateParameterCombinations(DRange::Vector{Float64}, νRange::Vector{Float64}, γRange::Vector{Float64})
    println("Generating parameter combinations")
    parameterCombinations = Vector{Parameters}()
    for D in DRange
        for γ in γRange
            for ν in νRange
                parameters = Parameters(M = 400, L = 200, D = D, ν = ν, γ = γ)
                push!(parameterCombinations, parameters)
            end
        end
    end
    return parameterCombinations
end
#---------------------------------------------------------------------------------------------------

#----- Sensitivity analysis -----#
function SensitivityAnalysis(empericalLogReturns::DataFrame, empericalMoments::Dict, parameterCombinations::Vector{Parameters})
    open("Data/SensitivityAnalysis/SensitivityAnalysisResults.csv", "w") do file
        println(file, "Type,M,L,D,Nu,Gamma,Seed,Mean,Std,Kurtosis,KS,Hurst,GPH,ADF,GARCH,Hill")
        for (i, parameters) in enumerate(parameterCombinations)
            try
                seed = 1
                lob_model¹ = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,mySourceTerm, myCouplingTerm, myRLPusher1, myRandomnessTerm);
                lob_model² = SLOB(num_paths, T, p₀, M, L, D, ν, α, γ,mySourceTerm, myCouplingTerm, myRLPusher2, myRandomnessTerm);
                Data = InteractOrderBooks([lob_model¹,lob_model²], -1, true);
                microPrices = Data[1][1].raw_price_paths[1:s]
                if isnothing(microPrices)
                    println("\nParameter Set: $(i-1) finished\n")
                    break
                end
                println("\nParameter Set: $(i)\n")
                # println(run(`free -m`))
                filter!(x -> !ismissing(x) && !(isnan(x)), microPrices)
                microPriceLogReturns = diff(log.(microPrices))
                simulatedMicroPriceMoments = Moments(microPriceLogReturns, empericalLogReturns.MicroPriceLogReturns)
                println(file, "MicroPrice,", parameters.M, ",", parameters.L, ",", parameters.D, ",", parameters.ν, ",", parameters.γ, ",", seed, ",", simulatedMicroPriceMoments.μ, ",", simulatedMicroPriceMoments.σ, ",", simulatedMicroPriceMoments.κ, ",", simulatedMicroPriceMoments.ks, ",", simulatedMicroPriceMoments.hurst, ",", simulatedMicroPriceMoments.gph, ",", simulatedMicroPriceMoments.adf, ",", simulatedMicroPriceMoments.garch, ",", simulatedMicroPriceMoments.hill)
                GC.gc()             # perform garbage collection

            catch e
                println(e)
                println(file, "MicroPrice,", parameters.M, ",", parameters.L, ",", parameters.D, ",", parameters.ν, ",", parameters.γ, ",", seed, ",", NaN, ",", NaN, ",", NaN, ",", NaN, ",", NaN, ",", NaN, ",", NaN, ",", NaN, ",", NaN)
            end
        end
    end
end
#---------------------------------------------------------------------------------------------------

## make sure these are the same for the stylized facts and Calibration
date = DateTime("2019-07-08")
startTime = date + Hour(9) + Minute(1)
endTime = date + Hour(16) + Minute(50) # Hour(17)

empericalLogReturns, empericalMoments = GenerateEmpericalReturnsAndMoments(startTime, endTime)

# NᴸₜRange = [3,6,9,12]
# NᴸᵥRange = [3,6,9,12]
# δRange = collect(range(0.01, 0.2, length = 4))
DRange = collect(range(0.1, 1, length = 4))
νRange = [1.0,3.0,6.0,9.0]
γRange = collect(range(0.0025, 0.025, length = 4))

parameterCombinations = GenerateParameterCombinations(DRange, νRange, γRange)

@time SensitivityAnalysis(empericalLogReturns, empericalMoments, parameterCombinations)

#---------------------------------------------------------------------------------------------------

#----- Visualizations -----#

#----- Compute Objective Function -----#
function ComputeObjective(empericalMoments::Dict)
    W = load("Data/Calibration/W.jld")["W"]
    sr = CSV.File(string("Data/SensitivityAnalysis/SensitivityAnalysisResults.csv")) |> DataFrame
    objs = Vector{Float64}()
    for i in 1:nrow(sr)
        if i % 2 == 0
            errors = [sr[i,:Mean]-empericalMoments["empericalMicroPriceMoments"].μ sr[i,:Std]-empericalMoments["empericalMicroPriceMoments"].σ sr[i,:KS]-empericalMoments["empericalMicroPriceMoments"].ks sr[i,:Hurst]-empericalMoments["empericalMicroPriceMoments"].hurst sr[i,:GPH]-empericalMoments["empericalMicroPriceMoments"].gph sr[i,:ADF]-empericalMoments["empericalMicroPriceMoments"].adf sr[i,:GARCH]-empericalMoments["empericalMicroPriceMoments"].garch sr[i,:Hill]-empericalMoments["empericalMicroPriceMoments"].hill]
        else
            errors = [sr[i,:Mean]-empericalMoments["empericalMidPriceMoments"].μ sr[i,:Std]-empericalMoments["empericalMidPriceMoments"].σ sr[i,:KS]-empericalMoments["empericalMidPriceMoments"].ks sr[i,:Hurst]-empericalMoments["empericalMidPriceMoments"].hurst sr[i,:GPH]-empericalMoments["empericalMidPriceMoments"].gph sr[i,:ADF]-empericalMoments["empericalMidPriceMoments"].adf sr[i,:GARCH]-empericalMoments["empericalMidPriceMoments"].garch sr[i,:Hill]-empericalMoments["empericalMidPriceMoments"].hill]
        end
        obj = errors * W * transpose(errors)
        push!(objs, obj[1])
    end
    insertcols!(sr, ncol(sr) + 1, "Objective" => objs)
    CSV.write("Data/SensitivityAnalysis/SensitivityAnalysisResultsObj.csv", sr)
end

ComputeObjective(empericalMoments)
#---------------------------------------------------------------------------------------------------

#----- Moment Values For Parameter Marginals -----#
function Winsorize(paramvalues, momentvalues)
    df = DataFrame(ParamValues = paramvalues, MomentValues = momentvalues)
    upper = quantile(df.MomentValues, 0.99)
    lower = quantile(df.MomentValues, 0.01)
    df_winsor = df[findall(x -> lower < x && x < upper, df.MomentValues),:]
    return df_winsor.ParamValues, df_winsor.MomentValues
end

function MomentViolinPlots(midmicro::String, winsorize::Bool)
    sr = CSV.File(string("Data/SensitivityAnalysis/SensitivityAnalysisResultsObj.csv")) |> DataFrame
    sr = sr[findall(x -> x == midmicro, sr.Type),:]
    colors = ["blue", "red", "green", "magenta", "orange", "purple"]
    parameters = [("Nt", raw"$N^{c}_{_{\mathrm{LT}}}$"), ("Nv", raw"$N^{f}_{_{\mathrm{LT}}}$"), ("Delta",raw"$\delta$"), ("Kappa", raw"$\kappa$"), ("Nu", raw"$\nu$"), ("SigmaV", raw"$\sigma_{f}$")]
    moments = [("Kurtosis", "Kurtosis"), ("KS", "KS"), ("Hurst", "Hurst Exponent"), ("GPH", "GPH Statistic"), ("ADF", "ADF Statistic"), ("GARCH", "GARCH"), ("Hill", "Hill Estimator"), ("Objective", "Objective Function")]
    for (i, (paramcol, paramlabel)) in enumerate(parameters)
        col = colors[i]
        for (momentcol, momentlabel) in moments
            params_sr = sr[:,paramcol]
            moments_sr = sr[:,momentcol]
            if winsorize
                params_sr, moments_sr = Winsorize(params_sr, moments_sr)
            end
            if paramcol == "Delta" || paramcol == "SigmaV"
                p = violin(string.(round.(params_sr, digits = 4)), moments_sr, quantiles = [0.025, 0.975], trim = true, show_median = true, tick_direction = :out, fillcolor = col, legend = false, xrotation = 30, yrotation = 30, tickfontsize = 12, guidefontsize = 22, xlabel = paramlabel, ylabel = momentlabel, left_margin = 5Plots.mm, bottom_margin = 5Plots.mm, fontfamily = "Computer Modern")
            else
                p = violin(round.(params_sr, digits = 4), moments_sr, quantiles = [0.025, 0.975], trim = true, show_median = true, tick_direction = :out, fillcolor = col, legend = false, xrotation = 30, yrotation = 30, tickfontsize = 12, guidefontsize = 22, xlabel = paramlabel, ylabel = momentlabel, left_margin = 5Plots.mm, bottom_margin = 5Plots.mm, fontfamily = "Computer Modern")
            end
            savefig(p, "../Images/SensitivityAnalysis/Violin/NoKurtosis/" * midmicro * "Images/" * paramcol * momentcol * ".pdf")
        end
    end
end
# MomentViolinPlots("MicroPrice", true)
# MomentViolinPlots("MidPrice", true)
#---------------------------------------------------------------------------------------------------

#----- Moment Surfaces For Parameter Interactions -----#
function MomentInteractionSurfaces(midmicro::String, winsorize::Bool)
    sr = CSV.File(string("../Data/SensitivityAnalysis/SensitivityAnalysisResultsObj.csv")) |> DataFrame
    sr = sr[findall(x -> x == midmicro, sr.Type),:]
    colors = ["blue", "red", "green", "magenta", "orange", "purple"]
    parameters = [("Nt", raw"$N^{c}_{_{\mathrm{LT}}}$"), ("Nv", raw"$N^{f}_{_{\mathrm{LT}}}$"), ("Delta",raw"$\delta$"), ("Kappa", raw"$\kappa$"), ("Nu", raw"$\nu$"), ("SigmaV", raw"$\sigma_{f}$")]
    pairwise_combinations = collect(combinations(parameters, 2))
    moments = [("Kurtosis", "Kurtosis"), ("KS", "Kolmogorov-Smirnov"), ("Hurst", "Hurst Exponent"), ("GPH", "GPH Statistic"), ("ADF", "ADF Statistic"), ("GARCH", "GARCH"), ("Hill", "Hill Estimator")]
    for params in pairwise_combinations
        for (momentcol, momentlabel) in moments
            sr_grouped = groupby(sr, [params[1][1], params[2][1]]) |> gdf -> combine(gdf, Symbol(momentcol) => mean)
            surface = plot(unique(sr_grouped[:, params[1][1]]), unique(sr_grouped[:, params[2][1]]), reshape(sr_grouped[:, momentcol * "_mean"], (4,4)), seriestype = :surface, xlabel = params[1][2], ylabel = params[2][2], zlabel = momentlabel, colorbar = false, camera=(45,60), seriesalpha = 0.8, left_margin = 5Plots.mm, right_margin = 15Plots.mm, colorscale = "Viridis", fontfamily = "Computer Modern") # cgrad(ColorScheme((colorant"green", colorant"red", length=10))), color = cgrad([:green, :springgreen4, :firebrick2, :red]),
            savefig(surface, "../Images/SensitivityAnalysis/MomentInteractionSurfaces/" * midmicro * "Images/" * params[1][1] * params[2][1] * momentcol * ".pdf")
        end
    end
end

# MomentInteractionSurfaces("MicroPrice", false)
# MomentInteractionSurfaces("MidPrice", false)
# ---------------------------------------------------------------------------------------------------

#----- Objective Function Surfaces For Parameter Interactions -----#
function ObjectiveInteractionSurfaces(midmicro::String, winsorize::Bool)
    sr = CSV.File(string("../Data/SensitivityAnalysis/SensitivityAnalysisResultsObj.csv")) |> DataFrame
    sr = sr[findall(x -> x == midmicro, sr.Type),:]
    colors = ["blue", "red", "green", "magenta", "orange", "purple"]
    parameters = [("Nt", raw"$N^{c}_{_{\mathrm{LT}}}$"), ("Nv", raw"$N^{f}_{_{\mathrm{LT}}}$"), ("Delta",raw"$\delta$"), ("Kappa", raw"$\kappa$"), ("Nu", raw"$\nu$"), ("SigmaV", raw"$\sigma_{f}$")]
    pairwise_combinations = collect(combinations(parameters, 2))
    for params in pairwise_combinations
        sr_grouped = groupby(sr, [params[1][1], params[2][1]]) |> gdf -> combine(gdf, :Objective => mean)
        surface = plot(unique(sr_grouped[:, params[1][1]]), unique(sr_grouped[:, params[2][1]]), reshape(sr_grouped[:, "Objective_mean"], (4,4)), seriestype = :surface, xlabel = params[1][2], ylabel = params[2][2], zlabel = "Objective", colorbar = false, camera=(45,60), seriesalpha = 0.8, left_margin = 5Plots.mm, right_margin = 15Plots.mm, colorscale = "Viridis", fontfamily = "Computer Modern") # cgrad(ColorScheme((colorant"green", colorant"red", length=10))), color = cgrad([:green, :springgreen4, :firebrick2, :red]),
        savefig(surface, "../Images/SensitivityAnalysis/ObjectiveInteractionSurfaces/NoKurtosis/" * midmicro * "Images/" * params[1][1] * params[2][1] * "Objective.pdf")
    end
end

# ObjectiveInteractionSurfaces("MicroPrice", false)
# ObjectiveInteractionSurfaces("MidPrice", false)
#---------------------------------------------------------------------------------------------------

#----- Correlations Matrix -----#
function ParameterMomentCorrelationMatrix(midmicro::String, winsorize::Bool)
    sr = CSV.File(string("Data/SensitivityAnalysis/SensitivityAnalysisResultsObj.csv")) |> DataFrame
    sr = sr[findall(x -> x == midmicro, sr.Type),:]
    # variables = [("Nt", "Nᴸₜ"), ("Nv", "Nᴸᵥ"), ("Delta","δ"), ("Kappa", "κ"), ("Nu", "ν"), ("SigmaV", "σᵥ"), ("Mean", "Mean"), ("Std", "Std"), ("Kurtosis", "Kurtosis"), ("KS", "KS"), ("Hurst", "Hurst"), ("GPH", "GPH"), ("ADF", "ADF"), ("GARCH", "GARCH"), ("Hill", "Hill"), ("Objective", "Objective")]
    variables = [("Nt", raw"$N^{c}_{_{\mathrm{LT}}}$"), ("Nv", raw"$N^{f}_{_{\mathrm{LT}}}$"), ("Delta",raw"$\delta$"), ("Kappa", raw"$\kappa$"), ("Nu", raw"$\nu$"), ("SigmaV", raw"$\sigma_{f}$"), ("Mean", "Mean"), ("Std", "Std"), ("KS", "KS"), ("Hurst", "Hurst"), ("GPH", "GPH"), ("ADF", "ADF"), ("GARCH", "GARCH"), ("Hill", "Hill"), ("Objective", "Objective")]
    sr = sr[:,first.(variables)]
    C = cor(Matrix(sr))
    (n,m) = size(C)
    H = heatmap(last.(variables), last.(variables), C, c = cgrad(:seismic, [0, 0.28, 0.56, 1]), xticks = (0.5:1:length(variables), last.(variables)), yticks = (0.5:1:length(variables), last.(variables)), xrotation = 45, yrotation = 0,  yflip=true, tickfontsize = 5, tick_direction = :out, alpha = 0.8, fontfamily = "Computer Modern")
    annotate!(H, [(j - 0.5, i - 0.5, text(round(C[i,j],digits=3), 5,:black, :center)) for i in 1:n for j in 1:m])
    savefig(H, "../Images/SensitivityAnalysis/CorrelationMatrix/CorrelationMatrix" * midmicro * ".pdf")
end

# ParameterMomentCorrelationMatrix("MicroPrice", false)
# ParameterMomentCorrelationMatrix("MidPrice", false)
#---------------------------------------------------------------------------------------------------

#----- Exposure Matrix -----#
function ExposureMatrix(midmicro::String)
    sr = CSV.File(string("Data/SensitivityAnalysis/SensitivityAnalysisResultsObj.csv")) |> DataFrame
    sr = sr[findall(x -> x == midmicro, sr.Type),:]
    moments = [("Mean", "Mean"), ("Std", "Std"), ("KS", "KS"), ("Hurst", "Hurst"), ("GPH", "GPH"), ("ADF", "ADF"), ("GARCH", "GARCH"), ("Hill", "Hill")]
    parameters = [("M", "M"), ("L", "L"), ("D","D"), ("Nu", "ν"), ("Gamma", "γ")]
    B = fill(0.0, length(parameters), length(moments))
    for (i, param) in enumerate(parameters)
        for (j, moment) in enumerate(moments)
            B[i,j] = cov(sr[:,first(param)], sr[:,first(moment)]) / var(sr[:,first(moment)])
        end
    end
    save("Data/SensitivityAnalysis/B.jld", "B", B)
end

ExposureMatrix("MicroPrice")
#---------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------


#----- Parameter Confidence Intervals -----#
function ParameterConfidenceIntervals(calibratedParams::Vector{Float64})
    W = load("Data/Calibration/W.jld")["W"]
    B = load("Data/SensitivityAnalysis/B.jld")["B"]
    sigmas = sqrt.(diag(B * inv(W) * transpose(B)))
    upper = calibratedParams .+ (1.96 .* sigmas)
    lower = calibratedParams .- (1.96 .* sigmas)
    parameters = [("M", "M"), ("L", "L"), ("D", "D"), ("Nu", "ν"), ("Gamma", "γ")]
    df = DataFrame(Parameters = first.(parameters), CalibratedParameters = calibratedParams, Lower = lower, Upper = upper)
    CSV.write("Data/Calibration/parameters.csv", df)
end

ParameterConfidenceIntervals([400, 200, 0.27, 12.55, 0.57])
