## Author: Patrick Chang
# Script file to investigate the correction of the Epps effect arising
# from asynchronous sampling on 40 days of JSE data

using LinearAlgebra, Plots, LaTeXStrings, StatsBase, Intervals, JLD2, ProgressMeter, Distributions, Optim, ArgCheck, FINUFFT, FFTW
#---------------------------------------------------------------------------

include("../Functions/Correlation Estimators/Dirichlet/NUFFTcorrDK-FGG.jl")
include("../Functions/Correlation Estimators/HY/HYcorr.jl")

#---------------------------------------------------------------------------
## RV correction
#---------------------------------------------------------------------------

function zeroticks(P)
    dif = diff(P,dims=1)
    n = size(dif)[1]
    m = size(dif)[2]
    p = zeros(m, 1)

    for i in 1:m
        r = dif[:,i]
        count = 0
        for j in 1:n
            if r[j] == 0
                count += 1
            end
        end
        p[i] = count/n
    end
    return p
end

function flattime(τ, t1, t2, T)
    t1 = [0; t1]
    t2 = [0; t2]
    syngrid = collect(0:τ:T)
    n = length(syngrid)
    γ1 = zeros(n,1)
    γ2 = zeros(n,1)
    for i in 1:n
        γ1[i] = maximum(filter(x-> x .<= syngrid[i], t1))
        γ2[i] = maximum(filter(x-> x .<= syngrid[i], t2))
    end
    ints = zeros(n-1,1)
    for i in 2:n
        a = γ1[i-1]..γ1[i]
        b = γ2[i-1]..γ2[i]
        c = intersect(a, b)
        ints[i-1] = c.last-c.first
    end
    return mean(ints) / (sqrt(mean(diff(γ1, dims = 1))*mean(diff(γ2, dims = 1))))
end

# Function to compute the correlations and corrections
function computecorrs(data, T = 28200, dt = collect(1:1:400))
    m = size(data)[1]
    measured = zeros(length(dt), m)
    prevtick = zeros(length(dt), m)
    overlap = zeros(length(dt), m)
    HY = zeros(m, 1)

    t1 = data[findall(!isnan, data[:,2]),1]
    t2 = data[findall(!isnan, data[:,3]),1]
    @showprogress "Computing..." for i in 1:length(dt)
        t = collect(0:dt[i]:T)
        n = length(t)
        p1 = zeros(n,1)
        p2 = zeros(n,1)
        for j in 1:n
            inds = findall(x -> x .<= t[j], data[:,1])
            p1[j] = filter(!isnan, data[inds,2])[end]
            p2[j] = filter(!isnan, data[inds,3])[end]
        end
        p = zeroticks([p1 p2])
        adj = flattime(dt[i], t1[2:end], t2[2:end], T)

        measured[i] = NUFFTcorrDKFGG([p1 p2], [t t])[1][1,2]
        prevtick[i] = measured[i] * ((1-p[1]*p[2]) / ((1-p[1])*(1-p[2])))
        overlap[i] = measured[i]/adj
    end

    P1 = filter(!isnan, data[:,2])
    P2 = filter(!isnan, data[:,3])
    HY = HYcorr(P1,P2,t1,t2)[1][1,2]

    return measured, prevtick, overlap, HY
end

# Streamline everything for computing correlations and corrections
function Empirical(Fulldata)
    res = computecorrs(Fulldata)
    return res
end

function generate_epps_plots_values(simulation_data)
    epps_value = [zeros(M, 1) for i = 1:num_paths]
    epps_mean = [zeros(M, 1) for i = 1:num_paths]
    total_epps_mean = zeros(size(epps_mean[1])[1])
    total_epps_value = zeros(size(epps_value[1])[1])
    m = 0

    for i in 1:num_paths
        path1 = simulation_data[i][1].raw_price_paths[1:s]
        path2 = simulation_data[i][2].raw_price_paths[1:s]
        index_vector = 0:1.0:(size(simulation_data[i][2].raw_price_paths[1:s])[1]-1)
        epps_data = hcat(index_vector, path1, path2)
        epps = Empirical(epps_data)
        m = size(epps_data)[1]
        # Save and Load
        # save("Computed Data/EppsCorrection/Empirical$i.jld", "epps$i", epps[i])
        # ComputedResults = load("Computed Data/EppsCorrection/Empirical$i.jld")
        # epps_result = ComputedResults["epps$i"]
        epps_value[i] = epps[1]
        epps_mean[i] = mean(epps[1], dims=2)
        total_epps_mean = total_epps_mean .+ epps_mean[i]
        total_epps_value = total_epps_value .+ epps_value[i]
    end

    average_epps_mean = total_epps_mean ./ num_paths
    average_epps_value = total_epps_value ./ num_paths

    return (average_epps_mean, average_epps_value, m)
end
