using Plots, StatsBase, Distributions, StatsPlots, Plots.PlotMeasures, LaTeXStrings

include("../order_classification.jl")

mutable struct StylizedFactsPlot
    price_path::Array{Float64,1}
    log_price_path::Array{Float64,1}
    log_returns::Array{Float64,1}
    order_signs::Array{Float64,1}
    order_flow_acf::Array{Float64,1}
    log_returns_acf::Array{Float64,1}
    abs_log_returns_acf::Array{Float64,1}
    N::Int
    L::Int
    l_conf::Float64
    u_conf::Float64
end

# +
function StylizedFactsPlot(price_path)
    log_price_path = log.(price_path) #Apply log onto price path
    log_returns = diff(log_price_path) # returns log_returns[i+1] - log_returns[i] for all i in 1 to end-1

    return StylizedFactsPlot(price_path, log_price_path, log_returns)

end
# -

function StylizedFactsPlot(price_path, log_price_path, log_returns)
    order_signs = tick_rule(price_path)
    order_flow_acf = autocor(order_signs)

    log_returns_acf = autocor(log_returns)
    abs_log_returns_acf = autocor(abs.(log_returns))

    N = size(price_path, 1)
    L = size(log_returns_acf, 1)
    l_conf = -1.96 / sqrt(N)
    u_conf = 1.96 / sqrt(N)

    return StylizedFactsPlot(price_path, log_price_path, log_returns,
        order_signs, order_flow_acf, log_returns_acf, abs_log_returns_acf,
        N, L, l_conf, u_conf)
end

function plot_acf(sf::StylizedFactsPlot)
    # Log returns
    p1 = plot(1:sf.L, sf.log_returns_acf, xlab="Lag",
        ylab="ACF Log Returns",
        dpi=300,
        legend=false, seriestype=:sticks)
    plot!(p1, 1:sf.L, x -> 0, linestyle=:solid, color="black")
    plot!(p1, 1:sf.L, x -> sf.l_conf, linestyle=:dot)
    plot!(p1, 1:sf.L, x -> sf.u_conf, linestyle=:dot)

    # Absolute Log returns
    plot!(1:sf.L,
        sf.abs_log_returns_acf,
        ylab="ACF Absolute Log Returns",
        legend=false,
        dpi=300,
        seriestype=:sticks,
        inset=(1, bbox(0.0, 0.0, 0.8, 0.8, :top, :right)), subplot=2
    )
    plot!(subplot=2, 1:sf.L, x -> 0, linestyle=:solid, color="black")
    plot!(subplot=2, 1:sf.L, x -> sf.l_conf, linestyle=:dot)
    plot!(subplot=2, 1:sf.L, x -> sf.u_conf, linestyle=:dot)


    # Order Flow
    plot!(1:sf.L,
        sf.order_flow_acf,
        dpi=300,
        ylab="ACF Order Flow",
        legend=false,
        seriestype=:sticks,
        inset=(1, bbox(0.0, 0.0, 0.6, 0.6, :top, :right)), subplot=3
    )
    plot!(subplot=3, 1:sf.L, x -> 0, linestyle=:solid, color="black")
    plot!(subplot=3, 1:sf.L, x -> sf.l_conf, linestyle=:dot)
    plot!(subplot=3, 1:sf.L, x -> sf.u_conf, linestyle=:dot)

    return p1
end

function plot_price_and_returns(sf::StylizedFactsPlot)
    p1 = plot(1:size(sf.price_path, 1),
        sf.price_path,
        legend=false,
        xlab="Time",
        dpi=300,
        # ylim=([229, 232]),
        ylab="Mid-Price"
    )
    plot!(1:size(sf.log_returns, 1),
        sf.log_returns,
        legend=false,
        xlab="Time",
        dpi=300,
        xlim=([0, 10 * 10^3]),
        xticks=([0, 10 * 10^3]),
        ylab="Log Returns",
        ylim=([-0.0003, 0.0003]),
        yticks=([-0.0003, 0, 0.0003]),
        inset=(1, bbox(0.0, 0.0, 0.4, 0.4, :top, :right)),
        subplot=2
    )
    return p1
end

function plot_hist_and_qq(sf::StylizedFactsPlot)
    pq = histogram(sf.log_returns,
        normalize=:pdf,
        ylab="Density",
        xlab="Log Returns",
        dpi=300,
    )
    f(x) = pdf(StatsBase.fit(Normal, sf.log_returns), x)
    plot!(pq, f, label="Normal Distribution", dpi=300)


    qqplot!(StatsBase.fit(Normal, sf.log_returns),
        sf.log_returns,
        xlab="Theoretical Quantiles",
        ylab="Sample Quantiles",
        markersize=2,
        dpi=300,
        markerstrokewidth=0.05,
        xticks=([-0.0002, 0, 0.0002]),
        xlims=([-0.0002, 0.00025]),
        inset=(1, bbox(0.0, 0.0, 0.4, 0.4, :top, :right)),
        subplot=2
    )
    return pq
end

function plot_all_stylized_facts(sf::StylizedFactsPlot)
    acf = plot_acf(sf)
    hist_qq = plot_hist_and_qq(sf)
    price_returns = plot_price_and_returns(sf)

    return acf, hist_qq, price_returns
end
