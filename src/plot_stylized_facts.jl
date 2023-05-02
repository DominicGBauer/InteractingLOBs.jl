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

    return StylizedFactsPlot(price_path,log_price_path, log_returns) 
    
end
# -

function StylizedFactsPlot(price_path,log_price_path,log_returns)
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

function hist_log_returns(sf::StylizedFactsPlot, title="Log Returns Histogram")
    pq = histogram(sf.log_returns, normalize=:pdf, label="observed returns",
        ylab="Density", xlab="Log Returns", title=title);
    f(x) = pdf(fit(Normal, sf.log_returns), x)
    plot!(pq, f, label="Normal Distribution");
    return pq
end

function plot_log_returns(sf::StylizedFactsPlot, title="Log Returns")
    p1 = plot(1:size(sf.log_returns, 1), sf.log_returns, legend=false,
    xlab="Time", ylab="Log Returns", title=title);
    return p1
end

function plot_qq_log_returns(sf::StylizedFactsPlot, title="Log Returns Normal Q-Q Plot")
    p1 = qqplot(fit(Normal, sf.log_returns), sf.log_returns,
    xlab="Theoretical Quantiles", ylab="Sample Quantiles", title=title, markersize=2,markerstrokewidth=0.05);
    return p1
end

function plot_acf_log_returns(sf::StylizedFactsPlot, title="Log Returns Autocorrelation")

    p1 = plot(1:sf.L, sf.log_returns_acf, xlab="Lag",
        ylab="ACF", title=title,
        legend=false, seriestype=:sticks);
    plot!(p1, 1:sf.L, x -> 0, linestyle=:solid, color="black")
    plot!(p1, 1:sf.L, x -> sf.l_conf, linestyle=:dot)
    plot!(p1, 1:sf.L, x -> sf.u_conf, linestyle=:dot)

    return p1
end





function plot_acf_order_flow(sf::StylizedFactsPlot, title="Order Flow Autocorrelation (Tick Rule)")

    p1 = plot(1:sf.L, sf.order_flow_acf, xlab="Lag",
        ylab="ACF", title=title,
        legend=false, seriestype=:sticks);
    plot!(p1, 1:sf.L, x -> 0, linestyle=:solid, color="black")
    plot!(p1, 1:sf.L, x -> sf.l_conf, linestyle=:dot)
    plot!(p1, 1:sf.L, x -> sf.u_conf, linestyle=:dot)
    return p1
end

function plot_acf_abs_log_returns(sf::StylizedFactsPlot, title="Absolute Log Returns Autocorrelation")

    p1 = plot(1:sf.L, sf.abs_log_returns_acf, xlab="Lag",
        ylab="ACF", title=title,
        legend=false, seriestype=:sticks);
    plot!(p1, 1:sf.L, x -> 0, linestyle=:solid, color="black")
    plot!(p1, 1:sf.L, x -> sf.l_conf, linestyle=:dot)
    plot!(p1, 1:sf.L, x -> sf.u_conf, linestyle=:dot)
    return p1
end

function plot_all_stylized_facts(sf::StylizedFactsPlot, plot_size=(1000, 800), titles_off=false)

    l = @layout [a b ; c d ; e f; g h]

    if titles_off
        p1 = plot(1:size(sf.price_path, 1), sf.price_path, legend=false,
        xlab="Time", ylab="Mid-Price");
        p2 = plot_log_returns(sf, title="")
        p3 = plot_qq_log_returns(sf, title="")
        p4 = plot_acf_order_flow(sf, title="")
        p5 = plot_acf_log_returns(sf, title="")
        p6 = plot_acf_abs_log_returns(sf, title="")
        p7 = hist_log_returns(sf, title="")
        p8 = plot(p1, p7, p2, p3, p4, p5, p6, layout=l, tickfontsize=6, guidefontsize=8,
            titlefontsize=10, right_margin=5mm, size=plot_size);
        return p8
    else
        p1 = plot(1:size(sf.price_path, 1), sf.price_path, legend=false,
        xlab="Time", ylab="Mid-Price", title="Mid-Price Path");
        p2 = plot_log_returns(sf)
        p3 = plot_qq_log_returns(sf)
        p4 = plot_acf_order_flow(sf)
        p5 = plot_acf_log_returns(sf)
        p6 = plot_acf_abs_log_returns(sf)
        p7 = hist_log_returns(sf)
        p8 = plot(p1, p7, p2, p3, p4, p5, p6, layout=l, tickfontsize=6, guidefontsize=8,
            titlefontsize=10, right_margin=5mm, size=plot_size);
        return p8
    end
end

