using Plots

function s(x, t, p, λ, μ)
    return λ * tanh(μ * (p(t) - x))
end

function plot_s(p, λ, μ)
    x = range(0, stop = 10, length = 100)  # Range of x values
    t = 0  # Time value
    y = s.(x, t, p, λ, μ)  # Evaluate s(1)(x, t) for each x

    plot(x, y, xlabel = "x", ylabel = "s(1)(x, t)", title = "Plot of s(1)(x, t)", legend = false)
end

function s2(x, t, p, λ, μ)
    return -λ * μ * (x - p(t)) * exp(μ * (x - p(t)))
end

function plot_s2(p, λ, μ)
    x = range(0, stop = 10, length = 100)  # Range of x values
    t = 0  # Time value
    y = s2.(x, t, p, λ, μ)  # Evaluate s(2)(x, t) for each x

    plot(x, y, xlabel = "x", ylabel = "s(2)(x, t)", title = "Plot of s(2)(x, t)", legend = false)
end

# Example usage
p(t) = 5 + t  # Define a function for p(t)
λ = 1
μ = 0.1  # Value for μ1

# plot_s(p,  λ, μ)
plot_s2(p, λ, μ)
