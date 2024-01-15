using Plots

dt = collect(1:1:400)

function generate_epps_plot(m, frequencies, power_spectrum, average_epps_mean, title=false)
  q = quantile.(TDist(m - 1), [0.975])

  p = plot(dt,
    average_epps_mean,
    legend=false,
    dpi=300,
    xlims=[0, Inf],
    xticks=([0, 100, 200, 300, 400]),
    xlabel=L"\Delta t\textrm{[sec]}",
    ylims=[0, 20 * 10^-5],
    yticks=([4 * 10^-5, 8 * 10^-5, 12 * 10^-5, 16 * 10^-5, 20 * 10^-5]),
    ylabel=L"\rho_{\Delta t}^{ij}",
    ribbon=(q .* std(average_epps_value, dims=2) * 0.001),
    fillalpha=0.15,
  )
  plot!(frequencies,
    power_spectrum,
    label=false,
    xlabel="Frequency (Hz)",
    xguidefontsize=6,
    xtickfontsize=6,
    ylabel="Power",
    yguidefontsize=6,
    ytickfontsize=6,
    xlims=[0, 3 * 10^-5],
    xticks=([3 * 10^-5]),
    ylims=[0, 1.5 * 10^-7],
    legend=false,
    inset=(1, bbox(0.15, 0, 0.3, 0.3, :top, :left)),
    subplot=2
  )

  return p
end

function generate_epps_plot_bottom_inset(m, frequencies, power_spectrum, average_epps_mean)
  q = quantile.(TDist(m - 1), [0.975])

  p = plot(dt,
    average_epps_mean,
    legend=false,
    dpi=300,
    xlims=[0, Inf],
    xticks=([0, 100, 200, 300, 400]),
    xlabel=L"\Delta t\textrm{[sec]}",
    ylims=[0, 8 * 10^-5],
    yticks=([4 * 10^-5, 8 * 10^-5]),
    ylabel=L"\rho_{\Delta t}^{ij}",
    ribbon=(q .* std(average_epps_value, dims=2) * 0.001),
    fillalpha=0.15,
  )
  plot!(frequencies,
    power_spectrum,
    label=false,
    xlabel="Frequency (Hz)",
    xguidefontsize=6,
    xtickfontsize=6,
    ylabel="Power",
    yguidefontsize=6,
    ytickfontsize=6,
    xlims=[0, 3 * 10^-5],
    xticks=([3 * 10^-5]),
    ylims=[0, 1.75 * 10^-7],
    legend=false,
    inset=(1, bbox(0.05, 0.1, 0.3, 0.3, :bottom, :right)),
    subplot=2
  )

  return p
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
    epps_value[i] = epps[1]
    epps_mean[i] = mean(epps[1], dims=2)
    total_epps_mean = total_epps_mean .+ epps_mean[i]
    total_epps_value = total_epps_value .+ epps_value[i]
  end

  average_epps_mean = total_epps_mean ./ num_paths
  average_epps_value = total_epps_value ./ num_paths

  return (average_epps_mean, average_epps_value, m)
end
