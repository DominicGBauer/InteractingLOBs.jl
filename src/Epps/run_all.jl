# include("./epps_plot_fixed_t_alternating_delta_x_1.jl")
# include("./epps_plot_fixed_t_alternating_delta_x_2.jl")
# include("./epps_plot_fixed_t_alternating_delta_x_3.jl")
# include("./epps_plot_fixed_t_alternating_delta_x_4.jl")
# include("./epps_plot_fixed_x_alternating_delta_t_1.jl")
# include("./epps_plot_fixed_x_alternating_delta_t_2.jl")
# include("./epps_plot_fixed_x_alternating_delta_t_3.jl")
# include("./epps_plot_fixed_x_alternating_delta_t_4.jl")

# ordering of plots is based on size of delta_x
final_plot_x = plot(p1_x, p2_x, p4_x, p3_x, layout = 4)
savefig(final_plot_x, "Plots/Epps/Epps_all_delta_x_changes.png")

# ordering of plots is based on size of delta_t
final_plot_t = plot(p1_t, p2_t, p3_t, p4_t, layout = 4)
savefig(final_plot_t, "Plots/Epps/Epps_all_delta_t_changes.png")
