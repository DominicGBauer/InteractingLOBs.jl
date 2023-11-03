using InteractingLOBs, LaTeXStrings, DataFrames, CSV, Plots, StatsBase, Distributions, StatsPlots, Plots.PlotMeasures
include("plot_stylized_facts.jl")
include("order_classification.jl")


real_data = vec(Matrix(CSV.read("Data/Original_Price_Bars_2300.csv", DataFrame)))

data_stylized_facts = StylizedFactsPlot(real_data);

p1 = plot_all_stylized_facts(data_stylized_facts,(1000,1200))

savefig(p1, "Plots/StylizedFacts/Empirical")
