using InteractingLOBs, CSV, DataFrames
include("plot_stylized_facts.jl")


real_data = vec(Matrix(CSV.read("Data/Original_Price_Bars_2300.csv", DataFrame)))

data_stylized_facts = StylizedFactsPlot(real_data);

(acf, hist_qq, price_returns) = plot_all_stylized_facts(data_stylized_facts)

savefig(acf, "Plots/StylizedFacts/Empirical_ACF")
savefig(hist_qq, "Plots/StylizedFacts/Empirical_Hist_QQ")
savefig(price_returns, "Plots/StylizedFacts/Empirical_Price_Returns")
