## Author: Patrick Chang
# Script file to investigate the Epps effect in calendar time
# for 4 banking stocks on the JSE.

## Preamble

using CSV, Plots, Random, ProgressMeter, StatsBase, JLD, LaTeXStrings, DataFrames, Dates, Distributions

# Estimators
include("../Functions/Correlation Estimators/Dirichlet/NUFFTcorrDK-FGG.jl")

# Read in the data
prices = CSV.read("Data/JSE_prices_2019-06-24_2019-06-28.csv", DataFrame)
times = CSV.read("Data/JSE_times_2019-06-24_2019-06-28.csv", DataFrame)
volume = CSV.read("Data/JSE_volume_2019-06-24_2019-06-28.csv", DataFrame)

# Remove extra column
prices = prices[:, 2:end];
times = times[:, 2:end];
volume = volume[:, 2:end];
# Pull out banking stocks
prices = prices[:, [6:8; 10]];
times = times[:, [6:8; 10]];
volume = volume[:, [6:8; 10]];
# Get the names of tickers
tickers = names(prices)

#---------------------------------------------------------------------------
## Supporting functions to streamline the process
# Function to split the data into the 5 days
# Results in dictionary of data, with a vector (of the two assets) of matricies
# where the matrix is price, times, volume
function DataSplit(A1::String, A2::String, P::DataFrame, t::DataFrame, V::DataFrame)
  # Filter out the pair of interest
  p1 = filter(!isnan, P[:, A1])
  p2 = filter(!isnan, P[:, A2])
  t1 = filter(!isnan, t[:, A1])
  t2 = filter(!isnan, t[:, A2])
  V1 = filter(!isnan, V[:, A1])
  V2 = filter(!isnan, V[:, A2])
  # Convert the times to dates for index extraction later on
  A1dates = Date.(unix2datetime.(t1))
  A2dates = Date.(unix2datetime.(t2))
  dates_unique = unique(A1dates)
  # Initialize storage
  data = Dict()
  # Loop through each day
  for i in 1:length(dates_unique)
    # Extract the data
    date_indsA1 = findall(x -> x == dates_unique[i], A1dates)
    date_indsA2 = findall(x -> x == dates_unique[i], A2dates)
    # Data for the day
    day_data = []
    push!(day_data, [p1[date_indsA1] t1[date_indsA1] V1[date_indsA1]])
    push!(day_data, [p2[date_indsA2] t2[date_indsA2] V2[date_indsA2]])
    # Add to dictionary
    push!(data, i => day_data)
  end
  return data
end
# Function to make the calendar time data into MM and HY format
# pad the data and pull the prices and times together
function MakeMMHYData(data)
  # Pull out the data
  A1 = data[1]
  A2 = data[2]
  # Get the dimensions
  n1 = size(A1, 1)
  n2 = size(A2, 1)
  n = max(n1, n2)
  # Initialize price and time matrix
  P = fill(NaN, n, 2)
  t = fill(NaN, n, 2)
  # Populate the price and time matrix
  P[1:n1, 1] = A1[:, 1]
  P[1:n2, 2] = A2[:, 1]
  t[1:n1, 1] = A1[:, 2]
  t[1:n2, 2] = A2[:, 2]
  return P, t
end
# Function to get the calendar time correlations
function getCTcorrs(tickers; P=prices, t=times, V=volume)
  # Compute the number of pairwise comparisons
  npairs = Int(factorial(length(tickers)) / (factorial(2) * factorial(length(tickers) - 2)))
  # Time scale of investigation
  dt = collect(1:1:300)
  # Initialize the estimates
  MM = zeros(length(dt), npairs)
  # Set up ind for storage
  ind = 1
  # Loop through the pairs
  @showprogress "Computing..." for i in 1:(length(tickers)-1)
    for j in (i+1):length(tickers)
      # Split the data into separate days
      data = DataSplit(tickers[i], tickers[j], P, t, V)
      # Initialize temporary storage matricies for the estimates
      MMtemp = zeros(length(dt), 5)
      # Loop through the different days
      for k in 1:length(data)
        # Extract data for the day
        day_data = data[k]
        # Create the MM and HY dataset
        MMHYData = MakeMMHYData(day_data)
        # Loop through the different time scales
        for l in 1:length(dt)
          # Get N for MM
          N = Int(floor((28200 / dt[l] - 1) / 2))
          # Compute correlations
          MMtemp[l, k] = NUFFTcorrDKFGG(MMHYData[1], MMHYData[2], N=N)[1][1, 2]
        end
      end
      # Store the 5 day average into the matrix
      MM[:, ind] = mean(MMtemp, dims=2)
      # Update the ind
      ind += 1
    end
  end
  return MM
end


## Obtain the results
res = getCTcorrs(tickers)

# Save results
save("Computed Data/Epps_Emperical.jld", "res", res)
# Load results
res = load("Computed Data/Epps_Emperical.jld")
res = res["res"]

## Extract label pairs

function test()
  indexs = 1
  pairnames = Matrix{Union{Nothing,String}}(nothing, 1, 6)
  for i in 1:(length(tickers)-1)
    for j in (i+1):length(tickers)
      # push!(pairnames, "$(tickers[i])"*"/"*"$(tickers[j])")
      pairnames[indexs] = "$(tickers[i])" * "/" * "$(tickers[j])"
      indexs += 1
    end
  end

  return pairnames
end

pairnames = test()

## Plot the results
dt = collect(1:1:300)

plot(dt, res[:, 1], legend=:bottomright, label=pairnames[1], ylims=(0, 0.82), dpi=300)
xlabel!(L"\Delta t\textrm{[sec]}")
ylabel!(L"\rho_{\Delta t}^{ij}")
savefig("Plots/Epps/Epps_Emperical.png")
