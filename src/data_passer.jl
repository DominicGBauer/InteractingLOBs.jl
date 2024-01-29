# -*- coding: utf-8 -*-
mutable struct DataPasser
    slob::Any
    lob_densities::Matrix{Float64}
    lob_densities_L::Matrix{Float64}
    lob_densities_R::Matrix{Float64}
    sources::Matrix{Float64}
    couplings::Matrix{Float64}
    rl_pushes::Matrix{Float64}
    raw_price_paths::Vector{Float64}
    obs_price_paths::Vector{Float64}
    P⁺s::Vector{Float64}
    P⁻s::Vector{Float64}
    Ps::Vector{Float64}
    V::Vector{Float64}
    x_shifts::Vector{Float64}
end
