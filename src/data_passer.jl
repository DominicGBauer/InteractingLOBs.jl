# -*- coding: utf-8 -*-
mutable struct DataPasser
    #test::Int64
    slob::Any
    lob_densities::Array{Float64,2}
    sources::Array{Float64,2}
    couplings::Array{Float64,2}  
    rl_pushes::Array{Float64,2}
    raw_price_paths::Array{Float64,1}
    obs_price_paths::Array{Float64,1}
    P⁺s::Array{Float64,1}
    P⁻s::Array{Float64,1}
    Ps::Array{Float64,1} 
    V::Array{Float64,1}
end
