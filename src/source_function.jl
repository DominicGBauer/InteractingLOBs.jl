# -*- coding: utf-8 -*-
mutable struct SourceTerm
    λ::Float64
    μ::Float64
end


# +
function (st::SourceTerm)(x::Float64,p¹::Float64,p²::Float64,t::Int64)
    f(y)=-st.λ*(st.μ*(y))*exp(-(st.μ*(y))^2) #y is a temporary variable
    return f(x-p¹)
    
end