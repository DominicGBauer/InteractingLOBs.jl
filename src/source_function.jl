# -*- coding: utf-8 -*-
mutable struct SourceTerm
    λ::Float64
    μ::Float64
    do_source::Bool
end


# +
function (st::SourceTerm)(D, slob_num, t)
    if !(st.do_source)
        temp = [0 for _ in D[slob_num].slob.x]
        return temp
    end

    p¹ = D[slob_num].raw_price_paths[t-1]

    # Use this for Lit Order Book
    # f(y)=-st.λ*st.μ*y*exp(-(st.μ*(y))^2) #y is a temporary variable
    # Use this for Latent Order Book
    f(y)=st.λ*tanh(st.μ*(y))

    return [f(xᵢ¹-p¹) for xᵢ¹ in D[slob_num].slob.x]
end
