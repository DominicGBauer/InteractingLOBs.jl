# -*- coding: utf-8 -*-
mutable struct SourceTerm
    λ::Float64
    μ::Float64
    do_source::Bool
end


# +
function (st::SourceTerm)(slob¹, φ_list¹, p_list¹, 
                          slob², φ_list², p_list², 
                          t)
    
    # extract most recent prices
    p¹ = p_list¹[t-1]
    #p² = p_list²[t-1]
    if !(st.do_source)
        temp = [0 for _ in slob¹.x]
        return temp
    end
    
    f(y)=-st.λ*(st.μ*(y))*exp(-(st.μ*(y))^2) #y is a temporary variable
    
    return [f(xᵢ¹-p¹) for xᵢ¹ in slob¹.x]
    
    
    #f(y)=-st.λ*(st.μ*(y))*exp(-(st.μ*(y))^2) #y is a temporary variable
    #return f(x-p¹)
    
end
