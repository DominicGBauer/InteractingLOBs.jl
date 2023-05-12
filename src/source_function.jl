# -*- coding: utf-8 -*-
mutable struct SourceTerm
    λ::Float64
    μ::Float64
    do_source::Bool
end


# +
function (st::SourceTerm)(D,
                                 slob_num, t)
    
    #slob¹ = slobs[slob_num]
    #p_list¹ = raw_price_paths[slob_num,:,path_num]
    #φ_list¹ = lob_densities[slob_num,:,:,path_num]
    #slob² = slobs[2+(1-slob_num)]
    #p_list² = raw_price_paths[2+(1-slob_num),:,path_num]
    #φ_list² = lob_densities[2+(1-slob_num),:,:,path_num]
    
    # extract most recent prices
    #p¹ = p_list¹[t-1]
    #p² = p_list²[t-1]
    
    if !(st.do_source)
        temp = [0 for _ in D[slob_num].slob.x]
        return temp
    end
    
    p¹ = D[slob_num].raw_price_paths[t-1]
    
    f(y)=-st.λ*st.μ*y*exp(-(st.μ*(y))^2) #y is a temporary variable
    
    return [f(xᵢ¹-p¹) for xᵢ¹ in D[slob_num].slob.x]
    
    #f(y)=-st.λ*(st.μ*(y))*exp(-(st.μ*(y))^2) #y is a temporary variable
    #return f(x-p¹)
    
end
