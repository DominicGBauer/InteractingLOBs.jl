# -*- coding: utf-8 -*-
mutable struct CouplingTerm
    μ::Float64
    a::Float64
    b::Float64
    c::Float64
    do_coupling::Bool
end

function CouplingTerm(μ::Float64, a::Float64, b::Float64, c::Float64)
    return CouplingTerm(μ, a, b, c, true)
end

# +
function (it::CouplingTerm)(D,slob_num, t)

    #slob¹ = slobs[slob_num]
    #p_list¹ = raw_price_paths[slob_num,:,path_num]
    #φ_list¹ = lob_densities[slob_num,:,:,path_num]
    #slob² = slobs[2+(1-slob_num)]
    #p_list² = raw_price_paths[2+(1-slob_num),:,path_num]
    #φ_list² = lob_densities[2+(1-slob_num),:,:,path_num]

    if (!(it.do_coupling))
        return [0 for xᵢ¹ in D[slob_num].slob.x] #for zero coupling
    end


    # extract most recent prices
    #p¹ = p_list¹[t-1]
    #p² = p_list²[t-1]

    p¹ = D[slob_num].raw_price_paths[t-1]
    p² = D[2+(1-slob_num)].raw_price_paths[t-1]

    function coupling_inner(x, p¹, p², t)
        f(y)=-(it.μ*(y))*exp(-(it.μ*(y))^2) #y is a temporary variable
        g = 1+tanh(abs(p²-p¹)/it.a)*it.c

        if (p²>p¹)
            if x>p¹
                return it.b* f(1/g*(x-p¹))
            else
                return it.b* g*f(x-p¹)
            end
        else
            if x<p¹
                return it.b* f(1/g*(x-p¹))
            else
                return it.b* g*f(x-p¹)
            end
        end
    end

    coupling¹ = [coupling_inner(xᵢ¹, p¹, p², t) for xᵢ¹ in D[slob_num].slob.x]

end
# -
