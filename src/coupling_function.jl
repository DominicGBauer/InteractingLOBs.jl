# -*- coding: utf-8 -*-
mutable struct CouplingTerm
    μ::Float64
    a::Float64
    b::Float64
    c::Float64
end

# +
function (it::CouplingTerm)(slob¹, φ_list¹, p_list¹, 
                            slob², φ_list², p_list², 
                            t)
    
    #return [0 for xᵢ¹ in slob¹.x] #for zero coupling
    
    
    # extract most recent prices
    p¹ = p_list¹[t-1]
    p² = p_list²[t-1]
    
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
    
    coupling¹ = [coupling_inner(xᵢ¹, p¹, p², t) for xᵢ¹ in slob¹.x]
    
end
