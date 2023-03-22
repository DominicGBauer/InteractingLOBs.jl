# -*- coding: utf-8 -*-
mutable struct CouplingTerm
    μ::Float64
    a::Float64
    b::Float64
    c::Float64
end

function (it::CouplingTerm)(x::Float64,p¹::Float64,p²::Float64,t::Int64)
    
    return 0 # for no coupling
    
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
