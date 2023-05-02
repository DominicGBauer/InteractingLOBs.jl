# -*- coding: utf-8 -*-
mutable struct RandomnessTerm
        σ::Float64  # Standard deviation in whatever distribution
        r::Float64 # Probability of self jump
        do_random::Bool
end

# +
function (rn::RandomnessTerm)(slob,V_list,t)
    if (rn.do_random)
        if (rand()<0.0)
            V_t = V_list[t-1]
        else
            ϵ = rand(Normal(0.0,rn.σ))
            V_t = sign(ϵ) * min(  abs(rn.σ * ϵ)  ,   slob.Δx / slob.Δt  ) # Ensures V_t¹ = ϵ¹σ ≤ Δx/Δt  ????
        end
        Z = (1/4) * (V_t * slob.Δx) / (slob.D)
        P = 1-rn.r
        P⁺ = (exp(Z))/(exp(Z) + exp(-Z))*rn.r
        P⁻ = 1 - P - P⁺
    else
        P = 1-rn.r
        P⁺ = (rn.r)/2
        P⁻ = P⁺¹
        V_t = 0
    end
    
    return P, P⁺, P⁻, V_t
    
end
