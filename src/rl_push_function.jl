# -*- coding: utf-8 -*-
mutable struct RLPushTerm
        StartTime::Int
        EndTime::Int
        Position::Int
        Amount::Float64
        Do::Bool
end

# +
function (rl::RLPushTerm)(slob¹, φ_list¹, p_list¹, 
                          slob², φ_list², p_list², 
                          t)
    
    #temp = φ_list¹[:,1].*0
    temp = [0 for xᵢ¹ in slob¹.x] 
    
    # return temp # return while its still 0 for no pushing 
    
    if rl.Do&&(t>=rl.StartTime)&&(t<rl.EndTime) #never run with t=1
        if (rl.Position>0)
            temp[rl.Position] = rl.Amount
        else 
            latest_φ = φ_list¹[:,t-1]
            
            my_p = extract_mid_price_index(slob¹,latest_φ)
            
            temp[my_p - rl.Position] = rl.Amount
        end
    end
    
    return temp
    
end
