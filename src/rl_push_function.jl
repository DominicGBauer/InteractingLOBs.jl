# -*- coding: utf-8 -*-
mutable struct RLPushTerm
        StartTime::Int
        EndTime::Int
        Position::Int
        Amount::Float64
        Do::Bool
end

# +
function (rl::RLPushTerm)(D,
                                 slob_num, t)
    
    #slob¹ = slobs[slob_num]
    #p_list¹ = raw_price_paths[slob_num,:,path_num]
    #φ_list¹ = lob_densities[slob_num,:,:,path_num]
    #slob² = slobs[2+(1-slob_num)]
    #p_list² = raw_price_paths[2+(1-slob_num),:,path_num]
    #φ_list² = lob_densities[2+(1-slob_num),:,:,path_num]
    
    temp = [0.0 for xᵢ¹ in D[slob_num].slob.x] 
    
    # return temp # return while its still 0 for no pushing 
    
    if rl.Do&&(t>=rl.StartTime)&&(t<rl.EndTime) #never run with t=1
        if (rl.Position>0)
            temp[rl.Position] = -rl.Amount
        else 
            latest_φ = D[slob_num].lob_densities[:,t-1]
            
            my_p = extract_mid_price_index(D[slob_num].slob,latest_φ)
            
            index_to_modify = my_p - rl.Position
            
            #left_side = [1:my_p;]
            
            # optional rescale of all points except index_to_modify (it is overwritten below) 
            # to have negative the area that index_to_modify has
            
            #area = sum(latest_φ[left_side])
            #rescaled_latest_φ = latest_φ[left_side] ./ area .* rl.Amount
            #temp[left_side] = rescaled_latest_φ
            
            temp[index_to_modify] = temp[index_to_modify] - rl.Amount
            
        end
    end
    
    
    
    return temp
    
end
