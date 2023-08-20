# -*- coding: utf-8 -*-
mutable struct RLPushTerm
        StartTime::Int
        EndTime::Int
        Position::Int
        Amount::Float64
        Do::Bool
end

# +
function (rl::RLPushTerm)(D,slob_num, t)

    temp = [0.0 for xᵢ¹ in D[slob_num].slob.x]

    # return temp # return while its still 0 for no pushing

    if rl.Do&&(t>=rl.StartTime)&&(t<rl.EndTime) #never run with t=1
        if (rl.Position>0)
            temp[rl.Position] = -rl.Amount
        else
            latest_φ = D[slob_num].lob_densities[:,t-1]

            my_p = extract_mid_price_index(D[slob_num].slob,latest_φ)

            index_to_modify = my_p - rl.Position

            temp[index_to_modify] = temp[index_to_modify] - rl.Amount

        end
    end

    return temp

end
