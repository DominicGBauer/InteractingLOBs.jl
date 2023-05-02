# -*- coding: utf-8 -*-
# +
function initial_conditions_numerical(D, 
                                        t, slob_num, V₀=0)
    slob¹ = D[slob_num].slob
    slob² = D[2+(1-slob_num)].slob
    
    if (t == 1)
        t = t + 1
    end# ?????????????????????????????????????????????????/
    
    # eq 36:
    ud = (-V₀/(2.0*slob¹.Δx) + slob¹.D/(slob¹.Δx^2)) * ones(Float64, slob¹.M)#upper diagonal
    md = ((-2.0*slob¹.D)/(slob¹.Δx^2) - slob¹.nu) * ones(Float64, slob¹.M+1) #middle diagonal
    ld = (V₀/(2.0*slob¹.Δx) + slob¹.D/(slob¹.Δx^2)) * ones(Float64, slob¹.M) #lower diagonal
    A = Tridiagonal(ld, md, ud)

    A[1,2] = 2*slob¹.D/(slob¹.Δx^2)
    A[end, end-1] = 2*slob¹.D/(slob¹.Δx^2)
    
    source   = slob¹.source_term(   D,
                                    slob_num, t) 

    coupling = slob¹.coupling_term( D,
                                    slob_num, t) 

    rl_push  = slob¹.rl_push_term(  D,
                                    slob_num, t) 
    
    B = .-(source.+coupling.+rl_push)#B is s in (37)
    
    #temp
    if slob¹.source_term.do_source
        φ = A \ B
    else 
        φ = source .* 0
    end
    
    D[slob_num].lob_densities[:,t] = φ
    D[slob_num].sources[:,t]       = source
    D[slob_num].couplings[:,t]     = coupling
    D[slob_num].rl_pushes[:,t]     = rl_push
    D[slob_num].V[t]               = V₀
    
    #give x such that A*x = B
    #return φ, source, coupling, rl_push, V₀
    
end
# -

function extract_mid_price_index(slob,lob_density)
    mid_price_ind = 2
    l = length(lob_density)
    #while (mid_price_ind<l)&&((lob_density[mid_price_ind] > 0) || (lob_density[mid_price_ind+1]>lob_density[mid_price_ind]))
    while ((mid_price_ind<l)&&(lob_density[mid_price_ind] > 0))
        mid_price_ind += 1
    end #scan array in x looking for cross over point of the mid price
    if ( mid_price_ind==l || mid_price_ind==2 ) 
        return -1
    end
    return mid_price_ind
end

function extract_mid_price(slob, lob_density)
    mid_price_ind = extract_mid_price_index(slob,lob_density)
    if mid_price_ind == -1
        return -1
    end
    
    y1 = lob_density[mid_price_ind-1]#y value (density) to the left
    y2 = lob_density[mid_price_ind]#y value (density) to the right
    x1 = slob.x[mid_price_ind-1]#x value value to the left

    mid_price = (-y1 * slob.Δx)/(y2 - y1) + x1 #solution to assuming straight approximation between 
                                               #left and right point (did the math) (see page 20)
    return mid_price
end


function SibuyaKernelDP(n,slob)
    return slob.SK_DP[n]
end

# +
#function intra_time_period_simulate_fractional(slobs, φ_list, p_list, V_list,
#                                                 t) 
function intra_time_period_simulate_fractional( D,
                                                t, slob_num) 
    slob = D[slob_num].slob
    myend = size(D[slob_num].lob_densities[:,t-1],1) #end and myend are now the same in all the below
    middle = 2:(myend-1) #all indices excluding the first and last

    P, P⁺, P⁻, V_t  = D[slob_num].slob.randomness_term(D[slob_num].slob, D[slob_num].V, t) 

    source   = slob.source_term(    D,
                                    slob_num, t) 

    coupling = slob.coupling_term(  D,
                                    slob_num, t) 

    rl_push  = slob.rl_push_term(   D,
                                    slob_num, t) 

    net_source = source .+ coupling .+ rl_push

    φ_next = zeros(Float64,myend)

    for t_m in max(1,t-D[slob_num].slob.cut_off):(t-1)

        front = SibuyaKernelDP(t-t_m,D[slob_num].slob)

        φ_next[1] += front * 
            (  P⁺ * D[slob_num].lob_densities[1        ,t_m] + 
               P⁻ * D[slob_num].lob_densities[2        ,t_m] + 
               P  * D[slob_num].lob_densities[1        ,t_m]  ) 

        φ_next[end] += front * 
        (  P⁻ * D[slob_num].lob_densities[end          ,t_m] + 
               P⁺ * D[slob_num].lob_densities[end-1    ,t_m] + 
               P  * D[slob_num].lob_densities[end      ,t_m]  )

        φ_next[middle] +=  front * 
            (  P⁺ * D[slob_num].lob_densities[middle.-1,t_m] + 
               P⁻ * D[slob_num].lob_densities[middle.+1,t_m] + 
               P  * D[slob_num].lob_densities[middle   ,t_m]  )   
    end

    φ_next +=  exp( - D[slob_num].slob.nu * slob.Δt) * D[slob_num].lob_densities[:,t-1] + slob.Δt * net_source

    # store all those intermediate values
    D[slob_num].lob_densities[:,t] = φ_next
    D[slob_num].sources[:,t]       = source
    D[slob_num].couplings[:,t]     = coupling
    D[slob_num].rl_pushes[:,t]     = rl_push
    D[slob_num].Ps[t]              = P
    D[slob_num].P⁺s[t]             = P⁺
    D[slob_num].P⁻s[t]             = P⁻
    D[slob_num].V[t]               = V_t

end


# +
function dtrw_solver_fractional(D) #handles change over time logic
    broke_point = -1
    slob = D[1].slob
    
    time_steps = to_simulation_time(slob.T, slob.Δt) 
    
    for slob_num in 1:length(D)
         D[slob_num].raw_price_paths[1] = D[slob_num].slob.p₀
         D[slob_num].obs_price_paths[1] = D[slob_num].slob.p₀
    end
    
    t = 1 #initial conditions take current t to read values
    
    for slob_num in 1:length(D)
        initial_conditions_numerical( D, t, slob_num) 
    end
    
    t = 2
    while (t <= time_steps) && (!(D[1].raw_price_paths[t-1]==-1) || (D[2].raw_price_paths[t-1]==-1)) 
        for slob_num in 1:length(D)
            intra_time_period_simulate_fractional( D, t, slob_num) 
            
            D[slob_num].raw_price_paths[t] = 
                                    extract_mid_price(D[slob_num].slob, D[slob_num].lob_densities[:,t])
        end
        
        t += 1 
    end
    
    #check if broken
    if (D[1].raw_price_paths[t-1]==-1) || (D[2].raw_price_paths[t-1]==-1)
        print("BE")
        broke_point = t-1
    end
    
    for slob_num in 1:length(D)
        D[slob_num].obs_price_paths = 
                    sample_mid_price_path(D[slob_num].slob, D[slob_num].raw_price_paths) 
    end
    
end
# -

# # New Junk

# +
# the following four functions are all equations (24) to (26)

# +
# function calculate_left_jump_probability(Z)
#     return (exp(-Z))/(exp(Z) + exp(-Z) + 1.0)
#     #return (exp(-Z))/(exp(Z) + exp(-Z))
# end

# +
# function calculate_right_jump_probability(Z)
#     return (exp(Z))/(exp(Z) + exp(-Z) + 1.0)
# end

# +
# function calculate_self_jump_probability(Z)
#     return (1.0)/(exp(Z) + exp(-Z) + 1.0)
# end

# +
# function calculate_jump_probabilities(slob, V_t)
#     Z = (3/4) * (V_t * slob.Δx) / (slob.D) #  equation 23 with modification for insertion into 24, 25 and 26
#                                             #  ≈ V_t * 3/8
#     p⁻ = calculate_left_jump_probability(Z)
#     p⁺ = calculate_right_jump_probability(Z)
#     p = calculate_self_jump_probability(Z)
    
    
#     #Z = (V_t * slob.Δx) / (2*slob.D) 
    
#     #p⁻ = calculate_left_jump_probability(Z)
#     #p⁺ = 1-p⁻
#     #p = 0
#     #p = calculate_self_jump_probability(Z)
#     return p⁺, p⁻, p
#
# -


#function initial_conditions_numerical(slob::SLOB, p_0¹, p_0², t)
    #ϵ = rand(Normal(0.0, 1.0))
    #V₀ = sign(ϵ) * min(abs(slob.σ * ϵ), slob.Δx / slob.Δt)
    #return initial_conditions_numerical(slob, p_0¹, p_0², V₀)
#end

# +
#function fractional_choose(n,k)
#    if (n==0)
#        return 0
#    end
#    #return gamma(n+1)/(gamma(n-k+1)*gamma(k+1)) 
#    return 1/(n+1) * 1/beta(k+1,n-k+1)
#end

# +
#function SibuyaKernelModified(n,γ,nu,Δt)
#    return (fractional_choose( 1-γ , n ) * (-1)^n + (n==1 ? 1 : 0)) * exp(-n*nu*Δt)
#end

# +
#function get_sub_period_time(slob, t, time_steps)
#    if slob.α <= 0.0 #if using LOB
#        return 0.0, time_steps - t + 1
#    end
#    τ = rand(Exponential(slob.α)) #random exponential time until recalculation
#    remaining_time = time_steps - t + 1
#    τ_periods = min(floor(Int, τ/slob.Δt), remaining_time)
#    return τ, τ_periods #return true time and number of simulation steps
#                        #until next solve of the initial conditions when using SLOB
#end


# +
# function dtrw_solver_no_recalc(slob¹::SLOB, slob²::SLOB)
#     time_steps = to_simulation_time(slob¹.T, slob¹.Δt) 
    
#     φ¹ = ones(Float64, slob¹.M + 1, time_steps + 1) #stores φ¹ for all time
#     φ² = ones(Float64, slob².M + 1, time_steps + 1) #stores φ² for all time
    
#     source¹ = ones(Float64, slob¹.M + 1, time_steps + 1) #stores sources for all time
#     source² = ones(Float64, slob².M + 1, time_steps + 1) #stores sources for all time
    
#     coupling¹ = ones(Float64, slob¹.M + 1, time_steps + 1) #stores couplings for all time
#     coupling² = ones(Float64, slob².M + 1, time_steps + 1) #stores couplings for all time
    
#     rl_push¹ = ones(Float64, slob¹.M + 1, time_steps + 1) #stores rl_pushes for all time
#     rl_push² = ones(Float64, slob².M + 1, time_steps + 1) #stores rl_pushes for all time

#     p¹ = ones(Float64, time_steps + 1)
#     p² = ones(Float64, time_steps + 1)
    
#     mid_prices¹ = ones(Float64, slob¹.T + 1) #stores mid prices for all time
#     mid_prices² = ones(Float64, slob².T + 1) #stores mid prices for all time

#     p¹[1] = slob¹.p₀
#     mid_prices¹[1] = slob¹.p₀
    
#     p²[1] = slob².p₀
#     mid_prices²[1] = slob².p₀

#     P⁺s¹ = fill(1/3, time_steps) #stores right jump prob for all time
#     P⁻s¹ = fill(1/3, time_steps) #store left jump prob for all time
#     Ps¹ = fill(1/3, time_steps) #store self jump prob for all time
#     V¹ = fill(0.0, time_steps+1)
    
#     P⁺s² = fill(1/3, time_steps) #stores right jump prob for all time
#     P⁻s² = fill(1/3, time_steps) #store left jump prob for all time
#     Ps² = fill(1/3, time_steps) #store self jump prob for all time
#     V² = fill(0.0, time_steps+1)
    
    
#     t = 1 #initial conditions take current t to read values
#     φ¹[:, t], source¹[:,t], coupling¹[:,t], rl_push¹[:,t], V¹[t] = initial_conditions_numerical(slob¹, φ¹, p¹,
#                                                                           slob², φ², p² ,
#                                                                           t, 0.0) #get initial shape given at t=1
    
#     φ²[:, t], source²[:,t], coupling²[:,t], rl_push²[:,t], V²[t] = initial_conditions_numerical(slob², φ², p² ,
#                                                                           slob¹, φ¹, p¹,
#                                                                           t, 0.0) #get initial shape given at t=1
#     t = 2
#     while t <= time_steps
#         φ¹[:, t], source¹[:,t], coupling¹[:,t], rl_push¹[:,t], P⁺s¹[t-1], P⁻s¹[t-1], Ps¹[t-1], V¹[t],  
#         φ²[:, t], source²[:,t], coupling²[:,t], rl_push²[:,t], P⁺s²[t-1], P⁻s²[t-1], Ps²[t-1], V²[t] =
#             intra_time_period_simulate(     slob¹, φ¹, p¹, V¹,
#                                             slob², φ², p², V²,
#                                             t                       )    #get φ at next time step
                                                        
#         p¹[t] = extract_mid_price(slob¹, φ¹[:, t])#get mid price of new φ. Raw mid prices
#         p²[t] = extract_mid_price(slob², φ²[:, t])#get mid price of new φ. Raw mid prices
        
#         t += 1 
        
#     end
    
#     sample_mid_prices¹ = sample_mid_price_path(slob¹, p¹) #calculate observed mid prices
#     sample_mid_prices² = sample_mid_price_path(slob², p²) #calculate observed mid prices
    
#     return   φ¹, source¹, coupling¹, rl_push¹, p¹, sample_mid_prices¹, P⁺s¹, P⁻s¹, Ps¹, V¹,
# #             φ², source², coupling², rl_push², p², sample_mid_prices², P⁺s², P⁻s², Ps², V²
# ##end
# -



# +
# function intra_time_period_simulate(slob¹, φ_list¹, p_list¹, V_list¹, slob², φ_list², p_list², V_list²,
#                                     t) 
#     #p_listⁱ is a list of all prices until the most recent
#     #φⁱ is the most recent order book only
    
#     ### corrections
#     φ¹ = φ_list¹[:, t-1]
#     φ² = φ_list²[:, t-1]
    
#     #p¹ = p_list¹[t-1]
#     #p² = p_list²[t-1]
#     ## end of corrections
    
    
    
#     ϵ¹ = rand(slob¹.dist)
#     ϵ² = rand(slob².dist)
    
#     V_t¹ = sign(ϵ¹) * min(  abs(slob¹.σ * ϵ¹)  ,   slob¹.Δx / slob¹.Δt  ) # Ensures V_t¹ = ϵ¹σ ≤ Δx/Δt  ????
#     V_t² = sign(ϵ²) * min(  abs(slob².σ * ϵ²)  ,   slob².Δx / slob².Δt  ) # Ensures V_t² = ϵ²σ ≤ Δx/Δt  ????
    
#     ##try self correlation
#     #if (rand()<0.3)
#     #   prev_V¹ = V_list¹[t-1]
#     #   V_t¹ = abs(V_t¹) * sign(prev_V¹)  
#     #end
#     #if (rand()<0.3)
#     #   prev_V² = V_list²[t-1]
#     #   V_t² = abs(V_t²) * sign(prev_V²)  
#     #end
#     ##
    

#     P⁺¹, P⁻¹, P¹ = calculate_jump_probabilities(slob¹, V_t¹)
#     P⁺², P⁻², P² = calculate_jump_probabilities(slob², V_t²)
    
#     myend = size(φ¹,1) #end and myend are now the same in all the below
#     middle = 2:(myend-1) #all indices excluding the first and last

#     φ₋₁¹ = φ¹[1]   #use left most value as extra ghost point
#     φ₊₁¹ = φ¹[end] #use right most value as extra ghost point
#     φ_next¹ = zeros(Float64, myend)
    
#     source¹ = zeros(Float64, myend)
#     coupling¹ = zeros(Float64, myend)
#     rl_push¹ = zeros(Float64, myend)
    
#     φ₋₁² = φ²[1]   #use left most value as extra ghost point
#     φ₊₁² = φ²[end] #use right most value as extra ghost point
#     φ_next² = zeros(Float64, myend)
    
#     source² = zeros(Float64, myend)
#     coupling² = zeros(Float64, myend)
#     rl_push² = zeros(Float64, myend)
    
#     #first calculate source terms:
#     source¹   = slob¹.source_term(    slob¹, φ_list¹, p_list¹, 
#                                       slob², φ_list², p_list², 
#                                       t) 
#     coupling¹ = slob¹.coupling_term(  slob¹, φ_list¹, p_list¹, 
#                                       slob², φ_list², p_list², 
#                                       t) 
#     rl_push¹  = slob¹.rl_push_term(   slob¹, φ_list¹, p_list¹, 
#                                       slob², φ_list², p_list², 
#                                       t) 
    
#     source²   = slob².source_term(    slob², φ_list², p_list², 
#                                       slob¹, φ_list¹, p_list¹, 
#                                       t)
#     coupling² = slob².coupling_term(  slob¹, φ_list¹, p_list¹, 
#                                       slob², φ_list², p_list², 
#                                       t) 
#     rl_push²  = slob².rl_push_term(   slob¹, φ_list¹, p_list¹, 
#                                       slob², φ_list², p_list², 
#                                       t)
    
    
#     net_source¹ = source¹ .+ coupling¹ .+ rl_push¹
#     net_source² = source² .+ coupling² .+ rl_push²
    
#     # the below 3 equations implements equation 9 
    
#     # order book 1
#     φ_next¹[1] = P⁺¹ * φ₋₁¹ + P⁻¹ * φ¹[2] + P¹ * φ¹[1] - 
#         slob¹.nu * slob¹.Δt * φ¹[1] + 
#         slob¹.Δt * net_source¹[1]

#     φ_next¹[end] = P⁻¹ * φ₊₁¹ + P⁺¹ * φ¹[end-1] + P¹ * φ¹[end] -
#         slob¹.nu * slob¹.Δt * φ¹[end] + 
#         slob¹.Δt * net_source¹[end]
     
#     φ_next¹[middle] = P⁺¹ * φ¹[middle.-1] + P⁻¹ * φ¹[middle.+1] + P¹ * φ¹[middle] -
#         slob¹.nu * slob¹.Δt * φ¹[middle] + 
#         slob¹.Δt * net_source¹[middle]
    
#     # order book 2
#     φ_next²[1] = P⁺² * φ₋₁² + P⁻² * φ²[2] + P² * φ²[1] - 
#         slob².nu * slob².Δt * φ²[1] + 
#         slob².Δt * net_source²[1]

#     φ_next²[end] = P⁻² * φ₊₁² + P⁺² * φ²[end-1] + P² * φ²[end] -

#         slob².nu * slob².Δt * φ²[end] + 
#         slob².Δt * net_source²[end]
   
#     φ_next²[middle] = P⁺² * φ²[middle.-1] + P⁻² * φ²[middle.+1] + P² * φ²[middle] -
#         slob².nu * slob².Δt * φ²[middle] + 
#         slob².Δt * net_source²[middle]

#     return   φ_next¹, source¹, coupling¹, rl_push¹, P⁺¹, P⁻¹, P¹, V_t¹,
#              φ_next², source², coupling², rl_push², P⁺², P⁻², P², V_t²
# end


# +
    # p_listⁱ is a list of all prices until the most recent
    # φⁱ is the most recent order book only
    
    
    ####for when force function is included#######
#     if (slob¹.do_random)
#         ϵ¹ = rand(slob¹.dist)
#         V_t¹ = sign(ϵ¹) * min(  abs(slob¹.σ * ϵ¹)  ,   slob¹.Δx / slob¹.Δt  ) # Ensures V_t¹ = ϵ¹σ ≤ Δx/Δt  ????
#         Z¹ = (1/4) * (V_t¹ * slob¹.Δx) / (slob¹.D)
#         P¹ = 1-slob¹.r
#         P⁺¹ = (exp(Z¹))/(exp(Z¹) + exp(-Z¹))*slob¹.r
#         P⁻¹ = 1 - P¹ - P⁺¹
#     else
#         P¹ = 1-slob¹.r
#         P⁺¹ = (slob¹.r)/2
#         P⁻¹ = P⁺¹
#         V_t¹ = 0
#     end
    
#     if (slob².do_random)
#         ϵ² = rand(slob².dist)
#         V_t² = sign(ϵ²) * min(  abs(slob².σ * ϵ²)  ,   slob².Δx / slob².Δt  ) # Ensures V_t² = ϵ²σ ≤ Δx/Δt  ????
#         Z² = (1/4) * (V_t² * slob².Δx) / (slob².D)
#         P² = 1-slob².r
#         P⁺² = (exp(Z²))/(exp(Z²) + exp(-Z²))*slob².r
#         P⁻² = 1 - P² - P⁺²
#     else
#         P² = 1-slob².r
#         P⁺² = (slob².r)/2
#         P⁻² = P⁺²
#         V_t² = 0
#     end
# -

# # Old junk

# +
# function dtrw_solver_with_recalc(slob¹::SLOB, slob²::SLOB,kick::Kicker)
#     time_steps = get_time_steps(slob¹.T, slob¹.Δt) #let's assume T,L,M,D,num_paths,α is the same for both 
    
#     φ¹ = ones(Float64, slob¹.M + 1, time_steps + 1) #stores φ¹ for all time
#     φ² = ones(Float64, slob².M + 1, time_steps + 1) #stores φ² for all time
    
#     source¹ = ones(Float64, slob¹.M + 1, time_steps + 1) #stores sources for all time
#     source² = ones(Float64, slob².M + 1, time_steps + 1) #stores sources for all time
    
#     coupling¹ = ones(Float64, slob¹.M + 1, time_steps + 1) #stores couplings for all time
#     coupling² = ones(Float64, slob².M + 1, time_steps + 1) #stores couplings for all time

#     p¹ = ones(Float64, time_steps + 1)
#     p² = ones(Float64, time_steps + 1)
    
#     mid_prices¹ = ones(Float64, slob¹.T + 1) #stores mid prices for all time
#     mid_prices² = ones(Float64, slob².T + 1) #stores mid prices for all time

#     p¹[1] = slob¹.p₀
#     mid_prices¹[1] = slob¹.p₀
    
#     p²[1] = slob².p₀
#     mid_prices²[1] = slob².p₀

#     P⁺s¹ = fill(1/3, time_steps) #stores right jump prob for all time
#     P⁻s¹ = fill(1/3, time_steps) #store left jump prob for all time
#     Ps¹ = fill(1/3, time_steps) #store self jump prob for all time
    
#     P⁺s² = fill(1/3, time_steps) #stores right jump prob for all time
#     P⁻s² = fill(1/3, time_steps) #store left jump prob for all time
#     Ps² = fill(1/3, time_steps) #store self jump prob for all time
    
#     t = 1
    
#     φ¹[:, t], source¹[:,t], coupling¹[:,t] = initial_conditions_numerical(slob¹, p¹[t], p²[t], 0.0) #get initial shape given at t=1
#                                                             #the source term and mid price position
#     φ²[:, t], source²[:,t], coupling²[:,t] = initial_conditions_numerical(slob², p²[t], p¹[t], 0.0) #get initial shape given at t=1
#                                                             #the source term and mid price position

#     t = 2
    
#     τ, τ_periods = get_sub_period_time(slob¹, t, time_steps) #get time periods till next recalc
    
#     time_since = 0
    
#     while t <= time_steps
        
#         kick.DoNow =  kick.Do&&(t>=kick.StartTime)&&(t<kick.EndTime) 
        
#         if (time_since <= τ_periods)
#             φ¹[:, t], source¹[:,t], coupling¹[:,t], P⁺s¹[t-1], P⁻s¹[t-1], Ps¹[t-1],  
#             φ²[:, t], source²[:,t], coupling²[:,t], P⁺s²[t-1], P⁻s²[t-1], Ps²[t-1] =
#             intra_time_period_simulate(     slob¹, φ¹[:, t-1], p¹[t-1],
#                                             slob², φ²[:, t-1], p²[t-1], 
#                                             kick                        )#get φ at next time step
#             time_since += 1
#         else 
#             φ¹[:, t], source¹[:,t], coupling¹[:,t] = initial_conditions_numerical(slob¹, p¹[t-1], p²[t-1])
#             φ²[:, t], source²[:,t], coupling²[:,t] = initial_conditions_numerical(slob², p²[t-1], p¹[t-1])
#             τ, τ_periods = get_sub_period_time(slob¹, t, time_steps) #get time periods till next recalc
#             time_since = 0
#         end

#         t += 1 
       
#         p¹[t] = extract_mid_price(slob¹, φ¹[:, t])#get mid price of new φ. Raw mid prices
#         p²[t] = extract_mid_price(slob², φ²[:, t])#get mid price of new φ. Raw mid prices

        
#     end
    
    
#     sample_mid_prices¹ = sample_mid_price_path(slob¹, p¹) #calculate observed mid prices
#     sample_mid_prices² = sample_mid_price_path(slob², p²)
#     return   φ¹, source¹, p¹, sample_mid_prices¹, P⁺s¹, P⁻s¹, Ps¹,
#                          φ², source², p², sample_mid_prices², P⁺s², P⁻s², Ps²
# end

# +
# function dtrw_solver_with_recalc(slob¹::SLOB, slob²::SLOB,kick::Kicker)
#     time_steps = get_time_steps(slob¹.T, slob¹.Δt) #let's assume T,L,M,D,num_paths,α is the same for both 
    
#     φ¹ = ones(Float64, slob¹.M + 1, time_steps + 1) #stores φ¹ for all time
#     φ² = ones(Float64, slob².M + 1, time_steps + 1) #stores φ² for all time
    
#     source¹ = ones(Float64, slob¹.M + 1, time_steps + 1) #stores sources for all time
#     source² = ones(Float64, slob².M + 1, time_steps + 1) #stores sources for all time
    
#     coupling¹ = ones(Float64, slob¹.M + 1, time_steps + 1) #stores couplings for all time
#     coupling² = ones(Float64, slob².M + 1, time_steps + 1) #stores couplings for all time

#     p¹ = ones(Float64, time_steps + 1)
#     p² = ones(Float64, time_steps + 1)
    
#     mid_prices¹ = ones(Float64, slob¹.T + 1) #stores mid prices for all time
#     mid_prices² = ones(Float64, slob².T + 1) #stores mid prices for all time

#     p¹[1] = slob¹.p₀
#     mid_prices¹[1] = slob¹.p₀
    
#     p²[1] = slob².p₀
#     mid_prices²[1] = slob².p₀

#     P⁺s¹ = fill(1/3, time_steps) #stores right jump prob for all time
#     P⁻s¹ = fill(1/3, time_steps) #store left jump prob for all time
#     Ps¹ = fill(1/3, time_steps) #store self jump prob for all time
    
#     P⁺s² = fill(1/3, time_steps) #stores right jump prob for all time
#     P⁻s² = fill(1/3, time_steps) #store left jump prob for all time
#     Ps² = fill(1/3, time_steps) #store self jump prob for all time
    
#     t = 1
    
#     φ¹[:, t], source¹[:,t], coupling¹[:,t] = initial_conditions_numerical(slob¹, p¹[t], p²[t], 0.0) #get initial shape given at t=1
#                                                             #the source term and mid price position
#     φ²[:, t], source²[:,t], coupling²[:,t] = initial_conditions_numerical(slob², p²[t], p¹[t], 0.0) #get initial shape given at t=1
#                                                             #the source term and mid price position

#     while t <= time_steps
#         τ, τ_periods = get_sub_period_time(slob¹, t, time_steps) #get time periods till next recalc

#         for τₖ = 1:τ_periods
#             t += 1 
            
#             kick.DoNow =  kick.Do&&(t>=kick.StartTime)&&(t<kick.EndTime) 
            
#             φ¹[:, t], source¹[:,t], coupling¹[:,t], P⁺s¹[t-1], P⁻s¹[t-1], Ps¹[t-1],  
#             φ²[:, t], source²[:,t], coupling²[:,t], P⁺s²[t-1], P⁻s²[t-1], Ps²[t-1] =
#             intra_time_period_simulate(     slob¹, φ¹[:, t-1], p¹[t-1],
#                                             slob², φ²[:, t-1], p²[t-1], 
#                                             kick                        )#get φ at next time step
            
                                                        
#             p¹[t] = extract_mid_price(slob¹, φ¹[:, t])#get mid price of new φ. Raw mid prices
#             p²[t] = extract_mid_price(slob², φ²[:, t])#get mid price of new φ. Raw mid prices

#         end
        
        
#         t += 1
        
#         if t > time_steps   #the loop ends if we happen to have exceeded the time
#             sample_mid_prices¹ = sample_mid_price_path(slob¹, p¹) #calculate observed mid prices
#             sample_mid_prices² = sample_mid_price_path(slob², p²)
#             return   φ¹, source¹, p¹, sample_mid_prices¹, P⁺s¹, P⁻s¹, Ps¹,
#                          φ², source², p², sample_mid_prices², P⁺s², P⁻s², Ps²
#         end
        
        
#         if slob¹.α > 0.0 || (t==kick.StartTime&&kick.Do==true)#SLOB resolves initial conditions at times ~ Exp(α) if α>0
#                         #otherwise this is just LOB which only solves initial conditions once
#             φ¹[:, t], source¹[:,t], coupling¹[:,t] = initial_conditions_numerical(slob¹, p¹[t-1], p²[t-1])
#             φ²[:, t], source²[:,t], coupling²[:,t] = initial_conditions_numerical(slob², p²[t-1], p¹[t-1])
            
#         end #move down 4 lines????????????????????????/

#         p¹[t] = extract_mid_price(slob¹, φ¹[:, t])#get mid price of new φ. Raw mid prices
#         p²[t] = extract_mid_price(slob², φ²[:, t])#get mid price of new φ. Raw mid prices

        
#     end
    
    
#     sample_mid_prices¹ = sample_mid_price_path(slob¹, p¹) #calculate observed mid prices
#     sample_mid_prices² = sample_mid_price_path(slob², p²)
#     return   φ¹, source¹, p¹, sample_mid_prices¹, P⁺s¹, P⁻s¹, Ps¹,
#                          φ², source², p², sample_mid_prices², P⁺s², P⁻s², Ps²
# end

# +
# Insertion 1   # to use these insertions, line up the this line with the corresponding line
                # corresponding line in the original code but include everything below
#             try
#                 p¹[t] = extract_mid_price(slob¹, φ¹[:, t])#get mid price of new φ. Raw mid prices
#                 p²[t] = extract_mid_price(slob², φ²[:, t])#get mid price of new φ. Raw mid prices
                
#             catch e
#                 println("Bounds Error at t=$t")
                
#                 mid_prices¹ = sample_mid_price_path(slob¹, p¹)
#                 mid_prices² = sample_mid_price_path(slob², p²)
                
#                 return   φ¹, p¹, mid_prices¹, P⁺s¹P⁻s¹, Ps¹,
#                          φ², p², mid_prices², P⁺s², P⁻s², Ps²
#             end

#             @info "Intra-period simulation. tick price = R$(p¹[t]) @t=$t"
#             @info "Intra-period simulation. tick price = R$(p²[t]) @t=$t"

# +
# Insertion 2   # to use these insertions, line up the this line with the corresponding line
                # corresponding line in the original code but include everything below
#         try 
#             p¹[t] = extract_mid_price(slob¹, φ¹[:, t])#get mid price of new φ. Raw mid prices
#             p²[t] = extract_mid_price(slob², φ²[:, t])#get mid price of new φ. Raw mid prices
#         catch e
#             println("Initial Conditions Bounds Error at t=$t")
#             mid_prices¹ = sample_mid_price_path(slob¹, p¹) #calculate observed mid prices
#             mid_prices² = sample_mid_price_path(slob², p²)
#             return   φ¹, p¹, mid_prices¹, P⁺s¹, P⁻s¹, Ps¹,
#                          φ², p², mid_prices², P⁺s², P⁻s², Ps²
#         end


#         @info "LOB Density recalculated. tick price = R$(p¹[t]) @t=$t"
#         @info "LOB Density recalculated. tick price = R$(p²[t]) @t=$t"

# +
# not being used
# function dtrw_solver(slob::SLOB)
#     time_steps = get_time_steps(slob.T, slob.Δt)
#     p = ones(Float64, time_steps + 1)
#     p[1] = slob.p₀

#     t = 1
#     φ = initial_conditions_numerical(slob, p[t], 0.0)   #get initial shape given
#                                                         #the source term and mid price position

#     while t <= time_steps
#         τ, τ_periods = get_sub_period_time(slob, t, time_steps)

#         for τₖ = 1:τ_periods
#             t += 1
#             φ, _, _, _  = intra_time_period_simulate(slob, φ, p[t-1]) #get φ at next time step
#             p[t] = extract_mid_price(slob, φ) #get mid price of new φ. Raw mid prices

#         end
        
#         if t > time_steps  #the loop ends if we happen to have exceeded the time
#             mid_prices = sample_mid_price_path(slob, p)
#             return mid_prices #calculate and return observed mid prices
#         end
        
#         t += 1
        
#         if slob.α > 0.0 #SLOB resolves initial conditions at times ~ Exp(α) if α>0
#                         #otherwise this is just LOB which only solves initial conditions once
#             φ = initial_conditions_numerical(slob, p[t-1])
#         end
        
#         p[t] = extract_mid_price(slob, φ)
        
#     end
    
#     mid_prices = sample_mid_price_path(slob, p)
#     return mid_prices #calculate and return observed mid prices
# end

# -

a = zeros(Float64,2).+1

a += 2*[1,2]

1:4-1


