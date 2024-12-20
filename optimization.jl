using Optimization
using OptimizationBBO
using Optimization, OptimizationOptimJL
using .syn_maturation_functions



############
# Practice
############



function optimise_synapticmaturation(x, p)
    total_pool_size = Int(p[1])
    total_time = p[2] 
    kesten_timestep = p[3]
    ε = p[4]
    η = p[5]
    σ_ε = p[6]
    σ_η = p[7]
    c, m, e, i = x
    rates = (c, m, e, i)

    # Run your simulation
    sol, synapse_sizes_diffeq, synapses_diffeq = syn_maturation_functions.run_simulation_diffeq(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, kesten_timestep)

    # Extract the mature population at the end of the simulation
    mature_population_diffeq = sol[2, :]
    
    target = 100
    # Calculate the error
    error = mature_population_diffeq[end] - target

    # Return the squared error
    return error^2
end


function optim(x,p)
    total_pool_size = Int(p[1])
    a1,a2,b1,b2,m,i = x
    final_I_value = (i*total_pool_size*(a1+b1))/(m*(a1+b1)+i*(a1+a2+b1+b2))
    final_M_value = (total_pool_size*(a1+b1))/(a1+b1+(i/m)*(a1+a2+b1+b2))

    target_M = 50
    target_I = 40

    errorM = final_M_value-target_M
    errorI = final_I_value-target_I

    return errorM^2 + errorI^2
end



# Define initial guesses for parameters
x0 = [0.5,0.5,0.5, 0.5, 0.1, 0.1]

# Define the parameters you want to pass to the objective function
total_pool_size = 100  # Example values
total_time = 100.0
kesten_timestep = 0.01
ε = 1.0
η = 0.
σ_ε = 0.01
σ_η = 0.01

p = [total_pool_size, total_time, kesten_timestep, ε, η, σ_ε, σ_η]


# Define bounds for the parameters
lower_bounds = [0.0, 0.0, 0.0, 0.0,0.0,0.0]
upper_bounds = [2.0, 2.0, 2.0, 2.0,2.0,2.0]

# Set up the optimization problem with bounds
opt_function = OptimizationFunction(optim, Optimization.AutoForwardDiff())
prob = Optimization.OptimizationProblem(opt_function, x0, p, lb=lower_bounds, ub=upper_bounds)

# Solve the optimization problem using LBFGS with bounds
result = Optimization.solve(prob, LBFGS())

# Output the results
println("Optimal parameters (c, m, e, i): ", result.minimizer)
println("Objective function value: ", result.minimum)





###########
# OPTIMISE for PARAMETERS
###########


using Optimization
using OptimizationBBO
using Optimization, OptimizationOptimJL
using OptimizationEvolutionary




function track_times_variable_rates_optimise(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, kesten_time_step)

    pool = total_pool_size
    steps = trunc(Int, total_time / kesten_time_step)
    immature = 0
    mature = 0
    pool_history = []
    immature_history = []
    mature_history = []

    state_records = zeros(total_pool_size, trunc(Int,total_time/kesten_time_step))

    I_pop = []
    M_pop = []
    pool_pop = [i for i in 1:total_pool_size]
    cr, m, el, A, lambda  = rates

    push!(pool_history, pool)
    push!(immature_history, immature)
    push!(mature_history, mature)
    
    # Synapse sizes for mature population
    synapse_sizes = Float64[]
    synapse_size_history = []
    
    
    # Simulation
    for t in 1:steps
        # 1 Transitions from pool to immature
        pool_to_immature = 0
        transition_prob1 = cr[t] * kesten_timestep
        for i in 1:pool
            if rand() < transition_prob1
                pool_to_immature += 1
            end
        end
        pool -= pool_to_immature
        immature += pool_to_immature
    
        # 2 Transitions from immature to mature
        immature_to_mature = 0
        transition_prob2 = m * kesten_timestep  # Probability for each immature synapse to become mature

        for i in 1:immature
            if rand() < transition_prob2
                immature_to_mature += 1
            end
        end

        # Update the counts
        immature -= immature_to_mature
        mature += immature_to_mature
    
        # Initialize new mature synapse sizes
        for i in 1:immature_to_mature
            push!(synapse_sizes, 0.0)  # Initial size of a new mature synapse
        end
    
        synapse_sizes = sort(synapse_sizes, rev=true)
    
        # 3 Transitions from mature to immature
        mature_to_immature_indices = []
        for (id, size) in enumerate(synapse_sizes)
            prob = A * exp(-size / lambda) * kesten_timestep #i*kesten_time_step
            if rand() < prob
                push!(mature_to_immature_indices, id)
            end
        end

        
    
        # Update states based on calculated probabilities
        mature_to_immature = length(mature_to_immature_indices)
        mature_to_immature = round(Int, mature_to_immature)
        mature -= mature_to_immature
        immature += mature_to_immature
    
        # Remove synapse sizes for synapses that became immature
        # Sort indices in reverse order
        sorted_indices = sort(mature_to_immature_indices, rev=true)
    
        # Delete elements at the specified indices
        for idx in sorted_indices
            deleteat!(synapse_sizes, idx)
        end
    
        # 4 Transitions from immature to pool
        immature_to_pool = 0
        transition_prob3 = el[t] * kesten_timestep  # Probability for each immature synapse to transition to the pool

        for i in 1:immature
            if rand() < transition_prob3
                immature_to_pool += 1
            end
        end

        # Update the counts
        immature -= immature_to_pool
        pool += immature_to_pool

        
        push!(pool_history, pool)
        push!(immature_history, immature)
        push!(mature_history, mature)


        # 1. Pool to immature
        pool_to_immature_count = abs(trunc(Int, pool_to_immature))
        pool_to_immature_indxs = safe_sample(pool_pop, pool_to_immature_count, replace=false)
        filter!(x -> !in(x, pool_to_immature_indxs), pool_pop)
        append!(I_pop, pool_to_immature_indxs)
        
        # 2. Immature to mature
        immature_to_mature_count = trunc(Int, immature_to_mature)
        immature_to_mature_indxs = safe_sample(I_pop, immature_to_mature_count, replace=false)
        filter!(x -> !in(x, immature_to_mature_indxs), I_pop)
        append!(M_pop, immature_to_mature_indxs)
        
        # 3. Mature to immature
        mature_to_immature_count = abs(trunc(Int, mature_to_immature))
        mature_to_immature_indxs = safe_sample(M_pop, mature_to_immature_count, replace=false)
        filter!(x -> !in(x, mature_to_immature_indxs), M_pop)
        append!(I_pop, mature_to_immature_indxs)
        
        # 4. Immature to pool
        immature_to_pool_count = trunc(Int, immature_to_pool)
        immature_to_pool_indxs = safe_sample(I_pop, immature_to_pool_count, replace=false)
        filter!(x -> !in(x, immature_to_pool_indxs), I_pop)
        append!(pool_pop, immature_to_pool_indxs)

        synapse_sizes = kesten_update_new(synapse_sizes, ε, η, σ_ε, σ_η)
    

        push!(synapse_size_history, synapse_sizes)

        # Now recording the current state of the synapses in the state record matrix
        for j in 1:total_pool_size
            if j in pool_pop
                state_records[j,t] = 0
            elseif j in I_pop
                state_records[j,t] = 1
            elseif j in M_pop
                state_records[j,t] = 2
            end
        end
    
    end

    return immature_history, mature_history, state_records, synapse_sizes

end


function find_optimal_parameters(x,p)
    total_pool_size = Int(p[1])
    total_time = Int(p[2])

    num_trials = 5
    
    a1,k1,b1,a2,k2,b2,m,A,lambda = x #,ε,η,σ_ε, σ_η = x
    kesten_timestep = 0.01


    creation_func(t) = a1 * exp(-t * k1) + b1
    elimination_func(t) = a2 * exp(-t * k2) + b2

    elim = elimination_func.(0:kesten_timestep:total_time)
    creat = creation_func.(0:kesten_timestep:total_time)


    rates_var = creat, m, elim, A, lambda


    state_recs_var_multiple = []

    for i in 1:num_trials
        ih_var, mh_var, state_record_var, syn_sizes_var = track_times_variable_rates_optimise(total_time, total_pool_size, rates_var, ε, η, σ_ε, σ_η, kesten_timestep);
        push!(state_recs_var_multiple, state_record_var)
    end

    develop_survival_multiplee = []
    adult_survival_multiplee = []

    for state_recs in state_recs_var_multiple
        developmental_period_16 = round(Int, (16/100)*size(state_recs,2))
        developmental_period_26 = round(Int, (26/100)*size(state_recs,2))
        
        adult_period = round(Int, (70/100)*size(state_recs,2))
        adult_period2 = round(Int, (88/100)*size(state_recs,2))

        developmental_survival_fraction1 = compute_survival_fraction(state_recs[:,developmental_period_16:developmental_period_26])
        adulthood_survival_fraction1 = compute_survival_fraction(state_recs[:,adult_period:adult_period2])
        push!(develop_survival_multiplee, developmental_survival_fraction1)
        push!(adult_survival_multiplee, adulthood_survival_fraction1)
    end

    dev_ids = collect(0:1:10)
    dev_ids = [round(Int, id/kesten_timestep) for id in dev_ids]

    adult_ids = [0,1,2,3,4,5,6,8,10,12,14,16,17,18]
    adult_ids = [round(Int, id/kesten_timestep) for id in adult_ids]
    adult_ids3 = [0,1,2,3,4,5,6,8,10,12,14,16,17,18]

    development_points_to_match_sim = [mean(develop_survival_multiplee)[id+1] for id in dev_ids]
    development_points_to_match_data = [1.0, 0.661896208, 0.52522361,0.468246877, 0.421466905, 0.397137735, 0.376028593, 0.364221812, 0.344543843, 0.348389962, 0.340339859]
    
    adulthood_points_to_match_sim = [mean(adult_survival_multiplee)[id+1] for id in adult_ids]
    adulthood_points_to_match_data = [1.0, 0.870199702, 0.82058372, 0.788018458, 0.775729644, 0.755248343, 0.7490909229625357, 0.7400000138716264, 0.7290909507057883, 0.7163636641068893, 0.7054545315829192, 0.694545468417081, 0.688556071, 0.681643617]

    development_survival_error = development_points_to_match_sim - development_points_to_match_data
    adulthood_survival_error = adulthood_points_to_match_sim - adulthood_points_to_match_data





    # work out the peak value of the combined populations and check if it occurs !at beginning and !end


    
    # Penalising if creation starts off higher than elimination
    total_error = sum(development_survival_error.^2) + sum(adulthood_survival_error.^2)


    if elim[1] < creat[1]
        total_error *= 2
    end

    # Penalising if elimination is higher than creation the whole time
    if all(elim .> creat)
        total_error *= 2
    end

    return total_error
end


# x = a1,k1,b1,a2,k2,b2,m,A,lambda,ε,η,σ_ε, σ_η
total_pool_size = 10
total_time = 100
kesten_timestep = 0.01
ε, η = .985, 0.015
σ_ε, σ_η = .05, .05

# Initial starting point x0
# x = [a1,k1,b1,a2,k2,b2,m,A,lambda] ,ε,η,σ_ε, σ_η
x0 = [0.4,1/30,0.2,0.8,1/10,0.2,0.01,0.1,1.0] #, ε, η, σ_ε, σ_η]

p = [total_pool_size, total_time]

# Define bounds for the parameters
lower_bounds = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] #, 0.0, 0.0, 0.0, 0.0]
upper_bounds = [10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0] #, 1.0, 10.0, 0.2, 0.2]

# Set up the optimization function with automatic differentiation
opt_function = OptimizationFunction(find_optimal_parameters, Optimization.AutoForwardDiff())

# Define the optimization problem with bounds
prob = OptimizationProblem(
    opt_function, 
    x0, 
    p, 
    lb = lower_bounds, 
    ub = upper_bounds
)

# Run the optimization with NelderMead
sol = solve(prob, NelderMead(), maxiters=100, show_trace=true)
