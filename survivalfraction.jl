# Function to compute the survival fraction from state_records
function compute_survival_fraction(state_records)
    total_time_steps = size(state_records, 2)

    # Identify the initial population
    initial_existing_synapses = findall(x -> x != 0, state_records[:, 1])

    # Initialize array to store the survival fraction over time
    survival_fraction = zeros(total_time_steps)

    # Loop over each time step and compute the fraction of surviving synapses
    for t in 1:total_time_steps
        # Find how many of the initial synapses are still in the existing population at time t
        surviving_synapses = count(x -> state_records[x, t] != 0, initial_existing_synapses)
        
        # Compute survival fraction as the ratio of surviving synapses to the initial population size
        survival_fraction[t] = surviving_synapses / length(initial_existing_synapses)
    end

    return survival_fraction
end

# function compute_survival_fraction(state_records)
#     total_time_steps = size(state_records, 2)

#     # Identify the initial population
#     initial_existing_synapses = findall(x -> x != 0, state_records[:, 1])

#     # Initialize array to store the survival fraction over time
#     survival_fraction = zeros(total_time_steps)

#     # Create an array to track if a synapse has ever left the population
#     has_left_population = falses(length(initial_existing_synapses))

#     # Loop over each time step and compute the fraction of surviving synapses
#     for t in 1:total_time_steps
#         surviving_synapses = 0

#         # Check each synapse in the initial population
#         for (index, synapse) in enumerate(initial_existing_synapses)
#             # If the synapse hasn't left and is still present, count it as surviving
#             if !has_left_population[index] && state_records[synapse, t] != 0
#                 surviving_synapses += 1
#             elseif state_records[synapse, t] == 0
#                 # Mark the synapse as having left the population permanently
#                 has_left_population[index] = true
#             end
#         end
        
#         # Compute survival fraction as the ratio of surviving synapses to the initial population size
#         survival_fraction[t] = surviving_synapses / length(initial_existing_synapses)
#     end

#     return survival_fraction
# end


total_pool_size = 1000
elim = elimination_func.(0:kesten_timestep:total_time)
creat = creation_func.(0:kesten_timestep:total_time)
rates_var = creat, m, elim, i
ih_var, mh_var, state_records_var, syn_sizes_var = track_times_variable_rates(total_time, total_pool_size, rates_var, ε, η, σ_ε, σ_η, kesten_timestep);


plot(0:0.01:total_time, ih_var)
plot!(0:0.01:total_time, mh_var)
plot!(0:0.01:total_time, ih_var+mh_var)
vline!([16,26], label="Developmental period")
vline!([70], label = "Adulthood")



developmental_period_16 = round(Int, 16/kesten_timestep)
developmental_period_26 = round(Int, 26/kesten_timestep)
adult_period = round(Int, 70/kesten_timestep)

# Compute survival fraction
developmental_survival_fraction = compute_survival_fraction(state_records_var[:,developmental_period_16:developmental_period_26])
adulthood_survival_fraction = compute_survival_fraction(state_records_var[:,adult_period:round(Int, 18/kesten_timestep)+adult_period])

# Plot survival fraction over time
developmental_survival_plot = plot(16:kesten_timestep:26, developmental_survival_fraction, xlabel="Postnatal Day", ylabel="Survival Fraction", xticks=16:1:26,
    title="Synapse Survival Fraction (Early Development)", lw=2, legend=false, ylim=(0,1))
     
adult_survival_plot = plot(0:kesten_timestep:18, adulthood_survival_fraction, xlabel="Days", ylabel="Survival Fraction",
title="Synapse Survival Fraction (Adulthood)", lw=2, legend=false, ylim=(0,1), xticks=0:1:18)
 
survival_fraction_plot = plot(developmental_survival_plot, adult_survival_plot, layout=(2,1))


# Multiple trials

num_trials = 2

state_recs_var_multiple = []

for i in 1:num_trials
    ih_var, mh_var, state_record_var, syn_sizes_var = track_times_variable_rates(total_time, total_pool_size, rates_var, ε, η, σ_ε, σ_η, kesten_timestep);
    push!(state_recs_var_multiple, state_record_var)
end

develop_survival_multiple = []
adult_survival_multiple = []

for state_recs in state_recs_var_multiple
    developmental_survival_fraction1 = compute_survival_fraction(state_recs[:,developmental_period_16:developmental_period_26])
    adulthood_survival_fraction1 = compute_survival_fraction(state_recs[:,adult_period:round(Int, 18/kesten_timestep)+adult_period])
    push!(develop_survival_multiple, developmental_survival_fraction1)
    push!(adult_survival_multiple, adulthood_survival_fraction1)
end


developmental_survival_plot_trials = plot(16:kesten_timestep:26, mean(develop_survival_multiple), ribbon=std(develop_survival_multiple), xlabel="Postnatal Day", ylabel="Survival Fraction", xticks=16:1:26,
    title="Synapse Survival Fraction (Early Development)", lw=2, legend=false, ylim=(0,1.05))
     
adult_survival_plot_trials = plot(0:kesten_timestep:18, mean(adult_survival_multiple), ribbon=std(adult_survival_multiple), xlabel="Days", ylabel="Survival Fraction",
title="Synapse Survival Fraction (Adulthood)", lw=2, legend=false, ylim=(0,1.05), xticks=0:1:18)
 
survival_fraction_plot_trials = plot(developmental_survival_plot_trials, adult_survival_plot_trials, layout=(2,1))

dev_ids = collect(0:1:10)
dev_ids = [round(Int, id/kesten_timestep) for id in dev_ids]

adult_ids = [0,1,2,3,4,5,17,18]
adult_ids = [round(Int, id/kesten_timestep) for id in adult_ids]

adult_ids1 = collect(0:1:18)
adult_ids2 = collect(0:1:18)
adult_ids2 = [round(Int, id/kesten_timestep) for id in adult_ids2]
adult_ids3 = [0,1,2,3,4,5,17,18]

development_points_to_match_sim = [mean(develop_survival_multiple)[id+1] for id in dev_ids]
development_points_to_match_data = [1.0, 0.661896208, 0.52522361,0.468246877, 0.421466905, 0.397137735, 0.376028593, 0.364221812, 0.344543843, 0.348389962, 0.340339859]
    
development_survival_error = development_points_to_match_sim - development_points_to_match_data

adulthood_points_to_match_sim = [mean(adult_survival_multiple)[id+1] for id in adult_ids]
adulthood_points_to_match_data = [1.0, 0.870199702, 0.82058372, 0.788018458, 0.775729644, 0.755248343, 0.688556071, 0.681643617]

adulthood_survival_error = adulthood_points_to_match_sim - adulthood_points_to_match_data

dev_scatter = scatter(16:1:26, development_points_to_match_sim, ylim=(0,1.05), label="Model", xticks=16:1:26)
scatter!(16:1:26, development_points_to_match_data, label="Data",title="Survival Fraction (Early Development)", xlabel="Postnatal day")
plot!(16:1:26, development_survival_error, color="red", label="Error")

adult_scatter = scatter(adult_ids1, [mean(adult_survival_multiple)[id+1] for id in adult_ids2], ylim=(0,1.05), label="Model", xticks=0:1:18)
scatter!(adult_ids3, adulthood_points_to_match_data, label="Data",title="Survival Fraction (Adulthood)", xlabel="Days")
plot!(adult_ids3, adulthood_survival_error, color="red", label="Error")


plot(dev_scatter, adult_scatter, layout=(2,1))

# Calculate loss/error
total_error = sum(development_survival_error.^2) + sum(adulthood_survival_error.^2)


###########
# OPTIMISE
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
    cr, m, el, i = rates

    push!(pool_history, pool)
    push!(immature_history, immature)
    push!(mature_history, mature)
    
    # Synapse sizes for mature population
    synapse_sizes = Float64[]
    synapse_size_history = []
    
    
    # Simulation
    for t in 1:steps
        # 1 Transitions from pool to immature
        # pool_to_immature = rand(Binomial(pool, cr[t] * kesten_timestep))
        # pool -= pool_to_immature
        # immature += pool_to_immature

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
        # immature_to_mature = rand(Binomial(immature, m * kesten_timestep))
        # immature -= immature_to_mature
        # mature += immature_to_mature
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
        # Calculate the probability (using exponential) for each mature synapse to become immature
        mature_to_immature_indices = []
        A = i
        for (id, size) in enumerate(synapse_sizes)
            prob = i*kesten_time_step #A * exp(-size / lambda) * kesten_timestep
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
        # immature_to_pool = rand(Binomial(immature, el[t] * kesten_timestep))
        # immature -= immature_to_pool
        # pool += immature_to_pool
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

        syn_maturation_functions.kesten_update!(synapse_sizes, ε, η, σ_ε, σ_η)
    

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
    total_time = p[2] 
    kesten_timestep = p[3]
    ε = p[4]
    η = p[5]
    σ_ε = p[6]
    σ_η = p[7]

    num_trials = 3
    b1,b2 = 0.2,0.2
    a1,k1,a2,k2,m,i = x

    creation_func(t) = a1 * exp(-t * k1) + b1
    elimination_func(t) = a2 * exp(-t * k2) + b2

    elim = elimination_func.(0:kesten_timestep:total_time)
    creat = creation_func.(0:kesten_timestep:total_time)


    rates_var = creat, m, elim, i


    state_recs_var_multiple = []

    for i in 1:num_trials
        ih_var, mh_var, state_record_var, syn_sizes_var = track_times_variable_rates_optimise(total_time, total_pool_size, rates_var, ε, η, σ_ε, σ_η, kesten_timestep);
        push!(state_recs_var_multiple, state_record_var)
    end

    develop_survival_multiple = []
    adult_survival_multiple = []

    for state_recs in state_recs_var_multiple
        developmental_survival_fraction1 = compute_survival_fraction(state_recs[:,developmental_period_16:developmental_period_26])
        adulthood_survival_fraction1 = compute_survival_fraction(state_recs[:,adult_period:round(Int, 18/kesten_timestep)+adult_period])
        push!(develop_survival_multiple, developmental_survival_fraction1)
        push!(adult_survival_multiple, adulthood_survival_fraction1)
    end

    dev_ids = collect(0:1:10)
    dev_ids = [round(Int, id/kesten_timestep) for id in dev_ids]

    adult_ids = [0,1,2,3,4,5,17,18]
    adult_ids = [round(Int, id/kesten_timestep) for id in adult_ids]


    development_points_to_match_sim = [mean(develop_survival_multiple)[id+1] for id in dev_ids]
    development_points_to_match_data = [1.0, 0.661896208, 0.52522361,0.468246877, 0.421466905, 0.397137735, 0.376028593, 0.364221812, 0.344543843, 0.348389962, 0.340339859]
    
    adulthood_points_to_match_sim = [mean(adult_survival_multiple)[id+1] for id in adult_ids]
    adulthood_points_to_match_data = [1.0, 0.870199702, 0.82058372, 0.788018458, 0.775729644, 0.755248343, 0.688556071, 0.681643617]

    development_survival_error = development_points_to_match_sim - development_points_to_match_data
    adulthood_survival_error = adulthood_points_to_match_sim - adulthood_points_to_match_data

    total_error = sum(development_survival_error.^2) + sum(adulthood_survival_error.^2)
    return total_error
end


# x = a1,k1,a2,k2,m,i
total_pool_size = 100
total_time = 100
kesten_timestep = 0.01
ε, η = .985, 0.015
σ_ε, σ_η = .05, .05

# Initial starting point
x0 = [0.4,1/30,0.8,1/10,0.2,0.1]
p = [total_pool_size, total_time, kesten_timestep, ε, η, σ_ε, σ_η]

# Define bounds for the parameters
lower_bounds = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
upper_bounds = [10.0, 10.0, 10.0, 10.0, 10.0, 10.0]

opt_function = OptimizationFunction(find_optimal_parameters, Optimization.AutoForwardDiff())
prob = Optimization.OptimizationProblem(opt_function, x0, p, lb=lower_bounds, ub=upper_bounds)
result = Optimization.solve(prob, LBFGS(); maxiters=10, g_tol=1e-3, f_tol=1e-3)


# Using CMAES
sigma = 0.5  # Initial standard deviation of the search distribution
population_size = 50  # Optional, population size
stop_tol = 1e-3  # Stop when the change in fitness is less than this

# Solve using CMAES optimizer
prob = Optimization.OptimizationProblem(opt_function, x0, p)
optimizer = CMAES()
result = Optimization.solve(prob, optimizer)

# Output the results
println("Optimal parameters: ", result.minimizer)
println("Objective function value: ", result.minimum)