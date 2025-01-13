using Pkg


Pkg.add("Random")
Pkg.add("StatsBase")
Pkg.add("DifferentialEquations")
Pkg.add("Distributions")
Pkg.add("Statistics")
Pkg.add("Optimization")
Pkg.add("OptimizationBBO")
Pkg.add("OptimizationOptimJL")
Pkg.add("OptimizationEvolutionary")
Pkg.add("Tables")
Pkg.add("DataFrames")
Pkg.add("CSV")

using Random, StatsBase, DifferentialEquations, Distributions, Statistics
using Optimization, OptimizationBBO, OptimizationOptimJL, OptimizationEvolutionary
using Tables, DataFrames, CSV


function safe_sample(list, num_samples; replace=false)
    return sample(list, min(num_samples, length(list)), replace=false)
end

function kesten_update_new(sizes, ε, η, σ_ε, σ_η)
    sizes=deepcopy(sizes)
    for i in 1:length(sizes)
        ε_i = rand(Normal(ε, σ_ε))
        η_i = rand(Normal(η, σ_η))
        new_size = ε_i * sizes[i] + η_i
        if sizes[i] < 0
            sizes[i] = 0.0
        end
        # sizes[i] = max(new_size, 0.0)  # Ensure size is non-negative
        if sizes[i] < 0.0
            sizes[i]=0.0
        end
    end

    return sizes
end

function time_average(data, window) where T
    # Check if the window size is valid
    if window < 1 || window > length(data)
        throw(ArgumentError("Window size must be between 1 and the length of the data."))
    end

    # Compute the moving average
    averaged_data = [mean(data[max(1, i - window + 1):i]) for i in 1:length(data)]
    return averaged_data
end


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
    ihs = []
    mhs = []

    for i in 1:num_trials
        ih_var, mh_var, state_record_var, syn_sizes_var = track_times_variable_rates_optimise(total_time, total_pool_size, rates_var, ε, η, σ_ε, σ_η, kesten_timestep);
        push!(state_recs_var_multiple, state_record_var)
        push!(ihs, ih_var)
        push!(mhs, mh_var)
    end

    develop_survival_multiplee = []
    adult_survival_multiplee = []

    for state_recs in state_recs_var_multiple
        developmental_period_16 = round(Int, (16/100)*size(state_recs,2))
        developmental_period_26 = round(Int, (26/100)*size(state_recs,2))
        
        adult_period = round(Int, (70/100)*size(state_recs,2))
        adult_period2 = round(Int, (88/100)*size(state_recs,2))

        developmental_survival_fraction1 = new_compute_survival_fraction(state_recs[:,developmental_period_16:developmental_period_26])
        adulthood_survival_fraction1 = new_compute_survival_fraction(state_recs[:,adult_period:adult_period2])
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

    total_error = sum(development_survival_error.^2) + sum(adulthood_survival_error.^2)

    # work out the peak value of the combined populations and check if it occurs !at beginning and !end
    smoothed_avg = time_average(mean(ihs)+mean(mhs),100)
    max_val = maximum(smoothed_avg)
    end_val = smoothed_avg[end]
    id_max = argmax(smoothed_avg)
    

    if id_max < 2500 || id_max > 8000
        total_error *= 2
    end
    
    if max_val - end_val < 0.1*max_val
        total_error *= 2
    end
    

    
    # # Penalising if creation starts off higher than elimination

    # if elim[1] < creat[1]
    #     total_error *= 2
    # end

    # # Penalising if elimination is higher than creation the whole time
    # if all(elim .> creat)
    #     total_error *= 2
    # end

    return total_error
end


# x = a1,k1,b1,a2,k2,b2,m,A,lambda,ε,η,σ_ε, σ_η
total_pool_size = 10
total_time = 100
kesten_timestep = 0.01


a1 = 0.9
k1 = 1/30
b1 = 0.2
a2 = 1.8
k2 = 1/10
b2 = 0.2

ε, η = .985, 0.015
σ_ε, σ_η = .05, .05

A = 0.05
lambda = 0.5

m=0.05



# Initial starting point x0
# x = [a1,k1,b1,a2,k2,b2,m,A,lambda] ,ε,η,σ_ε, σ_η
x0 = [a1,k1,b1,a2,k2,b2,m,A,lambda] #, ε, η, σ_ε, σ_η]

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
sol = solve(prob, NelderMead(), maxiters=2, show_trace=true)

df = DataFrame(Number = sol)
CSV.write("/users/jmcallister/SynapticMaturation_DataFiles/optimal_parameters.csv", df)


