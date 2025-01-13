function new_compute_survival_fraction007(state_records)
    total_time_steps = size(state_records, 2)

    # Identify the initial population
    initial_existing_synapses = findall(x -> x != 0, state_records[:, 1])
    
    # Create a set to track permanently eliminated synapses
    permanently_eliminated = Set{Int}()

    # Initialize array to store the survival fraction over time
    survival_fraction = zeros(total_time_steps)

    # Loop over each time step and compute the fraction of surviving synapses
    for t in 1:total_time_steps
        # Update the set of permanently eliminated synapses
        for idx in initial_existing_synapses
            if state_records[idx, t] == 0
                push!(permanently_eliminated, idx)
            end
        end
        
        # Find how many of the initial synapses are still in the existing population and not permanently eliminated
        surviving_synapses = count(x -> !(x in permanently_eliminated) && state_records[x, t] != 0, initial_existing_synapses)
        
        # Compute survival fraction as the ratio of surviving synapses to the initial population size
        survival_fraction[t] = surviving_synapses / length(initial_existing_synapses)
    end

    return survival_fraction
end


function safe_sample007(list, num_samples; replace=false)
    return sample(list, min(num_samples, length(list)), replace=false)
end

function kesten_update_new007(sizes, ε, η, σ_ε, σ_η)
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

function time_average007(data, window)
    # Check if the window size is valid
    if window < 1 || window > length(data)
        throw(ArgumentError("Window size must be between 1 and the length of the data."))
    end

    # Compute the moving average
    averaged_data = [mean(data[max(1, i - window + 1):i]) for i in 1:length(data)]
    return averaged_data
end


function track_times_variable_rates_007(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, kesten_time_step)

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

        synapse_sizes = kesten_update_new007(synapse_sizes, ε, η, σ_ε, σ_η)
    

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


total_pool_size = 100
total_time = 100
kesten_timestep = 0.01


# a1 = 0.9
# k1 = 1/30
# b1 = 0.2
# a2 = 1.8
# k2 = 1/10
# b2 = 0.2

# a1,k1,b1,a2,k2,b2,m,A,lambda = 0.9, 0.03333333333333333, 0.2, 1.8, 0.17500000000000002, 0.2, 0.05, 0.05, 0.5

# ε, η = .985, 0.015
# σ_ε, σ_η = .05, .05

# A = 0.05
# lambda = 0.5

# m=0.05
num_trials = 5


creation_func(t) = a1 * exp(-t * k1) + b1
elimination_func(t) = a2 * exp(-t * k2) + b2

elim = elimination_func.(0:kesten_timestep:total_time)
creat = creation_func.(0:kesten_timestep:total_time)


rates_var = creat, m, elim, A, lambda


state_recs_var_multiple = []
ihs = []
mhs = []

for i in 1:num_trials
    ih_var, mh_var, state_record_var, syn_sizes_var = track_times_variable_rates_007(total_time, total_pool_size, rates_var, ε, η, σ_ε, σ_η, kesten_timestep);
    push!(state_recs_var_multiple, state_record_var)
    push!(ihs, ih_var)
    push!(mhs, mh_var)
end

develop_survival_multiplee = []
adult_survival_multiplee = []

for state_recs in state_recs_var_multiple
    developmental_period_16 = round(Int, (16/total_time)*size(state_recs,2))
    developmental_period_26 = round(Int, (26/total_time)*size(state_recs,2))
    
    adult_period = round(Int, (70/total_time)*size(state_recs,2))
    adult_period2 = round(Int, (88/total_time)*size(state_recs,2))

    developmental_survival_fraction1 = new_compute_survival_fraction007(state_recs[:,developmental_period_16:developmental_period_26])
    adulthood_survival_fraction1 = new_compute_survival_fraction007(state_recs[:,adult_period:adult_period2])
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



developmentperiodplot = 16:10/length(develop_survival_multiplee[1]):26
adulthoodperiodplot = 0:18/length(adult_survival_multiplee[1]):18

l1 = length(developmentperiodplot) - length(develop_survival_multiplee[1])
l2 = length(adulthoodperiodplot) - length(adult_survival_multiplee[1])

length(developmentperiodplot[1:end-l1])
# Plot survival fraction over time
developmental_survival_plot = plot(developmentperiodplot[1:end-l1], mean(develop_survival_multiplee), xlabel="Postnatal Day", ylabel="Survival Fraction", xticks=16:1:26,
    title="Synapse Survival Fraction (Early Development)", lw=2, legend=false, ylim=(0,1.05))
    scatter!(16:1:26, development_points_to_match_data, label="Data",title="Survival Fraction (Early Development)", xlabel="Postnatal day")

adult_survival_plot = plot(adulthoodperiodplot[1:end-l2], mean(adult_survival_multiplee), xlabel="Days", ylabel="Survival Fraction",
    title="Synapse Survival Fraction (Adulthood)", lw=2, legend=false, ylim=(0,1.05), xticks=0:1:18,label="Model")
    scatter!(adult_ids3, adulthood_points_to_match_data, label="Data",title="Survival Fraction (Adulthood)", xlabel="Days", legend=:bottomleft)

survival_fraction_plot = plot(developmental_survival_plot, adult_survival_plot, layout=(2,1))













# work out the peak value of the combined populations and check if it occurs !at beginning and !end
smoothed_avg = time_average(mean(ihs)+mean(mhs),100)
max_val = maximum(smoothed_avg)
end_val = smoothed_avg[end]
id_max = argmax(smoothed_avg)
    

ihs[1]

plot(mean(ihs))
plot!(mean(mhs))
plot!(mean(ihs).+mean(mhs))
