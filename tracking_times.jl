using DifferentialEquations, Random

total_pool_size = 1000

function compute_durations(row)
    durations_of_0 = Int[]  # To store lengths of consecutive 0s
    durations_of_1 = Int[]  # To store lengths of consecutive 1s
    durations_of_2 = Int[]  # To store lengths of consecutive 2s

    current_value = row[1]  # Initialise with the first element
    current_length = 1      # Initialise length counter

    for i in 2:length(row)
        if row[i] == current_value
            current_length += 1  # If same as previous, increment the length
        else
            # Save the completed run for 0, 1, or 2
            if current_value == 0
                push!(durations_of_0, current_length)
            elseif current_value == 1
                push!(durations_of_1, current_length)
            elseif current_value == 2
                push!(durations_of_2, current_length)
            end
            # Update the current value and reset the length counter
            current_value = row[i]
            current_length = 1
        end
    end

    # Save the last run
    if current_value == 0
        push!(durations_of_0, current_length)
    elseif current_value == 1
        push!(durations_of_1, current_length)
    elseif current_value == 2
        push!(durations_of_2, current_length)
    end

    return durations_of_0, durations_of_1, durations_of_2
end


function get_durations_from_matrix(M)
    no_rows = size(M)[1]
    dur0s = []
    dur1s = []
    dur2s = []

    for i in 1:no_rows # run through all the matrix rows and store the durations of 0s, 1s, 2s
        d0s, d1s, d2s = compute_durations(M[i,:])
        append!(dur0s,d0s)
        append!(dur1s, d1s)
        append!(dur2s,d2s)
    end

    return dur0s, dur1s, dur2s
end

function track_times_constant_rates(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, kesten_time_step)

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
    c, m, e, i = rates
    A = i
    push!(pool_history, pool)
    push!(immature_history, immature)
    push!(mature_history, mature)
    
    # Synapse sizes for mature population
    synapse_sizes = Float64[]
    synapse_size_history = []
    
    
    # Simulation
    for t in 1:steps
        # 1 Transitions from pool to immature
        pool_to_immature = rand(Binomial(pool, c * kesten_time_step))
        pool -= pool_to_immature
        immature += pool_to_immature
    
        # 2 Transitions from immature to mature
        immature_to_mature = rand(Binomial(immature, m * kesten_time_step))
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
        for (i, size) in enumerate(synapse_sizes)
            prob = A * exp(-size / λ) * kesten_time_step
            if rand() < prob
                push!(mature_to_immature_indices, i)
            end
        end
    
        # Update states based on calculated probabilities
        mature_to_immature = length(mature_to_immature_indices) # if simulating it stochastically
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
        immature_to_pool = rand(Binomial(immature, e * kesten_time_step))
        immature -= immature_to_pool
        pool += immature_to_pool

        ######

        function safe_sample(list, num_samples; replace=false)
            return sample(list, min(num_samples, length(list)), replace=false)
        end

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
    
        push!(pool_history, pool)
        push!(immature_history, immature)
        push!(mature_history, mature)
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








# Define transition rates
c, m, e, i = 0.2, 0.2, 0.1, 0.05
rates = c, m, e, i


ih, mh, state_records, syn_sizes = track_times_constant_rates(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, kesten_timestep);


final_I_value = total_pool_size / (1 + m/i + e/c)
final_M_value = total_pool_size / (1 + i/m + (e*i)/(c*m))

plot(0:0.01:100, ih)
plot!(0:0.01:100, mh)
plot!(0:0.01:100, ih+mh)
hline!([final_I_value,final_M_value],label="Steady state solutions", linestyle= :dash,lw=3)


plot(state_records[1,:])

d0, d1, d2 = get_durations_from_matrix(state_records)

h0 = histogram(d0 .* kesten_time_step)
h1 = histogram(d1 .* kesten_time_step)
h2 = histogram(d2 .* kesten_time_step)


plot(h0,h1,h2, layout=(3,1))




# Count the number of 1s in each column
ones_per_column = sum(state_records .== 1, dims=1)
twos_per_column = sum(state_records .== 2, dims=1)

# Plot the result
plot(0:0.01:99.99, ones_per_column[:], xlabel="Column", ylabel="Number of 1s", label="1s in each column", lw=2)
plot!(0:0.01:99.99, twos_per_column[:])





#######
#######
# Variables rates
#######
#######

el(t) = 0.1 * exp(-t / 10) + 0.2
cr(t) = 0.2 * exp(-t / 30) + 0.2

elim = el.(1:100)
creat = cr.(1:100)

varrates = plot(elim, label="Elimination rate", ylim=(0,0.5), lw=2)
plot!(creat, label="Creation rate", lw=2, xlabel="time")

function synapse_dynamics_var!(du, u, p, t)
    m, i, λ, synapse_sizes = p 
    N_I, N_M, N_P = u
    A = i

    # Compute the rate of dematuration using the exponential probability distribution
    # Sum over all mature synapses' probabilities of transitioning to immature state
    if !isempty(synapse_sizes)
        dematuration_rate = A * sum(exp(- size / λ) for size in synapse_sizes) / length(synapse_sizes)
    else
        dematuration_rate = A
    end
    
    # Apply time-dependent e(t) and c(t)
    e_t = el(t)
    c_t = cr(t)
    # m = dematuration_rate

    du[1] = c_t * N_P - (m + e_t) * N_I + (dematuration_rate) * N_M  # dN_I/dt
    du[2] = m * N_I - (dematuration_rate) * N_M  # dN_M/dt
    du[3] = - du[1] - du[2]  # dN_P/dt
end


