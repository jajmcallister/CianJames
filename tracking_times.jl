using DifferentialEquations, Random, Distributions, Plots
using .syn_maturation_functions

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
        pool_to_immature = rand(Binomial(pool, c * kesten_timestep))
        pool -= pool_to_immature
        immature += pool_to_immature
    
        # 2 Transitions from immature to mature
        immature_to_mature = rand(Binomial(immature, m * kesten_timestep))
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
            prob = A * exp(-size / lambda) * kesten_timestep
            if rand() < prob
                push!(mature_to_immature_indices, i)
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
        immature_to_pool = rand(Binomial(immature, e * kesten_timestep))
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
c, m, e, i = 0.2, 0.2, 0.01, 0.05
rates = c, m, e, i
lambda = 2


ih, mh, state_records, syn_sizes = track_times_constant_rates(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, kesten_timestep);


final_I_value = total_pool_size / (1 + m/i + e/c)
final_M_value = total_pool_size / (1 + i/m + (e*i)/(c*m))

plot(0:0.01:total_time, ih)
plot!(0:0.01:total_time, mh)
plot!(0:0.01:total_time, ih+mh)
# hline!([final_I_value,final_M_value],label="Steady state solutions", linestyle= :dash,lw=3)

using Plots.PlotMeasures

durationplot = plot(0:0.01:99.99,state_records[1,:],yticks=([0,1,2],["Resource pool", "Immature", "Mature"]), xlabel="Time", bottommargin=5mm)

# savefig(durationplot, "C://Users/B00955735/OneDrive - Ulster University/Desktop/duration.png")

d0, d1, d2 = get_durations_from_matrix(state_records)

mean(d0 .* kesten_timestep)
mean(d1 .* kesten_timestep)
mean(d2 .* kesten_timestep)


bins = 0:2:100
h0 = histogram(d0 .* kesten_timestep,label=false, bins=bins,title="Resource pool")
h1 = histogram(d1 .* kesten_timestep,label=false, bins=bins,title="Immature")
h2 = histogram(d2 .* kesten_timestep,label=false, bins=bins,title="Mature", xlabel="Time in state")


duration_hist = plot(h0,h1,h2, layout=(3,1))

savefig(duration_hist, "C://Users/B00955735/OneDrive - Ulster University/Desktop/duration_hist.png")



# # Count the number of 1s and 2s in each column
# ones_per_column = sum(state_records .== 1, dims=1)
# twos_per_column = sum(state_records .== 2, dims=1)

# # Plot the result
# plot(0:0.01:99.99, ones_per_column[:], xlabel="Column", ylabel="Number of 1s", label="1s in each column", lw=2)
# plot!(0:0.01:99.99, twos_per_column[:])





# Multiple trials

ih_trials = []
mh_trials = []

for i in 1:20
    ih, mh, state_records, syn_sizes = track_times_constant_rates(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, kesten_timestep);
    push!(ih_trials,ih)
    push!(mh_trials, mh)
end

plot!(0:0.01:total_time, mean(ih_trials), ribbon=std(ih_trials))
plot!(0:0.01:total_time, mean(mh_trials), ribbon=std(mh_trials))
plot!(0:0.01:total_time, mean(ih_trials) .+ mean(mh_trials))

hline!([final_I_value,final_M_value],label="Steady state solutions", linestyle= :dash,lw=3)














#######
#######
# Variables rates
#######
#######



elimination_func(t) = 0.1 * exp(-t / 10) + 0.2
creation_func(t) = 0.2 * exp(-t / 30) + 0.2


function track_times_variable_rates(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, kesten_time_step)

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
        pool_to_immature = rand(Binomial(pool, cr[t] * kesten_timestep))
        pool -= pool_to_immature
        immature += pool_to_immature
    
        # 2 Transitions from immature to mature
        immature_to_mature = rand(Binomial(immature, m * kesten_timestep))
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
        for (id, size) in enumerate(synapse_sizes)
            prob = A * exp(-size / lambda) * kesten_timestep
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
        immature_to_pool = rand(Binomial(immature, el[t] * kesten_timestep))
        immature -= immature_to_pool
        pool += immature_to_pool

        
        push!(pool_history, pool)
        push!(immature_history, immature)
        push!(mature_history, mature)

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

total_time = 100

# Define transition rates
m, i = 0.2, 0.05
elim = elimination_func.(0:kesten_timestep:total_time)
creat = creation_func.(0:kesten_timestep:total_time)


rates_var = creat, m, elim, i
lambda = 2



ih_var, mh_var, state_records_var, syn_sizes_var = track_times_variable_rates(total_time, total_pool_size, rates_var, ε, η, σ_ε, σ_η, kesten_timestep);

final_I_value = total_pool_size / (1 + m/i + elim[end]/creat[end])
final_M_value = total_pool_size / (1 + i/m + (elim[end]*i)/(creat[end]*m))

plot(0:0.01:total_time, ih_var)
plot!(0:0.01:total_time, mh_var)
plot!(0:0.01:total_time, ih_var+mh_var)
hline!([final_I_value,final_M_value],label="Steady state solutions", linestyle= :dash,lw=3)


plot(0:0.01:99.99, state_records_var[1,:],yticks=([0,1,2],["Pool", "Immature", "Mature"]))

d0, d1, d2 = get_durations_from_matrix(state_records_var)

h0 = histogram(d0 .* kesten_timestep)
h1 = histogram(d1 .* kesten_timestep)
h2 = histogram(d2 .* kesten_timestep)


plot(h0,h1,h2, layout=(3,1))




# Count the number of 1s in each column
ones_per_column = sum(state_records_var .== 1, dims=1)
twos_per_column = sum(state_records_var .== 2, dims=1)

# Plot the result (from the state records)
plot!(0:0.01:99.99, ones_per_column[:], xlabel="Column", ylabel="Number of 1s", label="1s in each column", lw=2)
plot!(0:0.01:99.99, twos_per_column[:])



# Multiple trials
ih_trials_var = []
mh_trials_var = []

for i in 1:100
    ih_var_tr, mh_var_tr, state_records_var, syn_sizes_var = track_times_variable_rates(total_time, total_pool_size, rates_var, ε, η, σ_ε, σ_η, kesten_timestep);
    push!(ih_trials_var, ih_var_tr)
    push!(mh_trials_var, mh_var_tr)
end

plot!(0:0.01:total_time, mean(ih_trials_var), ribbon=std(ih_trials_var))
plot!(0:0.01:total_time, mean(mh_trials_var), ribbon=std(mh_trials_var))
plot!(0:0.01:total_time, mean(ih_trials_var) .+ mean(mh_trials_var), legend=false)



########
# Comparing constant v variable rates
########

ih, mh, state_records, syn_sizes = track_times_constant_rates(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, kesten_timestep);
ih_var, mh_var, state_records_var, syn_sizes_var = track_times_variable_rates(total_time, total_pool_size, rates_var, ε, η, σ_ε, σ_η, kesten_timestep);

d0, d1, d2 = get_durations_from_matrix(state_records)
d0_var, d1_var, d2_var = get_durations_from_matrix(state_records_var)

h0 = histogram(d0 .* kesten_timestep,label=false, title="Constant E&C Rates \n Pool", xlim=(0,100))
h1 = histogram(d1 .* kesten_timestep,label=false, title="Immature", xlim=(0,100))
h2 = histogram(d2 .* kesten_timestep,label=false, title="Mature", xlim=(0,100))

h0_var = histogram(d0_var .* kesten_timestep,label=false, title="Variable E&C rates \n Pool", xlim=(0,100))
h1_var = histogram(d1_var .* kesten_timestep,label=false, title="Immature", xlim=(0,100))
h2_var = histogram(d2_var .* kesten_timestep,label=false,title="Mature", xlim=(0,100))


plot(h0,h0_var,h1,h1_var,h2,h2_var,layout=(3,2))
