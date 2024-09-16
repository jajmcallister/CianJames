# Parameters
total_time = 100.0
total_pool_size = 100
m, i, λ = 0.1,0.05,2


elimination_func(t) = 0.1 * exp(-t / 10) + 0.2
creation_func(t) = 0.2 * exp(-t / 30) + 0.2


ε, η = 1, 0
σ_ε, σ_η = 0.5, 0.5
rates_var = creat, m, elim, i
kesten_time_step = 0.01


# Exponential parameter λ
A = i
λ = 2


function rand_walks_variable_rates(total_pool_size, total_time, rates, kesten_timestep)
    steps = Int(total_time / kesten_timestep)   
    pool = total_pool_size
    immature = 0
    mature = 0
    pool_history = []
    immature_history = []
    mature_history = []
    cre,m,eli,i = rates

    push!(pool_history, pool)
    push!(immature_history, immature)
    push!(mature_history, mature)

    # Synapse sizes for mature population
    synapse_sizes = Float64[]
    synapse_size_history = []

    # Simulation
    for t in 1:steps
        # Transitions from pool to immature
        pool_to_immature = rand(Binomial(pool, cre[t] * kesten_timestep))
        pool -= pool_to_immature
        immature += pool_to_immature

        # Transitions from immature to mature
        immature_to_mature = rand(Binomial(immature, m * kesten_timestep))
        immature -= immature_to_mature
        mature += immature_to_mature

        # Initialize new mature synapse sizes
        for i in 1:immature_to_mature
            push!(synapse_sizes, 0.0)  # Initial size of a new mature synapse
        end

        synapse_sizes = sort(synapse_sizes, rev=true)

        # Transitions from mature to immature
        # Calculate the probability (using exponential) for each mature synapse to become immature
        mature_to_immature_indices = []
        for (id, size) in enumerate(synapse_sizes)
            prob = A * exp(-size / λ) * kesten_timestep
            if rand() < prob
                push!(mature_to_immature_indices, id)
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


        # Transitions from immature to pool
        immature_to_pool = rand(Binomial(immature, eli[t] * kesten_timestep))
        immature -= immature_to_pool
        pool += immature_to_pool

        syn_maturation_functions.kesten_update!(synapse_sizes, ε, η, σ_ε, σ_η)

        push!(pool_history, pool)
        push!(immature_history, immature)
        push!(mature_history, mature)
        push!(synapse_size_history, synapse_sizes)

    end

    return immature_history, mature_history, synapse_sizes
end













time_walks = collect(0:0.01:total_time)

immature_history, mature_history, syn_sizes_rand_var = rand_walks_variable_rates(total_pool_size, total_time, rates_var, kesten_timestep)

var_plot_diffeq = plot!(time_array_var, immature_population_var, label = "Immature Synapses", color="red", lw=3, legend=:bottomright)
plot!(time_array_var, mature_population_var, label = "Mature Synapses", color="blue", lw=3, xlabel="Time",ylabel="Population size")
plot!(time_array_var, immature_population_var+mature_population_var, lw=3, label="Mature+Immature")


var_plot_randwalks = plot!(time_walks, immature_history, lw=1, label="Immature population")
plot!(time_walks, mature_history, lw=1, label="Mature population", legend=:right)
plot!(time_walks, immature_history+mature_history)




#### Trials


im_hist = []
mat_hist = []

for i in 1:100
    immature_history, mature_history, syn_sizes_rand_var = rand_walks_variable_rates(total_pool_size, total_time, rates_var, kesten_timestep)
    push!(im_hist, immature_history)
    push!(mat_hist, mature_history)
end

plot(0:kesten_timestep:total_time, mean(im_hist), ribbon=std(im_hist))
plot!(0:kesten_timestep:total_time, mean(mat_hist), ribbon=std(mat_hist))
plot!(0:kesten_timestep:total_time, mean(im_hist) .+ mean(mat_hist))

plot!(time_array_var, immature_population_var, label = "Immature Synapses", color="red", lw=3, legend=:bottomright)
plot!(time_array_var, mature_population_var, label = "Mature Synapses", color="blue", lw=3, xlabel="Time",ylabel="Population size")
plot!(time_array_var, immature_population_var+mature_population_var, lw=3, label="Mature+Immature")
hline!([immature_population_var[end]+mature_population_var[end]])

