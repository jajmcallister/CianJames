using Random
using Distributions
using Plots, StatsBase

# Gillespie step function with pool and restriction
function gillespie_step!(pool, synapses, num_synapses, rates)
    c, m, e, i = rates

    # Calculating propensities
    a_c = c * length(pool)
    a_m = m * count(==(0), synapses)
    a_e = e * count(==(0), synapses)
    a_i = i * count(==(1), synapses)
    a_total = a_c + a_m + a_e + a_i

    if a_total == 0
        return num_synapses, Inf
    end

    # Time to the next reaction
    τ = rand(Exponential(1/a_total))

    # Determine which reaction occurs
    r = rand() * a_total

    if r < a_c
        # Creation: Move a synapse from the "resource pool" to the immature state
        append!(synapses, 0)
        pop!(pool)
        num_synapses += 1
    elseif r < a_c + a_m
        # Maturation: An immature synapse becomes mature
        index = findfirst(==(0), synapses) # This takes an immature synapse (i.e. a 0 in the synapses array) and turns it into a 1
        synapses[index] = 1 # ....and turns it into a 1
        push!(synapse_sizes, 0.0)  # Initialize the size of new mature synapse to 0
    elseif r < a_c + a_m + a_e
        # Elimination: Remove an immature synapse
        index = findfirst(==(0), synapses)
        deleteat!(synapses, index)
        # deleteat!(synapse_sizes, index)
        num_synapses -= 1
        # Add synapse back to the pool
        push!(pool, 1)
    else
        # Immature: A mature synapse becomes immature (with the smallest synaptic size)
        index = findfirst(==(1), synapses) # finds the first mature synapse in the synapses list
        synapses[index] = 0 # ... and sets it back to 0
        smallest_idx = findfirst(==(minimum(synapse_sizes)), synapse_sizes) # finds the smallest synapse size
        deleteat!(synapse_sizes, smallest_idx)  # and eliminates it from synapse sizes
    end

    return num_synapses, τ
end

# Kesten process update function
function kesten_update!(sizes, ε, η, σ_ε, σ_η)
    for i in 1:length(sizes)
        ε_i = rand(Normal(ε, σ_ε))
        η_i = rand(Normal(η, σ_η))
        new_size = ε_i * sizes[i] + η_i
        sizes[i] = max(new_size, 0.0)  # Ensure size is non-negative
    end
end

# Simulation function
function run_simulation(total_time, total_pool_size, rates,ε, η, σ_ε, σ_η)
    pool = [1 for _ in 1:total_pool_size]  # Initialize pool with synapses
    synapses = Int[]  # Array to hold states of synapses
    num_synapses = 0
    synapse_sizes = Float64[]  # Sizes of mature synapses

    # Initialise tracking variables
    time_array = [0.0]
    immature_population = [0]
    mature_population = [0]

    current_time = 0.0

    while current_time < total_time
        push!(time_array, current_time)
        push!(immature_population, count(==(0), synapses))
        push!(mature_population, count(==(1), synapses))

        num_synapses, τ = gillespie_step!(pool, synapses, num_synapses, rates)
        current_time += τ

        # Initialize sizes for new mature synapses
        for i in 1:length(synapses)
            if synapses[i] == 1 && i > length(synapse_sizes)
                push!(synapse_sizes, 0.0)  # Initialize the size of new mature synapses to 0.0
            end
        end

        # Update sizes of mature synapses
        # mature_synapse_indices = findall(==(1), synapses)
        # mature_synapse_sizes = [synapse_sizes[i] for i in mature_synapse_indices]
        kesten_update!(synapse_sizes,ε, η, σ_ε, σ_η)

        # # Update sizes in the original list
        # for (idx, new_size) in zip(mature_synapse_indices, new_sizes)
        #     synapse_sizes[idx] = new_size
        # end

    end

    return time_array, immature_population, mature_population, synapse_sizes
end

# Parameters
total_time = 100.0
total_pool_size = 100
c, m, e, i = 0.2,0.2,0.01,0.01
ε, η = 1.0, 0.0
σ_ε, σ_η = .5, .5
rates = (c, m, e, i)


# Initialisations of arrays which keep track of synapses: the "synapses" array is made up of 0s and 1s; a 0 for every immature synapse, and a 1 for a mature synapse.
synapses = Int[]  
synapse_sizes = Float64[]  # Array to hold sizes of mature synapses
num_synapses = 0


# Run simulation
time_array_walks, immature_population_walks, mature_population_walks, synapse_sizes_walks = run_simulation(total_time, total_pool_size, rates,ε, η, σ_ε, σ_η)

# Plot results
walkplot = plot(time_array_walks, immature_population_walks, label="Immature", xlabel="Time", ylabel="Population", legend=:topright)
plot!(time_array_walks, mature_population_walks, label="Mature")

hist_randwalks = histogram(synapse_sizes_walks, title="Distribution of Synapse Sizes (RandWalks)",legend=false)
maximum(synapse_sizes_walks)
