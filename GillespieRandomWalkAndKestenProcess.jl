using Random
using Distributions
using Plots, StatsBase
using .syn_maturation_functions


# Parameters
total_time = 100.0
total_pool_size = 1000
c, m, e, i = 0.2,0.2,0.01,0.01
ε, η = 1.0, 0.0
σ_ε, σ_η = .5, .5
rates = (c, m, e, i)


# Initialisations of arrays which keep track of synapses: the "synapses" array is made up of 0s and 1s; a 0 for every immature synapse, and a 1 for a mature synapse.
synapses = Int[]  
synapse_sizes = Float64[]  # Array to hold sizes of mature synapses
num_synapses = 0


# Run simulation
time_array_walks, immature_population_walks, mature_population_walks, synapse_sizes_walks = syn_maturation_functions.run_simulation_randwalks(total_time, total_pool_size, synapse_sizes, rates,ε, η, σ_ε, σ_η);

# Plot results
walkplot = plot(time_array_walks, immature_population_walks, label="Immature", xlabel="Time", ylabel="Population", legend=:topright)
plot!(time_array_walks, mature_population_walks, label="Mature")

hist_randwalks = histogram(synapse_sizes_walks, title="Distribution of Synapse Sizes (RandWalks)",legend=false)
maximum(synapse_sizes_walks)
