using DifferentialEquations, Distributions, Plots
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
num_trials = 10

# Run diiferential equations
sol, synapse_sizes_diffeq, synapses_diffeq = syn_maturation_functions.run_simulation_diffeq(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η);
time_array_diffeq = sol.t
immature_population_diffeq = sol[1, :]
mature_population_diffeq = sol[2, :]


# Run multiple random walks simulations
time_walks = []
immature_total = []
mature_total = []
total_synapse_sizes = []
for i in 1:num_trials
    synapses = Int[]  
    synapse_sizes = Float64[]  # Array to hold sizes of mature synapses
    num_synapses = 0

    time_array_walks, immature_population_walks, mature_population_walks, synapse_sizes_walks = syn_maturation_functions.run_simulation_randwalks(total_time, total_pool_size, synapse_sizes, rates,ε, η, σ_ε, σ_η);
    push!(immature_total, immature_population_walks)
    push!(mature_total, mature_population_walks)
    push!(time_walks, time_array_walks)
    push!(total_synapse_sizes, synapse_sizes_walks)
end


# Plotting populations
rand_walks_plot = plot(time_walks, immature_total, color=:pink, label=false,title="Random Walks solution", lw=3, legend=:right)
plot!(time_walks, mature_total, color=:lightblue, label=false)
plot!(time_walks[1], immature_total[1], label="Immature Synapses", lw=3, color=:pink)
plot!(time_walks[1], mature_total[1], label="Mature Synapses", lw=3, color=:lightblue)

diffeq_plot = plot(time_array_diffeq, immature_population_diffeq, title="Differential equations solution", label = "Immature Synapses", color="red", lw=3, legend=:right)
plot!(time_array_diffeq, mature_population_diffeq, label = "Mature Synapses", color="blue", lw=3, xlabel="Time",ylabel="Population")


combined_plot = plot(rand_walks_plot, diffeq_plot, layout=(2,1))

# Plotting synapse sizes histograms
bin_edges = 0:2:10000


h = fit(Histogram, vcat(total_synapse_sizes...),bin_edges)
adjusted_weights = h.weights ./ num_trials
randwalks_hist = bar(h.edges, adjusted_weights, 
        ylabel="Frequency", 
        title="Average Distribution of Synapse Sizes (Rand Walks)", legend=false, xlim=(0,50),bar_width=2.0)

diffeq_hist = histogram(synapse_sizes_diffeq, nbins=100, title="Synapse Sizes (Diff Eq)", label=false, xlim=(0,50), ylabel="Frequency")


combined_hists = plot(randwalks_hist, diffeq_hist, layout=(2,1))