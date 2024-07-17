using DifferentialEquations, Distributions, Plots
using Random
using Distributions
using Plots, StatsBase
using .syn_maturation_functions

# Parameters
total_time = 50.0
total_pool_size = 1000
c, m, e, i = 0.2,0.1,0.01,0.05
ε, η = 1.0, 0.0
σ_ε, σ_η = .5, .5
rates = (c, m, e, i)
num_trials = 10
kesten_time_step = 0.01

# Run differential equations
total_synapse_sizes_diffeq = []
sol, synapse_sizes_diffeq, synapses_diffeq = syn_maturation_functions.run_simulation_diffeq(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, kesten_time_step);
time_array_diffeq = sol.t
immature_population_diffeq = sol[1, :]
mature_population_diffeq = sol[2, :]
for i in 1:num_trials
    sol, synapse_sizes_diffeq, synapses_diffeq = syn_maturation_functions.run_simulation_diffeq(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, kesten_time_step);
    push!(total_synapse_sizes_diffeq, synapse_sizes_diffeq)
end



# Run multiple random walks simulations
time_walks = []
immature_total = []
mature_total = []
total_synapse_sizes_randwalks = []
for i in 1:num_trials
    synapses = Int[]  
    synapse_sizes = Float64[]  # Array to hold sizes of mature synapses
    num_synapses = 0

    time_array_walks, immature_population_walks, mature_population_walks, synapse_sizes_walks = syn_maturation_functions.run_simulation_randwalks(total_time, total_pool_size, synapse_sizes, rates,ε, η, σ_ε, σ_η, kesten_time_step);
    push!(immature_total, immature_population_walks)
    push!(mature_total, mature_population_walks)
    push!(time_walks, time_array_walks)
    push!(total_synapse_sizes_randwalks, synapse_sizes_walks)
end


# Plotting populations
rand_walks_plot = plot(time_walks, immature_total, color=:pink, label=false,title="Random Walks solution", legend=:bottomright)
plot!(time_walks, mature_total, color=:lightblue, label=false)
plot!(time_walks[1], immature_total[1], label="Immature Synapses", lw=3, color=:pink)
plot!(time_walks[1], mature_total[1], label="Mature Synapses", lw=3, color=:lightblue)

diffeq_plot = plot(time_array_diffeq, immature_population_diffeq, title="Differential equations solution", label = "Immature Synapses", color="red", lw=3, legend=:bottomright)
plot!(time_array_diffeq, mature_population_diffeq, label = "Mature Synapses", color="blue", lw=3, xlabel="Time",ylabel="Population")


combined_plot = plot(rand_walks_plot, diffeq_plot, layout=(2,1))

# Plotting synapse sizes histograms
bin_edges = 0:2:10000
yticks = 0:100:1000

h1 = fit(Histogram, vcat(total_synapse_sizes_randwalks...),bin_edges)
adjusted_weights1 = h1.weights ./ num_trials
randwalks_hist = bar(h1.edges, adjusted_weights1, 
        ylabel="Frequency", 
        title="Average Distribution of Synapse Sizes (Rand Walks)", legend=false, yticks=yticks, xlim=(0,20),bar_width=2.0, ylim=(0,650))

h2 = fit(Histogram, vcat(total_synapse_sizes_diffeq...),bin_edges)
adjusted_weights2 = h2.weights ./ num_trials
diffeq_hist = bar(h2.edges, adjusted_weights2, 
                ylabel="Frequency", 
                title="Average Distribution of Synapse Sizes (Diff Eq)", legend=false, yticks=yticks, xlim=(0,20),bar_width=2.0,ylim=(0,650))
        


combined_hists = plot(randwalks_hist, diffeq_hist, layout=(2,1))


total_plots = plot(rand_walks_plot, randwalks_hist, diffeq_plot, diffeq_hist, layout=(2,2), size=(1200,400))

# savefig(total_plots, "C://Users/B00955735/OneDrive - Ulster University/Desktop/syn_mat_plots.png")