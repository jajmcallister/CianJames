using Random
using Distributions
using Plots, StatsBase
using .syn_maturation_functions


# Parameters
total_time = 100.0
total_pool_size = 1000
c, m, e, i = 0.2,0.2,0.01,0.05
ε, η = 1.0, 0.0
σ_ε, σ_η = .5, .5
rates = (c, m, e, i)
kesten_timestep = 0.01

# Initialisations of arrays which keep track of synapses: the "synapses" array is made up of 0s and 1s; a 0 for every immature synapse, and a 1 for a mature synapse.
synapses = Int[]  
synapse_sizes = Float64[]  # Array to hold sizes of mature synapses
num_synapses = 0


# Run simulation
time_array_walks, immature_population_walks, mature_population_walks, synapse_sizes_walks, synapse_size_history = syn_maturation_functions.run_simulation_randwalks(total_time, total_pool_size, synapse_sizes, rates,ε, η, σ_ε, σ_η, kesten_timestep);

# Plot results
walkplot = plot(time_array_walks, immature_population_walks, lw=4, label="Immature", xlabel="Time", ylabel="Population size", legend=:right, title="Random Walk Model of Mature and Immature Populations")
plot!(time_array_walks, mature_population_walks, lw=4, label="Mature", size=(800,500))
# plot!(time_array_walks,immature_population_walks .+ mature_population_walks, label="Combined Population")

hist_randwalks = histogram(synapse_sizes_walks, title="Distribution of Synapse Sizes (RandWalks)",legend=false,xlim=(0,30))

xbins=0:0.5:1000
syn_history1 = histogram(synapse_size_history[1],bins=xbins, xlim=(0,30), xlabel="Synapse size", ylabel="Frequency", label=false,ylim=(0,400),title="t = 1")
syn_history2 = histogram(synapse_size_history[5000],bins=xbins, xlim=(0,30),xlabel="Synapse size",label=false,ylim=(0,400), title="t = 50")
syn_history3 = histogram(synapse_size_history[10000],bins=xbins, xlim=(0,30),xlabel="Synapse size",label=false,ylim=(0,400), title = "t = 100")

using Plots.PlotMeasures
hist_plot_time = plot(syn_history1,syn_history2,syn_history3, layout=(1,3), size=(1000,300),bottommargin=8mm, leftmargin=5mm)





savefig(hist_plot_time, "C://Users/B00955735/OneDrive - Ulster University/Desktop/hist_plot_time.svg")
####################
