using .syn_maturation_functions
using DifferentialEquations, Distributions, Plots
using Random
using Distributions
using Plots, StatsBase


# Parameters
total_time = 100.0
total_pool_size = 1000
c, m, e, i = 0.2, 0.2, 0.01, 0.05
ε, η = 1, 0
σ_ε, σ_η = .5, .5
rates = (c, m, e, i)
kesten_time_step = 0.01
A = 0.05
lambda = 2
num_trials = 10


# Run multiple random walks simulations
immature_total = []
mature_total = []
total_synapse_sizes_randwalks = []

for i in 1:num_trials
    synapses = Int[]  
    synapse_sizes = Float64[]  # Array to hold sizes of mature synapses
    num_synapses = 0

    immature, mature, synapse_sizes, synapse_size_history = syn_maturation_functions.run_simulation_randwalks_exp(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, kesten_timestep, A, lambda)
    push!(immature_total, immature)
    push!(mature_total, mature)
    push!(total_synapse_sizes_randwalks, synapse_sizes)
end


time = 0:kesten_timestep:total_time
exp_plots = plot(time, immature_total, color=:pink, label=false,title="Over 10 trials", legend=:right)
plot!(time, mature_total, color=:lightblue, label=false)
plot!(time, immature_total[1], label="Immature Synapses (RandWalks)", lw=3, color=:pink)
plot!(time, mature_total[1], label="Mature Synapses (RandWalks)", lw=3, color=:lightblue)
plot!(time_array_diffeq, immature_population_diffeq, label = "Immature Synapses (DiffEq)", color="red", lw=3, legend=:right)
plot!(time_array_diffeq, mature_population_diffeq, label = "Mature Synapses (DiffEq)", color="blue", lw=3, xlabel="Time",ylabel="Population size")


savefig(exp_plots, "C://Users/B00955735/OneDrive - Ulster University/Desktop/exp_final_plots.svg")