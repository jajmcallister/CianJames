using DifferentialEquations, Distributions, Plots
using .syn_maturation_functions

# Parameters
total_time = 100.0
total_pool_size = 1000
c, m, e, i = 0.2, 0.2, 0.01, 0.01
ε, η = 1.0, 0.0
σ_ε, σ_η = .5, .5
rates = (c, m, e, i)
kesten_time_step = 0.01
# Run simulation
sol, synapse_sizes_diffeq, synapses_diffeq = syn_maturation_functions.run_simulation_diffeq(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, kesten_time_step);

# Plot the distribution of synapse sizes
hist_diffeq = histogram(synapse_sizes_diffeq, label="Synapse Sizes", title="Distribution of Synapse Sizes (DiffEqs)", xlabel="Size", ylabel="Frequency",legend=false)

# Extract the solution
time_array_diffeq = sol.t
immature_population_diffeq = sol[1, :]
mature_population_diffeq = sol[2, :]


diffeqplot = plot(time_array_diffeq, immature_population_diffeq, title="Differential equations solution", label = "Immature Synapses (DiffEq)", color="red", lw=3, legend=:right)
plot!(time_array_diffeq, mature_population_diffeq, label = "Mature Synapses (DiffEq)", color="blue", lw=3, xlabel="Time",ylabel="Population")
plot!(time_array_diffeq, mature_population_diffeq .+ immature_population_diffeq, label="Combined Population (Diff Eq)", color=:green, lw=3)
plot!(time_array_walks, immature_population_walks, label="Immature Synapses (RandWalks)")
plot!(time_array_walks, mature_population_walks, label="Mature Syanpses (RandWalks)")
plot!(time_array_walks,immature_population_walks .+ mature_population_walks, label="Combined Population (Rand Walks)")


plot(diffeqplot, hist_diffeq, layout=(2,1))


histogram_plots = plot(hist_randwalks,hist_diffeq,layout=(2,1))
