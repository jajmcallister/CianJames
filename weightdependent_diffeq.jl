using DifferentialEquations, Distributions, Plots
using .syn_maturation_functions

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



nis, nms, synapse_sizes_d, synapses = syn_maturation_functions.run_simulation_diffeq_weightdependent(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, λ, A, kesten_time_step);

# Extract the solution
# time_array_diffeq = solution.t
# immature_population_diffeq = solution[1, :]
# mature_population_diffeq = solution[2, :]

# diffeqplot = plot(time_array_diffeq, immature_population_diffeq, label = "Immature Synapses (DiffEq)", color="red", lw=3, legend=:right)
# plot!(time_array_diffeq, mature_population_diffeq, label = "Mature Synapses (DiffEq)", color="blue", lw=3, xlabel="Time",ylabel="Population size")

time_diffs = collect(0:0.01:100)[1:end-1]

plot(time_diffs, nis, label = "Immature Synapses (DiffEq)", color="red", lw=3, legend=:right)
plot!(time_diffs, nms, label = "Mature Synapses (DiffEq)", color="blue", lw=3, xlabel="Time",ylabel="Population size")

nis[end]