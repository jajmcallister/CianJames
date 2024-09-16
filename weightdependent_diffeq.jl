using DifferentialEquations, Distributions, Plots
using .syn_maturation_functions

# Parameters
total_time = 100.0
total_pool_size = 100
c, m, e, i = 0.2, 0.2, 0.1, 0.05
ε, η = 1, 0
σ_ε, σ_η = .5, .5
rates = (c, m, e, i)
kesten_timestep = 0.01
A = i
lambda = 2



nis, nms, synapse_sizes_d, synapses, soll = syn_maturation_functions.run_simulation_diffeq_weightdependent(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, lambda, A, kesten_timestep);

# Extract the solution
time_array_diffeq = soll.t
immature_population_diffeq = soll[1, :]
mature_population_diffeq = soll[2, :]
ppool = soll[3,:]

diffeqplot = plot(time_array_diffeq, immature_population_diffeq, label = "Immature Synapses (DiffEq)", color="red", lw=3, legend=:right)
plot!(time_array_diffeq, mature_population_diffeq, label = "Mature Synapses (DiffEq)", color="blue", lw=3, xlabel="Time",ylabel="Population size")
plot!(time_array_diffeq, immature_population_diffeq + mature_population_diffeq)


plot(time_array_diffeq, (c/e) .* ppool)
plot!(time_array_diffeq, immature_population_diffeq)


