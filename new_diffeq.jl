using DifferentialEquations, Distributions, Plots
using .syn_maturation_functions

# Parameters
total_time = 100.0
total_pool_size = 1000
c, m, e, i = 0.2,0.2,0.01,0.0191
ε, η = 1.0, 0.0
σ_ε, σ_η = .5, .5
rates = (c, m, e, i)
kesten_timestep = 0.01
lambda=2




# Run simulation
sol, synapse_sizes_diffeq, synapses_diffeq = syn_maturation_functions.run_simulation_diffeq_exp(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η,lambda, kesten_time_step);

# Plot the distribution of synapse sizes
xbins = 0:0.5:1000
hist_diffeq = histogram(synapse_sizes_diffeq, label="Synapse Sizes", bins=xbins, xlabel="Size", ylabel="Frequency",legend=false, xlim=(0,30))

# Extract the solution
time_array_diffeq = sol.t
immature_population_diffeq = sol[1, :]
mature_population_diffeq = sol[2, :]


diffeqplot = plot(time_array_diffeq, immature_population_diffeq, label = "Immature Synapses (DiffEq)", color="red", lw=3, legend=:right)
plot!(time_array_diffeq, mature_population_diffeq, label = "Mature Synapses (DiffEq)", color="blue", lw=3, xlabel="Time",ylabel="Population size")

# savefig(hist_diffeq, "C://Users/B00955735/OneDrive - Ulster University/Desktop/exp_hist_diffeq.svg")


