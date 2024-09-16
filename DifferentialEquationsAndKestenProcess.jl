using DifferentialEquations, Distributions, Plots
using .syn_maturation_functions

# Parameters
total_time = 100.0
total_pool_size = 100
# c, m, e, i = 1.3333327016873797, 0.8912588716692468, 2.4384797141525715e-15, 1.089316398653236
c, m, e, i = 0.2,0.2,0.05,0.05
ε, η = 1.0, 0.0
σ_ε, σ_η = .5, .5
rates = (c, m, e, i)
kesten_timestep = 0.01
# Run simulation
sol, synapse_sizes_diffeq1, synapses_diffeq1 = syn_maturation_functions.run_simulation_diffeq(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, kesten_timestep);

# Plot the distribution of synapse sizes
hist_diffeq = histogram(synapse_sizes_diffeq1, label="Synapse Sizes", title="Distribution of Synapse Sizes (DiffEqs)", xlabel="Size", ylabel="Frequency",legend=false, xlim=(0,30))

# Extract the solution
time_array_diffeq = sol.t
immature_population_diffeq = sol[1, :]
mature_population_diffeq = sol[2, :]

final_I_value = total_pool_size / (1 + m/i + e/c)
final_M_value = total_pool_size / (1 + i/m + (e*i)/(c*m))


diffeqplot = plot!(time_array_diffeq, immature_population_diffeq, title="Differential equation version of model", label = "Immature Synapses (DiffEq)", color="red", lw=3, legend=:right)
plot!(time_array_diffeq, mature_population_diffeq, label = "Mature Synapses (DiffEq)", color="blue", lw=3, xlabel="Time",ylabel="Population size")
plot!(time_array_diffeq, immature_population_diffeq+mature_population_diffeq, lw=3, label="Mature+Immature")

hline!([final_I_value,final_M_value],label="Steady state solutions", linestyle= :dash,lw=3)

# plot!(time_array_diffeq, mature_population_diffeq .+ immature_population_diffeq, label="Combined Population (Diff Eq)", color=:green, lw=3)
# plot!(time_array_walks, immature_population_walks, label="Immature Synapses (RandWalks)")
# plot!(time_array_walks, mature_population_walks, label="Mature Syanpses (RandWalks)")
# plot!(time_array_walks,immature_population_walks .+ mature_population_walks, label="Combined Population (Rand Walks)")


plot(diffeqplot, hist_diffeq, layout=(2,1))


# histogram_plots = plot(hist_randwalks,hist_diffeq,layout=(2,1))

# savefig(diffeqplot, "C://Users/B00955735/OneDrive - Ulster University/Desktop/sumplot.png")

plot(immature_population_diffeq, mature_population_diffeq, xlabel="N_I", ylabel="N_M", title="Phase Plane: N_I vs N_M", legend=false)

