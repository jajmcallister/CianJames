using Random
using Distributions
using Plots
using .syn_maturation_functions

# Parameters
N = 1000                    # Number of potential synapses
dt = 0.01                   # Time step
T = 100.0                    # Total simulation time
steps = Int(T/dt)           # Number of steps

# Parameters
total_time = 100.0
total_pool_size = 1000
c, m, e, i = 0.2, 0.2, 0.01, 0.05
ε, η = 1, 0
σ_ε, σ_η = .5, .5
rates = (c, m, e, i)
kesten_time_step = 0.01

# Exponential parameters
A = i
λ = 1

num_trials = 10

#Random walks
time_walks = []
immature_total = []
mature_total = []
total_synapse_sizes_randwalks = []

for i in 1:num_trials
    # States
    pool = N
    immature = 0
    mature = 0
    pool_history = []
    immature_history = []
    mature_history = []

    push!(pool_history, pool)
    push!(immature_history, immature)
    push!(mature_history, mature)

    # Synapse sizes for mature population
    synapse_sizes = Float64[]
    synapse_size_history = []


    # Simulation
    for t in 1:steps
        # Transitions from pool to immature
        pool_to_immature = rand(Binomial(pool, c * dt))
        pool -= pool_to_immature
        immature += pool_to_immature

        # Transitions from immature to mature
        immature_to_mature = rand(Binomial(immature, m * dt))
        immature -= immature_to_mature
        mature += immature_to_mature
        
        # Initialize new mature synapse sizes
        for i in 1:immature_to_mature
            push!(synapse_sizes, 0.0)  # Initial size of a new mature synapse
        end
        
        synapse_sizes = sort(synapse_sizes,rev=true)

        # Transitions from mature to immature
        # Calculate the probability (using exponential) for each mature synapse to become immature
        mature_to_immature_indices = []
        for (i, size) in enumerate(synapse_sizes)
            prob = A*exp(-size / λ)*dt
            if rand() < prob
                push!(mature_to_immature_indices, i)
            end
        end
        
        # Update states based on calculated probabilities
        #mature_to_immature = 0.0191*dt*mature # if taking the average fraction
        mature_to_immature = length(mature_to_immature_indices) # if simulating it stochastically
        mature_to_immature = round(Int, mature_to_immature)
        mature -= mature_to_immature
        immature += mature_to_immature
        
        # Remove synapse sizes for synapses that became immature
        # Sort indices in reverse order
        sorted_indices = sort(mature_to_immature_indices, rev=true)

        # Delete elements at the specified indices
        for idx in sorted_indices
            deleteat!(synapse_sizes, idx)
        end

        
        # Transitions from immature to pool
        immature_to_pool = rand(Binomial(immature, e * dt))
        immature -= immature_to_pool
        pool += immature_to_pool

        syn_maturation_functions.kesten_update!(synapse_sizes,ε, η, σ_ε, σ_η)

        push!(pool_history, pool)
        push!(immature_history, immature)
        push!(mature_history, mature)
        push!(synapse_size_history,synapse_sizes)

    end

    push!(immature_total, immature_history)
    push!(mature_total, mature_history)
end

time_walks = collect(0:0.01:100)

# Differential equations
Nis, Nms, synapse_sizes, synapses, solution = syn_maturation_functions.run_simulation_diffeq_weightdependent(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, λ, A, kesten_time_step);
time_array_diffeq = solution.t
immature_population_diffeq = solution[1, :]
mature_population_diffeq = solution[2, :]


trialsplot = plot(time_walks, immature_total, color=:pink, label=false,title="Weight Dependent Multiple Trials", legend=:right)
plot!(time_walks, mature_total, color=:lightblue, label=false)
plot!(time_walks, immature_total[1], label="Immature Synapses", lw=3, color=:pink)
plot!(time_walks, mature_total[1], label="Mature Synapses", lw=3, color=:lightblue)
# plot!(time_diffs, Nis, label = "Immature Synapses (DiffEq)", color="red", lw=3, legend=:right)
# plot!(time_diffs, Nms, label = "Mature Synapses (DiffEq)", color="blue", lw=3, xlabel="Time",ylabel="Population size")
plot!(time_array_diffeq, immature_population_diffeq,label = "Immature Synapses (DiffEq)", color="red", lw=3, legend=:right)
plot!(time_array_diffeq, mature_population_diffeq,label = "Mature Synapses (DiffEq)", color="blue", lw=3, xlabel="Time",ylabel="Population size")


savefig(trialsplot, "C://Users/B00955735/OneDrive - Ulster University/Desktop/trials_plot.png")