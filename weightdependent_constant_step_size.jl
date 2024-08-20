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

# Exponential parameter λ
A = i
λ = 2

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

time_walks = collect(0:0.01:100)

weight_dep_plot = plot(time_walks, immature_history, lw=3, label="Immature population")
plot!(time_walks, mature_history,lw=3, label="Mature population", legend=:right)
plot!(time_array_diffeq, immature_population_diffeq, label = "Immature Synapses (DiffEq)", color="red", lw=3, legend=:right)
plot!(time_array_diffeq, mature_population_diffeq, label = "Mature Synapses (DiffEq)", color="blue", lw=3, xlabel="Time",ylabel="Population size")



weight_dep_hist = histogram(synapse_sizes,xlim=(0,30),label=false,title="Final distribution of synapse sizes", xlabel="Synapse size", ylabel="Frequency")


# savefig(weight_dep_plot, "C://Users/B00955735/OneDrive - Ulster University/Desktop/weight_dep_plot_avg.svg")








xbins = 0:0.5:100
# Create an animation object
anim = @animate for t in 1:10000
    histogram(synapse_size_history[t], bins=xbins, xlim=(0, 10), ylim=(0, 500), 
              title="Time = $t", xlabel="Value", ylabel="Frequency",legend=false)
end

# Save the animation as a GIF
hist_gif = gif(anim, "histograms.gif", fps=50)


savefig(hist_gif, "C://Users/B00955735/OneDrive - Ulster University/Desktop/hist_gif.gif")