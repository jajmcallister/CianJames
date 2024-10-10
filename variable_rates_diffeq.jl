using DifferentialEquations, Distributions, Plots
using .syn_maturation_functions

function synapse_dynamics_var!(du, u, p, t)
    c_t, m, e_t, i, λ, synapse_sizes = p 
    N_I, N_M, N_P = u
    A = i

    # Compute the rate of dematuration using the exponential probability distribution
    # Sum over all mature synapses' probabilities of transitioning to immature state
    # if !isempty(synapse_sizes)
    #     dematuration_rate = A * sum(exp(- size / λ) for size in synapse_sizes) / length(synapse_sizes)
    # else
    #     dematuration_rate = A
    # end
    dematuration_rate = i
    # Apply time-dependent e(t) and c(t)
    e_t = elimination_func(t)
    c_t = creation_func(t)

    du[1] = c_t * N_P - (m + e_t) * N_I + (dematuration_rate) * N_M  # dN_I/dt
    du[2] = m * N_I - (dematuration_rate) * N_M  # dN_M/dt
    du[3] = - du[1] - du[2]  # dN_P/dt
end

# DiffEq Simulation function with Kesten process
function run_simulation_diffeq_var(total_time, total_pool_size, params, ε, η, σ_ε, σ_η, kesten_time_step)
    pool = fill(1, total_pool_size);  # Initialize resource pool with synapses
    synapses = Int[]  # Array to hold states of synapses (0s and 1s)
    synapse_sizes = Float64[]  # Sizes of mature synapses
    synapse_sizes_history = []
    m, i, λ = params
    # Initial conditions
    u0 = [0.0, 0.0, total_pool_size];
    tspan = (0.0, total_time);
    p = (m, i, λ, synapse_sizes)
    Ihist = []
    Mhist = []

    # Define ODE problem
    # prob = ODEProblem(synapse_dynamics_var!, u0, tspan, p);

    current_time = 0.0;

    prams = creation_func(current_time),m,elimination_func(current_time),i, λ, synapse_sizes
    probb = ODEProblem(synapse_dynamics_var!, u0, tspan, prams);

    while current_time < total_time
        prams = creation_func(current_time),m,elimination_func(current_time),i, λ, synapse_sizes
        probb = ODEProblem(synapse_dynamics_var!, u0, tspan, prams);
        sol = solve(probb, Tsit5(), saveat=current_time:kesten_time_step:current_time + kesten_time_step);
        N_I, N_M, P = sol.u[end];
        push!(Ihist, N_I)
        push!(Mhist,  N_M)
        current_time += kesten_time_step;

        N_M = max(0, N_M)
        N_I = max(0, N_I)
        # Update populations:  0s for synapses in N_I and 1s for in N_M
        synapses = vcat(fill(0, round(Int, N_I)), fill(1, round(Int, N_M)));
        
        P = max(0,P)
        # 1s for synapses in the pool
        pool = fill(1, round(Int, P));

        # Apply Kesten process to mature synapses in N_M
        if N_M > length(synapse_sizes) # If synapses have been added to the mature population since last time
            new_matures = round(Int, N_M) - length(synapse_sizes);
            append!(synapse_sizes, fill(0.0, new_matures));  # Initialize new mature synapses with size 0.0
        elseif N_M < length(synapse_sizes) # If synapses have dematured out of the mature population
            num_delete_matures = length(synapse_sizes) - round(Int, N_M); #find how many need to be deleted
            synapse_sizes = sort(synapse_sizes); # sort the synapse size array
            synapse_sizes = synapse_sizes[num_delete_matures+1:end] # ... and delete the first num_delete_matures
        end

        kesten_update_var!(synapse_sizes,ε, η, σ_ε, σ_η)

        push!(synapse_sizes_history, synapse_sizes)

    end

    solution = solve(probb);

    return solution, synapse_sizes, synapse_sizes_history, synapses, Ihist, Mhist
end

function kesten_update_var!(sizes, ε, η, σ_ε, σ_η)
    for i in 1:length(sizes)
        ε_i = rand(Normal(ε, σ_ε))
        η_i = rand(Normal(η, σ_η))
        new_size = ε_i * sizes[i] + η_i
        sizes[i] = max(new_size, 0.0)  # Ensure size is non-negative
    end
end

# Parameters
total_time = 100.0
total_pool_size = 100

a1 = 0.4
k1 = 1/30
b1 = 0.2
a2 = 1.8
k2 = 1/10
b2 = 0.2

creation_func(t) = a1 * exp(-t * k1) + b1
elimination_func(t) = a2 * exp(-t * k2) + b2

plot(creation_func.(0:1:200), ylim=(0,1), label="creation")
plot!(elimination_func.(0:1:200), label="elimination")


c, m, e, i = 0.2,0.2,0.01,0.05
m, i, λ = 0.2,0.1,2
params=(m, i, λ)

# a1,a2,k1,k2,b1,b2,m,i =0.5, 0.5, 0.5, 0.5, 0.1, 0.1, 0.5, 0.5
params=(m, i, λ)
# Run simulation
sol, synapse_sizes_var, synapse_sizes_history_var, synapses_var, ih, mh = run_simulation_diffeq_var(total_time, total_pool_size, params, ε, η, σ_ε, σ_η, kesten_timestep);

time_array_var = sol.t
immature_population_var = sol[1, :]
mature_population_var = sol[2, :]
poold = sol[3,:]


final_I_value = (i*total_pool_size*(a1+b1))/(m*(a1+b1)+i*(a1+a2+b1+b2))
final_M_value = (total_pool_size*(a1+b1))/(a1+b1+(i/m)*(a1+a2+b1+b2))

var_plot = plot(time_array_var, immature_population_var, label = "Immature Synapses", color="red", lw=1, legend=:bottomright)
plot!(time_array_var, mature_population_var, label = "Mature Synapses", color="blue", lw=1, xlabel="Time",ylabel="Population size")

plot!(time_array_var, immature_population_var+mature_population_var, lw=3, color="green", label="Mature+Immature")

hline!([immature_population_var[end] + mature_population_var[end]], label=false)
hline!([final_I_value,final_M_value],label="Steady state solutions", linestyle= :dash,lw=3)


# savefig(var_plot, "C://Users/B00955735/OneDrive - Ulster University/Desktop/variablerates1.png")




# Work out average i value based on average distribution of synapse sizes
avg_dematuration_rate = i * sum(exp(- size / λ) for size in synapse_sizes_var) / length(synapse_sizes_var)





histogram(synapse_sizes_var,xlim=(0,30), label=false)



h1 = histogram(synapse_sizes_history_var[1],xlim=(0,30), label=false)
h2 = histogram(synapse_sizes_history_var[round(Int,trunc(Int,length(synapse_sizes_history_var)/4))],xlim=(0,30), label=false)
h3 = histogram(synapse_sizes_history_var[3*round(Int,trunc(Int,length(synapse_sizes_history_var)/4))],xlim=(0,30), label=false)
h4 = histogram(synapse_sizes_history_var[end],xlim=(0,30), label=false)
plot(h1,h2,h3,h4, layout=(1,4), size=(1000,500))





# xbins = 0:0.5:100
# # Create an animation object
# anim = @animate for t in 1:10:length(synapse_sizes_history_var)
#     histogram(synapse_sizes_history_var[t], bins=xbins, xlim=(0, 10), ylim=(0, 500), 
#               title="Time = $t", xlabel="Value", ylabel="Frequency",legend=false)
# end

# # Save the animation as a GIF
# hist_gif = gif(anim, "var_histograms.gif", fps=50)


# savefig(hist_gif, "C://Users/B00955735/OneDrive - Ulster University/Desktop/hist_gif.gif")