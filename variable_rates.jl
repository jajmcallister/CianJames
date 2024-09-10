
# Parameters
total_time = 200.0
total_pool_size = 1000



el(t) = 0.1 * exp(-t / 10) + 0.2
cr(t) = 0.2 * exp(-t / 30) + 0.2

elim = el.(1:100)
creat = cr.(1:100)

varrates = plot(elim, label="Elimination rate", ylim=(0,0.5), lw=2)
plot!(creat, label="Creation rate", lw=2, xlabel="time")

savefig(varrates, "C://Users/B00955735/OneDrive - Ulster University/Desktop/varrates.png")

m = 0.05
i = 0.03

function synapse_dynamics_var!(du, u, p, t)
    m, i, λ, synapse_sizes = p 
    N_I, N_M, N_P = u
    A = i

    # Compute the rate of dematuration using the exponential probability distribution
    # Sum over all mature synapses' probabilities of transitioning to immature state
    if !isempty(synapse_sizes)
        dematuration_rate = A * sum(exp(- size / λ) for size in synapse_sizes) / length(synapse_sizes)
    else
        dematuration_rate = A
    end
    
    # Apply time-dependent e(t) and c(t)
    e_t = el(t)
    c_t = cr(t)
    # m = dematuration_rate

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

    # Define ODE problem
    prob = ODEProblem(synapse_dynamics_var!, u0, tspan, p);

    current_time = 0.0;

    while current_time < total_time

        sol = solve(prob, Tsit5(), saveat=current_time:kesten_time_step:current_time + kesten_time_step);
        N_I, N_M, P = sol.u[end];
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

    solution = solve(prob);

    return solution, synapse_sizes, synapse_sizes_history, synapses
end

function kesten_update_var!(sizes, ε, η, σ_ε, σ_η)
    for i in 1:length(sizes)
        ε_i = rand(Normal(ε, σ_ε))
        η_i = rand(Normal(η, σ_η))
        new_size = ε_i * sizes[i] + η_i
        sizes[i] = max(new_size, 0.0)  # Ensure size is non-negative
    end
end


m, i, λ = 0.1,0.05,2
params=(m, i, λ)

# Run simulation
sol, synapse_sizes_var, synapse_sizes_history_var, synapses_var = run_simulation_diffeq_var(total_time, total_pool_size, params, ε, η, σ_ε, σ_η, kesten_timestep);

time_array_var = sol.t
immature_population_var = sol[1, :]
mature_population_var = sol[2, :]
poold = sol[3,:]

final_I_value = total_pool_size / (1 + m/i + el(total_time)/cr(total_time))
final_M_value = total_pool_size / (1 + i/m + (el(total_time)*i)/(cr(total_time)*m))

var_plot = plot(time_array_var, immature_population_var, label = "Immature Synapses", color="red", lw=3, legend=:bottomright)
plot!(time_array_var, mature_population_var, label = "Mature Synapses", color="blue", lw=3, xlabel="Time",ylabel="Population size")
plot!(time_array_var, immature_population_var+mature_population_var, lw=3, label="Mature+Immature")

hline!([final_I_value,final_M_value],label="Steady state solutions", linestyle= :dash,lw=3)
hline!([(immature_population_var+mature_population_var)[end]],label=false)

savefig(var_plot, "C://Users/B00955735/OneDrive - Ulster University/Desktop/variablerates1.png")


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