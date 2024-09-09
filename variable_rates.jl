el(t) = 0.1 * exp(-t / 10) + 0.2
cr(t) = 0.2 * exp(-t / 30) + 0.2

elim = el.(1:100)
creat = cr.(1:100)

plot(elim, label="Elimination rate", ylim=(0,0.5), lw=2)
plot!(creat, label="Creation rate", lw=2)

m = 0.05
i = 0.5


rates = (c, m, e, i)

function synapse_dynamics_var!(du, u, p, t)
    N_I, N_M, N_P = u
    m, i = p  # Parameters that are not time-dependent
    
    # Apply time-dependent e(t) and c(t)
    e_t = el(t)
    c_t = cr(t)
    
    du[1] = c_t * N_P - (e_t + m + c_t) * N_I + (i - c_t) * N_M  # dN_I/dt
    du[2] = m * N_I - i * N_M                                    # dN_M/dt
    du[3] = - du[1] - du[2]                               # dN_P/dt
end


# DiffEq Simulation function with Kesten process
function run_simulation_diffeq_var(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, kesten_time_step)
    pool = fill(1, total_pool_size);  # Initialize resource pool with synapses
    synapses = Int[]  # Array to hold states of synapses (0s and 1s)
    synapse_sizes = Float64[]  # Sizes of mature synapses

    # Initial conditions
    u0 = [0.0, 0.0, total_pool_size];
    p = (rates..., ε, η);
    tspan = (0.0, total_time);

    # Define ODE problem
    prob = ODEProblem(synapse_dynamics_var!, u0, tspan, rates);

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

    end

    solution = solve(prob);

    return solution, synapse_sizes, synapses
end

function kesten_update_var!(sizes, ε, η, σ_ε, σ_η)
    for i in 1:length(sizes)
        ε_i = rand(Normal(ε, σ_ε))
        η_i = rand(Normal(η, σ_η))
        new_size = ε_i * sizes[i] + η_i
        sizes[i] = max(new_size, 0.0)  # Ensure size is non-negative
    end
end

# Run simulation
sol, synapse_sizes_diffeq, synapses_diffeq = run_simulation_diffeq_var(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, kesten_timestep);

time_array_diffeq = sol.t
immature_population_diffeq = sol[1, :]
mature_population_diffeq = sol[2, :]
poold = sol[3,:]

final_I_value = total_pool_size / (1 + m/i + el(total_time)/cr(total_time))
final_M_value = total_pool_size / (1 + i/m + (el(total_time)*i)/(cr(total_time)*m))

diffeqplot = plot(time_array_diffeq, immature_population_diffeq, label = "Immature Synapses", color="red", lw=3, legend=:bottomright)
plot!(time_array_diffeq, mature_population_diffeq, label = "Mature Synapses", color="blue", lw=3, xlabel="Time",ylabel="Population size")
plot!(time_array_diffeq, immature_population_diffeq+mature_population_diffeq, lw=3, label="Mature+Immature")

(immature_population_diffeq+mature_population_diffeq)[end]

# hline!([final_I_value,final_M_value],label="Steady state solutions", linestyle= :dash,lw=3)
hline!([(immature_population_diffeq+mature_population_diffeq)[end]],label=false)

# savefig(diffeqplot, "C://Users/B00955735/OneDrive - Ulster University/Desktop/variablerates.png")


histogram(synapse_sizes_diffeq,xlim=(0,30), label=false)