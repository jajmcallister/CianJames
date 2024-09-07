# M will have a bump if M/I > m/i

creation_rate(x) = 0.5*exp(-x) + 0.2

xvals = collect(0:0.1:10)

creation = creation_rate.(xvals)
elimination = creation_rate.(xvals)


total_time = 100.0
total_pool_size = 1000
c, m, e, i = 0.2,0.1,0.2,0.05
ε, η = 1.0, 0.0
σ_ε, σ_η = .5, .5
rates = (creation, m, elimination, i)
kesten_timestep = 0.01

function synapse_dynamics_var!(du, u, p, t)
    creat, m, elim, i = p
    N_I, N_M, P = u
    tt = trunc(Int,t)

    c = creat[tt+1]
    e = elim[tt+1]

    du[1] = c * P - (m + e) * N_I + i * N_M    # dN_I/dt
    du[2] = m * N_I - i * N_M                  # dN_M/dt
    du[3] = - du[1] - du[2]                  # dP/dt
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

        # Update populations:  0s for synapses in N_I and 1s for in N_M
        synapses = vcat(fill(0, round(Int, N_I)), fill(1, round(Int, N_M)));
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

final_I_value = total_pool_size / (1 + m/i + e/c)
final_M_value = total_pool_size / (1 + i/m + (e*i)/(c*m))

diffeqplot = plot(time_array_diffeq, immature_population_diffeq, title="Differential equation version of model", label = "Immature Synapses (DiffEq)", color="red", lw=3, legend=:right)
plot!(time_array_diffeq, mature_population_diffeq, label = "Mature Synapses (DiffEq)", color="blue", lw=3, xlabel="Time",ylabel="Population size")
plot!(time_array_diffeq, immature_population_diffeq+mature_population_diffeq, lw=3, label="Mature+Immature")

hline!([final_I_value,final_M_value],label="Steady state solutions", linestyle= :dash,lw=3)
