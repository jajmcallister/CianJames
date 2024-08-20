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


function synapse_dynamics_weightdependent!(du, u, p, t)
    c, m, e, i, λ, A = p
    N_I, N_M, P = u
    A = i

    dematuration_rate = A * mean(exp.(-synapse_sizes / λ))


    du[1] = c * P - (m + e) * N_I + dematuration_rate * N_M  # dN_I/dt
    du[2] = m * N_I - (dematuration_rate) * N_M  # dN_M/dt
    du[3] = - du[1] - du[2]  # dP/dt
end

function run_simulation_diffeq_weightdependent(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, λ, A, kesten_time_step)
    pool = fill(1, total_pool_size)  # Initialize resource pool with synapses
    synapses = Int[]  # Array to hold states of synapses (0s and 1s)
    synapse_sizes = Float64[]  # Sizes of mature synapses

    # Initial conditions
    u0 = [0.0, 0.0, total_pool_size]
    p = (rates..., ε, η, λ, A)  # Include λ for the exponential distribution parameter
    tspan = (0.0, total_time)

    # Define ODE problem
    prob = ODEProblem(synapse_dynamics_weightdependent!, u0, tspan, p)

    current_time = 0.0

    while current_time < total_time
        sol = solve(prob, Tsit5(), saveat=current_time:kesten_time_step:current_time + kesten_time_step)
        N_I, N_M, P = sol.u[end]
        current_time += kesten_time_step

        # Update populations:  0s for synapses in N_I and 1s for in N_M
        synapses = vcat(fill(0, round(Int, N_I)), fill(1, round(Int, N_M)))
        # 1s for synapses in the pool
        pool = fill(1, round(Int, P))

        # Apply Kesten process to mature synapses in N_M
        if N_M > length(synapse_sizes)  # If synapses have been added to the mature population since last time
            new_matures = round(Int, N_M) - length(synapse_sizes)
            append!(synapse_sizes, fill(0.0, new_matures))  # Initialize new mature synapses with size 0.0
        elseif N_M < length(synapse_sizes)  # If synapses have dematured out of the mature population
            num_delete_matures = length(synapse_sizes) - round(Int, N_M)  # Find how many need to be deleted
            # Remove smallest sizes (you could alternatively use other criteria)
            synapse_sizes = sort(synapse_sizes)[num_delete_matures + 1:end]
        end

        # Apply Kesten process
        syn_maturation_functions.kesten_update!(synapse_sizes, ε, η, σ_ε, σ_η)
    end

    solution = solve(prob)

    return solution, synapse_sizes, synapses
end

solution, synapse_sizes, synapses = run_simulation_diffeq_weightdependent(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, λ, A, kesten_time_step);

# Extract the solution
time_array_diffeq = solution.t
immature_population_diffeq = solution[1, :]
mature_population_diffeq = solution[2, :]

diffeqplot = plot(time_array_diffeq, immature_population_diffeq, label = "Immature Synapses (DiffEq)", color="red", lw=3, legend=:right)
plot!(time_array_diffeq, mature_population_diffeq, label = "Mature Synapses (DiffEq)", color="blue", lw=3, xlabel="Time",ylabel="Population size")
