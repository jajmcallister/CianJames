module syn_maturation_functions

using Random, Distributions, Plots, StatsBase
using DifferentialEquations, Distributions, Plots



export run_simulation_diffeq, run_simulation_randwalks

# Gillespie step function with pool and restriction
function gillespie_step!(pool, synapses, num_synapses, synapse_sizes, rates)
    c, m, e, i = rates

    # Calculating propensities
    a_c = c * length(pool)
    a_m = m * count(==(0), synapses)
    a_e = e * count(==(0), synapses)
    a_i = i * count(==(1), synapses)
    a_total = a_c + a_m + a_e + a_i

    if a_total == 0
        return num_synapses, Inf
    end

    # Time to the next reaction
    τ = rand(Exponential(1/a_total))

    # Determine which reaction occurs
    r = rand() * a_total

    if r < a_c
        # Creation: Move a synapse from the "resource pool" to the immature state
        append!(synapses, 0)
        pop!(pool)
        num_synapses += 1
    elseif r < a_c + a_m
        # Maturation: An immature synapse becomes mature
        index = findfirst(==(0), synapses) # This takes an immature synapse (i.e. a 0 in the synapses array) and turns it into a 1
        synapses[index] = 1 # ....and turns it into a 1
        push!(synapse_sizes, 0.0)  # Initialize the size of new mature synapse to 0
    elseif r < a_c + a_m + a_e
        # Elimination: Remove an immature synapse
        index = findfirst(==(0), synapses)
        deleteat!(synapses, index)
        # deleteat!(synapse_sizes, index)
        num_synapses -= 1
        # Add synapse back to the pool
        push!(pool, 1)
    else
        # Immature: A mature synapse becomes immature (with the smallest synaptic size)
        index = findfirst(==(1), synapses) # finds the first mature synapse in the synapses list
        synapses[index] = 0 # ... and sets it back to 0
        smallest_idx = findfirst(==(minimum(synapse_sizes)), synapse_sizes) # finds the smallest synapse size
        deleteat!(synapse_sizes, smallest_idx)  # and eliminates it from synapse sizes
    end

    return num_synapses, τ
end

# Kesten process update function
function kesten_update!(sizes, ε, η, σ_ε, σ_η)
    for i in 1:length(sizes)
        ε_i = rand(Normal(ε, σ_ε))
        η_i = rand(Normal(η, σ_η))
        new_size = ε_i * sizes[i] + η_i
        sizes[i] = max(new_size, 0.0)  # Ensure size is non-negative
    end
end

# RandWalks Simulation function
function run_simulation_randwalks(total_time, total_pool_size, synapse_sizes, rates,ε, η, σ_ε, σ_η, kesten_timestep)
    pool = [1 for _ in 1:total_pool_size]  # Initialize pool with synapses
    synapses = Int[]  # Array to hold states of synapses
    num_synapses = 0

    # Initialise tracking variables
    time_array = [0.0]
    immature_population = [0]
    mature_population = [0]

    current_time = 0.0

    while current_time < total_time
        push!(time_array, current_time)
        push!(immature_population, count(==(0), synapses))
        push!(mature_population, count(==(1), synapses))

        num_synapses, τ = gillespie_step!(pool, synapses, num_synapses, synapse_sizes, rates)
        current_time += τ

        # Initialize sizes for new mature synapses
        for i in 1:length(synapses)
            if synapses[i] == 1 && i > length(synapse_sizes)
                push!(synapse_sizes, 0.0)  # Initialize the size of new mature synapses to 0.0
            end
        end

        # Update sizes of mature synapses
        # mature_synapse_indices = findall(==(1), synapses)
        # mature_synapse_sizes = [synapse_sizes[i] for i in mature_synapse_indices]
            
        for j in 0:kesten_timestep:τ
            kesten_update!(synapse_sizes,ε, η, σ_ε, σ_η)
        end

        # # Update sizes in the original list
        # for (idx, new_size) in zip(mature_synapse_indices, new_sizes)
        #     synapse_sizes[idx] = new_size
        # end

    end

    return time_array, immature_population, mature_population, synapse_sizes
end

function synapse_dynamics!(du, u, p, t)
    c, m, e, i = p
    N_I, N_M, P = u

    # N_I is number of immature
    # N_M is number of mature
    # P is number in resource pool
    du[1] = c * P - (m + e) * N_I + i * N_M    # dN_I/dt
    du[2] = m * N_I - i * N_M                  # dN_M/dt
    du[3] = - du[1] - du[2]                  # dP/dt
end


# DiffEq Simulation function with Kesten process
function run_simulation_diffeq(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, kesten_time_step)
    pool = fill(1, total_pool_size);  # Initialize resource pool with synapses
    synapses = Int[]  # Array to hold states of synapses (0s and 1s)
    synapse_sizes = Float64[]  # Sizes of mature synapses

    # Initial conditions
    u0 = [0.0, 0.0, total_pool_size];
    p = (rates..., ε, η);
    tspan = (0.0, total_time);

    # Define ODE problem
    prob = ODEProblem(synapse_dynamics!, u0, tspan, rates);

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

        kesten_update!(synapse_sizes,ε, η, σ_ε, σ_η)

    end

    solution = solve(prob);

    return solution, synapse_sizes, synapses
end


end