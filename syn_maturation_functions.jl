module syn_maturation_functions

using Random, Distributions, Plots, StatsBase
using DifferentialEquations, Distributions, Plots



export run_simulation_diffeq, run_simulation_randwalks, run_simulation_diffeq_exp, run_simulation_randwalks_exp, synapse_dynamics_weightdependent, run_simulation_diffeq_weightdependent

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
    synapse_size_history = []

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
            syn_sizes = copy(synapse_sizes)
            push!(synapse_size_history, syn_sizes)
        end


    end

    return time_array, immature_population, mature_population, synapse_sizes, synapse_size_history
end

function run_simulation_randwalks_exp(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, kesten_timestep, A, lambda)
    c,m,e,i = rates
    dt = kesten_timestep
    steps = Int(total_time/kesten_timestep)  
    pool = total_pool_size
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
            prob = A*exp(-size / lambda)*dt
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

        # # Delete elements at the specified indices
        for idx in sorted_indices
            deleteat!(synapse_sizes, idx)
        end
        # for i in 1:mature_to_immature
        #     pop!(synapse_sizes)
        # end
        
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
    return immature_history, mature_history, synapse_sizes, synapse_size_history



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


function synapse_dynamics_exp!(du, u, p, t)
    c, m, e, i, λ, synapse_sizes, = p 
    N_I, N_M, P = u
    A = i

    # Compute the rate of dematuration using the exponential probability distribution
    # Sum over all mature synapses' probabilities of transitioning to immature state

    dematuration_rate = A * sum(exp(- size / λ) for size in synapse_sizes) / length(synapse_sizes)

    du[1] = c * P - (m + e) * N_I + (dematuration_rate) * N_M  # dN_I/dt
    du[2] = m * N_I - (dematuration_rate) * N_M  # dN_M/dt
    du[3] = - du[1] - du[2]  # dP/dt
end

function run_simulation_diffeq_exp(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, λ, kesten_time_step)
    pool = fill(1, total_pool_size)  # Initialize resource pool with synapses
    synapses = Int[]  # Array to hold states of synapses (0s and 1s)
    synapse_sizes = Float64[]  # Sizes of mature synapses

    # Initial conditions
    u0 = [0.0, 0.0, total_pool_size]

    p = (rates..., ε, η, λ, copy(synapse_sizes))
    tspan = (0.0, total_time)

    # Define ODE problem
    prob = ODEProblem(synapse_dynamics_exp!, u0, tspan, p)

    current_time = 0.0

    while current_time < total_time

        # p = (rates..., ε, η, λ, copy(synapse_sizes))
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


function synapse_dynamics_weightdependent!(du, u, p, t)
    c, m, e, i, dematuration_rate = p
    N_I, N_M, P = u

    # local num_mature_to_immature = 0
    # for (i, size) in enumerate(synapse_sizes)
    #     prob = A*exp(-size / λ)*0.01
    #     if rand() < prob
    #         num_mature_to_immature+=1
    #     end
    # end

    # dematuration_rate =  num_mature_to_immature/length(synapse_sizes) # mean(A*exp.(-synapse_sizes / λ))

    du[1] = c * P - (m + e) * N_I + dematuration_rate * N_M  # dN_I/dt
    du[2] = m * N_I - (dematuration_rate) * N_M  # dN_M/dt
    du[3] = - du[1] - du[2]  # dP/dt
end



function run_simulation_diffeq_weightdependent(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, λ, A, kesten_time_step)
    pool = fill(1, total_pool_size)  # Initialise resource pool with synapses
    synapses = Int[]  # Array to hold states of synapses (0s and 1s)
    synapse_sizes = Float64[]  # Sizes of mature synapses

    u0 = [0.0, 0.0, total_pool_size]
    tspan = (0.0, total_time)

    current_time = 0.0
    NIs = []
    NMs = []

    while current_time < total_time
        # local num_mature_to_immature = 0
        # for (i, size) in enumerate(synapse_sizes)
        #     prob = A*exp(-size / λ)
        #     if rand() < prob
        #         num_mature_to_immature+=1
        #     end
        # end

        # dematuration_rate =  num_mature_to_immature/length(synapse_sizes)
        if isempty(synapse_sizes)
            dematuration_rate = A
        else
            dematuration_rate = A * sum(exp(- size / λ) for size in synapse_sizes) / length(synapse_sizes)
        end

        c, m, e, i = rates
        p = (c, m, e, i , dematuration_rate)

        prob = ODEProblem(synapse_dynamics_weightdependent!, u0, tspan, p)

        sol = solve(prob, Tsit5(), saveat=current_time:kesten_time_step:current_time + kesten_time_step)
        N_I, N_M, P = sol.u[end]
        push!(NIs, N_I)
        push!(NMs, N_M)
        current_time += kesten_time_step

        # Update populations:  0s for synapses in N_I and 1s for in N_M
        synapses = vcat(fill(0, round(Int, N_I)), fill(1, round(Int, N_M)))
        # 1s for synapses in the pool
        pool = fill(1, round(Int, P))

        # Apply Kesten process to mature synapses in N_M
        if N_M > length(synapse_sizes)  # If synapses have been added to the mature population since last time
            new_matures = round(Int, N_M) - length(synapse_sizes)
            append!(synapse_sizes, fill(0.0, new_matures))  # Initialise new mature synapses with size 0.0
        elseif N_M < length(synapse_sizes)  # If synapses have dematured out of the mature population
            num_delete_matures = length(synapse_sizes) - round(Int, N_M)  # Find how many need to be deleted
            # Remove smallest
            synapse_sizes = sort(synapse_sizes)[num_delete_matures + 1:end]
        end

        # Apply Kesten process
        syn_maturation_functions.kesten_update!(synapse_sizes, ε, η, σ_ε, σ_η)
    end

    if isempty(synapse_sizes)
        dematuration_rate = A
    else
        dematuration_rate = A * sum(exp(- size / λ) for size in synapse_sizes) / length(synapse_sizes)
    end

    c, m, e, i = rates
    p = (c, m, e, i , dematuration_rate)

    prob = ODEProblem(synapse_dynamics_weightdependent!, u0, tspan, p)

    solution = solve(prob)
    return NIs, NMs, synapse_sizes, synapses, solution
end

end