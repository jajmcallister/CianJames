using DifferentialEquations, Random

total_pool_size = 100

function run_simulation_diffeq_tracktime(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, kesten_time_step)
    # Initialize the state of the pool, immature, and mature populations
    synapse_sizes = Float64[]  # Sizes of mature synapses
    state_records = Vector{Vector{Int}}(undef, total_pool_size)  # State history records for each synapse

    # Initialize state records for each synapse
    for i in 1:total_pool_size
        state_records[i] = []
    end

    # Initial conditions
    u0 = [0.0, 0.0, total_pool_size]
    p = (rates..., ε, η)
    tspan = (0.0, total_time)

    # Define ODE problem
    prob = ODEProblem(syn_maturation_functions.synapse_dynamics!, u0, tspan, rates)

    current_time = 0.0

    while current_time < total_time
        sol = solve(prob, Tsit5(), saveat=current_time:kesten_time_step:current_time + kesten_time_step)
        N_I, N_M, P = sol.u[end]
        current_time += kesten_time_step

        N_I = max(0, round(Int, N_I))
        N_M = max(0, round(Int, N_M))
        P = max(0, round(Int, P))

        # Create a new state vector for this time step
        new_states = vcat(fill(0, P),fill(1, N_I), fill(2, N_M))
        shuffle!(new_states)  # Shuffle to simulate randomness in state distribution

        # Update the state records for each synapse
        for i in 1:total_pool_size
            if i <= length(new_states)
                state_records[i] = vcat(state_records[i], new_states[i])
            else
                # If the number of synapses in the pool is less than total_pool_size, pad with zeros
                state_records[i] = vcat(state_records[i], 0)
            end
        end

        # Apply Kesten process to mature synapses in N_M
        if N_M > length(synapse_sizes)  # If synapses have been added to the mature population since last time
            new_matures = round(Int, N_M) - length(synapse_sizes)
            append!(synapse_sizes, fill(0.0, new_matures))  # Initialize new mature synapses with size 0.0
        elseif N_M < length(synapse_sizes)  # If synapses have dematured out of the mature population
            num_delete_matures = length(synapse_sizes) - round(Int, N_M)  # Find how many need to be deleted
            synapse_sizes = sort(synapse_sizes)  # Sort the synapse size array
            synapse_sizes = synapse_sizes[num_delete_matures + 1:end]  # Delete the first num_delete_matures
        end

        syn_maturation_functions.kesten_update!(synapse_sizes, ε, η, σ_ε, σ_η)
    end

    # Handle any remaining synapse states at the end of the simulation
    for i in (total_pool_size+1):length(state_records)
        state_records[i] = vcat(state_records[i], 0)  # Pad remaining records with 0
    end

    solution = solve(prob);

    return solution, synapse_sizes, state_records
end


sol, synapse_sizes_diffeq2, state_records = run_simulation_diffeq_tracktime(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, kesten_timestep);

tt = sol.t
immature_population_diff = sol[1, :]
mature_population_diff = sol[2, :]

final_I_value = total_pool_size / (1 + m/i + e/c)
final_M_value = total_pool_size / (1 + i/m + (e*i)/(c*m))

plot(tt, immature_population_diff)
plot!(tt, mature_population_diff)
plot!(tt, immature_population_diff+mature_population_diff)
hline!([final_I_value,final_M_value],label="Steady state solutions", linestyle= :dash,lw=3)



function calculate_state_durations(state_record::Vector{Int})
    state_durations = Dict{Int, Vector{Int}}()  # Dictionary to hold state and a vector of its durations

    current_state = state_record[1]
    current_duration = 0

    for state in state_record
        if state == current_state
            current_duration += 1
        else
            # Store the duration for the previous state
            if !haskey(state_durations, current_state)
                state_durations[current_state] = []
            end
            push!(state_durations[current_state], current_duration)

            # Reset for the new state
            current_state = state
            current_duration = 1
        end
    end

    # Store the duration for the last state
    if !haskey(state_durations, current_state)
        state_durations[current_state] = []
    end
    push!(state_durations[current_state], current_duration)

    return state_durations
end


# Initialize vectors to store durations for each state
durations_pool = Int[]
durations_immature = Int[]
durations_mature = Int[]

# Function to calculate durations for a given state record
function extract_durations(state_record::Vector{Int}, state_durations::Dict{Int, Vector{Int}})
    current_state = state_record[1]
    current_duration = 0

    for state in state_record
        if state == current_state
            current_duration += 1
        else
            # Store the duration for the previous state
            if !haskey(state_durations, current_state)
                state_durations[current_state] = []
            end
            push!(state_durations[current_state], current_duration)

            # Reset for the new state
            current_state = state
            current_duration = 1
        end
    end

    # Store the duration for the last state
    if !haskey(state_durations, current_state)
        state_durations[current_state] = []
    end
    push!(state_durations[current_state], current_duration)
end

# Process each synapse's state record
for record in state_records
    state_durations = Dict{Int, Vector{Int}}()
    extract_durations(record, state_durations)

    # Append durations to the corresponding vectors
    for (state, durations) in state_durations
        if state == 0
            append!(durations_pool, durations)
        elseif state == 1
            append!(durations_immature, durations)
        elseif state == 2
            append!(durations_mature, durations)
        end
    end
end

# Plot the distributions
p0 = histogram(durations_pool, bins=100, alpha=0.5, label="Pool", title="State Durations Distribution", xlabel="Duration", ylabel="Frequency")
p1 = histogram(durations_immature, bins=100, alpha=0.5, label="Immature", title="State Durations Distribution", xlabel="Duration", ylabel="Frequency")
p2 = histogram(durations_mature, bins=20, alpha=0.5, label="Mature", title="State Durations Distribution", xlabel="Duration", ylabel="Frequency")

# Combine the plots into a single figure
plot(p0, p1, p2, layout=(3,1), legend=:bottomright)

maximum(durations_pool)
maximum(durations_immature)
maximum(durations_mature)
mean(durations_pool)
mean(durations_immature)
mean(durations_mature)


#######
#######
# Variables rates
#######
#######

el(t) = 0.1 * exp(-t / 10) + 0.2
cr(t) = 0.2 * exp(-t / 30) + 0.2

elim = el.(1:100)
creat = cr.(1:100)

varrates = plot(elim, label="Elimination rate", ylim=(0,0.5), lw=2)
plot!(creat, label="Creation rate", lw=2, xlabel="time")

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

function kesten_update_var!(sizes, ε, η, σ_ε, σ_η)
    for i in 1:length(sizes)
        ε_i = rand(Normal(ε, σ_ε))
        η_i = rand(Normal(η, σ_η))
        new_size = ε_i * sizes[i] + η_i
        sizes[i] = max(new_size, 0.0)  # Ensure size is non-negative
    end
end


function run_simulation_diffeq_tracktime_var_old(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, kesten_time_step)
    # Initialize the state of the pool, immature, and mature populations
    synapse_sizes = Float64[]  # Sizes of mature synapses
    state_records = Vector{Vector{Int}}(undef, total_pool_size)  # State history records for each synapse

    # Initialize state records for each synapse
    for i in 1:total_pool_size
        state_records[i] = []
    end

    # Initial conditions
    u0 = [0.0, 0.0, total_pool_size]
    p = (m, i, λ, synapse_sizes)
    tspan = (0.0, total_time)

    # Define ODE problem
    prob = ODEProblem(synapse_dynamics_var!, u0, tspan, p);

    current_time = 0.0

    while current_time < total_time
        sol = solve(prob, Tsit5(), saveat=current_time:kesten_time_step:current_time + kesten_time_step)
        N_I, N_M, P = sol.u[end]
        current_time += kesten_time_step

        N_I = max(0, round(Int, N_I))
        N_M = max(0, round(Int, N_M))
        P = max(0, round(Int, P))

        # Create a new state vector for this time step
        new_states = vcat(fill(0, P),fill(1, N_I), fill(2, N_M))
        shuffle!(new_states)  # Shuffle to simulate randomness in state distribution

        # Update the state records for each synapse
        for i in 1:total_pool_size
            if i <= length(new_states)
                state_records[i] = vcat(state_records[i], new_states[i])
            else
                # If the number of synapses in the pool is less than total_pool_size, pad with zeros
                state_records[i] = vcat(state_records[i], 0)
            end
        end

        # Apply Kesten process to mature synapses in N_M
        if N_M > length(synapse_sizes)  # If synapses have been added to the mature population since last time
            new_matures = round(Int, N_M) - length(synapse_sizes)
            append!(synapse_sizes, fill(0.0, new_matures))  # Initialize new mature synapses with size 0.0
        elseif N_M < length(synapse_sizes)  # If synapses have dematured out of the mature population
            num_delete_matures = length(synapse_sizes) - round(Int, N_M)  # Find how many need to be deleted
            synapse_sizes = sort(synapse_sizes)  # Sort the synapse size array
            synapse_sizes = synapse_sizes[num_delete_matures + 1:end]  # Delete the first num_delete_matures
        end

        kesten_update_var!(synapse_sizes, ε, η, σ_ε, σ_η)
    end

    # Handle any remaining synapse states at the end of the simulation
    for i in (total_pool_size+1):length(state_records)
        state_records[i] = vcat(state_records[i], 0)  # Pad remaining records with 0
    end

    solution = solve(prob);

    return solution, synapse_sizes, state_records
end

function run_simulation_diffeq_tracktime_var(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, kesten_time_step)
    # Initialize pool, immature, and mature synapses
    pool_synapses = collect(1:total_pool_size)  # Index of all synapses in the pool
    immature_synapses = Int[]  # Synapses in immature state
    mature_synapses = Int[]  # Synapses in mature state
    state_records = Vector{Vector{Int}}(undef, total_pool_size)  # History records for each synapse
    synapse_sizes = Float64[] 

    # Initialize state records for each synapse (initially all in the pool, state 0)
    for i in 1:total_pool_size
        state_records[i] = fill(0, 0)  # Fill with initial state 0
    end

    # Initial conditions
    u0 = [0.0, 0.0, total_pool_size]
    p = (rates..., ε, η)
    tspan = (0.0, total_time)

    # Define ODE problem
    prob = ODEProblem(synapse_dynamics_var!, u0, tspan, p)

    # Time evolution of the system
    num_steps = round(Int, total_time / kesten_time_step)  # Number of steps
    current_time = 0.0

    for step in 1:num_steps
        # Solve the ODE for this time step
        sol = solve(prob, Tsit5(), saveat=current_time + kesten_time_step)
        N_I, N_M, P = sol.u[end]
        current_time += kesten_time_step

        N_I = max(0, round(Int, N_I))
        N_M = max(0, round(Int, N_M))
        P = max(0, round(Int, P))

        # Track number of synapses transitioning between states
        delta_P_to_I = max(0, N_I - length(immature_synapses))  # Pool -> Immature
        delta_I_to_M = max(0, N_M - length(mature_synapses))    # Immature -> Mature
        delta_M_to_I = max(0, length(mature_synapses) - N_M)    # Mature -> Immature
        delta_I_to_P = max(0, length(immature_synapses) - N_I)  # Immature -> Pool

        # Transition synapses based on the dynamics
        # # Move synapses from pool to immature (state 0 -> state 1)
        # if delta_P_to_I > 0
        #     selected_synapses = sample(pool_synapses, min(delta_P_to_I, length(pool_synapses)))
        #     for idx in selected_synapses
        #         push!(immature_synapses, idx)
        #         deleteat!(pool_synapses, findfirst(==(idx), pool_synapses))
        #     end
        # end

        # # Move synapses from immature to mature (state 1 -> state 2)
        # if delta_I_to_M > 0
        #     selected_synapses = sample(immature_synapses, min(delta_I_to_M, length(immature_synapses)))
        #     for idx in selected_synapses
        #         push!(mature_synapses, idx)
        #         deleteat!(immature_synapses, findfirst(==(idx), immature_synapses))
        #     end
        # end

        # # Move synapses from mature to immature (state 2 -> state 1)
        # if delta_M_to_I > 0
        #     selected_synapses = sample(mature_synapses, min(delta_M_to_I, length(mature_synapses)))
        #     for idx in selected_synapses
        #         push!(immature_synapses, idx)
        #         deleteat!(mature_synapses, findfirst(==(idx), mature_synapses))
        #     end
        # end

        # # Move synapses from immature back to pool (state 1 -> state 0)
        # if delta_I_to_P > 0
        #     selected_synapses = sample(immature_synapses, min(delta_I_to_P, length(immature_synapses)))
        #     for idx in selected_synapses
        #         push!(pool_synapses, idx)
        #         deleteat!(immature_synapses, findfirst(==(idx), immature_synapses))
        #     end
        # end

        # Move synapses from pool to immature (state 0 -> state 1)
        if delta_P_to_I > 0
            selected_synapses = sample(pool_synapses, min(delta_P_to_I, length(pool_synapses)))
            for idx in selected_synapses
                push!(immature_synapses, idx)
                pos = findfirst(==(idx), pool_synapses)
                if pos !== nothing
                    deleteat!(pool_synapses, pos)
                end
                append!(state_records[idx], 1)  # Update state record for immature (1)
            end
        end

        # Move synapses from immature to mature (state 1 -> state 2)
        if delta_I_to_M > 0
            selected_synapses = sample(immature_synapses, min(delta_I_to_M, length(immature_synapses)))
            for idx in selected_synapses
                push!(mature_synapses, idx)
                pos = findfirst(==(idx), immature_synapses)
                if pos !== nothing
                    deleteat!(immature_synapses, pos)
                end
                append!(state_records[idx], 2)  # Update state record for mature (2)
            end
        end

        # Move synapses from mature to immature (state 2 -> state 1)
        if delta_M_to_I > 0
            selected_synapses = sample(mature_synapses, min(delta_M_to_I, length(mature_synapses)))
            for idx in selected_synapses
                push!(immature_synapses, idx)
                pos = findfirst(==(idx), mature_synapses)
                if pos !== nothing
                    deleteat!(mature_synapses, pos)
                end
                append!(state_records[idx], 1)  # Update state record for immature (1)
            end
        end

        # Move synapses from immature back to pool (state 1 -> state 0)
        if delta_I_to_P > 0
            selected_synapses = sample(immature_synapses, min(delta_I_to_P, length(immature_synapses)))
            for idx in selected_synapses
                push!(pool_synapses, idx)
                pos = findfirst(==(idx), immature_synapses)
                if pos !== nothing
                    deleteat!(immature_synapses, pos)
                end
                append!(state_records[idx], 0)  # Update state record for pool (0)
            end
        end

        # At every step, record the state of each synapse
        for idx in pool_synapses
            push!(state_records[idx], 0)  # State 0 (Pool)
        end
        for idx in immature_synapses
            push!(state_records[idx], 1)  # State 1 (Immature)
        end
        for idx in mature_synapses
            push!(state_records[idx], 2)  # State 2 (Mature)
        end

        # Apply Kesten process to mature synapses (update sizes)
        for idx in mature_synapses
            kesten_update_var!(synapse_sizes, ε, η, σ_ε, σ_η)
        end
    end

    solution = solve(prob)


    return solution, synapse_sizes, state_records
end

# Run simulation
sol_var_times, synapse_sizes_var_times, state_records_var = run_simulation_diffeq_tracktime_var(total_time, total_pool_size, params, ε, η, σ_ε, σ_η, kesten_timestep);

state_records_var[1]
plot(sum(state_records_var))
time_array_var_times = sol_var_times.t
immature_population_var_times = sol_var_times[1, :]
mature_population_var_times = sol_var_times[2, :]
pool_tracktimes = sol_var_times[3,:]

final_I_value = total_pool_size / (1 + m/i + el(total_time)/cr(total_time))
final_M_value = total_pool_size / (1 + i/m + (el(total_time)*i)/(cr(total_time)*m))

var_times_plot = plot(time_array_var_times, immature_population_var_times, label = "Immature Synapses", color="red", lw=3, legend=:bottomright)
plot!(time_array_var_times, mature_population_var_times, label = "Mature Synapses", color="blue", lw=3, xlabel="Time",ylabel="Population size")
plot!(time_array_var_times, immature_population_var_times + mature_population_var_times, lw=3, label="Mature+Immature")

hline!([final_I_value,final_M_value],label="Steady state solutions", linestyle= :dash,lw=3)
hline!([(immature_population_var_times+mature_population_var_times)[end]],label=false)