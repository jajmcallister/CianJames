

ε = .9923
η = 1-ε
σ_ε, σ_η = .05, .03


creation_func1(t,A1,lambda1) = 0.5
elimination_func1(t,A2,lambda2) = 0.2
de_maturation_func1(t,A3,lambda3) = 0.2


creat1 = creation_func1.(0:kesten_timestep:total_time, 0.2, 0.2)
elim1 = elimination_func1.(0:kesten_timestep:total_time, 0.2, 0.2)

A = 0.2
lambda31 = 0.2
m1 = 0.2

rates_var1 = creat1, m1, elim1, A, lambda31


state_recs_var_multiple = []
ihs = []
mhs = []
synapse_sizes_multiple = []
syn_size_heatmaps_trials = []
synapse_size_history_multiple = []

num_trials = 10

for i in 1:num_trials
    ih_var, mh_var, state_record_var, syn_sizes_var, syn_heatmap, syn = track_times_variable_rates_007(total_time, total_pool_size, rates_var1, ε, η, σ_ε, σ_η, kesten_timestep);
    push!(state_recs_var_multiple, state_record_var)
    push!(ihs, ih_var)
    push!(mhs, mh_var)
    push!(synapse_sizes_multiple, syn_sizes_var)
    push!(syn_size_heatmaps_trials,syn_heatmap)
    push!(synapse_size_history_multiple, syn)
end


smoothed_avg = time_average007(mean(ihs)+mean(mhs),100)
max_val = maximum(smoothed_avg)
end_val = smoothed_avg[end]
id_max = argmax(smoothed_avg)


p1 = plot!(0:kesten_timestep:total_time,mean(mhs), ribbon=std(mhs)/sqrt(num_trials), lw=3, c=:green, label="Mature synapses", xlabel="Postnatal Day",ylabel="Number")
plot!(0:kesten_timestep:total_time, mean(ihs), ribbon=std(ihs)/sqrt(num_trials), lw=3, c=:magenta, label="Immature synapses")
plot!(0:kesten_timestep:total_time, mean(ihs).+mean(mhs), ribbon=std(ihs)/sqrt(num_trials), lw=3, linealpha=0.7, c=:grey, label="Total synapses")
plot!(title="Population Dynamics", lw=3, c=:black, label="Total synapses",legend=:bottomright)
plot!(grid=false,legendfontsize=12,ylim=(0,100))

#########


function synapse_dynamics_const!(du, u, p, t)
    # Unpack parameters
    A1, lambda1, A2, lambda2, m, A3, lambda3, synapse_sizes = p 
    N_I, N_M, N_P = u

    # Set all rates constant
    c_t = 0.5          # creation rate
    e_t = 0.2          # elimination rate
    dematuration_rate = 0.2  # dematuration rate

    du[1] = c_t * N_P - (m + e_t) * N_I + dematuration_rate * N_M  # dN_I/dt
    du[2] = m * N_I - dematuration_rate * N_M                       # dN_M/dt
    du[3] = -du[1] - du[2]                                          # dN_P/dt
end


function run_simulation_diffeq_var0071(total_time, total_pool_size, paras, ε, η, σ_ε, σ_η, kesten_time_step)
    pool = fill(1, total_pool_size);  # Initialize resource pool with synapses
    synapses = Int[]  # Array to hold states of synapses (0s and 1s)
    synapse_sizes = Float64[]  # Sizes of mature synapses
    synapse_sizes_history = []
    A1, lambda1, A2, lambda2, m, A3, lambda3 = paras
    # A, m, i, λ = paramss
    # Initial conditions
    u0 = [0.0, 0.0, total_pool_size];
    # u0 = [total_pool_size*0.4, total_pool_size*0.3, total_pool_size*0.3];
    tspan = (0.0, total_time);
    # p = (m, i, λ, synapse_sizes)
    Ihist = []
    Mhist = []

    # Define ODE problem
    # prob = ODEProblem(synapse_dynamics_var!, u0, tspan, p);

    current_time = 0.0;
    
    prams = A1, lambda1, A2, lambda2, m, A3, lambda3, synapse_sizes
    probb = ODEProblem(synapse_dynamics_const!, u0, tspan, prams);

    while current_time < total_time
        prams = A1, lambda1, A2, lambda2, m, A3, lambda3, synapse_sizes
        probb = ODEProblem(synapse_dynamics_const!, u0, tspan, prams);
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

        
        if N_M > length(synapse_sizes) # If synapses have been added to the mature population since last time
            new_matures = round(Int, N_M) - length(synapse_sizes);
            append!(synapse_sizes, fill(0.01, new_matures));  # Initialize new mature synapses with size 0.01
        # elseif N_M < length(synapse_sizes) # If synapses have dematured out of the mature population
        #     num_delete_matures = length(synapse_sizes) - round(Int, N_M); #find how many need to be deleted
        #     synapse_sizes = sort(synapse_sizes); # sort the synapse size array
        #     synapse_sizes = synapse_sizes[num_delete_matures+1:end] # ... and delete the first num_delete_matures
        # end

        elseif N_M < length(synapse_sizes) # If synapses have dematured out of the mature population
            num_delete_matures = length(synapse_sizes) - round(Int, N_M) #find how many need to be deleted
            sizes = synapse_sizes
            weights = A3 .* sizes #* exp.(-sizes ./ lambda3)
            weights ./= sum(weights)  # Normalise to sum to 1
        
            # Sample indices to delete based on weights, without replacement
            delete_indices = sample(1:length(sizes), Weights(weights), num_delete_matures; replace=false)
        
            # Remove selected synapses from the synapse_sizes array before applying Kesten process
            synapse_sizes = deleteat!(synapse_sizes, sort(delete_indices))
        end

        # Apply Kesten process to mature synapses in N_M
        synapse_sizes = kesten_update_new007(synapse_sizes,ε, η, σ_ε, σ_η)

        push!(synapse_sizes_history, synapse_sizes)

    end

    solution = solve(probb);

    return solution, synapse_sizes, synapse_sizes_history, synapses, Ihist, Mhist
end


parameters = (0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2)

sol, synapse_sizes_var, synapse_sizes_history_var, synapses_var, ih, mh = run_simulation_diffeq_var0071(total_time, total_pool_size, parameters, ε, η, σ_ε, σ_η, kesten_timestep);

time_array_var = sol.t
immature_population_var = sol[1, :]
mature_population_var = sol[2, :]
poold = sol[3,:]



p2 = plot(time_array_var, mature_population_var, lw=5, c=:green, label="Mature synapses", xlabel="Postnatal Day",ylabel="Number")
plot!(time_array_var, immature_population_var, lw=5, c=:magenta, label="Immature synapses")
plot!(time_array_var, immature_population_var.+mature_population_var, lw=5, c=:grey, linealpha=0.7,label="Total synapses")
plot!(title="Population Dynamics",legend=:bottomright)
plot!(grid=false,ylim=(0,100),legendfontsize=12, dpi=600)

plot!(0:kesten_timestep:total_time,mean(mhs), ribbon=std(mhs)/sqrt(num_trials), lw=3, c=:green, label=false)
plot!(0:kesten_timestep:total_time, mean(ihs), ribbon=std(ihs)/sqrt(num_trials), lw=3, c=:magenta, label=false)
plot!(0:kesten_timestep:total_time, mean(ihs).+mean(mhs), ribbon=std(ihs)/sqrt(num_trials), lw=3, linealpha=0.7, c=:grey, label=false)
plot!(title="Population Dynamics", lw=3, c=:black, label="Total synapses",legend=:bottomright)
plot!(grid=false,legendfontsize=12,ylim=(0,100))

savefig(p2, "C://Users/B00955735/OneDrive - Ulster University/Desktop/combined.png")