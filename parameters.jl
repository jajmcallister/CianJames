using Distributions, ProgressLogging


a1,k1,a2,k2,m,A,lambda = 0.9, 0.03333333333333333, 2, 0.17500000000000002, 0.05, 0.05, 2.
ε = .985
η = 1-ε
σ_ε, σ_η = .1, .1

δ1 = 1.0*a1
δ2 = 1.0*k1
δ3 = 1.0*a2
δ4 = 1.0*k2
δ5 = 1.0*m
δ6 = 1.0*A
δ7 = 1.0*lambda

# Define parameter ranges (adjust as needed)
param_ranges = Dict(
    :a1 => Uniform(a1 - δ1, a1 + δ1),
    :k1 => Uniform(k1 - δ2, k1 + δ2),
    :a2 => Uniform(a2 - δ3, a2 + δ3),
    :k2 => Uniform(k2 - δ4, k2 + δ4),
    :m  => Uniform(m - δ5, m + δ5),
    :A  => Uniform(A - δ6, A + δ6),
    :λ  => Uniform(lambda - δ7, lambda + δ7)
)

# Sample n parameter sets
function sample_params(n)
    [Dict(k => rand(v) for (k, v) in param_ranges) for _ in 1:n]
end

function run_model(a1, k1, a2, k2, m, A, lambda)
    creation_func(t) = a1 * exp(-t * k1) + b1
    elimination_func(t) = a2 * exp(-t * k2) + b2
    
    elim = elimination_func.(0:kesten_timestep:total_time)
    creat = creation_func.(0:kesten_timestep:total_time)

    rates_var = creat, m, elim, A, lambda

    state_recs_var_multiple = []
    ihs = []
    mhs = []
    synapse_sizes_multiple = []
    syn_size_heatmaps_trials = []
    synapse_size_history_multiple = []

    num_trials = 5

    for i in 1:num_trials
        ih_var, mh_var, state_record_var, syn_sizes_var, syn_heatmap, syn = track_times_variable_rates_007(total_time, total_pool_size, rates_var, ε, η, σ_ε, σ_η, kesten_timestep);
        push!(state_recs_var_multiple, state_record_var)
        push!(ihs, ih_var)
        push!(mhs, mh_var)
        push!(synapse_sizes_multiple, syn_sizes_var)
        push!(syn_size_heatmaps_trials,syn_heatmap)
        push!(synapse_size_history_multiple, syn)
    end
    return mean(ihs), mean(mhs)
end



n_samples = 100
param_sets = sample_params(n_samples)
param_sets
results = []

@progress for params in param_sets
    # Update parameters
    a1, k1, a2, k2, m, A, lambda = params[:a1], params[:k1], params[:a2], params[:k2], params[:m], params[:A], params[:λ]
    
    # Run your model here and collect outputs
    traj_imm, traj_mat = run_model(a1, k1, a2, k2, m, A, lambda)  # Assuming you have a model function

    # Find the time point of the maximum for each trajectory
    time_max_imm = argmax(traj_imm)  # Time when immature reaches max
    time_max_mat = argmax(traj_mat)  # Time when mature reaches max
    
    # Combined trajectory and its maximum time
    combined_trajectory = traj_imm .+ traj_mat
    time_max_combined = argmax(combined_trajectory) 

    push!(results, (
        params = params,
        time_max_imm = time_max_imm,
        time_max_mat = time_max_mat,
        time_max_combined = time_max_combined,
        final_imm = traj_imm[end],
        final_mat = traj_mat[end],
    ))
end



# making plots of parameters vs outputs for all parameters and outputs
# Define parameters and outputs to analyze
param_keys = [:a1, :k1, :a2, :k2, :m, :A, :λ]
output_keys = [:time_max_imm, :time_max_mat, :time_max_combined, :final_imm, :final_mat]



# Compute correlation matrix
correlation_matrix = [cor([r.params[p] for r in results], [r[o] for r in results])
                      for o in output_keys, p in param_keys]

# Plot heatmap

ylabs = ["Time when immature \n population reaches max", "Time when mature \n population reaches max", "Time when combined \n population reaches max", "Final immature \n population value", "Final mature \n population value"]
xlabs = ["\n a1", "Creation", "\n k1", "\n a2","Elimination", "\n k2", "Maturation \n m", "\n A", "De-maturation", "\n λ"]
xpoints = [1,1.5,2,3,3.5,4,5,6,6.5,7]
h1 = heatmap(correlation_matrix,
    xticks = (xpoints, xlabs),
    yticks = (1:length(output_keys), ylabs),
    xlabel = "Parameters",
    ylabel = "Outputs",
    title = "Parameter-Output Correlations",
    colorbar_title = "Correlation (r)", c=:bam, size=(700,500),clim=(-0.7,0.7)
)

# savefig(h1, "C://Users/B00955735/OneDrive - Ulster University/Desktop/heatmap_correlation2.png")



ylabs1 = ["Time when immature \n population reaches max", "Time when mature \n population reaches max", "Time when combined \n population reaches max"]
xlabs = ["\n a1", "Creation", "\n k1", "\n a2","Elimination", "\n k2", "Maturation \n m", "\n A", "De-maturation", "\n λ"]
xpoints = [1,1.5,2,3,3.5,4,5,6,6.5,7]
h1 = heatmap(correlation_matrix[1:3,:],
    xticks = (xpoints, xlabs),
    yticks = (1:3, ylabs1),
    xlabel = "Parameters",
    ylabel = "Outputs",
    title = "Parameter-Output Correlations",
    colorbar_title = "Correlation (r)", c=:bam, size=(700,500),clim=(-0.5,0.5)
)











scatt_plots = []

# Loop through each parameter-output pair
for (i, p) in enumerate(param_keys)
    for (j, o) in enumerate(output_keys)
        # Fetch correlation value from the precomputed matrix
        corr_value = correlation_matrix[j, i]

        # Create scatter plot with color based on correlation
        push!(scatt_plots, scatter(
            [r.params[p] for r in results],  # X-axis: parameter value
            [r[o] for r in results],        # Y-axis: output measure
            marker_z = fill(corr_value, length(results)),  # Color by correlation
            color = :bam,               # Use viridis colormap
            xlabel = string(p),
            ylabel = string(o),
            title = "$p vs $o (r=$(round(corr_value, digits=2)))",
            legend = false
        ))
    end
end

using Plots.PlotMeasures
# Display all plots in a grid layout
plot(scatt_plots..., layout= (length(param_keys),length(output_keys)), margin=10mm, size=(3000, 3000),markerstrokewidth=0.01, clim = (-0.4, 0.3))
