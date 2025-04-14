using Distributions, Statistics, Plots, ProgressLogging


# Parameters
a1,k1,a2,k2,m,A,λ = 0.9, 0.03333333333333333, 2, 0.17500000000000002, 0.05, 0.05, 2.0
ε, η = 0.985, 1 - 0.985
σ_ε, σ_η = 0.1, 0.1

# Perturbation widths
deltas = 0.05 .* [a1, k1, a2, k2, m, A, λ]

# Sample ranges
param_means = [a1, k1, a2, k2, m, A, λ]
param_ranges = [Uniform(p - d, p + d) for (p, d) in zip(param_means, deltas)]

# Sample n parameter sets
function sample_params_array(n)
    [rand.(param_ranges) for _ in 1:n]
end

# Run model for given parameters
function run_model(a1, k1, a2, k2, m, A, λ)
    creation_func(t) = a1 * exp(-t * k1) + b1
    elimination_func(t) = a2 * exp(-t * k2) + b2

    elim = elimination_func.(0:kesten_timestep:total_time)
    creat = creation_func.(0:kesten_timestep:total_time)

    rates_var = creat, m, elim, A, λ
    ihs, mhs = [], []

    for _ in 1:5
        ih, mh, _, _, _, _ = track_times_variable_rates_007(total_time, total_pool_size, rates_var, ε, η, σ_ε, σ_η, kesten_timestep)
        push!(ihs, ih)
        push!(mhs, mh)
    end
    return mean(ihs), mean(mhs)
end

# Simulate
n_samples = 100
param_sets = sample_params_array(n_samples)
results = Matrix{Float64}(undef, n_samples, 5)  # 5 outputs

@progress for (i, param_vec) in enumerate(param_sets)
    a1, k1, a2, k2, m, A, λ = param_vec
    traj_imm, traj_mat = run_model(a1, k1, a2, k2, m, A, λ)

    time_max_imm = argmax(traj_imm)
    time_max_mat = argmax(traj_mat)
    combined = traj_imm .+ traj_mat
    time_max_combined = argmax(combined)
    final_imm = traj_imm[end]
    final_mat = traj_mat[end]

    results[i, :] .= [time_max_imm, time_max_mat, time_max_combined, final_imm, final_mat]
end

# Compute correlation matrix (5 outputs × 7 parameters)
params_array = hcat(param_sets...)'  # (100, 7)
correlation_matrix = [cor(params_array[:, j], results[:, i]) for i in 1:5, j in 1:7]

# Plot
ylabs = ["Time when immature \n population reaches max", "Time when mature \n population reaches max", "Time when combined \n population reaches max", "Final immature \n population value", "Final mature \n population value"]
xlabs = ["\n a1", "Creation", "\n k1", "\n a2","Elimination", "\n k2", "Maturation \n m", "\n A", "De-maturation", "\n λ"]
xpoints = [1,1.5,2,3,3.5,4,5,6,6.5,7]

h0 = heatmap(correlation_matrix,
    xticks = (xpoints, xlabs),
    yticks = (1:5, ylabs),
    xlabel = "Parameters",
    ylabel = "Outputs",
    title = "Parameter-Output Correlations",
    colorbar_title = "Correlation (r)",
    c=:bam, size=(700,500), 
)

