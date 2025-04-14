using Distributions, Statistics, Plots, ProgressLogging


# Parameters
a1,k1,a2,k2,m,A,λ = 0.9, 0.03333333333333333, 2, 0.17500000000000002, 0.05, 0.05, 2.0
ε, η = 0.985, 1 - 0.985
σ_ε, σ_η = 0.1, 0.1

# Perturbation widths
deltas = 1.0 .* [a1, k1, a2, k2, m, A, λ]

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
    c=:bam, size=(700,500), clim=(-0.7,0.7)
)







###########
# 1 # 1D Sensitivity Analysis
###########

# Ordered parameters and baseline values
param_keys = [:a1, :k1, :a2, :k2, :m, :A, :λ]
baseline_params = [a1, k1, a2, k2, m, A, lambda]
deltas = 0.05 .* baseline_params

deltas

# Compute baseline outputs
base_imm, base_mat = run_model(baseline_params...)
base_outputs = [
    base_mat[end],                # final_mat
    base_imm[end],                # final_imm
    argmax(base_imm .+ base_mat),# time_max_combined
    argmax(base_mat),            # time_max_mat
    argmax(base_imm)             # time_max_imm
]

normalized_base_outputs = base_outputs ./ abs.(base_outputs)


# Initialize sensitivity matrix
n_outputs = length(base_outputs)
n_params = length(param_keys)
sensitivity_matrix = zeros(n_outputs, n_params)

# Loop over parameters
@progress for (i, delta) in enumerate(deltas)
    for sign in [-1.0, 1.0]
        perturbed = copy(baseline_params)
        perturbed[i] += sign * delta
        imm, mat = run_model(perturbed...)
        outputs = [
            mat[end],
            imm[end],
            argmax(imm .+ mat),
            argmax(mat),
            argmax(imm)
        ]
        normalized_outputs = outputs ./ abs.(base_outputs)
        Δoutput = normalized_outputs .- normalized_base_outputs
        for j in 1:n_outputs
            # sensitivity_matrix[j, i] += (outputs[j] - base_outputs[j]) / (2 * delta)
            sensitivity_matrix[j, i] += Δoutput[j] / (2 * delta)
        end
    end
end


# Plot
h1 = heatmap(sensitivity_matrix,
    xticks = (xpoints, xlabs),
    yticks = (1:length(output_keys), ylabs),
    xlabel = "Parameters",
    ylabel = "Outputs",
    title="One-at-a-time Sensitivity", colorbar_title="∂Output / ∂Param",
    c=:bam,size=(700, 500)
)

normalized_matrix = sensitivity_matrix .* (baseline_params' ./ normalized_base_outputs)
h1 = heatmap(normalized_matrix,
    xticks = (xpoints, xlabs),
    yticks = (1:length(output_keys), ylabs),
    xlabel = "Parameters",
    ylabel = "Outputs",
    title="One-at-a-time Sensitivity (normalised)", colorbar_title="∂Output / ∂Param",
    c=:bam,size=(700, 500),clim=(-1,1)
)

plot(h0,h1,layout=(1,2), size=(1500,600),margin=10mm)

##
# multiple sensitivity delta values

# Define multiple percentage deltas
delta_factors = [0.01, 0.03, 0.05]
normalized_base_outputs = base_outputs ./ abs.(base_outputs)
n_outputs = length(base_outputs)
n_params = length(param_keys)
sensitivity_matrix = zeros(n_outputs, n_params)

@progress for α in delta_factors
    deltas = α .* baseline_params
    for (i, delta) in enumerate(deltas)
        for sign in [-1.0, 1.0]
            perturbed = copy(baseline_params)
            perturbed[i] += sign * delta
            imm, mat = run_model(perturbed...)
            outputs = [
                mat[end],
                imm[end],
                argmax(imm .+ mat),
                argmax(mat),
                argmax(imm)
            ]
            normalized_outputs = outputs ./ abs.(base_outputs)
            Δoutput = normalized_outputs .- normalized_base_outputs
            for j in 1:n_outputs
                # sensitivity_matrix[j, i] += (outputs[j] - base_outputs[j]) / (2 * delta)
                sensitivity_matrix[j, i] += Δoutput[j] / (2 * delta)
            end
        end
    end
end

# Average over number of delta_factors and signs (2 directions)
sensitivity_matrix ./= (2 * length(delta_factors))

# Normalize by baseline params and outputs
normalized_matrix = sensitivity_matrix .* (baseline_params' ./ normalized_base_outputs)
h1 = heatmap(normalized_matrix,
    xticks = (xpoints, xlabs),
    yticks = (1:length(output_keys), ylabs),
    xlabel = "Parameters",
    ylabel = "Outputs",
    title="One-at-a-time Sensitivity (normalised)", colorbar_title="∂Output / ∂Param",
    c=:bam,size=(700, 500),clim=(-1,1)
)

plot(h0,h1,layout=(1,2), size=(1500,600),margin=10mm)



##############
# 2 pairwise Hessian  
##############

# Ordered parameters and baseline values
param_keys = [:a1, :k1, :a2, :k2, :m, :A, :λ]
baseline_params = [a1, k1, a2, k2, m, A, lambda]
deltas = 0.25 .* baseline_params

# Output index to analyze
output_labels = ["final_mat", "final_imm", "time_max_combined", "time_max_mat", "time_max_imm"]
output_index = 1  # corresponds to "final_mat"

# Function to extract output
function get_output(imm, mat)
    return [
        mat[end],
        imm[end],
        argmax(imm .+ mat),
        argmax(mat),
        argmax(imm)
    ]
end

# Initialize matrix
n_params = length(param_keys)
interaction_matrix = zeros(n_params, n_params)

# Compute second partial derivatives
@progress for i in 1:n_params
    for j in 1:n_params
        di, dj = deltas[i], deltas[j]

        function run_shift(shift_i, shift_j)
            pert = copy(baseline_params)
            pert[i] += shift_i * di
            pert[j] += shift_j * dj
            imm, mat = run_model(pert...)
            return get_output(imm, mat)[output_index]
        end

        fpp = run_shift(+1, +1)
        fpm = run_shift(+1, -1)
        fmp = run_shift(-1, +1)
        fmm = run_shift(-1, -1)

        interaction_matrix[i, j] = (fpp - fpm - fmp + fmm) / (4 * di * dj)
    end
end

# Plot
x_labels = ["a1", "k1", "a2", "k2", "m", "A", "λ"]
h11 = heatmap(interaction_matrix,
    xticks = (xpoints, xlabs),
    yticks = xticks = (xpoints, xlabs),yflip=true,
    xlabel="Parameter i", ylabel="Parameter j",
    title="Pairwise Interaction on $(output_labels[output_index])",
    colorbar_title="∂²Output / ∂pi∂pj", c=:amp, clim=(-0.05, 0.05), size=(700, 600)
)
