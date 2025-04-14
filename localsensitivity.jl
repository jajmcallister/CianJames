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
            sensitivity_matrix[j, i] += (outputs[j] - base_outputs[j]) / (2 * delta)
            # sensitivity_matrix[j, i] += Δoutput[j] #/ (2 * delta)
        end
    end
end

sensitivity_matrix
# Plot
h1 = heatmap(sensitivity_matrix,
    xticks = (xpoints, xlabs),
    yticks = (1:length(xpoints), ylabs),
    xlabel = "Parameters",
    ylabel = "Outputs",
    title="One-at-a-time Sensitivity", colorbar_title="∂Output / ∂Param",
    c=:bam,size=(700, 500),clim=(-100,100)
)


~plot(h0,h1,layout=(1,2), size=(1500,600),margin=10mm)

























# 1D Sensitivity Analysis with Array
n_params = length(baseline_params)
n_points = 20
n_outputs = 5

param_ranges = [range(p - d, p + d, length=n_points) for (p, d) in zip(baseline_params, deltas)]
results_1d = Array{Float64}(undef, n_params, n_points, n_outputs)  # (param, sample, output)

@progress for i in 1:n_params
    for j in 1:n_points
        p = copy(baseline_params)
        p[i] = param_ranges[i][j]
        traj_imm, traj_mat = run_model(p...)

        time_max_imm = argmax(traj_imm)
        time_max_mat = argmax(traj_mat)
        time_max_combined = argmax(traj_imm .+ traj_mat)
        final_imm = traj_imm[end]
        final_mat = traj_mat[end]

        results_1d[i, j, :] .= [time_max_imm, time_max_mat, time_max_combined, final_imm, final_mat]
    end
end

param_labels = ["a1", "k1", "a2", "k2", "m", "A", "λ"]
output_labels = ["t_max_imm", "t_max_mat", "t_max_sum", "final_imm", "final_mat"]

n_params, n_points, n_outputs = size(results_1d)
slopes = zeros(n_params, n_outputs)

for i in 1:n_params
    x = collect(param_ranges[i])
    for j in 1:n_outputs
        y = results_1d[i, :, j]
        X = hcat(ones(n_points), x)
        β = X \ y  # Linear regression
        slopes[i, j] = β[2]  # slope
    end
end

# Plot heatmap
heatmap(
    slopes';
    xlabel="Output", ylabel="Parameter",
    title="1D Sensitivity Slopes", color=:balance, cbar=true
)
