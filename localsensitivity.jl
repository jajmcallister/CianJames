###########
# 1 # 1D Sensitivity Analysis
###########

# Ordered parameters and baseline values
param_keys = [:a1, :k1, :a2, :k2, :m, :A, :λ]
baseline_params = [a1, k1, a2, k2, m, A, lambda]
deltas = 0.9 .* baseline_params

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


plot(h0,h1,layout=(1,2), size=(1500,600),margin=10mm)
























