using Distributed
using GlobalSensitivity
using QuasiMonteCarlo
using Statistics
using Plots

# Add worker processes if not already added
# This will use all available cores on your machine
if nprocs() == 1
    addprocs(Sys.CPU_THREADS - 1)
    @everywhere using GlobalSensitivity, Statistics, Distributions
end

# Make sure all workers have the necessary functions and parameters
@everywhere begin
    # Parameters
    a1, k1, a2, k2, m, A, λ = 0.9, 0.03333333333333333, 2, 0.17500000000000002, 0.05, 0.05, 2.0
    ε, η = 0.985, 1 - 0.985
    σ_ε, σ_η = 0.1, 0.1
    
    # Constants required for the model
    b1 = 0.1
    b2 = 0.1
    total_time = 100.0
    total_pool_size = 100
    kesten_timestep = 0.1

    
    # Model function that will be evaluated in parallel
    function model_func(p)
        a1, k1, a2, k2, m, A, λ = p
        
        # Define creation and elimination functions
        creation_func(t) = a1 * exp(-t * k1) + b1
        elimination_func(t) = a2 * exp(-t * k2) + b2
    
        # Calculate rates over time
        elim = elimination_func.(0:kesten_timestep:total_time)
        creat = creation_func.(0:kesten_timestep:total_time)
    
        rates_var = creat, m, elim, A, λ
        
        # Run multiple simulations and average the results
        ihs, mhs = [], []
        for _ in 1:3  # Reduced for faster calculation
            ih, mh, _, _, _, _ = track_times_variable_rates_007(total_time, total_pool_size, rates_var, ε, η, σ_ε, σ_η, kesten_timestep)
            push!(ihs, ih)
            push!(mhs, mh)
        end
        
        ih_avg = mean(ihs)
        mh_avg = mean(mhs)
        
        # Calculate outputs
        time_max_imm = argmax(ih_avg) * kesten_timestep
        time_max_mat = argmax(mh_avg) * kesten_timestep
        combined = ih_avg .+ mh_avg
        time_max_combined = argmax(combined) * kesten_timestep
        final_imm = ih_avg[end]
        final_mat = mh_avg[end]
        
        return [time_max_imm, time_max_mat, time_max_combined, final_imm, final_mat]
    end
end

# Parameter bounds for sensitivity analysis
param_bounds = Matrix{Float64}(undef, 7, 2)
param_means = [a1, k1, a2, k2, m, A, λ]
deltas = 1.0 .* param_means

for i in 1:7
    param_bounds[i, 1] = max(0.001, param_means[i] - deltas[i])
    param_bounds[i, 2] = param_means[i] + deltas[i]
end

# Set up Sobol analysis with parallel processing
println("Running parallelized Sobol sensitivity analysis with $(nprocs()-1) worker processes...")
sampler = SobolSample()

# You can reduce this for faster results if needed
n_samples = 100  # Reduced from 1000 to 500 for faster computation

# Generate design matrices
A, B = QuasiMonteCarlo.generate_design_matrices(n_samples, param_bounds[:,1], param_bounds[:,2], sampler)

# Run the analysis with parallel=true
sobol_result = gsa(model_func, Sobol(), A, B, batch=false, parallel=true)

# The rest of your code for visualization remains the same

ss1 = sobol_result.S1
sst = sobol_result.ST
# Parameter and output names for plotting
param_names = ["a1", "k1", "a2", "k2", "m", "A", "λ"]
output_names = [
    "Time to max immature", 
    "Time to max mature", 
    "Time to max combined", 
    "Final immature", 
    "Final mature"
]

# Create heatmap for first-order indices
p1 = heatmap(ss1,
    title="First-Order Sobol Indices",
    xticks = (xpoints, xlabs),
    yticks = (1:5, ylabs),
    xlabel = "Parameters",
    colorbar_title="Sensitivity Index",
    c=:viridis,
    clim=(0, max(maximum(ss1), 0.01)),  # Set colorbar limits with minimum of 0.01 to avoid issues with all zeros
    size=(600, 400)
)

# Create heatmap for total-order indices
p2 = heatmap(sst,
    xticks=(1:7, param_names),
    yticks=(1:5, output_names),
    title="Total-Order Sobol Indices",
    xlabel="Parameters",
    ylabel="Outputs",
    colorbar_title="Sensitivity Index",
    c=:viridis,
    clim=(0, max(maximum(sst), 0.01)),
    size=(600, 400)
)

plot(h0,p1,size=(1500,600))


##################


# Standard regression GSA
bounds = [[0, 10], [0, 10], [0, 10], [0, 10], [0, 10], [0, 10], [0, 10]]
reg_sens = gsa(model_func, RegressionGSA(true), bounds, samples = 1000)

rs1 = reg_sens.standard_regression
cc = maximum(abs.(rs1))
rs1h = heatmap(rs1,
    xticks = (xpoints, xlabs),
    yticks = (1:5, ylabs),
    xlabel = "Parameters",
    ylabel = "Outputs",
    title = "Sensitivity Analysis - Standard Regression",
    colorbar_title = "",
    c=:bam, size=(700,500), clim=(-cc,cc)
)

plot(h0,p1,rs1h,layout=(1,3),size=(2500,600))


rs2 = reg_sens.pearson
cc2 = maximum(abs.(rs2))
rs2h = heatmap(rs1,
    xticks = (xpoints, xlabs),
    yticks = (1:5, ylabs),
    xlabel = "Parameters",
    ylabel = "Outputs",
    title = "Sensitivity Analysis - Pearson",
    colorbar_title = "Correlation",
    c=:bam, size=(700,500), clim=(-cc2,cc2)
)


pp = plot(rs2h, p1, layout=(1,2), size=(1700,600), margin=10mm)

savefig(pp, "C://Users/B00955735/OneDrive - Ulster University/Desktop/sensitivity_analysis.png")