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
    A1,lambda1,A2,lambda2,m,A3,lambda3 = 0.9, 30, 2, 5, 0.05, 0.05, 2.
    ε, η = 0.985, 1 - 0.985
    σ_ε, σ_η = 0.1, 0.1
    
    # Constants required for the model
    b1 = 0.2
    b2 = 0.2
    total_time = 100.0
    total_pool_size = 100
    kesten_timestep = 0.1

    
    # Model function that will be evaluated in parallel
    function model_func(p)
        A1,lambda1,A2,lambda2,m,A3,lambda3 = p
        
        # Define creation and elimination functions
        creation_func(t) = A1 * exp(-t / lambda1) + b1
        elimination_func(t) = A2 * exp(-t / lambda2) + b2
    
        # Calculate rates over time
        elim = elimination_func.(0:kesten_timestep:total_time)
        creat = creation_func.(0:kesten_timestep:total_time)
    
        rates_var = creat, m, elim, A3, lambda3
        
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
param_means = [A1, lambda1, A2, lambda2, m, A3, lambda3]
deltas = 10.0 .* param_means

for i in 1:7
    param_bounds[i, 1] = max(0.000001, param_means[i] - deltas[i])
    param_bounds[i, 2] = param_means[i] + deltas[i]
end

# Set up Sobol analysis with parallel processing
sampler = SobolSample()


n_samples = 10

# Generate design matrices
A, B = QuasiMonteCarlo.generate_design_matrices(n_samples, param_bounds[:,1], param_bounds[:,2], sampler)

# Run the analysis with parallel=true
sobol_result = gsa(model_func, Sobol(), A, B, batch=false, parallel=true)

ylabs = ["Time when immature \n population reaches max", "Time when mature \n population reaches max", "Time when combined \n population reaches max", "Final immature \n population value", "Final mature \n population value"]
xlabs = ["\n A1", "Creation", "\n λ1", "\n A2","Elimination", "\n λ2", "Maturation \n m", "\n A3", "De-maturation", "\n λ3"]
xpoints = [1,1.5,2,3,3.5,4,5,6,6.5,7]

ss1 = sobol_result.S1
sst = sobol_result.ST
# Parameter and output names for plotting
param_names = ["A1", "λ1", "A2", "λ2", "m", "A3", "λ3"]
output_names = [
    "Time to max immature", 
    "Time to max mature", 
    "Time to max combined", 
    "Final immature", 
    "Final mature"
]

using Plots.PlotMeasures
# Create heatmap for first-order indices
p1 = heatmap(ss1,
    title="First-Order Sobol Indices",
    xticks = (xpoints, xlabs),
    yticks = (1:5, ylabs),
    xlabel = "Parameters",
    ylabel="Outputs",
    colorbar_title="\n First-Order Sensitivity Index",
    c=:viridis,
    clim=(0, max(maximum(ss1), 0.01)),  # Set colorbar limits with minimum of 0.01 to avoid issues with all zeros
    size=(800, 600),rightmargin=5mm
)

# savefig(p1, "C://Users/B00955735/OneDrive - Ulster University/Desktop/sobol_indices.png")

# Create heatmap for total-order indices
# p2 = heatmap(sst,
#     xticks=(1:7, param_names),
#     yticks=(1:5, output_names),
#     title="Total-Order Sobol Indices",
#     xlabel="Parameters",
#     ylabel="Outputs",
#     colorbar_title="Sensitivity Index",
#     c=:viridis,
#     clim=(0, max(maximum(sst), 0.01)),
#     size=(600, 400)
# )

# plot(h0,p1,size=(1500,600))


##################


# Standard regression GSA
#           A1      lambda1   A2   lambda2     m       A3        lambda3
bounds = [[0, 10], [0, 50], [0, 10], [0, 10], [0, 10], [0, 10], [0, 10]]
reg_sens = gsa(model_func, RegressionGSA(true), bounds, samples = 100)

rs1 = reg_sens.standard_regression
cc = maximum(abs.(rs1))
rs1h = heatmap(rs1,
    xticks = (xpoints, xlabs),
    yticks = (1:5, ylabs),
    xlabel = "Parameters",
    ylabel = "Outputs",
    title = "Regression Global Sensitivity Analysis",
    colorbar_title = "\n Standardised Regression Coefficient",rightmargin=5mm,
    c=:bam, size=(800,600), clim=(-cc,cc)
)
# savefig(rs1h, "C://Users/B00955735/OneDrive - Ulster University/Desktop/regression_sensitivity.png")

# plot(rs1h,rs2h, layout=(1,2),size=(2000,600))


# rs2 = reg_sens.pearson
# cc2 = maximum(abs.(rs2))
# rs2h = heatmap(rs1,
#     xticks = (xpoints, xlabs),
#     yticks = (1:5, ylabs),
#     xlabel = "Parameters",
#     ylabel = "Outputs",
#     title = "Sensitivity Analysis - Pearson",
#     colorbar_title = "Correlation",
#     c=:bam, size=(700,500), clim=(-cc2,cc2),
#     fontfamily="Computer Modern",
# )


pp = plot(rs1h, p1, layout=(2,1), size=(800,1000), leftmargin=10mm)

# savefig(pp, "C://Users/B00955735/OneDrive - Ulster University/Desktop/sensitivity_analysis.png")