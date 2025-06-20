using Distributed
using GlobalSensitivity
using QuasiMonteCarlo
using Statistics
using Plots

if nprocs() == 1
    addprocs(Sys.CPU_THREADS - 1)
    @everywhere using GlobalSensitivity, Statistics, Distributions
end

# Parameters
A1,lambda1,A2,lambda2,m,A3,lambda3 = 0.9, 30, 2, 5, 0.05, 0.05, 2.
ε, η = 0.985, 1 - 0.985
σ_ε, σ_η = 0.1, 0.1

# Constants required for the model
b1 = 0.2
b2 = 0.2
total_time = 120.0
total_pool_size = 100
kesten_timestep = 0.1

# Make sure all workers have the necessary functions and parameters
@everywhere begin

    # Model function that will be evaluated in parallel
    function model_func(p)
        A1,lambda1,A2,lambda2,m,A3,lambda3 = p
        
        # # Define creation and elimination functions
        # creation_func(t) = A1 * exp(-t / lambda1) + b1
        # elimination_func(t) = A2 * exp(-t / lambda2) + b2
    
        # # Calculate rates over time
        # elim = elimination_func.(0:kesten_timestep:total_time)
        # creat = creation_func.(0:kesten_timestep:total_time)
    
        # rates_var = creat, m, elim, A3, lambda3
        
        # # Run multiple simulations and average the results
        # ihs, mhs = [], []
        # for _ in 1:50
        #     ih, mh, _, _, _, _ = track_times_variable_rates_007(total_time, total_pool_size, rates_var, ε, η, σ_ε, σ_η, kesten_timestep)
        #     push!(ihs, ih)
        #     push!(mhs, mh)
        # end
        # ih_avg = mean(ihs)
        # mh_avg = mean(mhs)
        
        # # Calculate outputs
        # time_max_imm = argmax(ih_avg) * kesten_timestep
        # time_max_mat = argmax(mh_avg) * kesten_timestep
        # combined = ih_avg .+ mh_avg
        # time_max_combined = argmax(combined) * kesten_timestep
        # final_imm = ih_avg[end]
        # final_mat = mh_avg[end]

        paramss = (A1, lambda1, A2, lambda2, m, A3, lambda3)
        
        sol, synapse_sizes_var, synapse_sizes_history_var, synapses_var, ih, mh = run_simulation_diffeq_var007(total_time, total_pool_size, paramss, ε, η, σ_ε, σ_η, kesten_timestep);
        immature_population_var = sol[1, :]
        mature_population_var = sol[2, :]

        time_max_imm = sol.t[argmax(immature_population_var)]*kesten_timestep
        time_max_mat = sol.t[argmax(mature_population_var)]*kesten_timestep
        combined_population_var = immature_population_var .+ mature_population_var
        time_max_combined = sol.t[argmax(combined_population_var)]*kesten_timestep
        final_imm = immature_population_var[end]
        final_mat = mature_population_var[end]

        return [time_max_imm, time_max_mat, time_max_combined, final_imm, final_mat]
    end
end

# Parameter bounds for sensitivity analysis
param_bounds = Matrix{Float64}(undef, 7, 2)
param_means = [A1, lambda1, A2, lambda2, m, A3, lambda3]
deltas = 1 .* param_means

for i in 1:7
    param_bounds[i, 1] = max(0.000001, param_means[i] - deltas[i])
    param_bounds[i, 2] = param_means[i] + deltas[i]
end

for i in 1:7
    for j in 1:2
        param_bounds[i,j]=bounds[i][j]
    end
end

param_bounds


# Set up Sobol analysis with parallel processing
sampler = SobolSample()

ylabs = ["Time when immature \n population reaches max", "Time when mature \n population reaches max", "Time when combined \n population reaches max", "Final immature \n population value", "Final mature \n population value"]
xlabs = ["\n A1", "Creation", "\n λ1", "\n A2","Elimination", "\n λ2", "Maturation \n m", "\n A3", "De-maturation", "\n λ3"]
xpoints = [1,1.5,2,3,3.5,4,5,6,6.5,7]

# Parameter and output names for plotting
param_names = ["A1", "λ1", "A2", "λ2", "m", "A3", "λ3"]
output_names = [
    "Time to max immature", 
    "Time to max mature", 
    "Time to max combined", 
    "Final immature", 
    "Final mature"
]


##################

# A1,lambda1,A2,lambda2,m,A3,lambda3 = 0.9, 30, 2, 5, 0.05, 0.05, 2.



bounds = [[0, 5],    #A1  
            [0.1, 30],      #lambda1
            [0, 5],    #A2
            [0.1, 30],   #lambda2
            [0, 5],   #m
            [0, 5],    #A3
            [0.1, 30]]   #lambda3

# this is good
# bounds = [[0, 10],    #A1  
#             [1,1.2*lambda1],      #lambda1
#             [0, 10],    #A2
#             [1, 2*lambda2],   #lambda2
#             [0, 10],   #m
#             [0, 10],    #A3
#             [1, 1.5*lambda1]]   #lambda3


# Standard regression GSA
reg_sens = gsa(model_func, RegressionGSA(true), bounds, samples = 1000)

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

# rs1

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


plot(rs1h, h1, h2, layout=(3,1), size=(400,1000), legend=false)


########################################################
# Sobol indices

using QuasiMonteCarlo

# Construct randomized Sobol sampler
sampler = SobolSample(Shift())
n_samples = 1000
d = size(param_bounds, 1)

# Generate QMC samples
A_raw = QuasiMonteCarlo.sample(n_samples, d, sampler)
B_raw = QuasiMonteCarlo.sample(n_samples, d, sampler)

# Scale samples to parameter bounds
function scale_samples(samples, a, b)
    return a .+ (b .- a) .* samples
end

a = param_bounds[:, 1]
b = param_bounds[:, 2]

A = scale_samples(A_raw, a, b)
B = scale_samples(B_raw, a, b)


# getting confidence intervals

# sobol_result = gsa(model_func, Sobol(nboot=100), A, B; 
#                    batch=false, 
#                    parallel=true, 
#                    resampling=true, 
#                    n_resamples=10)

                   
sobol_result = gsa(model_func, Sobol(nboot=10), A, B; 
                    batch=false, 
                    parallel=true, 
                    resampling=true,
                    n_resamples=10)


s1 = sobol_result.S1  # First-order indices
s2 = sobol_result.ST  # Total-effect indices
s1ci = sobol_result.S1_Conf_Int  # Confidence intervals for first-order indices
s2ci = sobol_result.ST_Conf_Int   # Confidence intervals for total-effect indices

s1
s1ci 
s2
s2ci
heatmap(s1,
    title="First-Order Sobol Indices",
    xticks = (xpoints, xlabs),
    yticks = (1:5, ylabs),
    xlabel = "Parameters",
    ylabel="Outputs",
    colorbar_title="\n Total-Order Sensitivity Index",
    c=:viridis,
    size=(800, 600),rightmargin=5mm
)

vals = round.(s1; digits=3)
ci_widths = round.(s1ci; digits=3)

ann = []
for i in 1:size(s1, 1), j in 1:size(s1, 2)
    label = string(vals[i, j], " [±", ci_widths[i, j], "]")
    push!(ann, (j, i, text(label, :white, :center)))
end

h1 = heatmap(s1,
    title="First-Order Sobol Indices",
    xticks=(1:size(s1, 2), xlabs),
    yticks=(1:size(s1, 1), ylabs),
    xlabel="Parameters",
    ylabel="Outputs",
    colorbar_title="\nFirst-Order Index",
    c=:viridis,
    annotations=ann,
    size=(1600, 1200),
    right_margin=5mm
)


h2 = heatmap(s2,
    title="Total-Order Sobol Indices",
    xticks = (xpoints, xlabs),
    yticks = (1:5, ylabs),
    xlabel = "Parameters",
    ylabel="Outputs",
    colorbar_title="\n Total-Order Sensitivity Index",
    c=:viridis,
    size=(800, 600),rightmargin=5mm
)



plot(heatmap(s1), heatmap(s2), clim=(0,1))



minimum(s2)