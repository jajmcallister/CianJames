using Optimization
using OptimizationBBO
using Optimization, OptimizationOptimJL
using .syn_maturation_functions




function optimise_synapticmaturation1(x, p)
    total_pool_size = Int(p[1])
    total_time = p[2] 
    kesten_timestep = p[3]
    ε = p[4]
    η = p[5]
    σ_ε = p[6]
    σ_η = p[7]
    c, m, e, i = x
    rates = (c, m, e, i)

    # Run your simulation
    sol, synapse_sizes_diffeq, synapses_diffeq = syn_maturation_functions.run_simulation_diffeq(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, kesten_timestep)

    # Extract the mature population at the end of the simulation
    mature_population_diffeq = sol[2, :]

    # Calculate the error (you want mature_population_diffeq[end] = 500)
    error = mature_population_diffeq[end] - 500

    # Return the squared error
    return error^2
end

function optimise_synapticmaturation2(x, p)
    total_pool_size = Int(p[1])
    total_time = p[2] 
    kesten_timestep = p[3]
    ε = p[4]
    η = p[5]
    σ_ε = p[6]
    σ_η = p[7]
    c, m, e, i = x
    rates = (c, m, e, i)

    # Run your simulation
    sol, synapse_sizes_diffeq, synapses_diffeq = syn_maturation_functions.run_simulation_diffeq(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, kesten_timestep)

    # Extract immature and mature populations at the end of the simulation
    immature_population_diffeq = sol[1, :]
    mature_population_diffeq = sol[2, :]

    # Define the target values for the populations
    target_N_M = 500
    target_N_I = 600

    # Calculate the error for both populations
    error_mature = mature_population_diffeq[end] - target_N_M
    error_immature = immature_population_diffeq[end] - target_N_I

    # Return the combined squared error for both populations
    return error_mature^2 + error_immature^2
end

function optim3(x,p)
    total_pool_size = Int(p[1])
    total_time = p[2] 
    kesten_timestep = p[3]
    ε = p[4]
    η = p[5]
    σ_ε = p[6]
    σ_η = p[7]
    c, m, e, i = x
    final_I_value = total_pool_size / (1 + m/i + e/c)
    final_M_value = total_pool_size / (1 + i/m + (e*i)/(c*m))
    target_M = 500
    target_I = 400

    errorM = final_M_value-target_M
    errorI = final_I_value-target_I

    return errorM^2 + errorI^2
end



# Define initial guesses for parameters c, m, e, i
x0 = [0.5, 0.5, 0.1, 0.1]  # Initial guess for [c, m, e, i]

# Define the parameters you want to pass to the objective function
total_pool_size = 1000  # Example values
total_time = 100.0
kesten_timestep = 0.01
ε = 1.0
η = 0.
σ_ε = 0.01
σ_η = 0.01

p = [total_pool_size, total_time, kesten_timestep, ε, η, σ_ε, σ_η]


# Define bounds for the parameters c, m, e, i
lower_bounds = [0.0, 0.0, 0.0, 0.0]
upper_bounds = [2.0, 2.0, 2.0, 2.0]

# Set up the optimization problem with bounds
opt_function = OptimizationFunction(optimise_synapticmaturation2, Optimization.AutoForwardDiff())
prob = Optimization.OptimizationProblem(opt_function, x0, p, lb=lower_bounds, ub=upper_bounds)

# Solve the optimization problem using LBFGS with bounds
result = Optimization.solve(prob, LBFGS())

# Output the results
println("Optimal parameters (c, m, e, i): ", result.minimizer)
println("Objective function value: ", result.minimum)