using Optimization
using OptimizationBBO
using Optimization, OptimizationOptimJL
using .syn_maturation_functions




function optimise_synapticmaturation(x, p)
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
    
    target = 100
    # Calculate the error
    error = mature_population_diffeq[end] - target

    # Return the squared error
    return error^2
end


function optim(x,p)
    total_pool_size = Int(p[1])
    a1,a2,b1,b2,m,i = x
    final_I_value = (i*total_pool_size*(a1+b1))/(m*(a1+b1)+i*(a1+a2+b1+b2))
    final_M_value = (total_pool_size*(a1+b1))/(a1+b1+(i/m)*(a1+a2+b1+b2))

    target_M = 50
    target_I = 40

    errorM = final_M_value-target_M
    errorI = final_I_value-target_I

    return errorM^2 + errorI^2
end



# Define initial guesses for parameters
x0 = [0.5,0.5,0.5, 0.5, 0.1, 0.1]

# Define the parameters you want to pass to the objective function
total_pool_size = 100  # Example values
total_time = 100.0
kesten_timestep = 0.01
ε = 1.0
η = 0.
σ_ε = 0.01
σ_η = 0.01

p = [total_pool_size, total_time, kesten_timestep, ε, η, σ_ε, σ_η]


# Define bounds for the parameters
lower_bounds = [0.0, 0.0, 0.0, 0.0,0.0,0.0]
upper_bounds = [2.0, 2.0, 2.0, 2.0,2.0,2.0]

# Set up the optimization problem with bounds
opt_function = OptimizationFunction(optim, Optimization.AutoForwardDiff())
prob = Optimization.OptimizationProblem(opt_function, x0, p, lb=lower_bounds, ub=upper_bounds)

# Solve the optimization problem using LBFGS with bounds
result = Optimization.solve(prob, LBFGS())

# Output the results
println("Optimal parameters (c, m, e, i): ", result.minimizer)
println("Objective function value: ", result.minimum)