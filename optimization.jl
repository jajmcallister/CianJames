using Optimization


function optimise_synapticmaturation(x, p)
    total_pool_size,total_time,kesten_timestep,ε, η, σ_ε, σ_η = p
    c, m, e, i = x
    rates = (c, m, e, i)

    sol, synapse_sizes_diffeq, synapses_diffeq = syn_maturation_functions.run_simulation_diffeq(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, kesten_timestep);

    time_array_diffeq = sol.t
    immature_population_diffeq = sol[1, :]
    mature_population_diffeq = sol[2, :]

    ratio_mature_to_immature = mature_population_diffeq[end]/immature_population_diffeq[end]

    desired_ratio = 2

    error = mature_population_diffeq[end] - 500

    return error

end


rosenbrock(x, p) = (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
x0 = zeros(2)
p = [1.0, 100.0]


ε, η = 1.0, 0.0
σ_ε, σ_η = .5, .5
rates = (c, m, e, i)
kesten_timestep = 0.01


p = total_pool_size,total_time,kesten_timestep,ε, η, σ_ε, σ_η
x0 = zeros(4)
prob = OptimizationProblem(optimise_synapticmaturation, x0, p)


using OptimizationBBO
prob = OptimizationProblem(optimise_synapticmaturation, x0, p, lb = [0.0, 0.0, 0.0, 0.0], ub = [1.0, 1.0, 1.0, 1.0])
sol = solve(prob, BBO_adaptive_de_rand_1_bin_radiuslimited())

z = sol[1], sol[2], sol[3], sol[4]

