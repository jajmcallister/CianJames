
##########
# simplified setup

# Parameters: a1, k1, b1, a2, k2, b2, total_pool
a1, k1, b1, a2, k2, b2, total_pool_size = rand(),rand(),rand(),rand(),rand(),rand(),total_pool_size

finaltime = 100
# simple_setup(t) = (b1*total_pool_size*exp(k2*t+k1*t)+a1*total_pool_size*exp(k2*t))/((b2+b1)*exp(k2*t+k1*t)+a1*exp(k2*t)+a2*exp(k1*t))

simple_setup2(t) = - ((a1 * exp(-t * k1) + b1)*total_pool_size)/((a1 * exp(-t * k1) + b1) + a2 * exp(-t * k2) + b2)*exp(-(((a1 * exp(-t * k1) + b1)+a2 * exp(-t * k2) + b2))*t) +  ((a1 * exp(-t * k1) + b1)*total_pool_size)/((a1 * exp(-t * k1) + b1) + a2 * exp(-t * k2) + b2)
simple_steady_sol = b1*total_pool_size/(b1+b2)

C = (-(total_pool_size*b2/((b1+b2)*(b1+b2-k1)))*(a1-a2) - total_pool_size*b1/(b1+b2))*exp(a1/k1 + a2/k2)
simple_setup3(t) = C*exp(-((a1/k1)*exp(-k1*t)+b1*t+(a2/k2)*exp(-k2*t)+b2*t)) + (total_pool_size/((b1+b2)*(b1+b2-k1)))*(a1*b2*exp(-k1*t)-a2*b2*exp(-k2*t))+total_pool_size*b1/(b1+b2)

plot(0:0.1:100, simple_setup2.(0:0.1:100))
hline!([simple_steady_sol])



# Define the ODE
function ode1!(du, u, p, t)
    a1, k1, b1, a2, k2, b2, total = p
    du[1] = (a1 * exp(-k1 * t) + b1) * total - (a1 * exp(-k1 * t) + b1 + a2 * exp(-k2 * t) + b2) * u[1]
end


p = [a1, k1, b1, a2, k2, b2, total_pool_size]  # Example parameter values

# Initial condition for y(0)
u0 = [0.0]

# Time span for the solution
tspan = (0.0, finaltime)

# Define the ODE problem
prob = ODEProblem(ode1!, u0, tspan, p)

# Solve the ODE
sol = solve(prob)
sol.u

# Plot the solution for y(t)
plot(sol.t, sol[1,:], label="y(t)", xlabel="t", ylabel="y(t)", title="Solution to the ODE")

# Survival fraction

# Get the maximum value and its index
max_value, max_index = findmax(sol.u)  # max_value is the maximum, max_index is its index

# Get the time corresponding to the maximum value
max_time = sol.t[max_index]

# Calculate the survival fraction starting from the maximum value
S = sol.u[max_index:end] ./ max_value  # survival fraction S(t) = u(t) / u(max)

# Extract the corresponding time values starting from max_index
time_after_max = sol.t[max_index:end]

# Plot survival fraction starting from the maximum time
plot(vcat(S...), ylim=(0,1), xlabel="Time", ylabel="Survival Fraction", title="Synapse Survival Over Time", lw=2, legend=false)


# using QuadGK  # For numerical integration
# using Plots

# # Define the parameters (hardcoded values)
# a1, k1, b1, a2, k2, b2, l = 1.0, 0.5, 0.5, 0.8, 0.2, 0.4, 1.0

# # Define A(t) and B(t) as functions
# A(t) = a1 * exp(-k1 * t) + b1 + a2 * exp(-k2 * t) + b2
# B(t) = (a1 * exp(-k1 * t) + b1) * l

# # Define the analytical solution based on integrals
# function u_analytical(t)
#     # Compute the exponential factor: exp(-∫ A(t) dt)
#     exp_neg_int_A, _ = QuadGK.quadgk(s -> A(s), 0, t)  # Integral of A(t) from 0 to t
#     exp_neg_A = exp(-exp_neg_int_A)
    
#     # Compute the integral for the particular solution
#     integral_B_expA, _ = QuadGK.quadgk(s -> B(s) * exp(QuadGK.quadgk(A, 0, s)[1]), 0, t)
    
#     # Final solution
#     return exp_neg_A * integral_B_expA
# end

# # Time span for plotting
# time_vals = 0.0:0.1:100.0
# u_vals = [u_analytical(t_val) for t_val in time_vals]

# # Plotting the analytical solution
# plot(time_vals, u_vals, xlabel="Time (t)", ylabel="u(t)", title="Analytical Solution to ODE", lw=2)



# Single function to compute the ODE solution at time t
function full_soln_f(t, m, a1, a2, b1, b2, k1, k2)
    # Solve for the coefficients of the particular solution
    A = m * a1 / (b1 + b2 - k1)
    B = -m * a2 / (b1 + b2 - k2)
    C = m * b1 / (b1 + b2)
    
    # Compute the constant C1 using y(0) = 0
    C1 = -(A + B + C)
    
    # Homogeneous solution
    y_h = C1 * exp(-(b1 + b2) * t)
    
    # Particular solution
    y_p = A * exp(-k1 * t) + B * exp(-k2 * t) + C
    
    # General solution is the sum of homogeneous and particular solutions
    return y_h + y_p
end


plot(0:100,full_soln_f.(0:100, total_pool_size, a1, a2, b1, b2, k1, k2))
plot!(sol.t, sol[1,:], label="y(t)", xlabel="t", ylabel="y(t)", title="Solution to the ODE")



function compute_C(total_pool_size, a1, a2, b1, b2, k1, k2)
    part_a = (total_pool_size * a1) / (k1 + b1 + b2)
    part_b = (total_pool_size * a2) / (k2 + b1 + b2)
    part_c = (total_pool_size * b1) / (b1 + b2)
    
    exp_term = exp(-(a2 / k2 + a1 / k1))
    
    return -(part_a + part_b + part_c) * exp_term
end

C=compute_C(total_pool_size, a1, a2, b1, b2, k1, k2)
y_homogeneous(t) = C* exp((a2 / k2) * exp(-k2 * t) + (a1 / k1) * exp(-k1 * t) - (b2 * t) - (b1 * t))
y_particular(t) = (total_pool_size * a1) / (k1 + b1 + b2) * exp(-k1 * t) + 
          (total_pool_size * a2) / (k2 + b1 + b2) * exp(-k2 * t) + 
          (total_pool_size * b1) / (b1 + b2)

y_general(t) = y_homogeneous(t) + y_particular(t) 

plot(0:0.1:100, full_soln_f.(0:0.1:100))
plot!(sol.t, sol[1,:], label="y(t)", xlabel="t", ylabel="y(t)", title="Solution to the ODE")



















#########
#########
########


# Function to track times with just a resource pool and existing population
function simplified_track_times_constant_rates(total_time, total_pool_size, rates, kesten_time_step)
    pool = total_pool_size  # Total synapse pool
    existing = 0            # Initially, no synapses in the existing population
    steps = trunc(Int, total_time / kesten_time_step)

    pool_history = []
    existing_history = []
    state_records = zeros(total_pool_size, steps)  # To store the states (0 = pool, 1 = existing)

    pool_pop = [i for i in 1:total_pool_size]  # Synapses in pool
    existing_pop = []

    c, e = rates  # Rates for creation (c) and elimination (e)

    # Simulation
    for t in 1:steps
        # 1. Transitions from pool to existing (creation of synapses)
        pool_to_existing = rand(Binomial(pool, c * kesten_time_step))
        pool -= pool_to_existing
        existing += pool_to_existing

        # 2. Transitions from existing to pool (elimination of synapses)
        existing_to_pool = rand(Binomial(existing, e * kesten_time_step))
        existing -= existing_to_pool
        pool += existing_to_pool

        # Randomly sample synapses for transitions
        pool_to_existing_indxs = sample(pool_pop, min(pool_to_existing, length(pool_pop)), replace=false)
        filter!(x -> !in(x, pool_to_existing_indxs), pool_pop)  # Remove from pool
        append!(existing_pop, pool_to_existing_indxs)

        existing_to_pool_indxs = sample(existing_pop, min(existing_to_pool, length(existing_pop)), replace=false)
        filter!(x -> !in(x, existing_to_pool_indxs), existing_pop)  # Remove from existing population
        append!(pool_pop, existing_to_pool_indxs)

        # Record the states of each synapse at time t
        for j in 1:total_pool_size
            if j in pool_pop
                state_records[j, t] = 0  # In the pool
            elseif j in existing_pop
                state_records[j, t] = 1  # In the existing population
            end
        end

        # Store the population histories
        push!(pool_history, pool)
        push!(existing_history, existing)
    end

    return existing_history, state_records
end

# Example usage
total_time = 100.0
total_pool_size = 1000
rates = (0.2, 0.2)  # Creation rate (c), elimination rate (e)

existing_history, state_records = simplified_track_times_constant_rates(total_time, total_pool_size, rates, kesten_time_step)

# Plot the number of synapses in the existing population over time
plot(1:length(existing_history), existing_history, xlabel="Time Steps", ylabel="Number of Existing Synapses",
     title="Synapse Population Over Time", lw=2, legend=false)















# Function to track times with just a resource pool and existing population
function simplified_track_times_variable_rates(total_time, total_pool_size, rates, kesten_time_step)
    pool = total_pool_size  # Total synapse pool
    existing = 0            # Initially, no synapses in the existing population
    steps = trunc(Int, total_time / kesten_time_step)

    pool_history = []
    existing_history = []
    state_records = zeros(total_pool_size, steps)  # To store the states (0 = pool, 1 = existing)

    pool_pop = [i for i in 1:total_pool_size]  # Synapses in pool
    existing_pop = []

    cr,m,el,_ = rates  # Rates for creation (c) and elimination (e)

    # Simulation
    for t in 1:steps
        # 1. Transitions from pool to existing (creation of synapses)
        pool_to_existing = rand(Binomial(pool, cr[t] * kesten_time_step))
        pool -= pool_to_existing
        existing += pool_to_existing

        # 2. Transitions from existing to pool (elimination of synapses)
        existing_to_pool = rand(Binomial(existing, el[t] * kesten_time_step))
        existing -= existing_to_pool
        pool += existing_to_pool

        # Randomly sample synapses for transitions
        pool_to_existing_indxs = sample(pool_pop, min(pool_to_existing, length(pool_pop)), replace=false)
        filter!(x -> !in(x, pool_to_existing_indxs), pool_pop)  # Remove from pool
        append!(existing_pop, pool_to_existing_indxs)

        existing_to_pool_indxs = sample(existing_pop, min(existing_to_pool, length(existing_pop)), replace=false)
        filter!(x -> !in(x, existing_to_pool_indxs), existing_pop)  # Remove from existing population
        append!(pool_pop, existing_to_pool_indxs)

        # Record the states of each synapse at time t
        for j in 1:total_pool_size
            if j in pool_pop
                state_records[j, t] = 0  # In the pool
            elseif j in existing_pop
                state_records[j, t] = 1  # In the existing population
            end
        end

        # Store the population histories
        push!(pool_history, pool)
        push!(existing_history, existing)
    end

    return existing_history, state_records
end

# Example usage
total_time = 500.0
total_pool_size = 1000

# Define transition rates
m, i = 0.2, 0.05
elim = elimination_func.(0:kesten_time_step:total_time)
creat = creation_func.(0:kesten_time_step:total_time)

rates_var = creat, m, elim, i
lambda = 2


existing_history_var, state_records_var = simplified_track_times_variable_rates(total_time, total_pool_size, rates_var,kesten_time_step)

# Plot the number of synapses in the existing population over time
plot(1:length(existing_history_var), existing_history_var, xlabel="Time Steps", ylabel="Number of Existing Synapses",
     title="Synapse Population Over Time", lw=2, legend=false)


survival_fraction_var = compute_survival_fraction(state_records_var[:,10:end])

# Plot survival fraction over time
plot(1:length(survival_fraction_var), survival_fraction_var, xlabel="Time Steps", ylabel="Survival Fraction",
    title="Synapse Survival Fraction Over Time", lw=2, legend=false, ylim=(0,1))
     
