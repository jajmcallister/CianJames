
# This is for constant rates m,c,e,i
eigenvalue1, eigenvalue2 = (-(e+m+c+i)+sqrt((e+m+c+i)^2 - 4*(i*e+i*c+m*c)))/2, (-(e+m+c+i)-sqrt((e+m+c+i)^2 - 4*(i*e+i*c+m*c)))/2
eigenvector1, eigenvector2 = [1 m/(i+eigenvalue1)], [1 m/(i+eigenvalue2)]



using LinearAlgebra




# function Solution()
time_for_plot=0:0.1:100


function analytical_solution(t,c,m,e,i)

    matrix = [-(e+m+c) (i-c) ; m -i]

    eigenvalues,eigenvectors = eigen(matrix)

    eigenvalues
    eigvec1 = eigenvectors[:,1]
    eigvec2 = eigenvectors[:,2]

    inhomogeneous_term = [c*total_pool_size 0]

    particular_soln = matrix \ - inhomogeneous_term'

    constants = eigenvectors \ - particular_soln

    imm = constants[1] * eigvec1[1] * exp(eigenvalues[1] * t) + constants[2] * eigvec2[1] .* exp(eigenvalues[2] * t) .+ particular_soln[1]
    mmm = constants[1] * eigvec1[2] * exp(eigenvalues[1] * t) + constants[2] * eigvec2[2] .* exp(eigenvalues[2] * t) .+ particular_soln[2]
    return imm, mmm
    
end



matrix = [-(e+m+c) (i-c) ; m -i]

eigenvalues,eigenvectors = eigen(matrix)

inhomogeneous_term = [c*total_pool_size 0]

particular_soln = matrix \ - inhomogeneous_term'

constants = eigenvectors \ - particular_soln


# constant1 = - (eigenvalues[2]* c*total_pool_size * (i+ eigenvalues[1]))/((eigenvalues[1]-eigenvalues[2])*(e*i+c*i+c*m))
# constant2 = - (eigenvalues[1]* c*total_pool_size * (i+ eigenvalues[2]))/((eigenvalues[2]-eigenvalues[1])*(e*i+c*i+c*m))


eigvalue1 = eigenvalues[1]
eigvalue2 = eigenvalues[2]
eigvec1 = eigenvectors[:,1]
eigvec2 = eigenvectors[:,2]

kappa = - eigvalue2*constants[2]*sum(eigvec2) / eigvalue1*constants[1]*sum(eigvec1)



# Compute imm and mmm over the time range
results = analytical_solution.(time_for_plot,c,m,e,i)

# Separate imm and mmm values
immature_values = map(r -> r[1], results)  # Extract imm values
mature_values = map(r -> r[2], results)  # Extract mmm values

# Plot the results
plot(time_for_plot, immature_values, label="imm", xlabel="Time", ylabel="Value", legend=:topright)
plot!(time_for_plot, mature_values, label="mat", ylim=(0,1000),legend=false)
plot!(time_for_plot, immature_values+mature_values, lw=3)



##########
# simplified setup

# Parameters: a1, k1, b1, a2, k2, b2, total_pool
a1, k1, b1, a2, k2, b2, total_pool_size = rand(),rand(),rand(),rand(),rand(),rand(),total_pool_size


simple_setup(t) = (b1*total_pool_size*exp(k2*t+k1*t)+a1*total_pool_size*exp(k2*t))/((b2+b1)*exp(k2*t+k1*t)+a1*exp(k2*t)+a2*exp(k1*t))

tt = ((1/(k1-k2))*log((a1*b2*k1)/(a2*b1*k2)))

plot(simple_setup.(0:1:100))

# Define the ODE
function ode1!(du, u, p, t)
    a1, k1, b1, a2, k2, b2, l = p
    du[1] = (a1 * exp(-k1 * t) + b1) * l - (a1 * exp(-k1 * t) + b1 + a2 * exp(-k2 * t) + b2) * u[1]
    # du[1] = (b1) * l - (b1 + b2) * u[1]
end


p = [a1, k1, b1, a2, k2, b2, total_pool_size]  # Example parameter values

# Initial condition for y(0)
u0 = [0.0]

# Time span for the solution
tspan = (0.0, 100.0)

# Define the ODE problem
prob = ODEProblem(ode1!, u0, tspan, p)

# Solve the ODE
sol = solve(prob)
sol.u

# Plot the solution for y(t)
plot(sol.t, sol[1,:], label="y(t)", xlabel="t", ylabel="y(t)", title="Solution to the ODE")

using QuadGK  # For numerical integration
using Plots

# Define the parameters (hardcoded values)
a1, k1, b1, a2, k2, b2, l = 1.0, 0.5, 0.5, 0.8, 0.2, 0.4, 1.0

# Define A(t) and B(t) as functions
A(t) = a1 * exp(-k1 * t) + b1 + a2 * exp(-k2 * t) + b2
B(t) = (a1 * exp(-k1 * t) + b1) * l

# Define the analytical solution based on integrals
function u_analytical(t)
    # Compute the exponential factor: exp(-∫ A(t) dt)
    exp_neg_int_A, _ = QuadGK.quadgk(s -> A(s), 0, t)  # Integral of A(t) from 0 to t
    exp_neg_A = exp(-exp_neg_int_A)
    
    # Compute the integral for the particular solution
    integral_B_expA, _ = QuadGK.quadgk(s -> B(s) * exp(QuadGK.quadgk(A, 0, s)[1]), 0, t)
    
    # Final solution
    return exp_neg_A * integral_B_expA
end

# Time span for plotting
time_vals = 0.0:0.1:100.0
u_vals = [u_analytical(t_val) for t_val in time_vals]

# Plotting the analytical solution
plot(time_vals, u_vals, xlabel="Time (t)", ylabel="u(t)", title="Analytical Solution to ODE", lw=2)
