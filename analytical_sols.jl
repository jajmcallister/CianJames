
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
plot!(time_for_plot, immature_values, label="imm", xlabel="Time", ylabel="Value", legend=:topright)
plot!(time_for_plot, mature_values, label="mat", ylim=(0,1000),legend=false)
plot!(time_for_plot, immature_values+mature_values, lw=3)



##########



simple_setup(t) = (b1*total_pool_size*exp(k2*t+k1*t)+a1*total_pool_size*exp(k2*t))/((b2+b1)*exp(k2*t+k1*t)+a1*exp(k2*t)+a2*exp(k1*t))
plot(simple_setup.(0:1:200))
vline!([tt])

# Define the ODE
function ode1!(du, u, p, t)
    a1, k1, b1, a2, k2, b2, l = p
    du[1] = (a1 * exp(-k1 * t) + b1) * l - (a1 * exp(-k1 * t) + b1 + a2 * exp(-k2 * t) + b2) * u[1]
end


# Parameters: a1, k1, b1, a2, k2, b2, total_pool
a1, k1, b1, a2, k2, b2, total_pool_size = rand(),rand(),rand(),rand(),rand(),rand(),total_pool_size
p = [a1, k1, b1, a2, k2, b2, total_pool_size]  # Example parameter values

# Initial condition for y(0)
u0 = [0.0]

# Time span for the solution
tspan = (0.0, 200.0)

# Define the ODE problem
prob = ODEProblem(ode1!, u0, tspan, p)

# Solve the ODE
sol = solve(prob)
sol.u

# Plot the solution for y(t)
plot!(sol.t, vcat(sol.u...), label="y(t)", xlabel="t", ylabel="y(t)", title="Solution to the ODE")

using Symbolics

# Define the variables and parameters
@variables t k1 k2 a1 a2 b1 b2 total_pool_size

# Define the numerator and denominator
u = b1 * total_pool_size * exp((k2 + k1) * t) + a1 * total_pool_size * exp(k2 * t)
v = (b2 + b1) * exp((k2 + k1) * t) + a1 * exp(k2 * t) + a2 * exp(k1 * t)

# Apply the quotient rule: diff(u(t)/v(t))
f = u / v
f_prime = expand_derivatives(Differential(t)(f))

f_prime_func(t) = (a1*k2*total_pool_size*exp(k2*t) + b1*(k1 + k2)*total_pool_size*exp((k1 + k2)*t)) / (a1*exp(k2*t) + a2*exp(k1*t) + (b1 + b2)*exp((k1 + k2)*t)) - (a1*k2*exp(k2*t) + a2*k1*exp(k1*t) + (b1 + b2)*(k1 + k2)*exp((k1 + k2)*t))*((a1*total_pool_size*exp(k2*t) + b1*total_pool_size*exp((k1 + k2)*t)) / ((a1*exp(k2*t) + a2*exp(k1*t) + (b1 + b2)*exp((k1 + k2)*t))^2))

plot(f_prime_func.(0:1:100))

tt = (1/(k1-k2))*log((a1*b2*k1)/(a2*b1*k2))