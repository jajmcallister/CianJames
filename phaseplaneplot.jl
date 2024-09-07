using DifferentialEquations
using Plots

# Define the system of equations (right-hand side of the ODEs)
function phase_synapse_dynamics!(du, u, p, t)
    c, m, e, i = p  # Correctly unpack the parameters
    N_I, N_M, N_P = u
    
    du[1] = c * N_P - (e + m) * N_I + i * N_M  # dN_I/dt
    du[2] = m * N_I - i * N_M                  # dN_M/dt
    du[3] = e * N_I - c * N_P                  # dN_P/dt
end

# Parameters and initial conditions
c, m, e, i = 0.2,0.1,0.2,0.3

paramss = (c, m, e, i)  # Example rates
u01 = [0.0, 0.0, 1000.0]  # Initial conditions for N_I, N_M, N_P
u02 = [100.0, 300.0, 600.0]
u03 = [100.0, 700.0, 200.0]
u04 = [900.0, 100.0, 0.0]
tspan = (0.0, 1000.0)

# Solve the ODE system
prob1 = ODEProblem(phase_synapse_dynamics!, u01, tspan, paramss)
prob2 = ODEProblem(phase_synapse_dynamics!, u02, tspan, paramss)
prob3 = ODEProblem(phase_synapse_dynamics!, u03, tspan, paramss)
prob4 = ODEProblem(phase_synapse_dynamics!, u04, tspan, paramss)
sol1 = solve(prob1)
sol2 = solve(prob2)
sol3 = solve(prob3)
sol4 = solve(prob4)

# Extract N_I and N_M from the solution
N_I1 = sol1[1, :]  # First variable (N_I) over time
N_M1 = sol1[2, :]  # Second variable (N_M) over time
N_I2 = sol2[1, :] 
N_M2 = sol2[2, :]
N_I3 = sol3[1, :] 
N_M3 = sol3[2, :]
N_I4 = sol4[1, :] 
N_M4 = sol4[2, :]

final_I_value = total_pool_size / (1 + m/i + e/c)
final_M_value = total_pool_size / (1 + i/m + (e*i)/(c*m))

p = plot(N_I1,N_M1,lw=2, label="Trajectory", title="Phase plot of I and M trajectory", xlabel="I", ylabel="M")

# savefig(p, "C://Users/B00955735/OneDrive - Ulster University/Desktop/phase.svg")


# Vector field

meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))

x, y = meshgrid(0:100:1000, 0:100:1000)
u = @. c * (1000-x-y) - (e + m) * x + i * y
v = @. m * x - i * y

q = quiver(x, y, quiver=(0.3*u, 0.3*v),title="Phase plane", xlabel="I", ylabel="M")
plot!(N_I1, N_M1, lw=2, label="Trajectory 1")
plot!(N_I2, N_M2, lw=2, label="Trajectory 2")
plot!(N_I3, N_M3, lw=2, label="Trajectory 3")
plot!(N_I4, N_M4, lw=2, label="Trajectory 4")
scatter!([final_I_value],[final_M_value],label="Steady state solution")

# savefig(q, "C://Users/B00955735/OneDrive - Ulster University/Desktop/vectorfield.png")

plot(N_M4)
plot!(N_I4)
plot!(N_M3+N_I3)


m/i
maximum(N_M4 ./ N_I4)

# nullclines



function nullcline_NI(N_I)
    # Nullcline for dN_I/dt = 0
    return (1000*c-N_I * c-(e+m)*N_I)/(c-i)
end

function nullcline_NM(N_I)
    # Nullcline for dN_M/dt = 0
    return (m / i) * N_I
end

# Range for N_I values
N_I_vals = collect(0:10:1000)

# Compute corresponding N_M values for each nullcline
N_M_vals_NI = [nullcline_NI(N_I) for N_I in N_I_vals]  # dN_I/dt = 0
N_M_vals_NM = [nullcline_NM(N_I) for N_I in N_I_vals]  # dN_M/dt = 0

# Plot the nullclines

nullcline = quiver(x, y, quiver=(0.3*u, 0.3*v), title="Nullclines")
plot!(N_I1, N_M1, lw=2, label="Trajectory 1")
plot!(N_I2, N_M2, lw=2, label="Trajectory 2")
plot!(N_I3, N_M3, lw=2, label="Trajectory 3")
plot!(N_I4, N_M4, lw=2, label="Trajectory 4")
plot!(N_I_vals, N_M_vals_NI, label="I nullcline", linewidth=2, xlim=(0,1000), ylim=(0,1000), linestyle=:dash)
plot!(N_I_vals, N_M_vals_NM, label="M nullcline", linewidth=2, xlim=(0,1000), ylim=(0,1000), linestyle=:dash)

# savefig(nullcline, "C://Users/B00955735/OneDrive - Ulster University/Desktop/nullcline.png")