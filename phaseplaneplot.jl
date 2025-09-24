using DifferentialEquations
using Plots



#################
# CONSTANT RATES
#################

# Reduced system for phase plane plotting
function phase_synapse_dynamics!(du, u, p, t)
    # Parameters
    c, e, m, d, total_pool_size = p
    N_I, N_M = u
    N_P = total_pool_size - N_I - N_M  # eliminate N_P

    du[1] = c * N_P - (m + e) * N_I + d * N_M  # dN_I/dt
    du[2] = m * N_I - d * N_M                  # dN_M/dt
end


# Parameters
c, e, m, d = 0.5, 0.2, 0.2, 0.2
total_pool_size = 100
paramss = (c, e, m, d, total_pool_size)

# Initial conditions
u01 = [0.0, 0.0]       # I, M
u02 = [10.0, 30.0]
u03 = [10.0, 70.0]
u04 = [90.0, 10.0]
tspan = (0.0, 120.0)

# Solve
prob1 = ODEProblem(phase_synapse_dynamics!, u01, tspan, paramss)
prob2 = ODEProblem(phase_synapse_dynamics!, u02, tspan, paramss)
prob3 = ODEProblem(phase_synapse_dynamics!, u03, tspan, paramss)
prob4 = ODEProblem(phase_synapse_dynamics!, u04, tspan, paramss)

sol1, sol2, sol3, sol4 = solve.(Ref(prob1)), solve.(Ref(prob2)), solve.(Ref(prob3)), solve.(Ref(prob4));

default(tickdirection=:out)
# Plot trajectories
plot(sol1[1,:], sol1[2,:], lw=3, label="Traj 1", xlabel="N_I", ylabel="N_M", title="Phase plane")
plot!(sol2[1,:], sol2[2,:], lw=3, label="Traj 2")
plot!(sol3[1,:], sol3[2,:], lw=3, label="Traj 3")
plot!(sol4[1,:], sol4[2,:], lw=3, label="Traj 4")
scatter!([total_pool_size / (1 + m/d + e/c)], [total_pool_size / (1 + d/m + (e*d)/(c*m))], label="Steady State",markershape =:star, color=:red, markersize=12)


meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))

x, y = meshgrid(0:10:100, 0:10:100)
u = @. c * (100-x-y) - (e + m) * x + i * y
v = @. m * x - i * y

q = quiver(x, y, quiver=(0.3*u, 0.3*v),title="Phase plane", xlabel="I",c=:grey, alpha=0.5, ylabel="M")
plot!(sol1[1,:], sol1[2,:], lw=3, label="Traj 1", xlabel="N_I", ylabel="N_M", title="Phase plane")
plot!(sol2[1,:], sol2[2,:], lw=3, label="Traj 2")
plot!(sol3[1,:], sol3[2,:], lw=3, label="Traj 3")
plot!(sol4[1,:], sol4[2,:], lw=3, label="Traj 4")
scatter!([u01[1]], [u01[2]], label="Initial point",markershape =:circle, color=:black, markersize=4)
scatter!([u02[1]], [u02[2]], label=false,markershape =:circle, color=:black, markersize=4)
scatter!([u03[1]], [u03[2]], label=false,markershape =:circle, color=:black, markersize=4)
scatter!([u04[1]], [u04[2]], label=false,markershape =:circle, color=:black, markersize=4)
scatter!([total_pool_size / (1 + m/d + e/c)], [total_pool_size / (1 + d/m + (e*d)/(c*m))],dpi=600, grid=false,label="Steady State",markershape =:star, color=:red, markersize=12)

# savefig(q, "C://Users/B00955735/OneDrive - Ulster University/Desktop/phaseplane.png")

# nullclines
function nullcline_NI(N_I)
    # Nullcline for dN_I/dt = 0
    return (100*c-N_I * c-(e+m)*N_I)/(c-d)
end

function nullcline_NI(N_I)
    # Nullcline for dN_I/dt = 0
    return (m/(d+c))*N_I + (c*100)/(d+c)
end

function nullcline_NM(N_I)
    # Nullcline for dN_M/dt = 0
    return (m / d) * N_I
end

# Range for N_I values
N_I_vals = collect(0:1:100)

# Compute corresponding N_M values for each nullcline
N_M_vals_NI = [nullcline_NI(N_I) for N_I in N_I_vals]  # dN_I/dt = 0
N_M_vals_NM = [nullcline_NM(N_I) for N_I in N_I_vals]  # dN_M/dt = 0

# Plot the nullclines

N_I1 = sol1[1, :]  # First variable (N_I) over time
N_M1 = sol1[2, :]  # Second variable (N_M) over time
N_I2 = sol2[1, :] 
N_M2 = sol2[2, :]
N_I3 = sol3[1, :] 
N_M3 = sol3[2, :]
N_I4 = sol4[1, :] 
N_M4 = sol4[2, :]

nullcline = quiver(x, y, quiver=(0.3*u, 0.3*v), title="Nullclines",c=:grey, allpha=0.5,)
plot!(N_I1, N_M1, lw=2, label="Trajectory 1")
plot!(N_I2, N_M2, lw=2, label="Trajectory 2")
plot!(N_I3, N_M3, lw=2, label="Trajectory 3")
plot!(N_I4, N_M4, lw=2, label="Trajectory 4")
plot!(N_I_vals, N_M_vals_NI, label="I nullcline", linewidth=4, linestyle=:dash)
plot!(N_I_vals, N_M_vals_NM, label="M nullcline", linewidth=4, linestyle=:dash)
scatter!([u01[1]], [u01[2]], label="Initial point",markershape =:circle, color=:black, markersize=4)
scatter!([u02[1]], [u02[2]], label=false,markershape =:circle, color=:black, markersize=4)
scatter!([u03[1]], [u03[2]], label=false,markershape =:circle, color=:black, markersize=4)
scatter!([u04[1]], [u04[2]], label=false,markershape =:circle, color=:black, markersize=4)
scatter!([total_pool_size / (1 + m/d + e/c)], [total_pool_size / (1 + d/m + (e*d)/(c*m))],dpi=600, xlim=(0,100), ylim=(0,100), grid=false,label="Steady State",markershape =:star, color=:red, markersize=12)
plot!(legend=:right,dpi=600)

# savefig(nullcline, "C://Users/B00955735/OneDrive - Ulster University/Desktop/nullcline.png")



################
# VARIABLE RATES
################

function phase_synapse_dynamics_var!(du, u, p, t)
    A1, λ1, A2, λ2, m, A3, λ3, Ntot = p
    N_I, N_M = u
    N_P = Ntot - N_I - N_M

    # time-dependent rates
    c_t = creation_func(t, A1, λ1)
    e_t = elimination_func(t, A2, λ2)

    # for phase plot: approximate mean dematuration rate
    dematuration_rate = A3 * exp(-mean(N_M > 0 ? [1.0] : [0.0]) / λ3)  # or pick a representative value

    du[1] = c_t * N_P - (m + e_t) * N_I + dematuration_rate * N_M
    du[2] = m * N_I - dematuration_rate * N_M
end

# Example parameters
A1, λ1 = 0.9, 30.0
A2, λ2 = 2, 5.0
m      = 0.05
A3, λ3 = 0.05, .5
Ntot   = 100

paramss_var = (A1, λ1, A2, λ2, m, A3, λ3, Ntot)

u01 = [0.0, 0.0]
u02 = [10.0, 30.0]
u03 = [10.0, 70.0]
u04 = [90.0, 10.0]
tspan = (0.0, 200.0)

prob1 = ODEProblem(phase_synapse_dynamics_var!, u01, tspan, paramss_var)
prob2 = ODEProblem(phase_synapse_dynamics_var!, u02, tspan, paramss_var)
prob3 = ODEProblem(phase_synapse_dynamics_var!, u03, tspan, paramss_var)
prob4 = ODEProblem(phase_synapse_dynamics_var!, u04, tspan, paramss_var)

sol1, sol2, sol3, sol4 = solve.(Ref(prob1)), solve.(Ref(prob2)), solve.(Ref(prob3)), solve.(Ref(prob4));

xs, ys = 0:5:Ntot, 0:5:Ntot
X = Float64[]; Y = Float64[]; U = Float64[]; V = Float64[]
for x in xs, y in ys
    if x + y <= Ntot
        du = zeros(2)
        phase_synapse_dynamics_var!(du, [x,y], paramss_var, 0.0)  # pick t=0 or any fixed t
        push!(X, x); push!(Y, y); push!(U, du[1]); push!(V, du[2])
    end
end
mag = sqrt.(U.^2 .+ V.^2) .+ eps()
Un, Vn = U ./ mag, V ./ mag

scale = .20   # increase arrow length
Un .*= scale
Vn .*= scale

quiver(X, Y, quiver=(Un,Vn), color=:gray, alpha=0.5, xlabel="N_I", ylabel="N_M", title="Variable rate phase plane")
plot!(sol1[1,:], sol1[2,:], lw=3, label="Traj 1")
plot!(sol2[1,:], sol2[2,:], lw=3, label="Traj 2")
plot!(sol3[1,:], sol3[2,:], lw=3, label="Traj 3")
plot!(sol4[1,:], sol4[2,:], lw=3, label="Traj 4")
plot!(xlim=(0,100), ylim=(0,100))
scatter!([u01[1]], [u01[2]], label="Initial point",markershape =:circle, color=:black, markersize=4)
scatter!([u02[1]], [u02[2]], label=false,markershape =:circle, color=:black, markersize=4)
scatter!([u03[1]], [u03[2]], label=false,markershape =:circle, color=:black, markersize=4)
scatter!([u04[1]], [u04[2]], label=false,markershape =:circle, color=:black, markersize=4)
scatter!([sol1[1,end]], [sol1[2,end]],dpi=600, grid=false,label="End point",markershape =:star, color=:red, markersize=12)

