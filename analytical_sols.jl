
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


using LinearAlgebra

# Define the function that computes kappa
function compute_kappa(m, c, e, i, total_pool_size)
    # Define the matrix
    matrix = [-(e + m + c) (i - c);
              m         -i]

    # Compute eigenvalues and eigenvectors
    eigenvalues, eigenvectors = eigen(matrix)

    # Inhomogeneous term
    inhomogeneous_term = [c * total_pool_size 0]

    # Solve for particular solution
    particular_soln = matrix \ - inhomogeneous_term'

    # Solve for constants
    constants = eigenvectors \ -particular_soln

    # Extract the eigenvalues and eigenvectors
    eigvalue1 = eigenvalues[1]
    eigvalue2 = eigenvalues[2]
    eigvec1 = eigenvectors[:, 1]
    eigvec2 = eigenvectors[:, 2]

    # Calculate kappa
    kappa = - eigvalue2 * constants[2] * sum(eigvec2) / (eigvalue1 * constants[1] * sum(eigvec1))

    return kappa
end

# Example parameters (you can adjust the ranges)
total_pool_size = 1000
m_values = 1:1:100.0
c_values = 1:1:100.0
e_values = 1:1:100.0
i_values = 1:1:100.0

# Search for values where kappa > 0
for m in m_values, c in c_values, e in e_values, i in i_values
    kappa = compute_kappa(m, c, e, i, total_pool_size)
    if kappa > 0
        println("Found positive kappa: m=$m, c=$c, e=$e, i=$i, kappa=$kappa")
    end
end




