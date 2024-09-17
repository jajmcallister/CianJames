
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



# Compute imm and mmm over the time range
results = solution_total.(time_for_plot,c,m,e,i)

# Separate imm and mmm values
immature_values = map(r -> r[1], results)  # Extract imm values
mature_values = map(r -> r[2], results)  # Extract mmm values

# Plot the results
plot(time_for_plot, immature_values, label="imm", xlabel="Time", ylabel="Value", legend=:topright)
plot!(time_for_plot, mature_values, label="mat", ylim=(0,1000),legend=false)
