using QuadGK
using Plots
using Statistics

heaviside(x,threshold) = x >= threshold ? 1 : 0
function pareto_dist(x,alpha,x_m,threshold)
    return (alpha * x_m ^ alpha / x ^ (alpha + 1)) * heaviside(x, threshold)
end

function exp_dist(x,A,lambda)
    return A*exp(-x/lambda)
end


alpha = 0.01
x_m = 10
lambda = 0.1
A = 0.5
threshold = 0.0001

xvals=0:0.01:2



integral, error = quadgk(x -> pareto_dist(x,alpha,x_m,threshold) * exp_dist(x,A,lambda), 0, Inf)

plot(xvals, pareto_dist.(xvals,alpha,x_m,threshold))
plot!(xvals, exp_dist.(xvals,A,lambda))
plot!(xvals, pareto_dist.(xvals,alpha,x_m,threshold) .* exp_dist.(xvals,A,lambda))
