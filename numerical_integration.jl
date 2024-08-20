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

function log_normal(x,sigma,mu)
    return (1/(x*sigma*sqrt(2*pi)))*(exp(-((log(x)-mu)^2)/(2*sigma^2)))
end


alpha = 0.05
x_m = 100
lambda = 2
A = 0.05
threshold = 0.000
mu = 0.0
sigma = 1

xvals=0.01:0.01:10
v = rand(0:10,10)

mean(exp_dist.(v,A,1))
plot(xvals, exp_dist.(xvals,A,1))
hline!([mean(exp_dist.(xvals,A,1))])


comb_plot = plot(xvals, log_normal.(xvals,sigma,mu),lw=3, label="Log-normal limiting distribution",color=:orange)
plot!(xvals, exp_dist.(xvals,A,lambda),lw=3, label="Exponential pdf for weight-dependent dematuring", color=:green, xlabel="Synapse size", ylabel="Probability")
plot!(xvals, log_normal.(xvals,sigma,mu) .* exp_dist.(xvals,A,lambda),label="Product of other two functions to be integrated", lw=3, color=:purple,xlim=(0,5))

# savefig(comb_plot, "C://Users/B00955735/OneDrive - Ulster University/Desktop/comb_plot.svg")

integral, error = quadgk(x -> (1/(x*sigma*sqrt(2*pi)))*(exp(-((log(x)-mu)^2)/(2*sigma^2)))* exp_dist(x,A,lambda), 0.0, Inf)


integral2, error2 = quadgk(x -> pareto_dist(x,alpha,x_m,threshold)*(exp(-((log(x)-mu)^2)/(2*sigma^2)))* exp_dist(x,A,lambda), 0.0, Inf)
