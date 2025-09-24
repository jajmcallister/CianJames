




# function mutual_info_from_weights(weights)
#     if isempty(weights)
#         return 0.0
#     end
#     # closed-form formula
#     return 0.5 * log(1 + sum(abs2, weights))
# end

# function mutual_info_over_time(weight_dists)
#     return [mutual_info_from_weights(w) for w in weight_dists]
# end


# mi_vals_trials = [mutual_info_over_time(weight_dists) for weight_dists in trials_synapse_sizes]

# plot(0:0.2:120-0.2, mean(mi_vals_trials), ribbon=std(mi_vals_trials)/sqrt(length(mi_vals_trials)), xlabel="Days", xticks=0:20:120, ylabel="Mutual Information", lw=2, legend=false)


function differential_entropy_over_time(synapse_sizes)
    entropies = Float64[]
    for w in synapse_sizes
        σ² = sum(abs2, w)
        push!(entropies, 0.5 * log(2π * exp(1) * σ²))
    end
    return entropies
end

ent_vals_trials = [differential_entropy_over_time(weight_dists) for weight_dists in trials_synapse_sizes]

p = plot(0:0.2:120-0.2, mean(ent_vals_trials), ribbon=std(ent_vals_trials)/sqrt(length(ent_vals_trials)), c=:red, 
        lw=4, xlabel="Days", xticks=0:20:120, ylabel="Differential Entropy",legend=false, ylim=(-3,4), grid=false, title="Differential entropy of output y",dpi=600)

# savefig(p, "C://Users/B00955735/OneDrive - Ulster University/Desktop/diff_entropy_output.png")

using LinearAlgebra, StatsBase
function simulated_entropy(weights; N=10_000, nbins=60)
    if isempty(weights)
        return NaN
    end
    # sample x_i ~ N(0,1), compute y = sum(w_i * x_i)
    ys = [dot(weights, randn(length(weights))) for _ in 1:N]
    h = fit(Histogram, ys, nbins=nbins)
    p = h.weights ./ sum(h.weights)
    return -sum(p .* log.(p .+ eps()))
end

function simulated_entropy_mean(weights; N=10_000, nbins=60, reps=20)
    if isempty(weights)
        return NaN
    end
    vals = Float64[]
    for _ in 1:reps
        ys = [dot(weights, randn(length(weights))) for _ in 1:N]
        h = fit(Histogram, ys; nbins=nbins)
        p = h.weights ./ sum(h.weights)
        push!(vals, -sum(p .* log.(p .+ eps())))
    end
    return mean(vals)
end



simulated_entropies = [simulated_entropy_mean.(trials_synapse_sizes[i]) for i in 1:length(trials_synapse_sizes)]

plot(mean(simulated_entropies))


var_values = [var.(w) for w in trials_synapse_sizes]
plot(0:0.2:120-0.2, mean(var_values), ribbon=std(var_values)/sqrt(length(var_values)))



function lognormal_y_mean_approx(weights,mu,sigma)
    if isempty(weights)
        return NaN
    end
    sum_w = sum(weights)
    return sum_w * exp(mu + 0.5*sigma^2)
end

function lognormal_y_variance_approx(weights,mu,sigma)
    if isempty(weights)
        return NaN
    end
    sum_w2 = sum(abs2, weights)
    return sum_w2 * (exp(sigma^2) - 1) * exp(2*mu + sigma^2)
end

function lognormal_y_simulate_mean(weights,mu,sigma; N=10_000, nbins=60)
    if isempty(weights)
        return NaN
    end
    vals = []
    for _ in 1:N
        xs = rand(LogNormal(mu,sigma), length(weights))
        y = dot(weights, xs)
        push!(vals, y)
    end
    return mean(vals)
end

function lognormal_y_variance(weights,mu,sigma; N=10_000, nbins=60)
    if isempty(weights)
        return NaN
    end
    vals = []
    for _ in 1:N
        xs = rand(LogNormal(mu,sigma), length(weights))
        y = dot(weights, xs)
        push!(vals, y)
    end
    return var(vals)
end

function entropy_lognormal_approx(weights, mean, var)
    if isempty(weights)
        return NaN
    end
    
    return 0.5 * log(2π * exp(1) * var)
end

lognormal_means_approx = [lognormal_y_mean_approx.(trials_synapse_sizes[i],0.0,1.0) for i in 1:length(trials_synapse_sizes)]
lognormal_variances_approx = [lognormal_y_variance_approx.(trials_synapse_sizes[i],0.0,1.0) for i in 1:length(trials_synapse_sizes)]
lognormal_means_simulated = [lognormal_y_simulate_mean.(trials_synapse_sizes[i],0.0,1.0) for i in 1:length(trials_synapse_sizes)]
lognormal_variances_simulated = [lognormal_y_variance.(trials_synapse_sizes[i],0.0,1.0) for i in 1:length(trials_synapse_sizes)]


p1 = plot(0:0.2:120-0.2, mean(lognormal_means_approx), ribbon=std(lognormal_means_approx)/sqrt(length(lognormal_means_approx)), xlabel="Days", xticks=0:20:120, ylabel="Mean of y", lw=2, legend=false, title="Mean of y with lognormal inputs (μ=0, σ=1)")
plot!(p1, 0:0.2:120-0.2, mean(lognormal_means_simulated), ribbon=std(lognormal_means_simulated)/sqrt(length(lognormal_means_simulated)), c=:red, lw=2, label="Simulated")

p2 = plot(0:0.2:120-0.2, mean(lognormal_variances_approx), ribbon=std(lognormal_variances_approx)/sqrt(length(lognormal_variances_approx)), xlabel="Days", xticks=0:20:120, ylabel="Variance of y", lw=2, legend=false, title="Variance of y with lognormal inputs (μ=0, σ=1)")
plot!(p2, 0:0.2:120-0.2, mean(lognormal_variances_simulated), ribbon=std(lognormal_variances_simulated)/sqrt(length(lognormal_variances_simulated)), c=:red, lw=2, label="Simulated")


plot(p1, p2, layout=(1,2), dpi=600, size=(1000,400),margin=5mm,grid=false)

lognormal_entropies_approx = [entropy_lognormal_approx.(trials_synapse_sizes[i], lognormal_means_approx[i], lognormal_variances_approx[i]) for i in 1:length(trials_synapse_sizes)]

p = plot(0:0.2:120-0.2, mean(lognormal_entropies_approx), ribbon=std(lognormal_entropies_approx)/sqrt(length(lognormal_entropies_approx)), lw=4, c=:purple, xlabel="Days", 
xticks=0:20:120, ylabel="Differential Entropy",label="Lognormal input", grid=false, title="Differential Entropy of output y", ylim=(-3,5), dpi=600)

savefig(p, "C://Users/B00955735/OneDrive - Ulster University/Desktop/diff_entropy_approx.png")

plot!(0:0.2:120-0.2, mean(ent_vals_trials), ribbon=std(ent_vals_trials)/sqrt(length(ent_vals_trials)), c=:red, 
        lw=4, xlabel="Days", xticks=0:20:120, ylabel="Differential Entropy",label="Normal input", ylim=(-3,5), grid=false, title="Differential entropy of output y",dpi=600)


sim_entropies_trials = []
for k in 1:length(trials_synapse_sizes)

    sim_dists_normal = []
    sim_entropies_normal = []

    Nn = 1000
    for i in 1:length(trials_synapse_sizes[1])
        if isempty(trials_synapse_sizes[1][i])
            push!(sim_entropies_normal, NaN)
        else
            ys = Float64[]
        
            for j in 1:Nn
                x = randn(length(trials_synapse_sizes[1][i]))
                y = dot(trials_synapse_sizes[1][i], x)
                push!(ys, y)
            end

            h = histogram_entropy(ys)
            push!(sim_entropies_normal, h)
            push!(sim_dists_normal, ys)
        end
        
    end
    push!(sim_entropies_trials, sim_entropies_normal)
end

plot(mean(sim_entropies_trials))









# paul data
function expand_histogram(x::Vector{Float64}, y::Vector{Float64})
    # round counts to integers
    counts = round.(Int, y)
    # repeat each x[i] count[i] times
    return vcat([fill(x[i], counts[i]) for i in eachindex(x)]...)
end

# Convert dict into weight_dists array ordered by keys
keys_order = ["P15", "P35", "P55", "P120"]
weight_dists = [expand_histogram(data[k][1], data[k][2]) for k in keys_order]

mi_vals_paul = mutual_info_over_time(weight_dists)
plot(mi_vals_paul)