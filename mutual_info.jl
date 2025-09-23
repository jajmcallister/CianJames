




function mutual_info_from_weights(weights)
    if isempty(weights)
        return 0.0
    end
    # closed-form formula
    return 0.5 * log(1 + sum(abs2, weights))
end

function mutual_info_over_time(weight_dists)
    return [mutual_info_from_weights(w) for w in weight_dists]
end


mi_vals_trials = [mutual_info_over_time(weight_dists) for weight_dists in trials_synapse_sizes]

plot(0:0.2:120-0.2, mean(mi_vals_trials), ribbon=std(mi_vals_trials)/sqrt(length(mi_vals_trials)), xlabel="Days", xticks=0:20:120, ylabel="Mutual Information", lw=2, legend=false)



function entropy_over_time(synapse_sizes)
    entropies = Float64[]
    for w in synapse_sizes
        σ² = sum(abs2, w)
        push!(entropies, 0.5 * log(2π * exp(1) * σ²))
    end
    return entropies
end

ent_vals_trials = [entropy_over_time(weight_dists) for weight_dists in trials_synapse_sizes]

p = plot(0:0.2:120-0.2, mean(ent_vals_trials), ribbon=std(ent_vals_trials)/sqrt(length(ent_vals_trials)), c=:red, 
        lw=4, xlabel="Days", xticks=0:20:120, ylabel="Differential Entropy",legend=false, ylim=(-3,4), grid=false, title="Differential entropy of output y",dpi=600)

savefig(p, "C://Users/B00955735/OneDrive - Ulster University/Desktop/diff_entropy_output.png")






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