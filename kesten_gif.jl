using Random, Distributions, StatsPlots

# Parameters
ε = 0.9923
η = 1 - ε
σ_ε, σ_η = 0.05, 0.03
N = 1000                # number of synapses
T = 1000                 # number of timesteps

function kesten_stationary_moments(ε, ση, σε)
    # parameters
    η = 1 - ε
    a = ε
    α = ε^2 + σε^2
    b = η
    β = η^2 + ση^2

    # mean (if |a|<1)
    m = b / (1 - a)

    # second moment (if α<1)
    q = (2*a*b^2/(1 - a) + β) / (1 - α)

    # variance
    var = q - m^2

    return (mean=m, var=var, second_moment=q)
end

ε = 0.9923
ση, σε = 0.03, 0.05

stats = kesten_stationary_moments(ε, ση, σε)
println(stats)


# Initialize weights
weights = zeros(N)

# To store snapshots
snapshots = []

# Evolution loop
for t in 1:T
    # Draw noise
    ε_noise = rand.(Normal(0, σ_ε), N)
    η_noise = rand.(Normal(0, σ_η), N)

    # Update rule: w_{t+1} = (ε + noise) * w_t + (η + noise)
    weights .= (ε .+ ε_noise) .* weights .+ (η .+ η_noise)

    # Clip at zero (no negative weights)
    weights .= max.(weights, 0.0)

    # Store snapshot every few steps
    if t % 10 == 0
        push!(snapshots, copy(weights))
    end
end

mean(weights)
var(weights)

default(tickdirection=:out)
# Animate histogram evolution
@gif for (i, ws) in enumerate(snapshots)
    histogram(ws, bins=0:0.1:6, ylim=(0,200), grid=false, lw=0.2,
              xlabel="Synapse weight", legend=false, dpi=600,
              title="Kesten Process - Synapse Weight Distribution",
              xlim=(0, 6))
end every 1

hh = histogram(weights, bins=0:0.2:6, normalize=true,legend=false,c=:steelblue, title="Synaptic Weights", lw=0.3, grid=false, dpi=600)
# savefig(hh, "C://Users/B00955735/OneDrive - Ulster University/Desktop/kesten_weights.png")




# Precompute entropy for each snapshot
wss = [ws for ws in snapshots]
entropies = [histogram_entropy(filter(>(1), ws), nbinss=60) for ws in wss]
plot(entropies)

time_snap = 0:5:(length(snapshots)-1)*5  # match the snapshot interval you used

@gif for i in 1:length(snapshots)
    plot(layout=(2,1), size=(600,800))

    # Top: histogram
    histogram!(snapshots[i], bins=0:0.1:6, ylim=(0,200), grid=false, lw=0.2,
               xlabel="Synapse weight", legend=false,
               title="Kesten Process - Synapse Weight Distribution",
               xlim=(0,6), subplot=1)

    # Bottom: entropy line up to current snapshot
    plot!(time_snap[1:i], entropies[1:i],
          title="Entropy of synaptic weight distribution over time",
          c=:seagreen, lw=4,
          xlabel="Days", ylabel="Entropy (bits)",
          legend=false, grid=false,
          xticks=false,
          subplot=2)
end every 1


















# --- 2) Estimate alpha from E[A^α]=1 ---
M = 200000
A_samples = ε .+ σ_ε .* randn(M)

function g(α)
    mean(abs.(A_samples) .^ α) - 1.0
end

function find_alpha(lo::Float64, hi::Float64; tol=1e-6, maxiter=100)
    glo, ghi = g(lo), g(hi)
    if glo * ghi > 0
        error("Root not bracketed: try enlarging interval")
    end
    for _ in 1:maxiter
        mid = 0.5 * (lo + hi)
        gmid = g(mid)
        if abs(gmid) < tol
            return mid
        elseif glo * gmid <= 0
            hi, ghi = mid, gmid
        else
            lo, glo = mid, gmid
        end
    end
    return 0.5 * (lo + hi)
end

# Quadratic approximation of Kesten tail exponent
function kesten_alpha_quad(eps, sigma_eps)
    return 1 - (2 * eps^2 * log(eps)) / (sigma_eps^2)
end

alpha_est = kesten_alpha_quad(ε, σ_ε)

α = find_alpha(1e-3, 10.0)


x = 1:0.1:3
# power-law decay ~ x^(-α)
y = x .^ (-α)

# plot on log–log axes

histogram(weights,normalize=true,bins=0:0.1:8)
plot!(x, y, lw=4)


###############################


using Plots, StatsBase

# Histogram of your weights
h = fit(Histogram, weights, 0:0.1:8)
bin_centers = (h.edges[1][1:end-1] .+ h.edges[1][2:end]) ./ 2
pdf_vals = h.weights ./ sum(h.weights) ./ step(h.edges[1])  # normalized density

# Power-law reference
α = find_alpha(1e-3, 10.0)   # your root-finding estimate
x = bin_centers
y = x .^ (-α)

# Normalize the power law so it matches histogram scale
scale = pdf_vals[1] / y[1]
y_scaled = scale .* y

# Plot
plot(bin_centers, pdf_vals, seriestype=:line, label="Histogram (PDF)")
plot!(x, y_scaled, lw=3, label="Power law ~ x^-$(round(α,digits=2))")
plot!(xaxis=:log, yaxis=:log)
