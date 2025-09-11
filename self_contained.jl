using Random, StatsBase, DifferentialEquations, Distributions, Statistics, Plots, KernelDensity


#0.9, 0.03333333333333333, 0.2, 3, 0.17500000000000002, 0.2, 0.05, 0.05, 3.
#1.055102040816327, 0.046938776, 2.1020408163265314, 0.12448979591836733, 0.066326531, 0.066326531, 0.5897959183673469


function new_compute_survival_fraction007(state_records)
    total_time_steps = size(state_records, 2)

    # Identify the initial population
    initial_existing_synapses = findall(x -> x != 0, state_records[:, 1])
    
    # Create a set to track permanently eliminated synapses
    permanently_eliminated = Set{Int}()

    # Initialize array to store the survival fraction over time
    survival_fraction = zeros(total_time_steps)

    # Loop over each time step and compute the fraction of surviving synapses
    for t in 1:total_time_steps
        # Update the set of permanently eliminated synapses
        for idx in initial_existing_synapses
            if state_records[idx, t] == 0
                push!(permanently_eliminated, idx)
            end
        end
        
        # Find how many of the initial synapses are still in the existing population and not permanently eliminated
        surviving_synapses = count(x -> !(x in permanently_eliminated) && state_records[x, t] != 0, initial_existing_synapses)
        
        # Compute survival fraction as the ratio of surviving synapses to the initial population size
        survival_fraction[t] = surviving_synapses / length(initial_existing_synapses)
    end

    return survival_fraction
end

function safe_sample007(list, num_samples; replace=false)
    return sample(list, min(num_samples, length(list)), replace=false)
end

function kesten_update_new007(sizes, ε, η, σ_ε, σ_η)
    sizes=deepcopy(sizes)
    for i in 1:length(sizes)
        ε_i = rand(Normal(ε, σ_ε))
        η_i = rand(Normal(η, σ_η))
        new_size = ε_i * sizes[i] + η_i
        sizes[i]=new_size
        if sizes[i] < 0
            sizes[i] = 0.01
        end
        # sizes[i] = max(new_size, 0.0)  # Ensure size is non-negative
        if sizes[i] < 0.0
            sizes[i]=0.01
        end
    end

    return sizes
end

function time_average007(data, window)
    # # Check if the window size is valid
    # if window < 1 || window > length(data)
    #     throw(ArgumentError("Window size must be between 1 and the length of the data."))
    # end

    # Compute the moving average
    averaged_data = [mean(data[max(1, i - window + 1):i]) for i in 1:length(data)]
    return averaged_data
end

function track_times_variable_rates_007(total_time, total_pool_size, rates, ε, η, σ_ε, σ_η, kesten_time_step)
    pool = total_pool_size
    steps = trunc(Int, total_time / kesten_time_step)
    immature = 0
    mature = 0
    pool_history = []
    immature_history = []
    mature_history = []

    state_records = zeros(total_pool_size, trunc(Int,total_time/kesten_time_step))
    state_records_heatmap = zeros(total_pool_size, trunc(Int,total_time/kesten_time_step))

    I_pop = []
    M_pop = []
    pool_pop = [i for i in 1:total_pool_size]
    cr, m, el, A, lambda  = rates

    push!(pool_history, pool)
    push!(immature_history, immature)
    push!(mature_history, mature)
    
    # Synapse sizes for mature population
    synapse_sizes = Float64[]
    synapse_size_history = []
    
    
    # Simulation
    for t in 1:steps
        # 1 Transitions from pool to immature
        pool_to_immature = 0
        transition_prob1 = cr[t] * kesten_timestep
        for i in 1:pool
            if rand() < transition_prob1
                pool_to_immature += 1
            end
        end
        pool -= pool_to_immature
        immature += pool_to_immature
    
        # 2 Transitions from immature to mature
        immature_to_mature = 0
        transition_prob2 = m * kesten_timestep  # Probability for each immature synapse to become mature

        for i in 1:immature
            if rand() < transition_prob2
                immature_to_mature += 1
            end
        end

        # Update the counts
        immature -= immature_to_mature
        mature += immature_to_mature
    
        # Initialize new mature synapse sizes
        for i in 1:immature_to_mature
            push!(synapse_sizes, 0.01)  # Initial size of a new mature synapse
        end
    
        synapse_sizes = sort(synapse_sizes, rev=true)
    
        # 3 Transitions from mature to immature
        mature_to_immature_indices = []
        for (id, size) in enumerate(synapse_sizes)
            prob = A * exp(-size / lambda) * kesten_timestep
            if rand() < prob
                push!(mature_to_immature_indices, id)
            end
        end

        
    
        # Update states based on calculated probabilities
        mature_to_immature = length(mature_to_immature_indices)
        mature_to_immature = round(Int, mature_to_immature)
        mature -= mature_to_immature
        immature += mature_to_immature
    
        # Remove synapse sizes for synapses that became immature
        # Sort indices in reverse order
        sorted_indices = sort(mature_to_immature_indices, rev=true)
    
        # Delete elements at the specified indices
        for idx in sorted_indices
            deleteat!(synapse_sizes, idx)
        end
    
        # 4 Transitions from immature to pool
        immature_to_pool = 0
        transition_prob3 = el[t] * kesten_timestep  # Probability for each immature synapse to transition to the pool

        for i in 1:immature
            if rand() < transition_prob3
                immature_to_pool += 1
            end
        end

        # Update the counts
        immature -= immature_to_pool
        pool += immature_to_pool

        
        push!(pool_history, pool)
        push!(immature_history, immature)
        push!(mature_history, mature)


        # 1. Pool to immature
        pool_to_immature_count = abs(trunc(Int, pool_to_immature))
        pool_to_immature_indxs = safe_sample007(pool_pop, pool_to_immature_count, replace=false)
        filter!(x -> !in(x, pool_to_immature_indxs), pool_pop)
        append!(I_pop, pool_to_immature_indxs)
        
        # 2. Immature to mature
        immature_to_mature_count = trunc(Int, immature_to_mature)
        immature_to_mature_indxs = safe_sample007(I_pop, immature_to_mature_count, replace=false)
        filter!(x -> !in(x, immature_to_mature_indxs), I_pop)
        append!(M_pop, immature_to_mature_indxs)
        
        # 3. Mature to immature
        mature_to_immature_count = abs(trunc(Int, mature_to_immature))
        mature_to_immature_indxs = safe_sample007(M_pop, mature_to_immature_count, replace=false)
        filter!(x -> !in(x, mature_to_immature_indxs), M_pop)
        append!(I_pop, mature_to_immature_indxs)
        
        # 4. Immature to pool
        immature_to_pool_count = trunc(Int, immature_to_pool)
        immature_to_pool_indxs = safe_sample007(I_pop, immature_to_pool_count, replace=false)
        filter!(x -> !in(x, immature_to_pool_indxs), I_pop)
        append!(pool_pop, immature_to_pool_indxs)


        synapse_sizes = kesten_update_new007(synapse_sizes, ε, η, σ_ε, σ_η)
    

        push!(synapse_size_history, synapse_sizes)

        # Now recording the current state of the synapses in the state record matrix
        for j in 1:total_pool_size
            if j in pool_pop
                state_records[j,t] = 0
            elseif j in I_pop
                state_records[j,t] = 1
            elseif j in M_pop
                state_records[j,t] = 2
            end
        end

        # Record the current state in the state_records_heatmap matrix
        for j in 1:total_pool_size
            if j in pool_pop
                state_records_heatmap[j, t] = 0 # in the pool, size = 0
            elseif j in I_pop
                state_records_heatmap[j, t] = 0  # in immature population size is 0
            elseif j in M_pop
                idx = findfirst(==(j), M_pop)
                if idx !== nothing
                    state_records_heatmap[j, t] = synapse_sizes[idx]
                end
            end
        end
    
    end

    return immature_history, mature_history, state_records, synapse_sizes, state_records_heatmap, synapse_size_history
end

# function track_times_variable_rates_and_kesten_007(total_time, total_pool_size, rates, εs, ηs, σ_εs, σ_ηs, kesten_time_step)
#     pool = total_pool_size
#     steps = trunc(Int, total_time / kesten_time_step)
#     immature = 0
#     mature = 0
#     pool_history = []
#     immature_history = []
#     mature_history = []

#     state_records = zeros(total_pool_size, trunc(Int,total_time/kesten_time_step))
#     state_records_heatmap = zeros(total_pool_size, trunc(Int,total_time/kesten_time_step))

#     I_pop = []
#     M_pop = []
#     pool_pop = [i for i in 1:total_pool_size]
#     cr, m, el, A, lambda  = rates

#     push!(pool_history, pool)
#     push!(immature_history, immature)
#     push!(mature_history, mature)
    
#     # Synapse sizes for mature population
#     synapse_sizes = Float64[]
#     synapse_size_history = []
    
    
#     # Simulation
#     for t in 1:steps
#         # 1 Transitions from pool to immature
#         pool_to_immature = 0
#         transition_prob1 = cr[t] * kesten_timestep
#         for i in 1:pool
#             if rand() < transition_prob1
#                 pool_to_immature += 1
#             end
#         end
#         pool -= pool_to_immature
#         immature += pool_to_immature
    
#         # 2 Transitions from immature to mature
#         immature_to_mature = 0
#         transition_prob2 = m * kesten_timestep  # Probability for each immature synapse to become mature

#         for i in 1:immature
#             if rand() < transition_prob2
#                 immature_to_mature += 1
#             end
#         end

#         # Update the counts
#         immature -= immature_to_mature
#         mature += immature_to_mature
    
#         # Initialize new mature synapse sizes
#         for i in 1:immature_to_mature
#             push!(synapse_sizes, 0.0)  # Initial size of a new mature synapse
#         end
    
#         synapse_sizes = sort(synapse_sizes, rev=true)
    
#         # 3 Transitions from mature to immature
#         mature_to_immature_indices = []
#         for (id, size) in enumerate(synapse_sizes)
#             prob = A * exp(-size / lambda) * kesten_timestep #i*kesten_time_step
#             if rand() < prob
#                 push!(mature_to_immature_indices, id)
#             end
#         end

        
    
#         # Update states based on calculated probabilities
#         mature_to_immature = length(mature_to_immature_indices)
#         mature_to_immature = round(Int, mature_to_immature)
#         mature -= mature_to_immature
#         immature += mature_to_immature
    
#         # Remove synapse sizes for synapses that became immature
#         # Sort indices in reverse order
#         sorted_indices = sort(mature_to_immature_indices, rev=true)
    
#         # Delete elements at the specified indices
#         for idx in sorted_indices
#             deleteat!(synapse_sizes, idx)
#         end
    
#         # 4 Transitions from immature to pool
#         immature_to_pool = 0
#         transition_prob3 = el[t] * kesten_timestep  # Probability for each immature synapse to transition to the pool

#         for i in 1:immature
#             if rand() < transition_prob3
#                 immature_to_pool += 1
#             end
#         end

#         # Update the counts
#         immature -= immature_to_pool
#         pool += immature_to_pool

        
#         push!(pool_history, pool)
#         push!(immature_history, immature)
#         push!(mature_history, mature)


#         # 1. Pool to immature
#         pool_to_immature_count = abs(trunc(Int, pool_to_immature))
#         pool_to_immature_indxs = safe_sample007(pool_pop, pool_to_immature_count, replace=false)
#         filter!(x -> !in(x, pool_to_immature_indxs), pool_pop)
#         append!(I_pop, pool_to_immature_indxs)
        
#         # 2. Immature to mature
#         immature_to_mature_count = trunc(Int, immature_to_mature)
#         immature_to_mature_indxs = safe_sample007(I_pop, immature_to_mature_count, replace=false)
#         filter!(x -> !in(x, immature_to_mature_indxs), I_pop)
#         append!(M_pop, immature_to_mature_indxs)
        
#         # 3. Mature to immature
#         mature_to_immature_count = abs(trunc(Int, mature_to_immature))
#         mature_to_immature_indxs = safe_sample007(M_pop, mature_to_immature_count, replace=false)
#         filter!(x -> !in(x, mature_to_immature_indxs), M_pop)
#         append!(I_pop, mature_to_immature_indxs)
        
#         # 4. Immature to pool
#         immature_to_pool_count = trunc(Int, immature_to_pool)
#         immature_to_pool_indxs = safe_sample007(I_pop, immature_to_pool_count, replace=false)
#         filter!(x -> !in(x, immature_to_pool_indxs), I_pop)
#         append!(pool_pop, immature_to_pool_indxs)
        
#         ε, η, σ_ε, σ_η = εs[t], ηs[t], σ_εs[t], σ_ηs[t]

#         synapse_sizes = kesten_update_new007(synapse_sizes, ε, η, σ_ε, σ_η)
    

#         push!(synapse_size_history, synapse_sizes)

#         # Now recording the current state of the synapses in the state record matrix
#         for j in 1:total_pool_size
#             if j in pool_pop
#                 state_records[j,t] = 0
#             elseif j in I_pop
#                 state_records[j,t] = 1
#             elseif j in M_pop
#                 state_records[j,t] = 2
#             end
#         end

#         # Record the current state in the state_records_heatmap matrix
#         for j in 1:total_pool_size
#             if j in pool_pop
#                 state_records_heatmap[j, t] = 0 # in the pool, size = 0
#             elseif j in I_pop
#                 state_records_heatmap[j, t] = 0  # in immature population size is 0
#             elseif j in M_pop
#                 idx = findfirst(==(j), M_pop)
#                 if idx !== nothing
#                     state_records_heatmap[j, t] = synapse_sizes[idx]
#                 end
#             end
#         end
    
#     end

#     return immature_history, mature_history, state_records, synapse_sizes, state_records_heatmap
# end

total_pool_size = 100
total_time = 120
kesten_timestep = 0.1


# a1 = 0.9
# k1 = 1/30
# b1 = 0.2
# a2 = 1.8
# k2 = 1/10
# b2 = 0.2

b1 = 0.2
b2 = 0.2

# a1,k1,a2,k2,m,A,lambda = 0.51362973760933, 0.05301263362487851, 1.460204081632653, 0.13542274052478132, 0.07361516034985421, 0.0736151603498542, 0.629883381924198
A1,lambda1,A2,lambda2,m,A3,lambda3 = .9, 30, 2, 5, 0.05, 0.05, 2.
A1,lambda1,A2,lambda2,m,A3,lambda3 = .9, 30, 2, 5, 0.04557090317553659, 0.06622531490863738, 3.3026174908321067
ε = .985
η = 1-ε
σ_ε, σ_η = .1, .1


creation_func(t,A1,lambda1) = A1 * exp(-t / lambda1) + b1
elimination_func(t,A2,lambda2) = A2 * exp(-t / lambda2) + b2
de_maturation_func(t,A3,lambda3) = A3 * exp(-t / lambda3)


creat = creation_func.(0:kesten_timestep:total_time, A1, lambda1)
elim = elimination_func.(0:kesten_timestep:total_time, A2, lambda2)

ec_plot = plot(0:kesten_timestep:total_time, creat,lw=3,label="Creation rate")
plot!(0:kesten_timestep:total_time, elim, lw=3, label="Elimination rate", ylabel="Rate", xlabel="Days")

plot(de_maturation_func.(0:0.01:10, A3, lambda3), lw=3, label="De-maturation rate", ylabel="Rate", xlabel="Days", title="De-maturation rate over time")
# savefig(ec_plot, "C://Users/B00955735/OneDrive - Ulster University/Desktop/ec_rate.png")
# creat = sigmoid_creation.(0:kesten_timestep:total_time)


rates_var = creat, m, elim, A3, lambda3


state_recs_var_multiple = []
ihs = []
mhs = []
synapse_sizes_multiple = []
syn_size_heatmaps_trials = []
synapse_size_history_multiple = []

num_trials = 2

for i in 1:num_trials
    ih_var, mh_var, state_record_var, syn_sizes_var, syn_heatmap, syn = track_times_variable_rates_007(total_time, total_pool_size, rates_var, ε, η, σ_ε, σ_η, kesten_timestep);
    push!(state_recs_var_multiple, state_record_var)
    push!(ihs, ih_var)
    push!(mhs, mh_var)
    push!(synapse_sizes_multiple, syn_sizes_var)
    push!(syn_size_heatmaps_trials,syn_heatmap)
    push!(synapse_size_history_multiple, syn)
end

syns_ht = heatmap(syn_size_heatmaps_trials[1],clims=(0,3),xticks=(0:2000:12000,0:20:120),ylabel="Synapse ID", xlabel="Days",colorbar_title="Synapse size")
# savefig(syns_ht, "C://Users/B00955735/OneDrive - Ulster University/Desktop/syns_ht.png")


# Plots of synapse weights over time (P15, 25, 55 and 120)
v1 = vcat([synapse_size_history_multiple[i][trunc(Int,15/kesten_timestep)] for i in 1:num_trials]...)
v2 = vcat([synapse_size_history_multiple[i][trunc(Int,35/kesten_timestep)] for i in 1:num_trials]...)
v3 = vcat([synapse_size_history_multiple[i][trunc(Int,55/kesten_timestep)] for i in 1:num_trials]...)
v4 = vcat([synapse_size_history_multiple[i][trunc(Int,120/kesten_timestep)] for i in 1:num_trials]...) 

# Compute kernel density estimate
bw=0.1
kde_result1 = kde(v1,bandwidth=bw)
kde_result2 = kde(v2,bandwidth=bw)
kde_result3 = kde(v3,bandwidth=bw)
kde_result4 = kde(v4,bandwidth=bw)

color1 = RGB(176/255, 178/255, 240/255)
color2 = RGB(78/255,84/255,182/255)
color3 = RGBA(221/255, 221/255, 221/255, 255/255)
color4 = RGB(4/255, 4/255, 4/255)

distsplot = plot(kde_result1.x, kde_result1.density, linewidth=4, label="P15", color=color1)
plot!(kde_result2.x, kde_result2.density, linewidth=4, label="P35", color=color2)
plot!(kde_result3.x, kde_result3.density, linewidth=4, label="P55", color=color3)
plot!(kde_result4.x, kde_result4.density, linewidth=4, label="P120", color=color4,
        xlabel="Synaptic weight (a.u.)", grid=false, xlim=(-0.2,5), title="Distributions of synaptic weight over time",
        legendfontsize=12, ylabel="Density")

# savefig(distsplot, "C://Users/B00955735/OneDrive - Ulster University/Desktop/distplot1.png")
# savefig(distsplot, "C://Users/B00955735/OneDrive - Ulster University/Desktop/dist_plot.svg")

timepoints = [15,25,35,45,55,90,120]
synapticweightstime = [vcat([synapse_size_history_multiple[i][trunc(Int,tt/kesten_timestep)] for i in 1:num_trials]...) for tt in timepoints]

A_new = [[] for i in 1:7]
for i in 1:7
    for j in 1:size(synapticweightstime[i],1)
        if synapticweightstime[i][j] == 0.0
            push!(A_new[i], 0.01)
        else
            push!(A_new[i], synapticweightstime[i][j])
        end
    end
end

A_new[1]
for z in synapticweightstime
    filter!(x -> x != 0, z)
end 

function geometric_std(X)
    μg = geomean(X) #exp(mean(log.(X)))  # Geometric mean
    return exp(sqrt(sum((log.(X) .- log(μg)) .^ 2) / length(X)))
end

function top_25_mean(X)
    sorted_X = sort(X, rev=true)  # Sort in descending order
    top_25 = sorted_X[1:ceil(Int, 0.25 * length(X))]  # Select top 25%
    return mean(top_25)  # Compute mean
end

geommeanweights = geomean.(A_new)
meanweights = mean.(A_new)
geomstdweights = geometric_std.(A_new)
stdweights = std.(A_new)
top25meanweights = top_25_mean.(A_new)
skewnessweights = [skewness(Float64.(A_new[i])) for i in 1:7]

p1 = plot(timepoints,geommeanweights,xlabel="Postnatal Day",grid=false,c=:firebrick2,title="Geometric mean", xticks=timepoints,ylabel="Synaptic weight (a.u.)",label=false,lw=5)
p2 = plot(timepoints,geomstdweights,xlabel="Postnatal Day",grid=false,c=:coral2,title="Geometric standard deviation", xticks=timepoints,ylabel="Synaptic weight (a.u.)",label=false,lw=5)
p3 = plot(timepoints,top25meanweights,c=:lightslateblue,title="Mean of top 25% of synaptic weights",grid=false,xlabel="Postnatal Day", xticks=timepoints,ylabel="Synaptic weight (a.u.)",label=false,lw=5)

p4 = plot(timepoints, meanweights,xlabel="Postnatal Day",grid=false,ylim=(0,2),c=:firebrick2,title="Arithmetic mean", xticks=timepoints,ylabel="Synaptic weight (a.u.)",label=false,lw=5)
p5 = plot(timepoints,stdweights,xlabel="Postnatal Day",grid=false,c=:coral2,ylim=(0,1),title="Arithmetic standard deviation", xticks=timepoints,ylabel="Synaptic weight (a.u.)",label=false,lw=5)

p6 = plot(timepoints,skewnessweights,xlabel="Postnatal Day",c=:black,grid=false,ylim=(0,3),title="Skewness", xticks=timepoints,ylabel="Synaptic weight (a.u.)",label=false,lw=5)

# pp = plot(p1,p2,layout=(1,2),size=(1400,700), leftmargin=10mm,bottommargin=10mm,
# xlabelfontsize=18,ylabelfontsize=18,titlefontsize=20, xtickfontsize=18, ytickfontsize=18)

# savefig(p6,"C://Users/B00955735/OneDrive - Ulster University/Desktop/skew.png")



# checking total synaptic change over time
syn_sum = [sum.(synapse_size_history_multiple[i]) for i in 1:num_trials]
mean_syn_sum = mean(syn_sum, dims=1)
plot(mean_syn_sum)




hist_matrices = []
for i in 1:num_trials
    # Define histogram bins (consistent for all columns)
    bin_edges = 0:0.1:1.5  # Adjust bin size as needed
    bin_centers = (bin_edges[1:end-1] .+ bin_edges[2:end]) / 2  # Centers for plotting

    # Compute histograms for each time point (column)
    num_bins = length(bin_edges) - 1
    num_time_points = size(syn_size_heatmaps_trials[i], 2)
    hist_matrix = zeros(num_bins, num_time_points)

    for t in 1:num_time_points
        hist = fit(Histogram, syn_size_heatmaps_trials[i][:, t], bin_edges)
        hist_matrix[:, t] = hist.weights  # Store bin counts for this time point
    end
    push!(hist_matrices, hist_matrix)
end

bin_edges = 0:0.1:1.5  # Adjust bin size as needed
bin_centers = (bin_edges[1:end-1] .+ bin_edges[2:end]) / 2  # Centers for plotting
num_bins = length(bin_edges) - 1
num_time_points = size(syn_size_heatmaps_trials[1], 2)


# Plot the heatmap
hists = heatmap(
    1:num_time_points,               # Time points on x-axis
    bin_centers,                     # Bin centers on y-axis
    mean(hist_matrices),                     # Heatmap data
    xlabel = "Time",
    ylabel = "Synapse Size",
    colorbar_title = "Count",
    title = "Histograms of synapse size across time", clims=(0,20)
)



tt = 0.01:kesten_timestep:total_time
indi = plot(tt,syn_size_heatmaps_trials[1][60,:], label="Synapse 1")
plot!(tt,syn_size_heatmaps_trials[1][4,:], label="Synapse 2")
plot!(tt, syn_size_heatmaps_trials[1][1,:], label="Synapse 3",color="green", xticks=collect(0:20:120),xlabel="Days",ylabel="Synapse size")

# savefig(indi, "C://Users/B00955735/OneDrive - Ulster University/Desktop/indi.png")





develop_survival_multiplee = []
adult_survival_multiplee = []

for state_recs in state_recs_var_multiple
    developmental_period_16 = round(Int, (16/total_time)*size(state_recs,2))
    developmental_period_26 = round(Int, (26/total_time)*size(state_recs,2))
    
    adult_period = round(Int, (100/total_time)*size(state_recs,2))
    adult_period2 = round(Int, (118/total_time)*size(state_recs,2))

    developmental_survival_fraction1 = new_compute_survival_fraction007(state_recs[:,developmental_period_16:developmental_period_26])
    adulthood_survival_fraction1 = new_compute_survival_fraction007(state_recs[:,adult_period:adult_period2])
    push!(develop_survival_multiplee, developmental_survival_fraction1)
    push!(adult_survival_multiplee, adulthood_survival_fraction1)
end

dev_ids = collect(0:1:10)
dev_ids = [round(Int, id/kesten_timestep) for id in dev_ids]

adult_ids = [0,1,2,3,4,5,6,8,10,12,14,16,17,18]
adult_ids = [round(Int, id/kesten_timestep) for id in adult_ids]
adult_ids3 = [0,1,2,3,4,5,6,8,10,12,14,16,17,18]

development_points_to_match_sim = [mean(develop_survival_multiplee)[id+1] for id in dev_ids]
development_points_to_match_data = [1.0, 0.661896208, 0.52522361,0.468246877, 0.421466905, 0.397137735, 0.376028593, 0.364221812, 0.344543843, 0.348389962, 0.340339859]

adulthood_points_to_match_sim = [mean(adult_survival_multiplee)[id+1] for id in adult_ids]
adulthood_points_to_match_data = [1.0, 0.870199702, 0.82058372, 0.788018458, 0.775729644, 0.755248343, 0.7490909229625357, 0.7400000138716264, 0.7290909507057883, 0.7163636641068893, 0.7054545315829192, 0.694545468417081, 0.688556071, 0.681643617]

development_survival_error = development_points_to_match_sim - development_points_to_match_data
adulthood_survival_error = adulthood_points_to_match_sim - adulthood_points_to_match_data

total_error = sum(development_survival_error.^2) + sum(adulthood_survival_error.^2)

total_error = sum(development_survival_error.^2)/length(development_survival_error) + sum(adulthood_survival_error.^2)/length(adulthood_survival_error)


developmentperiodplot = 16:10/length(develop_survival_multiplee[1]):26
adulthoodperiodplot = 0:18/length(adult_survival_multiplee[1]):18

l1 = length(developmentperiodplot) - length(develop_survival_multiplee[1])
l2 = length(adulthoodperiodplot) - length(adult_survival_multiplee[1])


# Plot survival fraction over time
developmental_survival_plot = plot(developmentperiodplot[1:end-l1], mean(develop_survival_multiplee),lw=5, ribbon=std(develop_survival_multiplee)/sqrt(num_trials), fillalpha=0.2, xlabel="Postnatal Day", ylabel="Survival Fraction", xticks=16:1:26,
    title="Synapse Survival Fraction (Early Development)", legend=false, ylim=(0,1.05))
    scatter!(16:1:26, development_points_to_match_data, label="Data",title="Survival Fraction (Early Development)", xlabel="Postnatal day")

adult_survival_plot = plot(adulthoodperiodplot[1:end-l2], mean(adult_survival_multiplee), ribbon=std(adult_survival_multiplee)/sqrt(num_trials), fillalpha=0.2, xlabel="Days", ylabel="Survival Fraction",
    title="Synapse Survival Fraction (Adulthood)", lw=5, legend=false, ylim=(0,1.05), xticks=0:1:18,label="Model")
    scatter!(adult_ids3, adulthood_points_to_match_data, label="Data",title="Survival Fraction (Adulthood)", xlabel="Days", legend=:bottomleft)

survival_fraction_plot = plot(developmental_survival_plot, adult_survival_plot, layout=(2,1),markersize=6, grid=false)


# savefig(survival_fraction_plot, "C://Users/B00955735/OneDrive - Ulster University/Desktop/survivalfractionplot.svg")



# work out the peak value of the combined populations and check if it occurs !at beginning and !end
smoothed_avg = time_average007(mean(ihs)+mean(mhs),100)
max_val = maximum(smoothed_avg)
end_val = smoothed_avg[end]
id_max = argmax(smoothed_avg)


p1 = plot(0:kesten_timestep:total_time,mean(mhs), ribbon=std(mhs)/sqrt(num_trials), lw=5, c=:green, label="Mature synapses", xlabel="Postnatal Day",ylabel="Number")
plot!(0:kesten_timestep:total_time, mean(ihs), ribbon=std(ihs)/sqrt(num_trials), lw=5, c=:magenta, label="Immature synapses")
plot!(0:kesten_timestep:total_time, mean(ihs).+mean(mhs), ribbon=std(ihs)/sqrt(num_trials), lw=5, linealpha=0.7, c=:grey, label="Total synapses")
plot!(title="Population Dynamics (Random Walks Model)", lw=5, c=:black, label="Total synapses",legend=:bottomright)
plot!(grid=false,legendfontsize=12,ylim=(0,total_pool_size))






############
############


function synapse_dynamics_var007!(du, u, p, t)
    A1, lambda1, A2, lambda2, m, A3, lambda3, synapse_sizes = p 
    N_I, N_M, N_P = u

    # Compute the rate of dematuration using the exponential probability distribution
    # Sum over all mature synapses' probabilities of transitioning to immature state
    if !isempty(synapse_sizes)
        dematuration_rate = A3* sum(exp(- size / lambda3) for size in synapse_sizes) / length(synapse_sizes)
    else
        dematuration_rate = 0
    end
    # dematuration_rate = i
    # Apply time-dependent e(t) and c(t)
    c_t = creation_func(t, A1, lambda1)
    e_t = elimination_func(t, A2, lambda2)
    

    du[1] = c_t * N_P - (m + e_t) * N_I + (dematuration_rate) * N_M  # dN_I/dt
    du[2] = m * N_I - (dematuration_rate) * N_M  # dN_M/dt
    du[3] = - du[1] - du[2]  # dN_P/dt
end

function run_simulation_diffeq_var007(total_time, total_pool_size, paras, ε, η, σ_ε, σ_η, kesten_time_step)
    pool = fill(1, total_pool_size);  # Initialize resource pool with synapses
    synapses = Int[]  # Array to hold states of synapses (0s and 1s)
    synapse_sizes = Float64[]  # Sizes of mature synapses
    synapse_sizes_history = []
    A1, lambda1, A2, lambda2, m, A3, lambda3 = paras
    # A, m, i, λ = paramss
    # Initial conditions
    u0 = [0.0, 0.0, total_pool_size];
    # u0 = [total_pool_size*0.4, total_pool_size*0.3, total_pool_size*0.3];
    tspan = (0.0, total_time);
    # p = (m, i, λ, synapse_sizes)
    Ihist = []
    Mhist = []

    # Define ODE problem
    # prob = ODEProblem(synapse_dynamics_var!, u0, tspan, p);

    current_time = 0.0;
    
    prams = A1, lambda1, A2, lambda2, m, A3, lambda3, synapse_sizes
    probb = ODEProblem(synapse_dynamics_var007!, u0, tspan, prams);

    while current_time < total_time
        prams = A1, lambda1, A2, lambda2, m, A3, lambda3, synapse_sizes
        probb = ODEProblem(synapse_dynamics_var007!, u0, tspan, prams);
        sol = solve(probb, Tsit5(), saveat=current_time:kesten_time_step:current_time + kesten_time_step);
        N_I, N_M, P = sol.u[end];
        push!(Ihist, N_I)
        push!(Mhist,  N_M)
        current_time += kesten_time_step;

        N_M = max(0, N_M)
        N_I = max(0, N_I)
        # Update populations:  0s for synapses in N_I and 1s for in N_M
        synapses = vcat(fill(0, round(Int, N_I)), fill(1, round(Int, N_M)));
        
        P = max(0,P)
        # 1s for synapses in the pool
        pool = fill(1, round(Int, P));

        
        if N_M > length(synapse_sizes) # If synapses have been added to the mature population since last time
            new_matures = round(Int, N_M) - length(synapse_sizes);
            append!(synapse_sizes, fill(0.01, new_matures));  # Initialize new mature synapses with size 0.01
        # elseif N_M < length(synapse_sizes) # If synapses have dematured out of the mature population
        #     num_delete_matures = length(synapse_sizes) - round(Int, N_M); #find how many need to be deleted
        #     synapse_sizes = sort(synapse_sizes); # sort the synapse size array
        #     synapse_sizes = synapse_sizes[num_delete_matures+1:end] # ... and delete the first num_delete_matures
        # end

        elseif N_M < length(synapse_sizes) # If synapses have dematured out of the mature population
            num_delete_matures = length(synapse_sizes) - round(Int, N_M) #find how many need to be deleted
            sizes = synapse_sizes
            weights = A3 * exp.(-sizes ./ lambda3)
            weights ./= sum(weights)  # Normalise to sum to 1
        
            # Sample indices to delete based on weights, without replacement
            delete_indices = sample(1:length(sizes), Weights(weights), num_delete_matures; replace=false)
        
            # Remove selected synapses from the synapse_sizes array before applying Kesten process
            synapse_sizes = deleteat!(synapse_sizes, sort(delete_indices))
        end

        # Apply Kesten process to mature synapses in N_M
        synapse_sizes = kesten_update_new007(synapse_sizes,ε, η, σ_ε, σ_η)

        push!(synapse_sizes_history, synapse_sizes)

    end

    solution = solve(probb);

    return solution, synapse_sizes, synapse_sizes_history, synapses, Ihist, Mhist
end

total_time = 140
total_pool_size = 1000
ε = 0.985
η = 1-ε
σ_ε, σ_η = .1, .1
kesten_timestep = 0.001
i=A3
m = 0.05
parameters = (A1, lambda1, A2, lambda2, m, A3, lambda3)

sol, synapse_sizes_var, synapse_sizes_history_var, synapses_var, ih, mh = run_simulation_diffeq_var007(total_time, total_pool_size, parameters, ε, η, σ_ε, σ_η, kesten_timestep);

time_array_var = sol.t
immature_population_var = sol[1, :]
mature_population_var = sol[2, :]
poold = sol[3,:]
maximum(immature_population_var+mature_population_var)
id_of_bump = argmax(immature_population_var+mature_population_var)

xtickss = collect(0:20:total_time)
push!(xtickss, trunc(Int,time_array_var[id_of_bump]))
xtickss=sort(xtickss)

using Plots.PlotMeasures
var_plot = plot(0:kesten_timestep:total_time, mean(ihs), lw=2, fillalpha=0.2,  ribbon=std(ihs)/num_trials, label="Immature (Random walks)", color=:pink)
plot!(0:kesten_timestep:total_time, mean(mhs), lw=2, ribbon=std(mhs)/num_trials, color=:skyblue1, label="Mature (Random walks)")
plot!(0:kesten_timestep:total_time, mean(ihs)+mean(mhs), color=:palegreen, lw=2,ribbon=std(ihs+mhs)/num_trials, label="Combined (Random walks)")
plot!(time_array_var, immature_population_var, lw=3, label = "Immature (Diff Eq)", legend=:bottomright, color=:red)
plot!(time_array_var, mature_population_var, color=:blue, label = "Mature (Diff Eq)", lw=3, xlabel="Time (Days)",ylabel="Population size")
plot!(time_array_var, immature_population_var+mature_population_var, lcolor=:green, w=3,label="Combined (Diff Eq)", title="Populations")
vline!([time_array_var[id_of_bump]],label="Where bump happens")
vspan!([16,26],fillalpha=0.1,label="Developmental period")
vspan!([70,88], fillalpha=0.1,label="Adulthood period", xticks=xtickss,legend=:outerright,size=(1000,500),bottommargin=5mm, leftmargin=5mm)
# hline!([85.222],xlim=(35,50),ylim=(80,90))

# savefig(var_plot, "C://Users/B00955735/OneDrive - Ulster University/Desktop/varplot.png")

p2 = plot(time_array_var, mature_population_var, lw=5, c=:green, label="Mature synapses", xlabel="Postnatal Day",ylabel="Number")
plot!(time_array_var, immature_population_var, lw=5, c=:magenta, label="Immature synapses")
plot!(time_array_var, immature_population_var.+mature_population_var, lw=5, c=:grey, linealpha=0.7,label="Total synapses")
plot!(title="Population Dynamics (Differential Equations Model)", lw=5, c=:black, label="Total synapses",legend=:bottomright)
plot!(grid=false,ylim=(0,total_pool_size),legendfontsize=12)


plot(ih)
plot!(mh)
plot!(ih+mh)

# savefig(p1,"C://Users/B00955735/OneDrive - Ulster University/Desktop/populations_randwalks.png")
# savefig(p2,"C://Users/B00955735/OneDrive - Ulster University/Desktop/populations_diffeqs.png")
# savefig(p1,"C://Users/B00955735/OneDrive - Ulster University/Desktop/populations_randwalks.svg")
# savefig(p2,"C://Users/B00955735/OneDrive - Ulster University/Desktop/populations_diffeqs.svg")

# Plot elimination and creation rates
plot(0:kesten_timestep:total_time, creat,lw=3,label="Creation rate")
plot!(0:kesten_timestep:total_time, elim, lw=3, label="Elimination rate",ylim=(0,2.5))


using KernelDensity

bit_to_add = 20
synapse_sizes_history_var
synapse_sizes_history_var
timepoint1 = trunc(Int,(15+bit_to_add)/kesten_timestep)
timepoint2 = trunc(Int,(35+bit_to_add)/kesten_timestep)
timepoint3 = trunc(Int,(55+bit_to_add)/kesten_timestep)
timepoint4 = trunc(Int,(120+bit_to_add)/kesten_timestep)
v1 = synapse_sizes_history_var[timepoint1]
v2 = synapse_sizes_history_var[timepoint2]
v3 = synapse_sizes_history_var[timepoint3]
v4 = synapse_sizes_history_var[timepoint4]




# Compute kernel density estimate
bw=0.5
kde_result1 = kde(v1,bandwidth=bw)
kde_result2 = kde(v2,bandwidth=bw)
kde_result3 = kde(v3,bandwidth=bw)
kde_result4 = kde(v4,bandwidth=bw)

# Scale by total synaptic strength to get "unnormalised" curves
scaled_density1 = kde_result1.density .* sum(v1)
scaled_density2 = kde_result2.density .* sum(v2)
scaled_density3 = kde_result3.density .* sum(v3)
scaled_density4 = kde_result4.density .* sum(v4)

color1 = RGB(176/255, 178/255, 240/255)
color2 = RGB(78/255,84/255,182/255)
color3 = RGBA(221/255, 221/255, 221/255, 255/255)
color4 = RGB(4/255, 4/255, 4/255)

# Plot the smooth density curve
plot(kde_result1.x, kde_result1.density, linewidth=4, label="P15", color=color1)
plot!(kde_result2.x, kde_result2.density, linewidth=4, label="P35", color=color2)
plot!(kde_result3.x, kde_result3.density, linewidth=4, label="P55", color=color3)
plot!(kde_result4.x, kde_result4.density, linewidth=4, label="P120", color=color4, xlabel="Synaptic weight", grid=false)

plot(kde_result1.x, scaled_density1, linewidth=4, label="P15", color=color1)
plot!(kde_result2.x, scaled_density2, linewidth=4, label="P35", color=color2)
plot!(kde_result3.x, scaled_density3, linewidth=4, label="P55", color=color3)
plot!(kde_result4.x, scaled_density4, linewidth=4, label="P120", color=color4,
      xlabel="Synaptic weight", ylabel="Total synaptic strength", grid=false, xlim=(-0.2,4))

bins=0:.1:5
histogram(v1, fillalpha=0.5, bins=bins)
histogram!(v2, fillalpha=0.5, bins=bins)
histogram!(v3, fillalpha=0.5, bins=bins)
histogram!(v4, fillalpha=0.5, xlabel="Synaptic weight", ylabel="Count", bins=bins)

combined_synapse_sizes = []
trials_synapse_sizes = []

tttrials = 10
for i in 1:tttrials
    sol, synapse_sizes_var, synapse_sizes_history_var, synapses_var, ih, mh = run_simulation_diffeq_var007(total_time, total_pool_size, parameters, ε, η, σ_ε, σ_η, kesten_timestep);
    push!(combined_synapse_sizes, sum.(synapse_sizes_history_var))
    push!(trials_synapse_sizes, synapse_sizes_history_var)
end

xticks1 = collect(0:20:total_time)
p = plot(0:kesten_timestep:total_time, mean(combined_synapse_sizes), xticks=xticks1, ribbon=std(combined_synapse_sizes)/sqrt(tttrials), grid=false, lw=3, title="Total synaptic weight over time", xlabel="Days", ylabel="Total synaptic weight", legend=false)

# savefig(p, "C://Users/B00955735/OneDrive - Ulster University/Desktop/total_syn_weight.png")


########################
########################
# working out information stored in distribution of synaptic weights


using SpecialFunctions
function bartol_bits(sizes; CV=0.083, overlap=0.69, min_clip=1e-3)
    sizes = filter(>(0.01), sizes)
    if isempty(sizes)
        return NaN
    end
    
    smin = max(minimum(sizes), min_clip)
    smax = maximum(sizes)
    range_factor = smax / smin

    # z corresponding to 69% discrimination
    z = sqrt(2) * erfinv(overlap)
    spacing_factor = 1 + 2 * CV * z

    # number of distinguishable states
    N = log(range_factor) / log(spacing_factor)
    return log2(N)  # bits per synapse
end


function odonnell_bits(w, c)
    w = filter(>(0.01), w)
    log_w = log.(w)
    sigma_w_sq = var(log_w)

    sigma_eta_sq = log(1+c^2)
    sigma_x_sq = sigma_w_sq - sigma_eta_sq

    I_w = (1/(2*log(2)))*log2(1+(sigma_x_sq)/(sigma_eta_sq))

    return I_w
end

bits_timecourse = [[bartol_bits(syn_hist) for syn_hist in trials_synapse_sizes[i]] for i in 1:tttrials]
bits_timecourse2 = [[odonnell_bits(syn_hist, 0.083) for syn_hist in trials_synapse_sizes[i]] for i in 1:tttrials]

tt = collect(0:1:total_time)
plot(0:kesten_timestep:total_time, mean(bits_timecourse2), ribbon = std(bits_timecourse2)/sqrt(tttrials), lw=3,
     xlabel="Time", ylabel="Bits per synapse", 
     label="Information capacity (O'Donnell)", color=:purple,size=(700,500))
plot!(0:kesten_timestep:total_time, mean(bits_timecourse), lw=3, ribbon=std(bits_timecourse)/sqrt(tttrials), label="Information capacity (Bartol)", color=:darkgreen)


##############



using Distributions

normal_data1 = rand(Normal(10,1), 1000)
normal_data2 = rand(Normal(11,1), 1000)
normal_data3 = rand(Normal(12,1), 1000)
uniform_data1 = rand(Uniform(10,20), 1000)
uniform_data2 = rand(Uniform(10,30), 1000)
uniform_data3 = rand(Uniform(10,40), 1000)
lognormal_data1 = rand(LogNormal(2,0.3), 1000)
lognormal_data2 = rand(LogNormal(2.5,0.3), 1000)
lognormal_data3 = rand(LogNormal(3,0.3), 1000)

normal_data = [normal_data1, normal_data2, normal_data3]
uniform_data = [uniform_data1, uniform_data2, uniform_data3]
lognormal_data = [lognormal_data1, lognormal_data2, lognormal_data3]

normal_info = [bartol_bits(data) for data in normal_data]
normal_info2 = [odonnell_bits(data, 0.1) for data in normal_data]
uniform_info = [bartol_bits(data) for data in uniform_data]
uniform_info2 = [odonnell_bits(data, 0.1) for data in uniform_data]
lognormal_info = [bartol_bits(data) for data in lognormal_data]
lognormal_info2 = [odonnell_bits(data, 0.1) for data in lognormal_data]












### paul's distribution data


using Trapz  # for trapezoidal integration

# Data as dictionaries of vectors
data = Dict(
    "P15" => ([0.6732111463850229, 1.3043471980192245, 2.019635044217049, 2.776997985165697,
               3.997194993683278, 5.4698457808297505, 7.8681630980521104, 11.192146945350368,
               14.263673803895774, 16.998596694310645, 20.028051668229192],
              [15.690377445561595, 266.7364165745415, 1114.0167986348492, 1060.669491378325,
               571.1297629600457, 326.3598987509062, 131.79907477625085, 50.20915994256428,
               31.380754891122745, 15.690377445561595, 12.55234983968099]),

    "P35" => ([0.7994396407614467, 1.5147258818972922, 2.019635044217049, 2.650771095851251,
               3.5764376259271438, 4.417952361439413, 5.427770686078928, 6.732119489160131,
               8.07854178193018, 9.046285011818872, 11.82328299698457, 14.726509476526694,
               18.176718608077405, 21.54277755012648, 24.82468630267392, 27.81206297171768],
              [28.24272728524214, 269.87456388850165, 790.7950471979033, 743.7239148612194,
               527.1966582292422, 423.6401910301539, 304.3931069693449, 222.803431551818,
               197.69885158053575, 160.041802061493, 100.41831988512811, 65.89953738812565,
               40.79507712492357, 34.5187824970029, 25.10469967936198, 21.966432657321914]),

    "P55" => ([1.0098183246395143, 1.3884989925828477, 2.7349228904148744, 3.450210736612699,
               4.586255950566658, 6.016830037900329, 7.237027046417912, 8.499299149686312,
               9.382888979949403, 10.561010893716167, 13.169705289754615, 14.389902298272196,
               15.988781579795091, 17.166900283437894, 18.80785465971161, 21.16409527712117,
               23.772792883283582, 25.960729911565895, 29.957924905249172],
              [18.828405051441756, 56.48545457048472, 461.29688142495723, 524.0583912072021,
               451.8829183153963, 367.1546167515897, 285.564941334063, 244.76986420913968,
               219.66528423785786, 172.59415190117375, 119.24696435272955, 100.41831988512811,
               91.00423706748762, 72.17583201604587, 43.93310473080329, 34.5187824970029,
               25.10469967936198, 18.828405051441756, 18.828405051441756]),

    "P120" => ([0.8415147355122694, 1.3884989925828477, 2.692847795664053, 3.3660589420490776,
                4.628331045317479, 5.427770686078928, 6.185133627027572, 7.237027046417912,
                8.120616876681003, 10.350632209838098, 12.74894792199848, 15.063113444657219,
                17.798036335072098, 20.911641498492283, 24.109396851414107, 27.18092692008348,
                30.0],
               [15.690377445561595, 59.62348217636466, 643.3054752680116, 646.4436225819717,
                520.9206030174819, 411.08796089855286, 338.91212888250726, 266.7364165745415,
                235.35566168341902, 144.35142461593185, 87.86620946160703, 65.89953738812565,
                37.656810102883064, 25.10469967936198, 28.24272728524214, 18.828405051441756,
                9.41432223380083]))



# Desired order
plot_order = ["P15","P35","P55","P120"]
colors = [RGB(0.25,0.85,1.0), RGB(0.4,0.7,1.0), RGB(0.2,0.5,1.0), RGB(0.0,0.2,0.8)]

p_auc = plot()
for (i, label) in enumerate(plot_order)
    x, y = data[label]  # get the correct curve
    plot!(x, y, label=label, color=colors[i], lw=5)
end
xlabel!("Size * Intensity Mean")
ylabel!("Frequency")
title!("Distributions", size=(600,400), grid=false)

# Compute areas
areas = Dict()
auc_data = []
for (label, (x, y)) in data
    areas[label] = trapz(x, y)
    push!(auc_data, areas[label])
end


for label in plot_order
    x, y = data[label]
    areas[label] = trapz(x, y)
    push!(auc_data, trapz(x, y))
end

auc_data

xpos = [15, 35, 55, 120]
p2 = plot(xpos,auc_data, c=:grey, lw=2, ylim=(0,4500), xticks=([15,35,55,120], ["P15", "P35", "P55", "P120"]), title="Areas under curves", xlabel="Days", ylabel="Area under curve")
scatter!([xpos[1]],[auc_data[1]], m=:circle, ms=8, c=colors[1])
scatter!([xpos[2]],[auc_data[2]], m=:circle, ms=8, c=colors[2])
scatter!([xpos[3]],[auc_data[3]], m=:circle, ms=8, c=colors[3])
scatter!([xpos[4]],[auc_data[4]], m=:circle, ms=8, c=colors[4],legend=false)


p3 = plot(p_auc,p2, layout=(1,2), size=(1000,400), bottommargin=5mm, leftmargin=5mm)

# savefig(p3, "C://Users/B00955735/OneDrive - Ulster University/Desktop/paul_distributions_and_areas.png")


function sample_from_histogram(x, y)
    total = sum(y)
    probs = y ./ total
    samples = []
    for (xi, pi) in zip(x, probs)
        n = round(Int, pi * 1000)  # 1000 pseudo-samples per distribution
        append!(samples, fill(xi, n))
    end
    return Float64.(samples)
end

# Apply your synapse information function
function synapse_information_from_histogram(x, y, cv)
    samples = sample_from_histogram(x, y)
    return odonnell_bits(samples,cv), bartol_bits(samples; CV=cv)
end

# Compute bits for each distribution
bits_per_dist_odonnell = Dict()
bits_per_dist_bartol = Dict()
for label in ["P15", "P35", "P55", "P120"]
    x, y = data[label]
    bits_per_dist_odonnell[label], bits_per_dist_bartol[label] = synapse_information_from_histogram(x, y,0.083)
end


xpos = [15, 35, 55, 120]
dists_bits_ordered_odonnell = [bits_per_dist_odonnell[label] for label in plot_order]
dists_bits_ordered_bartol = [bits_per_dist_bartol[label] for label in plot_order]
p = plot(xpos,dists_bits_ordered_odonnell, label="O'Donnell Info", lw=5, xticks=([15,35,55,120], ["P15", "P35", "P55", "P120"]), 
xlabel="Distribution", ylabel="Bits per synapse", c=:black, title="Information from distributions",
marker=:circle,         # use :circle, :diamond, etc.
    markercolor=:white,     # white fill
    markerstrokecolor=:black,markersize=8,)
plot!(xpos,dists_bits_ordered_bartol, c=:grey, label="Bartol Info", lw=5, ylim=(0,8),marker=:diamond,         # use :circle, :diamond, etc.
    markercolor=:white,     # white fill
    markerstrokecolor=:black,markersize=8,)
scatter!([xpos[1]],[dists_bits_ordered_bartol[1]], m=:diamond, ms=8, c=colors[1],label=false)
scatter!([xpos[2]],[dists_bits_ordered_bartol[2]], m=:diamond, ms=8, c=colors[2],label=false)
scatter!([xpos[3]],[dists_bits_ordered_bartol[3]], m=:diamond, ms=8, c=colors[3],label=false)
scatter!([xpos[4]],[dists_bits_ordered_bartol[4]], m=:diamond, ms=8, c=colors[4],label=false)
scatter!([xpos[1]],[dists_bits_ordered_odonnell[1]], m=:circle, ms=8, c=colors[1],label=false)
scatter!([xpos[2]],[dists_bits_ordered_odonnell[2]], m=:circle, ms=8, c=colors[2],label=false)
scatter!([xpos[3]],[dists_bits_ordered_odonnell[3]], m=:circle, ms=8, c=colors[3],label=false)
scatter!([xpos[4]],[dists_bits_ordered_odonnell[4]], m=:circle, ms=8, c=colors[4],label=false,legend=:topright)
plot!(xlabel="Days")

# savefig(p, "C://Users/B00955735/OneDrive - Ulster University/Desktop/information_from_distributions.png")












########################
########################



# using Optimization
# rosenbrock(x, p) = (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
# x0 = zeros(2)
# p = [1.0, 100.0]

# prob = OptimizationProblem(rosenbrock, x0, p)

# using OptimizationOptimJL
# sol = solve(prob, NelderMead(), g_tol=1e-9)



#############################
###########################
########################
# Variable Kesten Process Values
########################
############################
#############################

q1,q2,q3,q4 = .985, 0.015, .05, .05
w1, w2, w3, w4 = 1/100,1/100,1/100,1/100
εs_func(t) = q1 #*exp(-t*w1)
ηs_func(t) = q2 #*exp(-t*w2)
σ_εs_func(t) = q3*exp(-t*w3)
σ_ηs_func(t) = q4*exp(-t*w4)


plot(ts,εs_func.(ts))

εs = εs_func.(0:kesten_timestep:total_time)
ηs = ηs_func.(0:kesten_timestep:total_time)
σ_εs = σ_εs_func.(0:kesten_timestep:total_time)
σ_ηs = σ_ηs_func.(0:kesten_timestep:total_time)

rates_var = creat, m, elim, A, lambda


state_recs_var_multiple = []
ihs = []
mhs = []
synapse_sizes_multiple = []
syn_size_heatmaps_trials = []

num_trials = 5

for i in 1:num_trials
    ih_var, mh_var, state_record_var, syn_sizes_var, syn_heatmap = track_times_variable_rates_and_kesten_007(total_time, total_pool_size, rates_var, εs, ηs, σ_εs, σ_ηs, kesten_timestep);
    push!(state_recs_var_multiple, state_record_var)
    push!(ihs, ih_var)
    push!(mhs, mh_var)
    push!(synapse_sizes_multiple, syn_sizes_var)
    push!(syn_size_heatmaps_trials,syn_heatmap)
end

syns_ht = heatmap(syn_size_heatmaps_trials[1],clims=(0,3),xticks=(0:2000:12000,0:20:120),ylabel="Synapse ID", xlabel="Days",colorbar_title="Synapse size")
# savefig(syns_ht, "C://Users/B00955735/OneDrive - Ulster University/Desktop/syns_ht.png")

hist_matrices = []
for i in 1:num_trials
    # Define histogram bins (consistent for all columns)
    bin_edges = 0:0.1:2.0  # Adjust bin size as needed
    bin_centers = (bin_edges[1:end-1] .+ bin_edges[2:end]) / 2  # Centers for plotting

    # Compute histograms for each time point (column)
    num_bins = length(bin_edges) - 1
    num_time_points = size(syn_size_heatmaps_trials[i], 2)
    hist_matrix = zeros(num_bins, num_time_points)

    for t in 1:num_time_points
        hist = fit(Histogram, syn_size_heatmaps_trials[i][:, t], bin_edges)
        hist_matrix[:, t] = hist.weights  # Store bin counts for this time point
    end
    push!(hist_matrices, hist_matrix)
end

bin_edges = 0:0.1:2.0  # Adjust bin size as needed
bin_centers = (bin_edges[1:end-1] .+ bin_edges[2:end]) / 2  # Centers for plotting
num_bins = length(bin_edges) - 1
num_time_points = size(syn_size_heatmaps_trials[1], 2)


# Plot the heatmap
hists = heatmap(
    1:num_time_points,               # Time points on x-axis
    bin_centers,                     # Bin centers on y-axis
    mean(hist_matrices),                     # Heatmap data
    xlabel = "Time",
    ylabel = "Synapse Size",
    colorbar_title = "Count",
    title = "Histograms of synapse size across time",clims=(1,2)
)

# savefig(hists, "C://Users/B00955735/OneDrive - Ulster University/Desktop/hists.png")



tt = 0.01:kesten_timestep:total_time
indi = plot(tt,syn_size_heatmaps_trials[1][60,:], label="Synapse 1")
plot!(tt,syn_size_heatmaps_trials[1][4,:], label="Synapse 2")
plot!(tt, syn_size_heatmaps_trials[1][1,:], label="Synapse 3",color="green", xticks=collect(0:20:120),xlabel="Days",ylabel="Synapse size")
