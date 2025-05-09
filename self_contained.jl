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
            sizes[i] = 0.0
        end
        # sizes[i] = max(new_size, 0.0)  # Ensure size is non-negative
        if sizes[i] < 0.0
            sizes[i]=0.0
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
            prob = A * exp(-size / lambda) * kesten_timestep #i*kesten_time_step
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
A1,lambda1,A2,lambda2,m,A3,lambda3 = 0.9, 30, 2, 5, 0.05, 0.05, 2.
ε = .985
η = 1-ε
σ_ε, σ_η = .1, .1


creation_func(t) = A1 * exp(-t / lambda1) + b1
elimination_func(t) = A2 * exp(-t / lambda2) + b2

elim = elimination_func.(0:kesten_timestep:total_time)
creat = creation_func.(0:kesten_timestep:total_time)

ec_plot = plot(0:kesten_timestep:total_time, creat,lw=3,label="Creation rate")
plot!(0:kesten_timestep:total_time, elim, lw=3, label="Elimination rate", ylabel="Rate", xlabel="Days")

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
kde_result1 = kde(v1,bandwidth=0.5)
kde_result2 = kde(v2,bandwidth=0.5)
kde_result3 = kde(v3,bandwidth=0.5)
kde_result4 = kde(v4,bandwidth=0.5)

color1 = RGB(176/255, 178/255, 240/255)
color2 = RGB(78/255,84/255,182/255)
color3 = RGBA(221/255, 221/255, 221/255, 255/255)
color4 = RGB(4/255, 4/255, 4/255)

distsplot = plot(kde_result1.x, kde_result1.density, linewidth=4, label="P15", color=color1)
plot!(kde_result2.x, kde_result2.density, linewidth=4, label="P35", color=color2)
plot!(kde_result3.x, kde_result3.density, linewidth=4, label="P55", color=color3)
plot!(kde_result4.x, kde_result4.density, linewidth=4, label="P120", color=color4,
        xlabel="Synaptic weight (a.u.)", grid=false, xlim=(-0.2,3), title="Distributions of synaptic weight over time",
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

p4 = plot(timepoints, meanweights,xlabel="Postnatal Day",grid=false,ylim=(0,1),c=:firebrick2,title="Arithmetic mean", xticks=timepoints,ylabel="Synaptic weight (a.u.)",label=false,lw=5)
p5 = plot(timepoints,stdweights,xlabel="Postnatal Day",grid=false,c=:coral2,ylim=(0,1),title="Arithmetic standard deviation", xticks=timepoints,ylabel="Synaptic weight (a.u.)",label=false,lw=5)

p6 = plot(timepoints,skewnessweights,xlabel="Postnatal Day",c=:black,grid=false,ylim=(0,3),title="Skewness", xticks=timepoints,ylabel="Synaptic weight (a.u.)",label=false,lw=5)

# pp = plot(p1,p2,layout=(1,2),size=(1400,700), leftmargin=10mm,bottommargin=10mm,
# xlabelfontsize=18,ylabelfontsize=18,titlefontsize=20, xtickfontsize=18, ytickfontsize=18)

# savefig(p6,"C://Users/B00955735/OneDrive - Ulster University/Desktop/skew.png")










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
    title = "Histograms of synapse size across time", clims=(0,2)
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
    c_t, m, e_t, i, λ, synapse_sizes = p 
    N_I, N_M, N_P = u
    A = i

    # Compute the rate of dematuration using the exponential probability distribution
    # Sum over all mature synapses' probabilities of transitioning to immature state
    if !isempty(synapse_sizes)
        dematuration_rate = A * sum(exp(- size / λ) for size in synapse_sizes) / length(synapse_sizes)
    else
        dematuration_rate = 0
    end
    # dematuration_rate = i
    # Apply time-dependent e(t) and c(t)
    e_t = elimination_func(t)
    c_t = creation_func(t)

    du[1] = c_t * N_P - (m + e_t) * N_I + (dematuration_rate) * N_M  # dN_I/dt
    du[2] = m * N_I - (dematuration_rate) * N_M  # dN_M/dt
    du[3] = - du[1] - du[2]  # dN_P/dt
end

function run_simulation_diffeq_var007(total_time, total_pool_size, paramss, ε, η, σ_ε, σ_η, kesten_time_step)
    pool = fill(1, total_pool_size);  # Initialize resource pool with synapses
    synapses = Int[]  # Array to hold states of synapses (0s and 1s)
    synapse_sizes = Float64[]  # Sizes of mature synapses
    synapse_sizes_history = []
    m, i, λ = paramss
    # Initial conditions
    u0 = [0.0, 0.0, total_pool_size];
    tspan = (0.0, total_time);
    # p = (m, i, λ, synapse_sizes)
    Ihist = []
    Mhist = []

    # Define ODE problem
    # prob = ODEProblem(synapse_dynamics_var!, u0, tspan, p);

    current_time = 0.0;

    prams = creation_func(current_time),m,elimination_func(current_time),i, λ, synapse_sizes
    probb = ODEProblem(synapse_dynamics_var007!, u0, tspan, prams);

    while current_time < total_time
        prams = creation_func(current_time),m,elimination_func(current_time),i, λ, synapse_sizes
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

        elseif N_M < length(synapse_sizes)
            num_delete_matures = length(synapse_sizes) - round(Int, N_M) #find how many need to be deleted
            sizes = synapse_sizes
            weights = exp.(-sizes ./ λ)
            weights ./= sum(weights)  # Normalise to sum to 1
        
            # Sample indices to delete based on weights, without replacement
            delete_indices = sample(1:length(sizes), Weights(weights), num_delete_matures; replace=false)
        
            # Remove selected synapses from the synapse_sizes array before applying Kesten process
            synapse_sizes = deleteat!(synapse_sizes, sort(delete_indices; rev=true))
        end

        # Apply Kesten process to mature synapses in N_M
        synapse_sizes = kesten_update_new007(synapse_sizes,ε, η, σ_ε, σ_η)

        push!(synapse_sizes_history, synapse_sizes)

    end

    solution = solve(probb);

    return solution, synapse_sizes, synapse_sizes_history, synapses, Ihist, Mhist
end

total_pool_size = 100
# ε, η = .985, 1-ε
# σ_ε, σ_η = .1, .1
kesten_timestep = .2
i=A
paramss = (m, i, lambda3)

sol, synapse_sizes_var, synapse_sizes_history_var, synapses_var, ih, mh = run_simulation_diffeq_var007(total_time, total_pool_size, paramss, ε, η, σ_ε, σ_η, kesten_timestep);

time_array_var = sol.t
immature_population_var = sol[1, :]
mature_population_var = sol[2, :]
poold = sol[3,:]
maximum(immature_population_var+mature_population_var)
id_of_bump = argmax(immature_population_var+mature_population_var)

xtickss = collect(0:20:120)
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

# savefig(p1,"C://Users/B00955735/OneDrive - Ulster University/Desktop/populations_randwalks.png")
# savefig(p2,"C://Users/B00955735/OneDrive - Ulster University/Desktop/populations_diffeqs.png")
# savefig(p1,"C://Users/B00955735/OneDrive - Ulster University/Desktop/populations_randwalks.svg")
# savefig(p2,"C://Users/B00955735/OneDrive - Ulster University/Desktop/populations_diffeqs.svg")

# Plot elimination and creation rates
plot(0:kesten_timestep:total_time, creat,lw=3,label="Creation rate")
plot!(0:kesten_timestep:total_time, elim, lw=3, label="Elimination rate",ylim=(0,2.5))


using KernelDensity
synapse_sizes_history_var
v1 = synapse_sizes_history_var[15]
v2 = synapse_sizes_history_var[35]
v3 = synapse_sizes_history_var[55]
v4 = synapse_sizes_history_var[120]

# Compute kernel density estimate
kde_result1 = kde(v1,bandwidth=0.3)
kde_result2 = kde(v2,bandwidth=0.3)
kde_result3 = kde(v3,bandwidth=0.3)
kde_result4 = kde(v4,bandwidth=0.3)

color1 = RGB(176/255, 178/255, 240/255)
color2 = RGB(78/255,84/255,182/255)
color3 = RGBA(221/255, 221/255, 221/255, 255/255)
color4 = RGB(4/255, 4/255, 4/255)

# Plot the smooth density curve
plot(kde_result1.x, kde_result1.density, linewidth=4, label="P15", color=color1)
plot!(kde_result2.x, kde_result2.density, linewidth=4, label="P35", color=color2)
plot!(kde_result3.x, kde_result3.density, linewidth=4, label="P55", color=color3)
plot!(kde_result4.x, kde_result4.density, linewidth=4, label="P120", color=color4, xlabel="Synaptic weight", grid=false, xlim=(-0.2,4))



























########################
########################



using Optimization
rosenbrock(x, p) = (p[1] - x[1])^2 + p[2] * (x[2] - x[1]^2)^2
x0 = zeros(2)
p = [1.0, 100.0]

prob = OptimizationProblem(rosenbrock, x0, p)

using OptimizationOptimJL
sol = solve(prob, NelderMead(), g_tol=1e-9)



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
