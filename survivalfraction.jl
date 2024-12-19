# Function to compute the survival fraction from state_records
function compute_survival_fraction(state_records)
    total_time_steps = size(state_records, 2)

    # Identify the initial population
    initial_existing_synapses = findall(x -> x != 0, state_records[:, 1])

    # Initialize array to store the survival fraction over time
    survival_fraction = zeros(total_time_steps)

    # Loop over each time step and compute the fraction of surviving synapses
    for t in 1:total_time_steps
        # Find how many of the initial synapses are still in the existing population at time t
        surviving_synapses = count(x -> state_records[x, t] != 0, initial_existing_synapses)
        
        # Compute survival fraction as the ratio of surviving synapses to the initial population size
        survival_fraction[t] = surviving_synapses / length(initial_existing_synapses)
    end

    return survival_fraction
end

# function compute_survival_fraction(state_records)
#     total_time_steps = size(state_records, 2)

#     # Identify the initial population
#     initial_existing_synapses = findall(x -> x != 0, state_records[:, 1])

#     # Initialize array to store the survival fraction over time
#     survival_fraction = zeros(total_time_steps)

#     # Create an array to track if a synapse has ever left the population
#     has_left_population = falses(length(initial_existing_synapses))

#     # Loop over each time step and compute the fraction of surviving synapses
#     for t in 1:total_time_steps
#         surviving_synapses = 0

#         # Check each synapse in the initial population
#         for (index, synapse) in enumerate(initial_existing_synapses)
#             # If the synapse hasn't left and is still present, count it as surviving
#             if !has_left_population[index] && state_records[synapse, t] != 0
#                 surviving_synapses += 1
#             elseif state_records[synapse, t] == 0
#                 # Mark the synapse as having left the population permanently
#                 has_left_population[index] = true
#             end
#         end
        
#         # Compute survival fraction as the ratio of surviving synapses to the initial population size
#         survival_fraction[t] = surviving_synapses / length(initial_existing_synapses)
#     end

#     return survival_fraction
# end



a1,k1,b1,a2,k2,b2,m,A,lambda,ε,η,σ_ε, σ_η = sol

total_pool_size = 100
total_time = 1000
kesten_timestep = 0.01
creation_func(t) = a1 * exp(-t * k1) + b1
elimination_func(t) = a2 * exp(-t * k2) + b2

# parameters with pool = 100: 0.2123729118369409, 0.09117010118164337, 1.1700516097047382, 0.04654061119446891, 0.11116109898864276, 0.19446676360218018
elim = elimination_func.(0:kesten_timestep:total_time)
creat = creation_func.(0:kesten_timestep:total_time)

plot(elim)
plot!(creat)

using Distributions
rates_var = creat, m, elim, i, A, lambda;
ih_var, mh_var, state_records_var, syn_sizes_var = track_times_variable_rates(total_time, total_pool_size, rates_var, ε, η, σ_ε, σ_η, kesten_timestep);


plot(0:kesten_timestep:total_time, ih_var)
plot!(0:kesten_timestep:total_time, mh_var)
plot!(0:kesten_timestep:total_time, ih_var+mh_var)
vline!([16,26], label="Developmental period")
vline!([70], label = "Adulthood")

developmental_period_16 = round(Int, (16/total_time)*100/kesten_timestep)
developmental_period_26 = round(Int, (26/total_time)*100/kesten_timestep)

adult_period = round(Int, (70/total_time)*100/kesten_timestep)

# Compute survival fraction
developmental_survival_fraction = compute_survival_fraction(state_records_var[:,developmental_period_16:developmental_period_26])
adulthood_survival_fraction = compute_survival_fraction(state_records_var[:,adult_period:round(Int, 18/kesten_timestep)+adult_period])

plot(developmental_survival_fraction)
plot(adulthood_survival_fraction)
collect(16:kesten_timestep:26)

# Plot survival fraction over time
developmental_survival_plot = plot(16:kesten_timestep:26, developmental_survival_fraction, xlabel="Postnatal Day", ylabel="Survival Fraction", xticks=16:1:26,
    title="Synapse Survival Fraction (Early Development)", lw=2, legend=false, ylim=(0,1.05))
     
adult_survival_plot = plot(0:kesten_timestep:18, adulthood_survival_fraction, xlabel="Days", ylabel="Survival Fraction",
title="Synapse Survival Fraction (Adulthood)", lw=2, legend=false, ylim=(0,1.05), xticks=0:1:18)
 
survival_fraction_plot = plot(developmental_survival_plot, adult_survival_plot, layout=(2,1))


# Multiple trials

num_trials = 10

state_recs_var_multiple = []

for i in 1:num_trials
    ih_var, mh_var, state_record_var, syn_sizes_var = track_times_variable_rates(total_time, total_pool_size, rates_var, ε, η, σ_ε, σ_η, kesten_timestep);
    push!(state_recs_var_multiple, state_record_var)
end

develop_survival_multiple = []
adult_survival_multiple = []

for state_recs in state_recs_var_multiple
    developmental_survival_fraction1 = compute_survival_fraction(state_recs[:,developmental_period_16:developmental_period_26])
    adulthood_survival_fraction1 = compute_survival_fraction(state_recs[:,adult_period:round(Int, 18/kesten_timestep)+adult_period])
    push!(develop_survival_multiple, developmental_survival_fraction1)
    push!(adult_survival_multiple, adulthood_survival_fraction1)
end


developmental_survival_plot_trials = plot(16:kesten_timestep:26, mean(develop_survival_multiple), ribbon=std(develop_survival_multiple), xlabel="Postnatal Day", ylabel="Survival Fraction", xticks=16:1:26,
    title="Synapse Survival Fraction (Early Development)", lw=2, legend=false, ylim=(0,1.05))
     
adult_survival_plot_trials = plot(0:kesten_timestep:18, mean(adult_survival_multiple), ribbon=std(adult_survival_multiple), xlabel="Days", ylabel="Survival Fraction",
title="Synapse Survival Fraction (Adulthood)", lw=2, legend=false, ylim=(0,1.05), xticks=0:1:18)
 
survival_fraction_plot_trials = plot(developmental_survival_plot_trials, adult_survival_plot_trials, layout=(2,1))

dev_ids = collect(0:1:10)
dev_ids = [round(Int, id/kesten_timestep) for id in dev_ids]

adult_ids = [0,1,2,3,4,5,6,8,10,12,14,16,17,18]
adult_ids = [round(Int, id/kesten_timestep) for id in adult_ids]

adult_ids1 = collect(0:1:18)
adult_ids2 = collect(0:1:18)
adult_ids2 = [round(Int, id/kesten_timestep) for id in adult_ids2]
adult_ids3 = [0,1,2,3,4,5,6,8,10,12,14,16,17,18]




development_points_to_match_sim = [mean(develop_survival_multiple)[id+1] for id in dev_ids]
development_points_to_match_data = [1.0, 0.661896208, 0.52522361,0.468246877, 0.421466905, 0.397137735, 0.376028593, 0.364221812, 0.344543843, 0.348389962, 0.340339859]
    
development_survival_error = development_points_to_match_sim - development_points_to_match_data

adulthood_points_to_match_sim = [mean(adult_survival_multiple)[id+1] for id in adult_ids]
# adulthood_points_to_match_data = [1.0, 0.870199702, 0.82058372, 0.788018458, 0.775729644, 0.755248343, 0.688556071, 0.681643617]
adulthood_points_to_match_data = [1.0, 0.870199702, 0.82058372, 0.788018458, 0.775729644, 0.755248343, 0.7490909229625357, 0.7400000138716264, 0.7290909507057883, 0.7163636641068893, 0.7054545315829192, 0.694545468417081, 0.688556071, 0.681643617]


adulthood_survival_error = adulthood_points_to_match_sim - adulthood_points_to_match_data

dev_scatter = scatter(16:1:26, development_points_to_match_sim, ylim=(0,1.05), label="Model", xticks=16:1:26)
scatter!(16:1:26, development_points_to_match_data, label="Data",title="Survival Fraction (Early Development)", xlabel="Postnatal day")
# plot!(16:1:26, development_survival_error.^2, color="red", label="Error^2")

adult_scatter = scatter(adult_ids1, [mean(adult_survival_multiple)[id+1] for id in adult_ids2], ylim=(0,1.05), label="Model", xticks=0:1:18)
scatter!(adult_ids3, adulthood_points_to_match_data, label="Data",title="Survival Fraction (Adulthood)", xlabel="Days", legend=:bottomleft)
# plot!(adult_ids3, adulthood_survival_error.^2, color="red", label="Error^2")


pp = plot(dev_scatter, adult_scatter, layout=(2,1))

# savefig(pp, "C://Users/B00955735/OneDrive - Ulster University/Desktop/parameter_fit.png")


developmental_survival_plot_trials = plot(16:kesten_timestep:26, mean(develop_survival_multiple), ribbon=std(develop_survival_multiple), xlabel="Postnatal Day", ylabel="Survival Fraction", xticks=16:1:26,
    title="Synapse Survival Fraction (Early Development)", lw=2, legend=false, ylim=(0,1.05),label="Model")
scatter!(16:1:26, development_points_to_match_data, label="Data",title="Survival Fraction (Early Development)", xlabel="Postnatal day")

adult_survival_plot_trials = plot(0:kesten_timestep:18, mean(adult_survival_multiple), ribbon=std(adult_survival_multiple), xlabel="Days", ylabel="Survival Fraction",
title="Synapse Survival Fraction (Adulthood)", lw=2, legend=false, ylim=(0,1.05), xticks=0:1:18, label="Model")
scatter!(adult_ids3, adulthood_points_to_match_data, label="Data",title="Survival Fraction (Adulthood)", xlabel="Days", legend=:bottomleft)

ppp = plot(developmental_survival_plot_trials, adult_survival_plot_trials, layout=(2,1))
# savefig(ppp, "C://Users/B00955735/OneDrive - Ulster University/Desktop/parameter_fit.png")


