using Optim

function objective(pars)
    A3, lambda3, m = pars
    # re-run your code with these values of A3, lambda3
    rates_var = creat, m, elim, A3, lambda3

    state_recs_var_multiple = []
    ihs, mhs = [], []
    develop_survival_multiplee, adult_survival_multiplee = [], []

    num_trials =2

    for i in 1:num_trials
        ih_var, mh_var, state_record_var, syn_sizes_var, syn_heatmap, syn =
            track_times_variable_rates_007(total_time, total_pool_size,
                                           rates_var, ε, η, σ_ε, σ_η, kesten_timestep)
        push!(state_recs_var_multiple, state_record_var)
        push!(ihs, ih_var)
        push!(mhs, mh_var)

        developmental_period_16 = round(Int, (16/total_time)*size(state_record_var,2))
        developmental_period_26 = round(Int, (26/total_time)*size(state_record_var,2))
        adult_period  = round(Int, (100/total_time)*size(state_record_var,2))
        adult_period2 = round(Int, (118/total_time)*size(state_record_var,2))

        developmental_survival_fraction1 =
            new_compute_survival_fraction007(state_record_var[:, developmental_period_16:developmental_period_26])
        adulthood_survival_fraction1 =
            new_compute_survival_fraction007(state_record_var[:, adult_period:adult_period2])

        push!(develop_survival_multiplee, developmental_survival_fraction1)
        push!(adult_survival_multiplee, adulthood_survival_fraction1)
    end

    # errors
    development_points_to_match_sim = [mean(develop_survival_multiplee)[id+1] for id in dev_ids]
    adulthood_points_to_match_sim   = [mean(adult_survival_multiplee)[id+1] for id in adult_ids]

    development_survival_error = development_points_to_match_sim - development_points_to_match_data
    adulthood_survival_error   = adulthood_points_to_match_sim - adulthood_points_to_match_data

    return sum(development_survival_error.^2) + sum(adulthood_survival_error.^2)
end

# Run the optimisation
result = optimize(objective, [0.05, 2.0, 0.05])  # initial guess [A3, lambda3]
println(result)
println("Optimal A3 = ", Optim.minimizer(result)[1])
println("Optimal λ3 = ", Optim.minimizer(result)[2])
println("Optimal m  = ", Optim.minimizer(result)[3])
println("Minimum error = ", Optim.minimum(result))
