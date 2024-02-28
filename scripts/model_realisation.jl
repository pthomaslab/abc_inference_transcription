include("model.jl")

function run_sim(θ::Vector{Float64}, iv::Vector{Float64}, agevec::Vector{Float64}, cycle::Float64, pulsevec::Vector{Float64}, chasevec::Vector{Float64}, t0::Float64, vary_map::Vector{Any},
    scaling::Int64, n_steps::Int64, downsampling::Bool, betas::Vector{Float64}, age::Vector{Int64}, pulse_idx::Vector{Int64}, chase_idx::Vector{Int64}, age_dist::Matrix{Float64})
    s = Matrix{Float64}(undef,(length(agevec),5))
    s_pulse = Matrix{Float64}(undef,(length(agevec),2))
    s_chase = Matrix{Float64}(undef,(length(agevec),2))
    s_ratios = Vector{Float64}(undef,length(pulsevec))
    s_mean_corr = Vector{Float64}(undef,length(pulsevec))    #mean covariance over age distribution
    s_corr_mean = Vector{Float64}(undef,length(pulsevec))   #covariance of age-dependent means   
    ss_iv::Vector{Float64} = get_steady_state_iv(θ,vary_map,scaling,n_steps,iv,cycle)            
    for j in 1:6
        s = get_synthetic_data(θ,vary_map,scaling,n_steps,ss_iv,agevec,cycle,pulsevec[j],chasevec[j],t0,downsampling,betas[pulse_idx],age[pulse_idx])
        s_ratios[j] = sum(age_dist[:,j] .* s[:,2]) ./ (sum(age_dist[:,j] .* s[:,1]) + sum(age_dist[:,j] .* s[:,2]))
        stds::Float64 = sqrt(abs((sum(age_dist[:,j] .* s[:,3]) + weighted_cov(s[:,1],s[:,1],age_dist[:,j])) * (sum(age_dist[:,j] .* s[:,5]) + weighted_cov(s[:,2],s[:,2],age_dist[:,j]))))
        s_mean_corr[j] = sum(age_dist[:,j] .* s[:,4]) / stds
        s_corr_mean[j] = weighted_cov(s[:,1],s[:,2],age_dist[:,j]) / stds
        if j == 6
            ε::Vector{Float64} = [0.0 * (>(x,0.0)) + 0.0001 *(==(x,0.0)) for x in s[:,1] + s[:,2]]
            s_pulse[:,1] = s[:,1] + s[:,2]
            s_pulse[:,2] = (s[:,3] + 2*s[:,4] + s[:,5]) ./ (s[:,1] + s[:,2] + ε)
        end
    end
    for j in 7:length(pulsevec)
        s = get_synthetic_data(θ,vary_map,scaling,n_steps,ss_iv,agevec,cycle,pulsevec[j],chasevec[j],t0,downsampling,betas[chase_idx],age[chase_idx])
        s_ratios[j] = sum(age_dist[:,j] .* s[:,2]) ./ (sum(age_dist[:,j] .* s[:,1]) + sum(age_dist[:,j] .* s[:,2]))
        stds::Float64 = sqrt(abs((sum(age_dist[:,j] .* s[:,3]) + weighted_cov(s[:,1],s[:,1],age_dist[:,j])) * (sum(age_dist[:,j] .* s[:,5]) + weighted_cov(s[:,2],s[:,2],age_dist[:,j]))))
        s_mean_corr[j] = sum(age_dist[:,j] .* s[:,4]) / stds
        s_corr_mean[j] = weighted_cov(s[:,1],s[:,2],age_dist[:,j]) / stds
        if j == 7
            ε::Vector{Float64} = [0.0 * (>(x,0.0)) + 0.0001 *(==(x,0.0)) for x in s[:,1] + s[:,2]]
            s_chase[:,1] = s[:,1] + s[:,2]
            s_chase[:,2] = (s[:,3] + 2*s[:,4] + s[:,5]) ./ (s[:,1] + s[:,2] + ε)
        end
    end
    return s_pulse, s_chase, s_ratios, s_mean_corr, s_corr_mean
end

function weighted_cov(x::Vector{Float64},y::Vector{Float64}, w::Vector{Float64})
    return sum(w .* ((x .- sum(w .* x)) .* (y .- sum(w .* y))))
end

function full_model(paramvals::Vector{Float64}, vary_map::Vector{Any}, scaling::Int64, t_steps::Vector{Float64}, iv::Vector{Float64}, tmin::Float64,
   tmax::Float64, cycle::Float64, texp::Float64, pulse::Float64)
   tspan::Tuple{Float64, Float64} = (tmin, tmax)
   _moment_odes(dμ, μ, θ, t) = momentodes(dμ, μ, θ, t, cycle, texp, pulse, vary_map, t_steps, scaling)
   odeprob = ODEProblem(_moment_odes, iv, tspan, paramvals)
   sol = solve(odeprob, CVODE_BDF())
   return sol.u, sol.u[end], sol.t
end

function full_trajectories(paramvals::Vector{Float64}, vary_map::Vector{Any}, scaling::Int64, t_steps::Vector{Float64}, iv::Vector{Float64}, age::Float64,
   cycle::Float64, pulse::Float64, t0::Float64, texp::Float64)
    local endpoint::Vector{Float64}
    local nsols::Int64
    gene_covidx::Vector{Int64} = [5,6]
    meanidx::Vector{Int64} = [2,3]
    varidx::Vector{Int64} = [7,9]
    covidx::Vector{Int64} = [8] 
    if age > 0 && age < cycle
        nsols = floor((age-t0)/cycle) + 1
    else
        nsols = floor((age-t0)/cycle)
    end
    sol = Vector{Vector{Vector{Float64}}}(undef, nsols)
    timespan = Vector{Vector{Float64}}(undef, nsols)
    τ = t0
    k = 1
    for k in 1:nsols
        if τ < 0
            tf = τ + cycle
        else
            tf = age
        end
        if τ == t0
            sol[k], endpoint, timespan[k] = full_model(paramvals, vary_map, scaling, t_steps, iv, τ, tf, cycle, texp, pulse)
        else
            #set new initial conditions for next cell cycle
            gene_coviv = endpoint[gene_covidx]
            meaniv = endpoint[meanidx]
            variv = endpoint[varidx]
            coviv = endpoint[covidx]
            #binomial partitioning
            endpoint[gene_covidx] = gene_coviv/2
            endpoint[covidx] = coviv/4
            endpoint[varidx] = variv/4 + meaniv/4
            endpoint[meanidx] = meaniv/2
            sol[k], endpoint, timespan[k] = full_model(paramvals, vary_map, scaling, t_steps, endpoint, τ, tf, cycle, texp, pulse)
        end
        τ = tf
    end
    return vcat(sol...), vcat(timespan...)
end

function get_solutions(paramvals::Vector{Float64}, vary_map::Vector{Any}, n_steps::Int64, scaling::Int64, iv::Vector{Float64}, age::Float64, cycle::Float64, pulse::Float64, chase::Float64, t0::Float64)
    texp::Float64 = age - pulse - chase
    meanidx::Vector{Int64} = [2,3]
    #varidx::Vector{Int64} = [12,14]
    #covidx::Vector{Int64} = [13]
    allcovidx::Vector{Int64} = [7:9;]
    ss_iv::Vector{Float64} = get_steady_state_iv(paramvals,vary_map,scaling,n_steps,iv,cycle)
    t_steps::Vector{Float64} = [0.0:cycle/n_steps:cycle;]
    sol_, time_ = full_trajectories(paramvals,vary_map,scaling,t_steps,ss_iv,age,cycle,pulse,t0,texp)
    mean_sol = [s_[meanidx] for s_ in sol_]
    cov_sol = [s_[allcovidx] for s_ in sol_]
    mean_mat = Matrix{Float64}(undef, (length(time_), length(meanidx)))
    cov_mat = Matrix{Float64}(undef, (length(time_), length(allcovidx)))
    for j in 1:size(mean_mat)[2]
       mean_mat[:,j] = [mean_sol[i][j] for i in 1:length(time_)]
    end
    for j in 1:size(cov_mat)[2]
       cov_mat[:,j] = [cov_sol[i][j] for i in 1:length(time_)]
    end
    return mean_mat, cov_mat, time_
end

function solution_plotter(paramvals::AbstractArray{Float64}, vary_map::Vector{Any}, n_steps::Int64, scaling::Int64, iv::Vector{Float64}, t0::Float64,
                cycle::Float64, age::Float64, pulse::Float64, chase::Float64)
    spec_labels = ["unlabelled" "labelled"]
    mean_sol, cov_sol, timespan = get_solutions(paramvals, vary_map, scaling, n_steps, iv, age, cycle, pulse, chase, t0)
    error_bars = Matrix{Float64}(undef, (length(timespan),4))
    error_bars[:,1] = mean_sol[:,1] .- sqrt.(abs.(cov_sol[:,1]))
    error_bars[:,2] = mean_sol[:,1] .+ sqrt.(abs.(cov_sol[:,1]))
    error_bars[:,3] = mean_sol[:,2] .- sqrt.(abs.(cov_sol[:,3]))
    error_bars[:,4] = mean_sol[:,2] .+ sqrt.(abs.(cov_sol[:,3]))
    colours = [:steelblue, :darkorange]
    pl = plot(timespan,
                hcat(error_bars[:,1], error_bars[:,2]),
                label = nothing,
                xlabel = "time",
                ylabel = "mean transcript levels",
                fillrange = mean_sol[:,1],
                fillalpha = 0.12,
                linealpha = [0.12 0.12],
                fg_legend = :transparent,
                background_color_legend = :transparent,
                x_ticks = false,
                fillcolor = [colours[1] colours[1]],
                linecolor = [colours[1] colours[1]],
                size = (400,300))
    pl = plot!(timespan,
                hcat(error_bars[:,3],error_bars[:,4]),
                label = nothing,
                xlabel = "time",
                ylabel = "mean transcript levels",
                fillrange = mean_sol[:,2],
                fillalpha = 0.12,
                linealpha = [0.12 0.12],
                x_ticks = false,
                fillcolor = [colours[2] colours[2]],
                linecolor = [colours[2] colours[2]])
    pl = plot!(timespan, mean_sol,
                label = spec_labels,
                xlabel = "time",
                ylabel = "mean transcript levels",
                linewidth = 3,
                x_ticks = false,
                linecolor = [colours[1] colours[2]])
    return pl
end

function plot_rate(θ::Vector{Float64}, vary_flag::Vector{Int64}, n_steps::Int64, scaling::Int64, cycle::Float64, m::Int64)
    col = [:skyblue4, :lightcyan4, :skyblue4, :seagreen4, :gold][m]
    rate_label = ["synthesis rate", "synthesis rate", 
                "burst frequency", "synthesis rate", "decay rate"][m]                
    tspan::Vector{Float64} = [0.0:0.1:cycle-0.1;]
    t_steps::Vector{Float64} = [0.0:cycle/n_steps:cycle;]
    vary_map = get_vary_map(vary_flag,n_steps)
    if m >= 3
        vary_idx::Vector{Int64} = findall(x->x>0,vary_flag)
    else
        vary_idx = [3]
    end
    local denom::Vector{Float64}
    local p::Plots.Plot{Plots.GRBackend}
    y = Matrix{Float64}(undef,(length(vary_flag),length(tspan)))
    for (k,t) in enumerate(tspan)
        y[:,k] = 10 .^(get_rate(θ[1:end-1], vary_map, scaling, t_steps, cycle, t))
    end
    if vary_idx == 3
        denom = y[2,:]
    else
        denom = ones(length(tspan))
    end
    p = plot(tspan ./ cycle, y[vary_idx[1],:] ./ denom,
                    linewidth = 5.0,
                    xlabel = "cell cycle progression",
                    ylabel = rate_label,
                    label = false,
                    size = (300, 200),
                    ylims = (==(m,2)) .* (0,2*maximum(y[vary_idx[1],:])) .+ (!=(m,2)) .* (0,maximum(y[vary_idx[1],:])+1.0),
                    left_margin = 4Plots.mm,
                    bottom_margin = 5Plots.mm,
                    linecolor = col,
                    legend = :topleft,
                    legendfontsize = 12,
                    fg_legend = :transparent)
    return p 
end


function plot_id_specific_stats(data::Vector{Float64}, m::Int64, which_data::Int64)
    data_label = ["mean labelled vs total ratio", "mean correlation", "correlation of means"][which_data]
    col = [:skyblue4, :lightcyan4, :skyblue4, :seagreen4, :gold][m]
    condition_id = hcat([1/4,1/2,3/4,1,2,3,22,22,22,22,22], [0,0,0,0,0,0,0,1,2,4,6])
    ticks = vcat(condition_id[1:6,1], condition_id[7:end,1] .+ condition_id[7:end,2] .- 17)
    x_tick_labels = vcat(["$j" for j in condition_id[1:3,1]], ["$(round(Int64,j))" for j in condition_id[4:6,1]], ["$(round(Int64,i))"* " + " *"$(round(Int64,j))" for i in condition_id[7:end,1] for j in condition_id[7:end,2]])
    idx::Vector{Int64} = [1,4,5,6,7,8,9,10,11]
    p = vspan([0,3], color = [:lightcyan3], alpha = 0.35, label = "pulse conditions")
    p = vspan!([5,11], color = [:mediumpurple4], alpha = 0.15, label = "chase conditions")
    p = plot!(ticks, data, xlabel = "labelling time (h)", ylabel = data_label, ylims = (-0.2,0.8), background_color_legend = :transparent, fg_legend = :transparent, legend = :topright,
             linewidth = 5, markersize = 5, label = false, linecolor = col, xrotation=30, xticks = (ticks[idx], x_tick_labels[idx]), margin = 2mm, size = (400,300))
    return p
end



function plot_moments(s_mean::Vector{Float64}, s_ff::Vector{Float64}, m::Int64)
    col = [:skyblue4, :lightcyan4, :skyblue4, :seagreen4, :gold][m]
    model_label = ["constant scaling", "constant non-scaling", 
                "varying burst frequency", "varying burst size", "varying decay rate"][m]
    agevec = [0.1:0.2:0.9;]        #maximum(s_ff,dims=1)[1,:]
    #lims = [(0.0,maximum(vcat(data[:,1] + data_se[:,1], s_mean_max)) + 1.0), (0.0,5.0)];
    local mean_plot::Plots.Plot{Plots.GRBackend}
    local ff_plot::Plots.Plot{Plots.GRBackend}
    mean_plot = plot(agevec,
            s_mean,
            ylims = (0,100),
            linecolor = col,
            label = false,
            linewidth = 5.0,
            size = (400,300),
            xlabel = "cell cycle progression",
            ylabel = "mean transcript levels")
    ff_plot = plot(agevec,
            s_ff,
            ylims = (0,20),
            linecolor = col,
            label = false,
            linewidth = 5.0,
            size = (400,300),
            xlabel = "cell cycle progression",
            ylabel = "Fano factor of transcript levels")   
    return mean_plot,ff_plot
end

function full_gene_realisation(m::Int64,θ::Vector{Float64},iv::Vector{Float64},agevec::Vector{Float64},cycle::Float64,pulsevec::Vector{Float64},chasevec::Vector{Float64},t0::Float64,
                            n_steps::Int64,downsampling::Bool,betas::Vector{Float64},age::Vector{Float64},pulse_idx::Vector{Float64},chase_idx::Vector{Float64},age_id_distribution::Vector{Float64})
    vary_flag = [[0,0,0,0],[0,0,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,1]][m]
    m_name = ["const_scaling","const_non_scaling","burst_freq","burst_size","decay_rate"][m]
    plot_rate_name = ["const_scaling_synthesis","const_non_scaling_synthesis","burst_freq","burst_size","decay_rate"][m]
    vary_map = get_vary_map(vary_flag,n_steps)
    scaling = 1 * (!=(m,2))
    #compute summary statistics
    s_pulse,s_chase,s_ratios,s_mean_corr,s_corr_mean = run_sim(θ,iv,agevec,cycle,pulsevec,chasevec,t0,vary_map,scaling,n_steps,downsampling,betas,age,pulse_idx,chase_idx,age_id_distribution);
    #plot example single-cell realisation 
    p = solution_plotter(θ,vary_map,scaling,n_steps,iv,t0,cycle,agevec[end],pulsevec[7],chasevec[7]);
    savefig(p,"data/paper_figures/supplement/"*m_name*"_realisation.pdf")
    #plot characteristic kinetic rate along the cell cycle
    p = plot_rate(θ,vary_flag,n_steps,scaling,cycle,m);
    savefig(p,"data/paper_figures/supplement/"*plot_rate_name*".pdf")
    #plot cell cycle-dependent mean and fano factor
    p1,p2 = plot_moments(s_pulse[:,1], s_pulse[:,2],m);
    savefig(p1,"data/paper_figures/supplement/"*m_name*"_mean.pdf")
    savefig(p2,"data/paper_figures/supplement/"*m_name*"_ff.pdf")
    #plot labelling condition-specific ratios
    p = plot_id_specific_stats(s_ratios, m, 1);
    savefig(p,"data/paper_figures/supplement/"*m_name*"_ratios.pdf")
    #plot labelling condition-specific mean correlations
    p = plot_id_specific_stats(s_mean_corr, m, 2);
    savefig(p,"data/paper_figures/supplement/"*m_name*"_mean_corr.pdf")
    #plot labelling condition-specific correlations of means
    p = plot_id_specific_stats(s_corr_mean, m, 3);
    savefig(p,"data/paper_figures/supplement/"*m_name*"_corr_mean.pdf"
end
    
condition_id = hcat([1/4,1/2,3/4,1,2,3,22,22,22,22,22],[0,0,0,0,0,0,0,1,2,4,6]);
cond_idx = [1:11;];
n_steps = length(τ_);
iv = zeros(9);  
iv[1] = 1/2;
cycle = 20.0;
t0 = -4 * cycle;
agevec = τ_ .* cycle;     
pulsevec = condition_id[cond_idx,1];
chasevec = condition_id[cond_idx,2];     
downsampling = false;

# simulate one gene for each model
################################################

# constant scaling gene
m = 1;
θ = log.(10,[0.5,1.0,20.0,0.1,0.7]);
full_gene_realisation(m,θ,iv,agevec,cycle,pulsevec,chasevec,t0,n_steps,downsampling,betas,age,pulse_idx,chase_idx,age_id_distribution)

# constant non-scaling gene
m = 2;
θ = log.(10,[0.5,1.0,20.0,0.1,0.7]);
full_gene_realisation(m,θ,iv,agevec,cycle,pulsevec,chasevec,t0,n_steps,downsampling,betas,age,pulse_idx,chase_idx,age_id_distribution)

# burst frequency gene
m = 3;
θ = log.(10,[0.1,0.1,30.0,0.1,0.1,1.0,15.0,0.1,0.7]);
full_gene_realisation(m,θ,iv,agevec,cycle,pulsevec,chasevec,t0,n_steps,downsampling,betas,age,pulse_idx,chase_idx,age_id_distribution)

# burst size gene
m = 4;
θ = log.(10,[1.0,1.0,1.0,1.0,40.0,1.0,1.0,0.1,0.7]);
full_gene_realisation(m,θ,iv,agevec,cycle,pulsevec,chasevec,t0,n_steps,downsampling,betas,age,pulse_idx,chase_idx,age_id_distribution)

# decay rate gene
m = 5;
θ = log.(10,[2.0,1.0,85.0,2.0,3.0,1.0,1.0,1.5,0.7]);
full_gene_realisation(m,θ,iv,agevec,cycle,pulsevec,chasevec,t0,n_steps,downsampling,betas,age,pulse_idx,chase_idx,age_id_distribution)

