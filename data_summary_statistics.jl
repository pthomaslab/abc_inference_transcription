using DelimitedFiles
using CSV
using DataFrames
using Statistics
using Distributions
using StatsBase

using("data_processing.jl")


########################  1st-2nd order moments  ########################
function get_mean_ccp_data(arr::Vector{Float64}, n_clusters::Int64, age_clusters::Vector{Int64})
    m = Vector{Float64}(undef,n_clusters)
    for c in sort(unique(age_clusters))
        cells = findall(x->x==c, age_clusters)
        m[c] = mean(arr[cells])
    end
    if length(unique(age_clusters)) < n_clusters
        missing_idx = findall(x->x ∉ sort(unique(age_clusters)), [1:n_clusters;])
        for idx in missing_idx
            splice!(m, idx:idx-1, 0.0)
        end
    end
    return m
end

function get_ff_ccp_data(arr::Vector{Float64}, n_clusters::Int64, age_clusters::Vector{Int64})
    ff = Vector{Float64}(undef,n_clusters)
    local ε::Float64
    for c in sort(unique(age_clusters))
        cells = findall(x->x==c, age_clusters)
        ε = 0.0 * (>(mean(arr[cells]),0)) + 0.0001 * (==(mean(arr[cells]),0))
        ff[c] = var(arr[cells]) / (mean(arr[cells]) + ε)
    end
    if length(unique(age_clusters)) < n_clusters
        missing_idx = findall(x->x ∉ sort(unique(age_clusters)), [1:n_clusters;])
        for idx in missing_idx
            splice!(ff, idx:idx-1, 0.0)
        end
    end
    return ff
end

function get_ccp_data(data::Vector{Float64}, n_clusters::Int64, age_clusters::Vector{Int64})
    mean_data::Vector{Float64} = get_mean_ccp_data(data, n_clusters, age_clusters)
    ff_data::Vector{Float64} = get_ff_ccp_data(data, n_clusters, age_clusters)
    return hcat(mean_data, ff_data)
end

########################  ratios  ########################
function get_ratios(u_data::Vector{Float64}, l_data::Vector{Float64}, experiment::Vector{Int64}, cond_vec::Vector{Int64})
    ratio_data = Vector{Float64}(undef, length(cond_vec))
    for j in 1:length(cond_vec)
        cond_idx::Vector{Int64} = findall(x->x==cond_vec[j], experiment)
        ratio_data[j] = get_cond_specific_ratios(u_data[cond_idx], l_data[cond_idx])
    end
    return ratio_data
end

function get_cond_specific_ratios(arr1::Vector{Float64}, arr2::Vector{Float64})
    if mean(arr1) + mean(arr2) > 0
        r = mean(arr2) / (mean(arr1) + mean(arr2))
    else
        r = 0.0
    end
    return r
end

########################  correlations  ########################
function get_correlations(u_data::Vector{Float64}, l_data::Vector{Float64}, n_clusters::Int64, age_clusters::Vector{Int64}, experiment::Vector{Int64}, cond_vec::Vector{Int64},age_id_dist::Matrix{Float64})
   mean_corr_data = Vector{Float64}(undef, length(cond_vec))
   corr_mean_data = Vector{Float64}(undef, length(cond_vec))
   for j in 1:length(cond_vec)
       cond_idx::Vector{Int64} = findall(x->x==cond_vec[j], experiment)
       mean_corr_data[j], corr_mean_data[j] = get_cond_specific_correlations(u_data[cond_idx],l_data[cond_idx],n_clusters,age_clusters[cond_idx],age_id_dist[:,j])
   end
   return mean_corr_data, corr_mean_data
end

function get_cond_specific_correlations(arr1::Vector{Float64}, arr2::Vector{Float64}, n_clusters::Int64, age_clusters::Vector{Int64}, age_dist::Vector{Float64})
    m1 = Vector{Float64}(undef,n_clusters)
    m2 = Vector{Float64}(undef,n_clusters)
    var1 = Vector{Float64}(undef,n_clusters)
    var2 = Vector{Float64}(undef,n_clusters)
    cov12 = Vector{Float64}(undef,n_clusters)
    for c in sort(unique(age_clusters))
        cells = findall(x->x==c, age_clusters)
        m1[c] = mean(arr1[cells])
        m2[c] = mean(arr2[cells])
        var1[c] = var(arr1[cells])
        var2[c] = var(arr2[cells])
        cov12[c] = cov(arr1[cells],arr2[cells])
    end
    if length(unique(age_clusters)) < n_clusters
        missing_idx = findall(x->x ∉ sort(unique(age_clusters)), [1:n_clusters;])
        for idx in missing_idx
            splice!(m1, idx:idx-1, 0.0)
            splice!(m2, idx:idx-1, 0.0)
            splice!(var1, idx:idx-1, 0.0)
            splice!(var2, idx:idx-1, 0.0)
            splice!(cov12, idx:idx-1, 0.0)
        end
    end
    total_var1::Float64 = sum(age_dist .* var1) + weighted_cov(m1,m1,age_dist)
    total_var2::Float64 = sum(age_dist .* var2) + weighted_cov(m2,m2,age_dist)
    if cov12 != zeros(n_clusters) && total_var1 != 0.0 && total_var2 != 0.0
        mean_corr::Float64 = sum(age_dist .* cov12) / sqrt(abs(total_var1 * total_var2))
    else
        mean_corr = 0.0
    end
    if total_var1 != 0.0 && total_var2 != 0.0
        corr_mean::Float64 = weighted_cov(m1,m2,age_dist) / sqrt(abs(total_var1 * total_var2))
    else
        corr_mean = 0.0
    end   
    return mean_corr, corr_mean
end


function weighted_cov(x::Vector{Float64},y::Vector{Float64}, w::Vector{Float64})
    return sum(w .* ((x .- sum(w .* x)) .* (y .- sum(w .* y))))
end


########################  bootstrap standard errors  ########################
function get_standard_errors(u_data::Vector{Float64},l_data::Vector{Float64},n_clusters::Int64,age_clusters::Vector{Int64},experiment::Vector{Int64})
    bootstraps::Int64 = 100
    stats_pulse_mean = Matrix{Float64}(undef, (n_clusters,bootstraps))
    stats_pulse_ff = Matrix{Float64}(undef, (n_clusters,bootstraps))
    stats_chase_mean = Matrix{Float64}(undef, (n_clusters,bootstraps))
    stats_chase_ff = Matrix{Float64}(undef, (n_clusters,bootstraps))
    stats_ratios = Matrix{Float64}(undef, (length(unique(experiment))-2,bootstraps))
    stats_mean_corr = Matrix{Float64}(undef, (length(unique(experiment))-2,bootstraps))
    stats_corr_mean = Matrix{Float64}(undef, (length(unique(experiment))-2,bootstraps))
    for j in 1:bootstraps
        sample_idx::Vector{Int64} = sample([1:1:length(u_data);], length(u_data); replace = true)
        u_::Vector{Float64} = u_data[sample_idx]
        l_ ::Vector{Float64} = l_data[sample_idx]
        age_::Vector{Int64} = age_clusters[sample_idx]
        experiment_ = experiment[sample_idx]
        cond_vec::Vector{Int64} = sort(unique(experiment))[2:end-1]
        pulse_idx::Vector{Int64} = findall(x->x<=6, experiment)
        chase_idx::Vector{Int64} = findall(x->x>=7, experiment)
        #compute new age distribution
        cells_per_id = [findall(x->x==e,experiment_) for e in sort(unique(experiment_))[2:end-1]]
        cells_per_age = [findall(x->x==τ,age_) for τ in sort(unique(age_))]
        cells_age_id = [[intersect(cells_age,cells_id) for cells_age in cells_per_age] for cells_id in cells_per_id]
        age_id_dist_ = [length.(cells) / sum(length.(cells)) for cells in cells_age_id]
        #compute new summary stats
        temp_pulse::Matrix{Float64} = get_ccp_data(u_[pulse_idx]+l_[pulse_idx],n_clusters,age_[pulse_idx])
        stats_pulse_mean[:,j] = temp_pulse[:,1]
        stats_pulse_ff[:,j] = temp_pulse[:,2]
        temp_chase::Matrix{Float64} = get_ccp_data(u_[chase_idx]+l_[chase_idx],n_clusters,age_[chase_idx])
        stats_chase_mean[:,j] = temp_chase[:,1]
        stats_chase_ff[:,j] = temp_chase[:,2]
        stats_ratios[:,j] = get_ratios(u_,l_,experiment,cond_vec)
        stats_mean_corr[:,j], stats_corr_mean[:,j] = get_correlations(u_,l_,n_clusters,age_,experiment_,cond_vec,age_id_dist_)
    end
    pulse_mean_se = Vector{Float64}(undef, n_clusters)
    pulse_ff_se = Vector{Float64}(undef, n_clusters)
    chase_mean_se = Vector{Float64}(undef, n_clusters)
    chase_ff_se = Vector{Float64}(undef, n_clusters)
    for i in 1:n_clusters
        pulse_mean_se[i] = sqrt((1 / (bootstraps - 1)) * sum((stats_pulse_mean[i,:] .- mean(stats_pulse_mean[i,:])).^2))
        pulse_ff_se[i] = sqrt((1 / (bootstraps - 1)) * sum((stats_pulse_ff[i,:] .- mean(stats_pulse_ff[i,:])).^2))
        chase_mean_se[i] = sqrt((1 / (bootstraps - 1)) * sum((stats_chase_mean[i,:] .- mean(stats_chase_mean[i,:])).^2))
        chase_ff_se[i] = sqrt((1 / (bootstraps - 1)) * sum((stats_chase_ff[i,:] .- mean(stats_chase_ff[i,:])).^2))
    end
    ratio_se = Vector{Float64}(undef, length(cond_vec))
    mean_corr_se = Vector{Float64}(undef, length(cond_vec))
    corr_mean_se = Vector{Float64}(undef, length(cond_vec))
    for i in 1:length(cond_vec)
        ratio_se[i] = sqrt((1 / (bootstraps - 1)) * sum((stats_ratios[i,:] .- mean(stats_ratios[i,:])).^2))
        mean_corr_se[i] = sqrt((1 / (bootstraps - 1)) * sum((stats_mean_corr[i,:] .- mean(stats_mean_corr[i,:])).^2))
        corr_mean_se[i] = sqrt((1 / (bootstraps - 1)) * sum((stats_corr_mean[i,:] .- mean(stats_corr_mean[i,:])).^2))
    end
    return hcat(pulse_mean_se,pulse_ff_se), hcat(chase_mean_se,chase_ff_se), ratio_se, mean_corr_se, corr_mean_se
end


function remove_bootstrap_nan(v::AbstractArray{Float64})
    v[isnan.(v)] .= 0.0
    return v
end

function get_summary_stats(u_data::Vector{Float64}, l_data::Vector{Float64}, n_clusters::Int64, age_clusters::Vector{Int64}, experiment::Vector{Int64}, 
            cond_vec::Vector{Int64}, pulse_idx::Vector{Int64}, chase_idx::Vector{Int64}, age_id_dist::Matrix{Float64})
    pulse_data::Matrix{Float64} = get_ccp_data(u_data[pulse_idx]+l_data[pulse_idx],n_clusters,age_clusters[pulse_idx])
    chase_data::Matrix{Float64} = get_ccp_data(u_data[chase_idx]+l_data[chase_idx],n_clusters,age_clusters[chase_idx])
    ratio_data::Vector{Float64} = get_ratios(u_data,l_data,experiment,cond_vec)
    mean_corr_data::Vector{Float64},corr_mean_data::Vector{Float64} = get_correlations(u_data,l_data,n_clusters,age_clusters,experiment,cond_vec,age_id_dist)
    pulse_se,chase_se,ratio_se,mean_corr_se,corr_mean_se = get_standard_errors(u_data,l_data,n_clusters,age_clusters,experiment)
    return pulse_data,pulse_se,chase_data,chase_se,ratio_data,ratio_se,mean_corr_data,mean_corr_se,corr_mean_data,corr_mean_se
end

all_genes = Int64.(readdlm("Julia/all_data/selected_genes.txt")[:,1])
gene_id = readdlm("Julia/all_data/all_gene_ids.txt")[:,1]
all_gene_names = readdlm("Julia/all_data/all_gene_names.txt")[:,1]
genes = all_genes[Int64.(readdlm("Julia/all_data/selected_genes_main_subset.txt")[:,1])]
ids = gene_id[genes];

gene_names = replace(all_gene_names[genes], NaN => "NaN")
uu_data, us_data, lu_data, ls_data, theta, rfp, gfp, experiment, gene_id = read_all_data("Julia/all_data/",".csv")
#total_data = uu_data + us_data + lu_data + ls_data
ncells,ngenes = size(us_data);

u_data = uu_data[:,genes] + us_data[:,genes];
l_data = lu_data[:,genes] + ls_data[:,genes];
#t_data = total_data[:,genes];

n_clusters = 5
age, age_idx, mean_age, τ_ = age_clusters(theta, n_clusters, "equidistant")
#age_distribution = length.(age_idx)/ncells

cond_vec = sort(unique(experiment))[2:end-1]
pulse_idx = findall(x->x<=6, experiment)
chase_idx = findall(x->x>=7, experiment)

cells_per_id = [findall(x->x==e,experiment) for e in sort(unique(experiment))[2:end-1]]
cells_per_age = [findall(x->x==τ,age) for τ in sort(unique(age))]
cells_age_id = [[intersect(cells_age,cells_id) for cells_age in cells_per_age] for cells_id in cells_per_id]

age_id_distribution = hcat([length.(cells) / sum(length.(cells)) for cells in cells_age_id]...)


for g in 1:length(genes)
    p_data,p_se,c_data,c_se,r_data,r_se,mc_data,mc_se,cm_data,cm_se = get_summary_stats(u_data[:,g],l_data[:,g],n_clusters,age,experiment,cond_vec,pulse_idx,chase_idx,age_id_distribution)
    open("Julia/all_data/summary_stats/pulse_data.txt", "a") do io
        writedlm(io, transpose(p_data))
    end 
    open("Julia/all_data/summary_stats/pulse_se.txt", "a") do io
        writedlm(io, transpose(p_se))
    end
    open("Julia/all_data/summary_stats/chase_data.txt", "a") do io
        writedlm(io, transpose(c_data))
    end
    open("Julia/all_data/summary_stats/chase_se.txt", "a") do io
        writedlm(io, transpose(c_se))
    end
    open("Julia/all_data/summary_stats/ratio_data.txt", "a") do io
        writedlm(io, transpose(r_data))
    end
    open("Julia/all_data/summary_stats/ratio_se.txt", "a") do io
        writedlm(io, transpose(r_se))
    end
    open("Julia/all_data/summary_stats/mean_corr_data.txt", "a") do io
        writedlm(io, transpose(mc_data))
    end
    open("Julia/all_data/summary_stats/mean_corr_se.txt", "a") do io
        writedlm(io, transpose(mc_se))
    end
    open("Julia/all_data/summary_stats/corr_mean_data.txt", "a") do io
        writedlm(io, transpose(cm_data))
    end
    open("Julia/all_data/summary_stats/corr_mean_se.txt", "a") do io
        writedlm(io, transpose(cm_se))
    end
end
