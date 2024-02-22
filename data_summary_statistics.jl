function get_mean_subset(data::Matrix{Float64})
    sub = Matrix{Float64}(undef,(Int64(size(data)[1]/2),size(data)[2]))
    for i in 1:size(sub)[1]
        sub[i,:] = data[2*i-1,:]
    end
    return sub
end

function get_ff_subset(data::Matrix{Float64})
    sub = Matrix{Float64}(undef,(Int64(size(data)[1]/2),size(data)[2]))
    for i in 1:size(sub)[1]
        sub[i,:] = data[2*i,:]
    end
    return sub
end

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

function get_ccp_se(data::Vector{Float64}, n_clusters::Int64, age_clusters::Vector{Int64})
    n_bootstraps::Int64 = 100
    stats_mean = Matrix{Float64}(undef, (n_clusters,n_bootstraps))
    stats_ff = Matrix{Float64}(undef, (n_clusters,n_bootstraps))
    for j in 1:n_bootstraps
        sample_idx = sample([1:1:length(data);], length(data); replace = true)
        data_ = data[sample_idx]
        age_ = age_clusters[sample_idx]
        stats_mean[:,j] = get_mean_ccp_data(data_, n_clusters, age_)
        stats_ff[:,j] = get_ff_ccp_data(data_, n_clusters, age_)
    end
    mean_se = Vector{Float64}(undef, n_clusters)
    ff_se = Vector{Float64}(undef, n_clusters)
    for i in 1:n_clusters
        mean_se[i] = sqrt((1 / (n_bootstraps - 1)) * sum((stats_mean[i,:] .- mean(stats_mean[i,:])).^2))
        ff_se[i] = sqrt((1 / (n_bootstraps - 1)) * sum((stats_ff[i,:] .- mean(stats_ff[i,:])).^2))
    end
    return hcat(mean_se, ff_se)
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

function get_ratio_se(u_data::Vector{Float64}, l_data::Vector{Float64}, experiment::Vector{Int64}, cond_vec::Vector{Int64})
    n_bootstraps::Int64 = 100
    stats_ratios = Matrix{Float64}(undef, (length(cond_vec), n_bootstraps))
    for k in 1:n_bootstraps
        sample_idx = sample([1:1:length(u_data);], length(u_data); replace = true)
        u_data_ = u_data[sample_idx]
        l_data_ = l_data[sample_idx]
        experiment_ = experiment[sample_idx]
        for j in 1:length(cond_vec)
            cond_idx::Vector{Int64} = findall(x->x==cond_vec[j], experiment_)
            stats_ratios[j,k] = get_cond_specific_ratios(u_data_[cond_idx], l_data_[cond_idx])
        end
    end
    ratio_se = Vector{Float64}(undef, length(cond_vec))
    for i in 1:length(cond_vec)
        ratio_se[i] = sqrt((1 / (n_bootstraps - 1)) * sum((stats_ratios[i,:] .- mean(stats_ratios[i,:])).^2))
    end
    return ratio_se
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

function get_correlation_se(u_data::Vector{Float64}, l_data::Vector{Float64}, n_clusters::Int64, age_clusters::Vector{Int64}, experiment::Vector{Int64}, cond_vec::Vector{Int64})
    n_bootstraps::Int64 = 100
    stats_mean_corr = Matrix{Float64}(undef, (length(cond_vec),n_bootstraps))
    stats_corr_mean = Matrix{Float64}(undef, (length(cond_vec),n_bootstraps))
    for k in 1:n_bootstraps
        sample_idx = sample([1:1:length(u_data);], length(u_data); replace = true)
        u_data_ = u_data[sample_idx]
        l_data_ = l_data[sample_idx]
        experiment_ = experiment[sample_idx]
        age_ = age_clusters[sample_idx]
        #cells_per_id = [findall(x->x==e,experiment_) for e in sort(unique(experiment_))[2:end-1]]
        cells_per_id = [findall(x->x==e,experiment_) for e in sort(unique(experiment_))]
        cells_per_age = [findall(x->x==τ,age_) for τ in sort(unique(age_))]
        cells_age_id = [[intersect(cells_age,cells_id) for cells_age in cells_per_age] for cells_id in cells_per_id]
        age_id_dist_ = hcat([length.(cells) / sum(length.(cells)) for cells in cells_age_id]...)
        age_id_dist_[isnan.(age_id_dist_)] .= 0.0
        for j in 1:length(cond_vec)
            cond_idx::Vector{Int64} = findall(x->x==cond_vec[j], experiment_)
            stats_mean_corr[j,k], stats_corr_mean[j,k] = get_cond_specific_correlations(u_data_[cond_idx],l_data_[cond_idx],n_clusters,age_[cond_idx],age_id_dist_[:,j])
        end
    end
    mean_corr_se = Vector{Float64}(undef, length(cond_vec))
    corr_mean_se = Vector{Float64}(undef, length(cond_vec))
    for i in 1:length(cond_vec)
        mean_corr_se[i] = sqrt((1 / (n_bootstraps - 1)) * sum((stats_mean_corr[i,:] .- mean(stats_mean_corr[i,:])).^2))
        corr_mean_se[i] = sqrt((1 / (n_bootstraps - 1)) * sum((stats_corr_mean[i,:] .- mean(stats_corr_mean[i,:])).^2))
    end
    return mean_corr_se, corr_mean_se
end
 
function weighted_cov(x::Vector{Float64},y::Vector{Float64}, w::Vector{Float64})
    return sum(w .* ((x .- sum(w .* x)) .* (y .- sum(w .* y))))
end

function get_summary_stats(u_data::Vector{Float64}, l_data::Vector{Float64}, n_clusters::Int64, age_clusters::Vector{Int64}, experiment::Vector{Int64}, 
            cond_vec::Vector{Int64}, pulse_idx::Vector{Int64}, chase_idx::Vector{Int64}, age_id_dist::Matrix{Float64})
    pulse_data::Matrix{Float64} = get_ccp_data(u_data[pulse_idx]+l_data[pulse_idx],n_clusters,age_clusters[pulse_idx])
    pulse_se::Matrix{Float64} = get_ccp_se(u_data[pulse_idx]+l_data[pulse_idx],n_clusters,age_clusters[pulse_idx])
    chase_data::Matrix{Float64} = get_ccp_data(u_data[chase_idx]+l_data[chase_idx],n_clusters,age_clusters[chase_idx])
    chase_se::Matrix{Float64} = get_ccp_se(u_data[chase_idx]+l_data[chase_idx],n_clusters,age_clusters[chase_idx])
    ratio_data::Vector{Float64} = get_ratios(u_data,l_data,experiment,cond_vec)
    ratio_se::Vector{Float64} = get_ratio_se(u_data,l_data,experiment,cond_vec)
    mean_corr_data::Vector{Float64},corr_mean_data::Vector{Float64} = get_correlations(u_data,l_data,n_clusters,age_clusters,experiment,cond_vec,age_id_dist)
    mean_corr_se::Vector{Float64},corr_mean_se::Vector{Float64} = get_correlation_se(u_data,l_data,n_clusters,age_clusters,experiment,cond_vec)
    return pulse_data,pulse_se,chase_data,chase_se,ratio_data,ratio_se,mean_corr_data,mean_corr_se,corr_mean_data,corr_mean_se
end

for g in 1:length(genes)
    p_data,p_se,c_data,c_se,r_data,r_se,mc_data,mc_se,cm_data,cm_se = get_summary_stats(u_data[:,g],l_data[:,g],n_clusters,age,experiment,cond_vec,pulse_idx,chase_idx,age_id_distribution)
    open("data/summary_stats/pulse_data.txt", "a") do io
        writedlm(io, transpose(p_data))
    end 
    open("data/summary_stats/pulse_se.txt", "a") do io
        writedlm(io, transpose(p_se))
    end
    open("data/summary_stats/chase_data.txt", "a") do io
        writedlm(io, transpose(c_data))
    end
    open("data/summary_stats/chase_se.txt", "a") do io
        writedlm(io, transpose(c_se))
    end
    open("data/summary_stats/ratio_data.txt", "a") do io
        writedlm(io, transpose(r_data))
    end
    open("data/summary_stats/ratio_se.txt", "a") do io
        writedlm(io, transpose(r_se))
    end
    open("data/summary_stats/mean_corr_data.txt", "a") do io
        writedlm(io, transpose(mc_data))
    end
    open("data/summary_stats/mean_corr_se.txt", "a") do io
        writedlm(io, transpose(mc_se))
    end
    open("data/summary_stats/corr_mean_data.txt", "a") do io
        writedlm(io, transpose(cm_data))
    end
    open("data/summary_stats/corr_mean_se.txt", "a") do io
        writedlm(io, transpose(cm_se))
    end
end
