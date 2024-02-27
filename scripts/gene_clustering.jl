function preprocess_nc_rates(sets::Matrix{Float64},m::Int64)
    vary_flag::Vector{Int64} = [[1,0,0,0],[0,0,1,0],[0,0,0,1]][m]
    vary_map::Vector{Any} = get_vary_map(vary_flag,5)
    vary_idx::Vector{Int64} = filter(x->length(x)>1,vary_map)[1]
    x::Matrix{Float64} = 10 .^sets[:,vary_idx]
    norm_x::Matrix{Float64} = (x .- mean(x,dims=2)) ./ std(x,dims=2)
    return norm_x
end

function get_cosine_similarity(sets::Matrix{Float64},m::Int64)
    norm_sets::Matrix{Float64} = preprocess_nc_rates(sets,m)
    cos_dist::Matrix{Float64} = pairwise(CosineDist(), norm_sets, dims=1)
    return norm_sets,cos_dist
end

function hcl_grid_search(dist::Matrix{Float64},link_method::Symbol,n_clust_max::Int64)
    k_range::Vector{Int64} = [1:n_clust_max;]
    wss::Vector{Float64} = zeros(length(k_range))
    hcl = Clustering.hclust(dist, linkage = link_method)
    for (i,k) in enumerate(k_range)
        assign = Clustering.cutree(hcl, k = k)
        cluster_idx = [findall(x->x==j,assign) for j in unique(assign)]
        for l in 1:k 
            wss[i] += sum(dist[cluster_idx[l],cluster_idx[l]] .^2) 
        end
    end
    wss = wss ./ 2
    tss::Float64 = sum(dist.^2) / 2
    bss::Vector{Float64} = tss .- wss
    pve::Vector{Float64} = bss ./ tss
    return wss,pve
end

function get_hierarchical_clustering(sets::Matrix{Float64},dist::Matrix{Float64},link_method::Symbol,gene_vec::Vector{Int64},n_clusters::Int64)
    hcl = Clustering.hclust(dist, linkage = link_method)
    assign::Vector{Int64} = Clustering.cutree(hcl, k = n_clusters)
    cluster_idx::Vector{Vector{Int64}} = [findall(x->x==i,assign) for i in unique(assign)]
    clustered_sets::Vector{Matrix{Float64}} = [sets[idx,:] for idx in cluster_idx]
    genes::Vector{Int64} = vcat(gene_vec...)
    gene_cluster::Vector{Vector{Int64}} = [genes[idx] for idx in cluster_idx]
    return sets[hcl.order,:],clustered_sets,gene_cluster,cluster_idx
end

id_label = ["pulse_15", "pulse_30", "pulse_45",
            "pulse_60", "pulse_120", "pulse_180", "chase_0",
            "chase_60", "chase_120", "chase_240", "chase_360"];
cycle_label = ["G1", "G1/S", "S", "S/G2", "G2/M"];

#observed mean expression
mean_data = hcat([mean(t_data[cells,:], dims = 1)[1,:] for cells in cells_per_age]...);
mean_ccp = [mean_data[sel_,:] for sel_ in sel_genes];

#model-recovered mean expression
s_u_mean = [[readdlm("Julia/recovered_statistics/"*mn*"/"*id*"/mean_u.txt") for id in id_label] for mn in model_names];
s_l_mean = [[readdlm("Julia/recovered_statistics/"*mn*"/"*id*"/mean_l.txt") for id in id_label] for mn in model_names];
s_mean = [s_u .+ s_l for (s_u,s_l) in zip(s_u_mean,s_l_mean)];
s_mean_ccp = [sm[7] for sm in s_mean];

#parameter MAP estimates
map_sets = [readdlm("data/posterior_estimates/map_sets_"*name*".txt") for name in model_names];


