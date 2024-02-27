gene_vec = vcat(sel_genes[1:2]...);

linkage_method = :ward;
n_clust_max = 20;

# clustering based on recovered mean expression
mn_ccp = vcat(s_mean_ccp[1:2]...);
norm_sets = (mn_ccp .- mean(mn_ccp,dims=2)) ./ std(mn_ccp,dims=2);
cos_dist = pairwise(CosineDist(), norm_sets, dims=1);

wss, pve = hcl_grid_search(cos_dist,linkage_method,n_clust_max);
p = plot([1:n_clust_max;], wss, linewidth = 4, xlabel = "number of clusters", ylabel = "within-cluster distance", color = :steelblue4, label = false)
p = plot([1:n_clust_max;], pve, linewidth = 4, xlabel = "number of clusters",ylabel = "proportion of variance explained", color = :steelblue4, label = false)

n_clusters = 2;
ordered_sets_0,clustered_sets_0,gene_cluster_0,cluster_idx_0 = get_hierarchical_clustering(norm_sets,cos_dist,linkage_method,sel_genes[1:2],n_clusters);

p = heatmap(vcat(clustered_sets_0...), xticks = ([1:5;],cycle_label),title = "recovered mean expression",ylabel = "gene count", xlabel = "cell cycle progression", 
            colorbar_title = "\n z-score(mean expression)", color = :viridis,clims = (-1.79,1.79), right_margin = 5mm, size = (450,250), dpi=200);
savefig(p,"data/paper_figures/figure_4/recov_mean_expression.pdf") 

# clustering based on observed mean expression
perm = vcat(cluster_idx_0...);
gene_vec = vcat(sel_genes[1:2]...)[perm];
mn_ccp = mean_data[gene_vec,:];
norm_sets = (mn_ccp .- mean(mn_ccp,dims=2)) ./ std(mn_ccp,dims=2);

p = heatmap(norm_sets, xticks = ([1:5;],cycle_label),title = "observed mean expression",xlabel = "cell cycle progression", ylabel = "gene count", colorbar_title = "\n z-score(mean expression)", 
        color = :viridis, clims = (-1.79,1.79), right_margin = 5mm, size = (450,250), dpi=200);
savefig(p,"data/paper_figures/figure_4/observed_mean_expression.pdf") 

