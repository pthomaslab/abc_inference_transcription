nc_map_sets = map_sets[3:5];
linkage_method = :ward;

# decay rate genes
m = 3;
n_clust_max = minimum([20,size(nc_map_sets[m])[1]]);
norm_sets, cos_dist = get_cosine_similarity(nc_map_sets[m],m);
wss, pve= hcl_grid_search(cos_dist,linkage_method,n_clust_max);
plot([1:n_clust_max;], wss, linewidth = 4, xlabel = "number of clusters", ylabel = "within-cluster distance", color = :steelblue4, label = false)
plot([1:n_clust_max;], pve, linewidth = 4, xlabel = "number of clusters",ylabel = "proportion of variance explained", color = :steelblue4, label = false)

n_clusters = 2;
ordered_sets,clustered_sets,gene_cluster,cluster_idx = get_hierarchical_clustering(norm_sets,cos_dist,linkage_method,sel_genes[m+2],n_clusters);

perm = [2,1];
perm_idx = vcat(cluster_idx[perm]...);
p = heatmap(norm_sets[perm_idx,:]', yticks = ([1:5;],cycle_label),xlabel = "gene count", ylabel = "cell cycle progression", 
        clim = (minimum(norm_sets), maximum(norm_sets)), colorbar_title = "\nz-score(rate)",color = :Blues,right_margin = 3mm, size = (450,300), dpi = 200);
savefig(p,"data/paper_figures/figure_5/gamma_clusters.svg")

#heatmap of mean expression for these genes
mn_ccp = mean_ccp[m+2];
x = (mn_ccp .- mean(mn_ccp,dims=2)) ./ std(mn_ccp,dims=2);
p = heatmap(x[perm_idx,:]', yticks = ([1:5;],cycle_label),xlabel = "gene count", ylabel = "cell cycle progression",
         colorbar_title = "\nz-score(mean expression)",color = :Blues,right_margin = 3mm, size = (450,300), dpi = 200);
savefig(p,"data/paper_figures/supplement/gamma_expr_clusters.svg")



# burst size genes
m = 2;
n_clust_max = minimum([20,size(nc_map_sets[m])[1]]);
norm_sets, cos_dist = get_cosine_similarity(nc_map_sets[m],m);
wss, pve = hcl_grid_search(cos_dist,linkage_method,n_clust_max);
plot([1:n_clust_max;], wss, linewidth = 4, xlabel = "number of clusters", ylabel = "within-cluster distance", color = :steelblue4, label = false)
plot([1:n_clust_max;], pve, linewidth = 4, xlabel = "number of clusters",ylabel = "proportion of variance explained", color = :steelblue4, label = false)

n_clusters = 1;
ordered_sets,clustered_sets,gene_cluster,cluster_idx = get_hierarchical_clustering(norm_sets,cos_dist,linkage_method,sel_genes[m+2],n_clusters);

hm = heatmap(hcat(ordered_sets'[:,2:end],ordered_sets'[:,1]), yticks = ([1:5;],cycle_label),ylabel = "cell cycle progression",xlabel = "gene count", colorbar_title = "\nz-score(rate)",
                color = :Blues, clim = (minimum(norm_sets), maximum(norm_sets)),right_margin = 3mm,size = (450,300), dpi = 200)
savefig(hm,"data/paper_figures/figure_5/alpha_clusters.svg")

mn_ccp = mean_ccp[m+2];
x = (mn_ccp .- mean(mn_ccp,dims=2)) ./ std(mn_ccp,dims=2);
perm_idx = vcat(cluster_idx...);
p = heatmap(x[perm_idx,:]', yticks = ([1:5;],cycle_label), xlabel = "gene count", ylabel = "cell cycle progression",
         colorbar_title = "\nz-score(mean expression)",color = :Blues, right_margin = 3mm, size = (450,300), dpi = 200);
savefig(p,"data/paper_figures/supplement/alpha_expr_clusters.svg")



# burst frequency genes
m = 1;
norm_sets, cos_dist = get_cosine_similarity(nc_map_sets[m],m);
wss, pve = hcl_grid_search(cos_dist,linkage_method,n_clust_max);
plot([1:n_clust_max;], wss, linewidth = 4, xlabel = "number of clusters", ylabel = "within-cluster distance", color = :steelblue4, label = false)
plot([1:n_clust_max;], pve, linewidth = 4, xlabel = "number of clusters",ylabel = "proportion of variance explained", color = :steelblue4, label = false)

n_clusters = 3;
ordered_sets,clustered_sets,gene_cluster,cluster_idx = get_hierarchical_clustering(norm_sets,cos_dist,linkage_method,sel_genes[m+2],n_clusters);

perm = [1,3,2];
perm_idx = vcat(cluster_idx[perm]...);
p = heatmap(norm_sets[perm_idx,:]', yticks = ([1:5;],cycle_label),xlabel = "gene count", ylabel = "cell cycle progression", tick_direction = 
        clim = (minimum(norm_sets), maximum(norm_sets)), colorbar_title = "\nz-score(rate)",color = :Blues,right_margin = 3mm, size = (450,300), dpi = 200);
savefig(p,"data/paper_figures/figure_5/kon_clusters.svg")


mn_ccp = mean_ccp[m+2];
x = (mn_ccp .- mean(mn_ccp,dims=2)) ./ std(mn_ccp,dims=2);
p = heatmap(x[perm_idx,:]', yticks = ([1:5;],cycle_label), xlabel = "gene count", ylabel = "cell cycle progression",
         colorbar_title = "\nz-score(mean expression)",color = :Blues,right_margin = 3mm, size = (450,300), dpi = 200)
savefig(p,"data/paper_figures/supplement/kon_expr_clusters.svg")
