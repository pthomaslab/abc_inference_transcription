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


# compare properties of constant scaling vs non-scaling genes
mat_1 = mean_data[sel_genes[1],:] ./ mean_data[sel_genes[1],1];
mat_2 = mean_data[sel_genes[2],:] ./ mean_data[sel_genes[2],1];
col = [:skyblue4, :lightcyan4];
h0 = boxplot(mat_1[:,end], linewidth = 1, color = col[1], label = false, legendfontsize = 11, xtickfontsize = 11, legend = :topright, 
            dpi = 200, size = (400,300), xticks = ([1,2],["scaling genes","non-scaling genes"]), ylabel = "mean expression fold change", outliers = false);
h0 = boxplot!(mat_2[:,end], linewidth = 1, color = col[2], label = false, outliers = false);
savefig(h0,"data/paper_figures/figure_4/fold_change_const_genes.svg")


mat_1 = mean_data[sel_genes[1],:];
mat_2 = mean_data[sel_genes[2],:];
der_1 = [mean(diff(mat_1[i,:])) for i in 1:size(mat_1)[1]];
der_2 = [mean(diff(mat_2[i,:])) for i in 1:size(mat_2)[1]];
p = boxplot(der_1, linewidth = 1, color = col[1], label = false, legendfontsize = 11, legend = :topright, size = (400,300),
            dpi = 200, xtickfontsize = 11,xticks = ([1,2],["scaling genes","non-scaling genes"]), ylabel = "mean rate of change",outliers=false);
p = boxplot!(der_2, linewidth = 1, color = col[2], label = false,outliers=false);
savefig(p, "data/paper_figures/figure_4/mean_deriv_const_genes.svg")


p = boxplot(10 .^map_sets[1][:,4], linewidth = 1, color = col[1], label = false, legendfontsize = 11,ylabel = "decay rate",xticks = ([1,2],["scaling genes","non-scaling genes"]),
     legend = :topright, dpi = 200, size = (400,300), xtickfontsize = 11, outliers = false);
p = boxplot!(10 .^map_sets[2][:,4], linewidth = 1,  color = col[2], label = false, outliers = false);
savefig(p, "data/paper_figures/figure_4/decay_rate.svg")
