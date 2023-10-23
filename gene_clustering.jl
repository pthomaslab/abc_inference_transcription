using DelimitedFiles
using DataFrames
using Statistics
using Distributions
using CSV
using Distances
using Clustering
using Plots
using Plots.PlotMeasures
using Colors
using StatsPlots
using StatsBase

include("load_data.jl")
include("model.jl")


function load_param_sets(path::String, model_idx::Int64, n_groups::Int64)
    n::Int64 = 10^6;
    model_name::String = ["const","const_const","kon","alpha","gamma"][model_idx]
    if model_idx <= 2
        sets = Matrix{Float64}(undef,(n_groups*n,5))
    else
        sets = Matrix{Float64}(undef,(n_groups*n,9))
    end
    for i in 1:n_groups
        sets[n*(i-1)+1:n*i,:] = readdlm(path*"/simulations_$i/"*model_name*"/sets_"*model_name*".txt")
    end
    return sets
end

function get_posterior_estimate(sets::Matrix{Float64},posterior_idx::Vector{Vector{Int64}},gene_vec::Vector{Int64},estimate::String)
    estimates = Matrix{Float64}(undef,(length(gene_vec),size(sets)[2]))
    if estimate == "map"
        for (i,g) in enumerate(gene_vec)
            estimates[i,:] = sets[posterior_idx[g][1],:]
        end
    elseif estimate == "mean"
        for (i,g) in enumerate(gene_vec)
            posterior_sets::Matrix{Float64} = sets[posterior_idx[g],:]
            estimates[i,:] = mean(posterior_sets,dims=1)[1,:]
        end
    elseif estimate == "median"
        for (i,g) in enumerate(gene_vec)
            posterior_sets::Matrix{Float64} = sets[posterior_idx[g],:]
            estimates[i,:] = median(posterior_sets,dims=1)[1,:]
        end
    end
    return estimates
end
#################################################################################################
#################################################################################################
#################################################################################################

all_genes = Int64.(readdlm("data/selected_genes.txt")[:,1])
gene_id = readdlm("data/all_gene_ids.txt")[:,1]
all_gene_names = readdlm("data/all_gene_names.txt")[:,1]
genes = all_genes[Int64.(readdlm("data/selected_genes_main_subset.txt")[:,1])]
ids = gene_id[genes];
gene_names = replace(all_gene_names[genes], NaN => "NaN")


uu_data, us_data, lu_data, ls_data, theta, rfp, gfp, experiment, gene_id = read_all_data("data/",".csv");
total_data = uu_data + us_data + lu_data + ls_data;
ncells,ngenes = size(ls_data);
pulse_idx = findall(x->x<=6, experiment)
chase_idx = findall(x->x>=7, experiment)
n_clusters = 5
age, cluster_idx, mean_age, τ_ = age_clusters(theta, n_clusters, "equidistant")
age_distribution = length.(cluster_idx)/ncells;

u_data = uu_data[:,genes] + us_data[:,genes];
l_data = lu_data[:,genes] + ls_data[:,genes];

t_data = u_data + l_data;

mean_data = hcat([mean(t_data[cells,:], dims = 1)[1,:] for cells in cells_per_age]...);


#=
mean_ccp = Matrix{Float64}(undef,(size(t_data)[2],n_clusters))
for j in 1:size(t_data)[2]
    mean_ccp[j,:] = get_cell_cycle_mean_data(t_data[:,j],n_clusters,clusters)
end
=#


ε = 4.8;
model_names = ["const","const_const","kon","alpha","gamma"];     
posterior = Vector{Vector{Vector{Int64}}}(undef,length(model_names))
for (i,name) in enumerate(model_names)
    file = open("data/posteriors/$ε/particles_"*name*".txt", "r")
    posterior[i] = Vector{Vector{Int64}}(undef,length(genes))
    for (j,line) in enumerate(eachline(file))
        posterior[i][j] = parse.(Int64,split(line,"\t")) ## or whatever else you want to do with the line.
    end
    close(file)
end  

n_particles_model = get_n_particles(posterior)
n_particles_mat = hcat(n_particles_model...)
n_particles = sum(n_particles_mat, dims=2)[:,1]


n_groups_sims = 1
sets = [load_param_sets("data/large_scale_simulations/",m,n_groups_sims) for m in 1:5]

maps = Int64.(readdlm("data/posteriors/map_gene_per_model.txt")[:,1:5])


sel_genes = [Int64.(readdlm("data/model_selection/"*mn*"_genes.txt")[:,1]) for mn in model_names]

#map_sets = [sets[m][maps[sel_genes[m],m],:] for m in 1:lastindex(model_names)];


mn = [mean_data[sel_,:] for sel_ in sel_genes]

######################## clustering based on cosine similarity of point estimates of varying kinetics ########################

function preprocess_rates(sets::Vector{Matrix{Float64}},m::Vector{Int64})
    nc_vary_flags::Vector{Vector{Int64}} = [[1,0,0,0],[0,0,1,0],[0,0,0,1]][m];
    vary_maps::Vector{Vector{Any}} = [get_vary_map(vary_flag,5) for vary_flag in nc_vary_flags];
    vary_idx::Vector{Vector{Int64}} = [filter(x->length(x)>1,vary_map)[1] for vary_map in vary_maps]
    vary_kinetics::Vector{Matrix{Float64}} = [sets[i][:,idx] for (i,idx) in zip(m,vary_idx)]
    x::Matrix{Float64} = vcat([10 .^vk for vk in vary_kinetics]...)
    norm_x::Matrix{Float64} = (x .- mean(x,dims=2)) ./ std(x,dims=2)
    return norm_x
end

function get_cosine_similarity(sets::Vector{Matrix{Float64}},m::Vector{Int64})
    norm_sets = preprocess_rates(sets,m)
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
            wss[i] += sum(dist[cluster_idx[l],cluster_idx[l]] .^2) #/ length(cluster_idx[l])
        end
    end
    wss = wss ./ 2
    tss::Float64 = sum(dist.^2) / 2
    bss::Vector{Float64} = tss .- wss
    pve::Vector{Float64} = bss ./ tss
    #k_best::Int64 = findfirst(x->x<0,diff(diff(wss))) - 1
    #k_alt::Int64 = findfirst(x->x>0,diff(diff(pve))) - 1
    return wss,pve
end


function get_hierarchical_clustering(sets::Matrix{Float64},dist::Matrix{Float64},link_method::Symbol,gene_vec::Vector{Vector{Int64}},n_clusters::Int64)
    hcl = Clustering.hclust(dist, linkage = link_method)
    assign::Vector{Int64} = Clustering.cutree(hcl, k = n_clusters)
    cluster_idx::Vector{Vector{Int64}} = [findall(x->x==i,assign) for i in unique(assign)]
    clustered_sets::Vector{Matrix{Float64}} = [sets[idx,:] for idx in cluster_idx]
    genes::Vector{Int64} = vcat(gene_vec...)
    gene_cluster::Vector{Vector{Int64}} = [genes[idx] for idx in cluster_idx]
    #hcl.order
    return sets[hcl.order,:],clustered_sets,gene_cluster,cluster_idx
end


cycle_label = ["G1", "G1/S", "S", "S/G2", "G2/M"];


############################      clustering mean expression of constant genes      ############################


gene_vec = vcat(sel_genes[1:2]...)

linkage_method = :ward;
n_clust_max = 20;

mean_ccp = vcat(s_mean_ccp[1:2]...)

norm_sets = (mean_ccp .- mean(mean_ccp,dims=2)) ./ std(mean_ccp,dims=2)

cos_dist = pairwise(CosineDist(), norm_sets, dims=1)
wss, pve = hcl_grid_search(cos_dist,linkage_method,n_clust_max)

p = plot([1:n_clust_max;], wss, linewidth = 4, xlabel = "number of clusters", ylabel = "within-cluster distance", color = :steelblue4, label = false)
#savefig(p,"wss_maps.png")
p = plot([1:n_clust_max;], pve, linewidth = 4, xlabel = "number of clusters",ylabel = "proportion of variance explained", color = :steelblue4, label = false)
#savefig(p,"pve_maps.png")

n_clusters = 2
ordered_sets_0,clustered_sets_0,gene_cluster_0,cluster_idx_0 = get_hierarchical_clustering(norm_sets,cos_dist,linkage_method,sel_genes[1:2],n_clusters)
hm = heatmap(ordered_sets_0, xticks = ([1:5;],cycle_label),xlabel = "cell cycle progression",ylabel = "normalised mean expression", title = "all constant genes", color = :viridis,
                    size = (450,300))
#savefig(hm,"hcl_mean_expr_const.pdf")

hs_0 = [];
for j in 1:n_clusters
    push!(hs_0,heatmap(clustered_sets_0[j], xticks = ([1:5;],cycle_label),title = "constant genes - cluster $j",ylabel = "normalised mean expression",
            xlabel = "cell cycle progression", color = :viridis,clims = (-1.7828,1.78502), size = (450,300), dpi=300))
end

p = heatmap(vcat(clustered_sets_0...), xticks = ([1:5;],cycle_label),title = "recovered mean expression",ylabel = "gene count", 
        xlabel = "cell cycle progression", colorbar_title = "\n z-score(mean expression)", color = :viridis,clims = (-1.7828,1.78502), right_margin = 5mm, size = (450,250), dpi=300)

savefig(p,"data/paper_figures/figure_4/recov_mean_expression.pdf") 


perm = vcat(cluster_idx_0...)
gene_vec = vcat(sel_genes[1:2]...)[perm]
mean_ccp = mean_data[gene_vec,:]    #mean_ccp = vcat(s_mean[1:2]...)
norm_sets = (mean_ccp .- mean(mean_ccp,dims=2)) ./ std(mean_ccp,dims=2)

p = heatmap(norm_sets, xticks = ([1:5;],cycle_label),title = "observed mean expression",xlabel = "cell cycle progression", ylabel = "gene count", colorbar_title = "\n z-score(mean expression)", 
        color = :viridis,clims = (minimum(norm_sets),maximum(norm_sets)), right_margin = 5mm, size = (450,250), dpi=300)

savefig(p,"data/paper_figures/figure_4/observed_mean_expression.pdf") 



##############################   kinetics and mean expression clustering of non-constant models   ############################## 

map_sets = [get_posterior_estimate(sets[m],posterior[m],sel_genes[m],"map") for m in 3:5]


linkage_method = :ward;
n_clust_max = 20;

# decay rate genes
m = [3]

norm_sets, cos_dist = get_cosine_similarity(map_sets,m)
wss, pve= hcl_grid_search(cos_dist,linkage_method,n_clust_max)

p = plot([1:n_clust_max;], wss, linewidth = 4, xlabel = "number of clusters", ylabel = "within-cluster distance", color = :steelblue4, label = false)
#savefig(p,"wss_maps.png")
p = plot([1:n_clust_max;], pve, linewidth = 4, xlabel = "number of clusters",ylabel = "proportion of variance explained", color = :steelblue4, label = false)
#savefig(p,"pve_maps.png")

n_clusters = 2
ordered_sets_1,clustered_sets_1,gene_cluster_1,cluster_idx_1 = get_hierarchical_clustering(norm_sets,cos_dist,linkage_method,sel_genes[m .+ 2],n_clusters)


hs_1 = [];
for j in 1:n_clusters
    push!(hs_1,heatmap(clustered_sets_1[j], xticks = ([1:5;],cycle_label),ylabel = "normalised MAP rates", xlabel = "cell cycle progression",
        clim = (minimum(norm_sets), maximum(norm_sets)), title = "decay rate genes - cluster $j", color = :viridis, size = (450,300), dpi = 300))
    #savefig(hs_1[j], "hcl_decay_cluster_$j")
end

perm = [2,1]
perm_idx = vcat(cluster_idx_1[perm]...)
p = heatmap(norm_sets[perm_idx,:]', yticks = ([1:5;],cycle_label),xlabel = "gene count", ylabel = "cell cycle progression", 
        clim = (minimum(norm_sets), maximum(norm_sets)), colorbar_title = "\nz-score(rate)",color = :Blues,right_margin = 3mm, size = (450,300), dpi = 300)


savefig(p,"data/paper_figures/figure_5/gamma_clusters.svg")

#heatmap of mean expression for these genes
mean_ccp = mn[5]
x = (mean_ccp .- mean(mean_ccp,dims=2)) ./ std(mean_ccp,dims=2)
p = heatmap(x[perm_idx,:]', yticks = ([1:5;],cycle_label),xlabel = "gene count", ylabel = "cell cycle progression",
         colorbar_title = "\nz-score(mean expression)",color = :Blues,right_margin = 3mm, size = (450,300), dpi = 300)
savefig(p,"data/paper_figures/supplement/gamma_expr_clusters.svg")




# burst size genes
m = [2]


norm_sets, cos_dist = get_cosine_similarity(map_sets,m)
wss, pve = hcl_grid_search(cos_dist,linkage_method,n_clust_max)

p = plot([1:n_clust_max;], wss, linewidth = 4, xlabel = "number of clusters", ylabel = "within-cluster distance", color = :steelblue4, label = false)
p = plot([1:n_clust_max;], pve, linewidth = 4, xlabel = "number of clusters",ylabel = "proportion of variance explained", color = :steelblue4, label = false)


n_clusters = 1
ordered_sets_1,clustered_sets_1,gene_cluster_1,cluster_idx_1 = get_hierarchical_clustering(norm_sets,cos_dist,linkage_method,sel_genes[m .+ 2],n_clusters)

hm = heatmap(hcat(ordered_sets_1'[:,2:end],ordered_sets_1'[:,1]), yticks = ([1:5;],cycle_label),ylabel = "cell cycle progression",xlabel = "gene count", colorbar_title = "\nz-score(rate)",
                     color = :Blues, clim = (minimum(norm_sets), maximum(norm_sets)),right_margin = 3mm,size = (450,300), dpi = 300)
savefig(hm,"data/paper_figures/figure_5/alpha_clusters.svg")

hs_1 = [];
for j in 1:n_clusters
    push!(hs_1,heatmap(clustered_sets_1[j], xticks = ([1:5;],cycle_label),ylabel = "normalised MAP rates", xlabel = "cell cycle progression",
        clim = (minimum(norm_sets), maximum(norm_sets)), title = "burst size genes - cluster $j", color = :viridis, size = (450,300), dpi = 300))
end

mean_ccp = mn[4]
x = (mean_ccp .- mean(mean_ccp,dims=2)) ./ std(mean_ccp,dims=2)
perm_idx = vcat(cluster_idx_1...)
p = heatmap(x[perm_idx,:]', yticks = ([1:5;],cycle_label), xlabel = "gene count", ylabel = "cell cycle progression",
         colorbar_title = "\nz-score(mean expression)",color = :Blues, right_margin = 3mm, size = (450,300), dpi = 300)
savefig(p,"data/paper_figures/supplement/alpha_expr_clusters.svg")






# burst frequency genes

m = [1]
norm_sets, cos_dist = get_cosine_similarity(map_sets,m)
wss, pve = hcl_grid_search(cos_dist,linkage_method,n_clust_max)

p = plot([1:n_clust_max;], wss, linewidth = 4, xlabel = "number of clusters", ylabel = "within-cluster distance", color = :steelblue4, label = false)
#savefig(p,"wss_maps.png")
p = plot([1:n_clust_max;], pve, linewidth = 4, xlabel = "number of clusters",ylabel = "proportion of variance explained", color = :steelblue4, label = false)
#savefig(p,"pve_maps.png")

n_clusters = 3
ordered_sets_1,clustered_sets_1,gene_cluster_1,cluster_idx_1 = get_hierarchical_clustering(norm_sets,cos_dist,linkage_method,sel_genes[m .+ 2],n_clusters)


hs_1 = [];
for j in 1:n_clusters
    push!(hs_1,heatmap(clustered_sets_1[j], xticks = ([1:5;],cycle_label),ylabel = "normalised MAP rates", xlabel = "cell cycle progression",
        clim = (minimum(norm_sets), maximum(norm_sets)), title = "burst frequency genes - cluster $j", color = :viridis, size = (450,300), dpi = 300))
    #savefig(hs_1[j], "hcl_burst_freq_cluster_$j")
end

perm = [1,3,2]
perm_idx = vcat(cluster_idx_1[perm]...)
p = heatmap(norm_sets[perm_idx,:]', yticks = ([1:5;],cycle_label),xlabel = "gene count", ylabel = "cell cycle progression", tick_direction = 
        clim = (minimum(norm_sets), maximum(norm_sets)), colorbar_title = "\nz-score(rate)",color = :Blues,right_margin = 3mm, size = (450,300), dpi = 300)


savefig(p,"data/paper_figures/figure_5/kon_clusters.svg")


mean_ccp = mn[3]
x = (mean_ccp .- mean(mean_ccp,dims=2)) ./ std(mean_ccp,dims=2)
p = heatmap(x[perm_idx,:]', yticks = ([1:5;],cycle_label), xlabel = "gene count", ylabel = "cell cycle progression",
         colorbar_title = "\nz-score(mean expression)",color = :Blues,right_margin = 3mm, size = (450,300), dpi = 300)
savefig(p,"data/paper_figures/supplement/kon_expr_clusters.svg")




