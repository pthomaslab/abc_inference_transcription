using DelimitedFiles
using DataFrames
using Statistics
using Distributions
using CSV
using StatsBase

include("load_data.jl")

uu_data, us_data, lu_data, ls_data, theta, rfp, gfp, experiment, gene_id = read_all_data("data/",".csv");
total_data = uu_data + us_data + lu_data + ls_data;

###############  1st filtering  ###############

#=
tc_l = sum(lu_data .+ ls_data,dims=1)[1,:]
all_genes = findall(x->x>=500.0, tc_l)
writedlm("data/selected_genes.txt", all_genes)
writedlm("data/all_gene_ids.txt", gene_id[all_genes,1])
writedlm("data/all_gene_names.txt", gene_id[all_genes,2])
=#

all_genes = Int64.(readdlm("data/selected_genes.txt")[:,1])
gene_id = readdlm("data/all_gene_ids.txt")[:,1]
all_gene_names = readdlm("data/all_gene_names.txt")[:,1]


pulse_idx = findall(x->x<=6, experiment)
chase_idx = findall(x->x>=7, experiment)

n_clusters = 5
age, age_idx, mean_age, τ_ = age_clusters(theta, n_clusters, "equidistant")
age_distribution = length.(age_idx)/ncells

u_data = uu_data[:,genes] + us_data[:,genes]
l_data = lu_data[:,genes] + ls_data[:,genes]

t_data = u_data + l_data


###############  2nd filtering  ###############

s = sum(total_data[:,all_genes],dims=1)[1,:]
min_tc = size(t_data)[1];
age_freq = length.(age_idx)
gene_tc_per_age = [sum(total_data[age,all_genes],dims=1)[1,:] for age in age_idx]


#select the genes that have on average at least *min_count* counts in at least one cell-cycle phase
sel_1 = [];
min_count = 1;
for g in 1:lastindex(all_genes)
    if sum([gene_tc_per_age[i][g] >= min_count*age_freq[i] for i in 1:n_clusters]) >= 1
        push!(sel_1,g)
    end
end

sel_2 = [];
for g in 1:lastindex(all_genes)
    if s[g] >= min_tc
        push!(sel_2,g)
    end
end

sel = Int64.(union(sel_1,sel_2))

#writedlm("data/selected_genes_main_subset.txt",sel)

