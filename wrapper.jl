############################################################################################################
########## Load required packages and raw data ########## 

using DifferentialEquations
using Sundials
using DelimitedFiles
using CSV
using DataFrames
using Statistics
using Distributions
using StatsBase
using Plots
using Plots.PlotMeasures
using Colors
using StatsPlots
using Catalyst
using LsqFit
using Distances

Plots.default(dpi = 200, size = (400,300), fg_legend = :transparent, background_color_legend = :transparent);

include("load_data.jl")

#load raw data

dir = "~/raw_data/";  #edit with your personal directory
uu_data, us_data, lu_data, ls_data, theta, rfp, gfp, experiment = read_all_data(dir,".csv");

#total counts matrix:
total_data = uu_data + us_data + lu_data + ls_data;


########## Pre-process data ##########

#cluster cells with respect to cell cycle position
n_clusters = 5
age, age_idx, τ_ = age_clusters(theta, n_clusters)

#empirical distributions of cells across the 5 cell cycle stages and 11 labelling conditions
cells_per_id = [findall(x->x==e,experiment) for e in sort(unique(experiment))[2:end-1]]
cells_per_age = [findall(x->x==τ,age) for τ in sort(unique(age))]
cells_age_id = [[intersect(cells_age,cells_id) for cells_age in cells_per_age] for cells_id in cells_per_id]

include("gene_selection.jl")
#gene list
all_genes = Int64.(readdlm("data/selected_genes.txt")[:,1])
gene_id = readdlm("data/all_gene_ids.txt")[:,1]
all_gene_names = readdlm("data/all_gene_names.txt")[:,1]

#unlabelled, labelled and total mRNA count datasets for selected genes
u_data = uu_data[:,genes] + us_data[:,genes];
l_data = lu_data[:,genes] + ls_data[:,genes];
t_data = u_data + l_data;

########## Estimate cell-specific capture efficiencies ##########
age_dist_pulse = length.([intersect(idx,pulse_idx) for idx in age_idx]) ./ length(pulse_idx)
age_dist_chase = length.([intersect(idx,chase_idx) for idx in age_idx]) ./ length(chase_idx)

mean_beta_pulse = 0.1
ratios_pulse = [tc_pulse[i] / atc_pulse_age[age[pulse_idx[i]]] for i in 1:lastindex(pulse_idx)]
betas_pulse = mean_beta_pulse .* ratios_pulse

mean_beta_chase = 0.2
ratios_chase = [tc_chase[i] / atc_chase_age[age[chase_idx[i]]] for i in 1:lastindex(chase_idx)]
betas_chase = mean_beta_chase .* ratios_chase

betas = Vector{Float64}(undef,ncells)
betas[pulse_idx] = betas_pulse
betas[chase_idx] = betas_chase

