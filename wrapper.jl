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
age, τ_ = age_clusters(theta, n_clusters)

#distinguishing between pulse-treated and chase-treated cells
pulse_idx = findall(x->x>=0 && x<=6, experiment)
chase_idx = findall(x->x>=7, experiment)

#empirical distributions of cells across the 5 cell cycle stages and 11 labelling conditions
cells_per_id = [findall(x->x==e,experiment) for e in sort(unique(experiment))[2:end-1]]
cells_per_age = [findall(x->x==τ,age) for τ in sort(unique(age))]
cells_age_id = [[intersect(cells_age,cells_id) for cells_age in cells_per_age] for cells_id in cells_per_id]

age_dist_pulse = length.([intersect(cells,pulse_idx) for cells in cells_per_age]) ./ length(pulse_idx)
age_dist_chase = length.([intersect(cells,chase_idx) for cells in cells_per_age]) ./ length(chase_idx)

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
include("total_count_capture_efficiency.jl")

betas = estimate_betas(total_data,age,pulse_idx,chase_idx,cells_per_age)


