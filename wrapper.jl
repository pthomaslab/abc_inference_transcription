############################################################################################################
# 1. Data

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
cells_per_id = [findall(x->x==e,experiment) for e in sort(unique(experiment))[2:end-1]];
cells_per_age = [findall(x->x==τ,age) for τ in sort(unique(age))];
cells_age_id = [[intersect(cells_age,cells_id) for cells_age in cells_per_age] for cells_id in cells_per_id];

n_cells_id = sum(length.(cells_per_id));
age_id_distribution = [length.(cells) ./ n_cells_id for cells in cells_age_id]
age_dist_pulse = length.([intersect(cells,pulse_idx) for cells in cells_per_age]) ./ length(pulse_idx);
age_dist_chase = length.([intersect(cells,chase_idx) for cells in cells_per_age]) ./ length(chase_idx);

#gene filtering
include("gene_selection.jl")

#load gene list
all_gene_ids = readdlm("data/all_gene_ids.txt")[:,1]
all_gene_names = readdlm("data/all_gene_names.txt")[:,1]
genes = Int64.(readdlm("data/selected_genes.txt")[:,1])
gene_names = all_gene_names[genes];
gene_ids = all_gene_ids[genes];

#define unlabelled, labelled and total mRNA count datasets for selected genes
u_data = uu_data[:,genes] + us_data[:,genes];
l_data = lu_data[:,genes] + ls_data[:,genes];
t_data = u_data + l_data;

########## Estimate cell-specific capture efficiencies ##########
include("total_count_capture_efficiency.jl")
betas = readdlm("data/capture_efficiencies.txt")[:,1];

########## Compute summary statistics ##########
include("data_summary_statistics.jl")

pulse_mean,pulse_ff,pulse_mean_se,pulse_ff_se,chase_mean,chase_ff,chase_mean_se,chase_ff_se,ratio_data,ratio_se,mean_corr_data,mean_corr_se,corr_mean_data,corr_mean_se = load_summary_stats("data/summary_stats/",".txt");

############################################################################################################
# 2. Modelling and ABC simulations

# m: modelling hypothesis index  --> {1: constant rates - scaling with size
#                                     2: constant rates - non-scaling
#                                     3: cell cycle-dependent burst frequency - scaling with size
#                                     4: cell cycle-dependent burst size - scaling with size
#                                     5: cell cycle-dependent decay rate - scaling with size }
#                                       
 
#choose model index to simulate
m = 1   #2, 3, 4, 5
#select number of simulations in current run
n_trials = 250000
#submit simulations locally or on high-performance computing cluster (can submit multiple simulation runs)
submit = 1

include("abc_simulation.jl")

############################################################################################################
# 3. Parameter inference and model selection




