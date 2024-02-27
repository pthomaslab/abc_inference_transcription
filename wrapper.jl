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
using Clustering

Plots.default(dpi = 200, size = (400,300), fg_legend = :transparent, background_color_legend = :transparent);

dir = "~/raw_data/";  #edit with your personal directory

#load raw data
include("scripts/load_data.jl")

#gene filtering
include("scripts/gene_selection.jl")

#load gene list
genes = Int64.(readdlm("data/selected_genes.txt")[:,1])
gene_names = all_gene_names[genes];
gene_ids = all_gene_ids[genes];

########## Estimate cell-specific capture efficiencies ##########
include("scripts/total_count_capture_efficiency.jl")
betas = readdlm("data/capture_efficiencies.txt")[:,1];

########## Compute summary statistics ##########
include("scripts/data_summary_statistics.jl")

#load summary statistics
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
n_trials = 250000;
#submit simulations locally or on high-performance computing cluster (can submit multiple simulation runs)
submit = 1;

#run simulations and generate summary statistics
include("scripts/abc_simulation.jl")

############################################################################################################
# 3. Parameter inference and model selection

#choose model index to compute errors against the data
m = 1  #2, 3, 4, 5
model_name = ["const","const_const","kon","alpha","gamma"][m]

#compute errors across models-parameters and genes
include("scripts/compute_errors.jl")
#convert error data to a compressed format
include("scripts/process_error_files.jl")

#obtain posterior distributions
include("scripts/accepted_particles.jl")

#compute constant and non-constant model probabilities
include("scripts/model_probs.jl")
#compute conditional model probabilities (for constant and non-constant models separately)
include("scripts/constant_model_probs.jl")
include("scripts/non_constant_model_probs.jl")

#classify genes wrt to the 5 models
include("scripts/model_selection.jl") 

#indices of genes across the 5 different models
sel_genes = [Int64.(readdlm("data/model_selection/"*mn*"_genes.txt")[:,1]) for mn in model_names];
#see model selection outcome
ms_df = CSV.read("data/model_selection/model_selection_results.txt", DataFrame);

############################################################################################################
# 4. Transcription kinetics
#select model class
m = 1
#pick a gene index
g = sample(sel_genes[m]);
#plot posterior densities of the gene's parameters and distributions of point estimates of parameters across genes
include("scripts/posterior_kinetics.jl")

############################################################################################################
# 5. Statistics recovery and noise decomposition
#select model class
m = 1  #2, 3, 4, 5
#compute all cell cycle-dependent and labelling-dependent mean and noise statistics (without technical noise)
include("scripts/recover_statistics.jl")

#compute and visualise noise-mean relationships across genes
include("scripts/noise_decomposition.jl")

############################################################################################################
# 6. Constant genes: Scaling properties of gene expression with cell size

#load functions for clustering genes with respect to mean expression or kinetic rate similarity
include("scripts/gene_clustering.jl")

#analyse properties of constant scaling and non-scaling genes
include("scripts/constant_genes.jl")

############################################################################################################
# 7. Non-constant genes: Properties of cell cycle-dependent gene expression

#cluster non-constant genes with respect to their cell cycle-dependent kinetic rates
include("scripts/non_constant_genes.jl")

