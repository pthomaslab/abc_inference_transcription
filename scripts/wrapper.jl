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
include("load_data.jl")

#gene filtering
include("gene_selection.jl")

#load gene list
genes = Int64.(readdlm("data/selected_genes.txt")[:,1])
gene_names = all_gene_names[genes];
gene_ids = all_gene_ids[genes];

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

#run simulations and generate summary statistics
include("abc_simulation.jl")

############################################################################################################
# 3. Parameter inference and model selection

#choose model index to compute errors against the data
m = 1  #2, 3, 4, 5
model_name = ["const","const_const","kon","alpha","gamma"][m]

#compute errors across models-parameters and genes
include("compute_errors.jl")
#convert error data to a compressed format
include("process_error_files.jl")

#obtain posterior distributions
include("accepted_particles.jl")

#compute constant and non-constant model probabilities
include("model_probs.jl")
#compute conditional model probabilities (for constant and non-constant models separately)
include("constant_model_probs.jl")
include("non_constant_model_probs.jl")

#classify genes wrt to the 5 models
include("model_selection.jl") 

#see model selection outcome
ms_df = CSV.read("data/model_selection/model_selection_results.txt", DataFrame)

############################################################################################################
# 4. Transcription kinetics analysis
include("posterior_kinetics.jl")

############################################################################################################
# 5. Noise decomposition


############################################################################################################
# 6. Constant genes: Scaling properties of gene expression with cell size


############################################################################################################
# 7. Non-constant genes: Properties of cell cycle-dependent gene expression

