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

#edit with your personal directory:
dir = "~/data/";
uu_data, us_data, lu_data, ls_data, theta, rfp, gfp, experiment = read_all_data(dir,".csv");

#total counts matrix:
total_data = uu_data + us_data + lu_data + ls_data;

#Gene list
all_genes = Int64.(readdlm("data/selected_genes.txt")[:,1])
gene_id = readdlm("data/all_gene_ids.txt")[:,1]
all_gene_names = readdlm("data/all_gene_names.txt")[:,1]
