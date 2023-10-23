#################################################################################################
##################### Model selection based on approximate MAP estimates ########################
#################################################################################################
using DelimitedFiles
using DataFrames
using Statistics
using Distributions
using JDF
using CSV
using Plots
using Plots.PlotMeasures
using Colors
using StatsPlots
using StatsBase
using Clustering

include("load_data.jl")
include("model.jl")



df = CSV.read("Julia/gene_enrichment/c_scaling_2.csv",DataFrame)
col_names = names(df)

categ = unique(df[:,1])

filter_1 = findall(x->x==categ[1], df[:,1])
func = df[filter_1,2]
func = [split(f,"~")[2] for f in func]
p_val = df[filter_1,5]
all_ticks = findall(x->x<0.05,p_val)

ticks = all_ticks[[4,10,13,16,17,19,25,26,49,50,65,78,95,96,134,199]]
p = scatter(df[ticks,3], func[ticks], zcolor = p_val[ticks], xticks = ([0:400:800;],[0:200:800;]),
    ylims = (0,length(ticks)+1),yticks = ([1:length(ticks);],func[ticks]),colorbar_title = "\n p-value", right_margin = 6mm, label = false, xlabel = "gene count", 
            markersize = 7, bottom_margin = 3mm,alpha = 0.8,color = :viridis,size = (550, 250))

savefig(p,"Julia/paper_figures/figure_4/c_scaling_go.svg")
savefig(p,"Julia/paper_figures/figure_4/c_scaling_go.pdf")
    
df = CSV.read("Julia/gene_enrichment/c_scaling_2.csv",DataFrame)

func = df[:,2]
func = [split(f,"~")[2] for f in func]
p_val = df[:,5]
all_ticks = findall(x->x<0.05,p_val)

perm = length(all_ticks) .- [0:1:length(all_ticks)-1;]
ticks = all_ticks[perm]

p = scatter(df[ticks,3], func[ticks], zcolor = p_val[ticks], titlefontsize = 11, xlims = (1,600),
ylims = (0,length(ticks)), yticks = ([0.5:1:length(ticks)-0.5;],func[ticks]),colorbar_title = "\np-value", right_margin = 4mm, label = false, xlabel = "gene count", 
markersize = 7, top_margin = 5mm,bottom_margin = 3mm,alpha = 0.8,color = :viridis,size = (550, 250))
#p = annotate!([29.5],[9.8], "x 10⁻³",annotationfontsize = 10)

savefig(p,"Julia/paper_figures/figure_4/const_scaling_go.svg")
savefig(p,"Julia/paper_figures/figure_4/const_scaling_go.pdf")
 

df = CSV.read("Julia/gene_enrichment/c_non_scaling_3.csv",DataFrame);

func = df[:,2]
func = [split(f,"~")[2] for f in func]
p_val = df[:,5]
all_ticks = findall(x->x<0.05,p_val)

perm = length(all_ticks) .- [0:1:length(all_ticks)-1;]
ticks = all_ticks[perm]

p = scatter(df[ticks,3], func[ticks], zcolor = p_val[ticks], titlefontsize = 11, xlims = (10,100),
ylims = (0,length(ticks)), yticks = ([0.5:1:length(ticks)-0.5;],func[ticks]),colorbar_title = "\np-value", right_margin = 4mm, label = false, xlabel = "gene count", 
markersize = 7, top_margin = 5mm,bottom_margin = 3mm,alpha = 0.8,color = :viridis,size = (550, 250))
#p = annotate!([29.5],[9.8], "x 10⁻³",annotationfontsize = 10)

savefig(p,"Julia/paper_figures/figure_4/const_non_scaling_go.svg")
savefig(p,"Julia/paper_figures/figure_4/const_non_scaling_go.pdf")
 




df = CSV.read("Julia/gene_enrichment/kon.csv",DataFrame)
col_names = names(df)

categ = unique(df[:,1])

filter_1 = findall(x->x==categ[2], df[:,1])
func = df[filter_1,2]
func = [split(f,"~")[2] for f in func]
p_val = df[filter_1,5]
all_ticks = findall(x->x<0.05,p_val)

ticks = all_ticks[vcat([1:10;],14,15,16)]
p = scatter(df[ticks,3], func[ticks], zcolor = p_val[ticks], titlefontsize = 11, xlims = (1,15),
            colorbar_title = "\n p-value", right_margin = 6mm, label = false, xlabel = "gene count", 
            markersize = 7, bottom_margin = 3mm,alpha = 0.8,color = :viridis,size = (550, 250))

savefig(p,"Julia/paper_figures/figure_5/burst_freq_go.svg")
savefig(p,"Julia/paper_figures/figure_5/burst_freq_go.pdf")

df = CSV.read("Julia/gene_enrichment/kon_2.csv",DataFrame)

func = df[:,2]
func = [split(f,"~")[2] for f in func]
p_val = df[:,5]
all_ticks = findall(x->x<0.05,p_val)

perm = length(all_ticks) .- [0:1:length(all_ticks)-1;]
ticks = all_ticks[perm]

p = scatter(df[ticks,3], func[ticks], zcolor = p_val[ticks] .* 10^3, titlefontsize = 11, xlims = (5,26),
ylims = (0,length(ticks)), colorbar_title = "\np-value", right_margin = 4mm, label = false, xlabel = "gene count", 
markersize = 7, top_margin = 5mm,bottom_margin = 3mm,alpha = 0.8,color = :Blues,size = (450, 200))
p = annotate!([29.5],[9.8], "x 10⁻³",annotationfontsize = 10)

savefig(p,"Julia/paper_figures/figure_5/burst_freq_go.svg")
savefig(p,"Julia/paper_figures/figure_5/burst_freq_go.pdf")
 

df = CSV.read("Julia/gene_enrichment/alpha.csv",DataFrame)
col_names = names(df)
categ = unique(df[:,1])

filter_1 = findall(x->x==categ[2], df[:,1])
func = df[filter_1,2]
func = [split(f,"~")[2] for f in func]
p_val = df[filter_1,5]
all_ticks = findall(x->x<0.05,p_val)

ticks = all_ticks
p = scatter(df[ticks,3], func[ticks], zcolor = p_val[ticks], titlefontsize = 11, xlims = (1,6),
            ylims = (0,length(ticks)+1),yticks = ([1:length(ticks);],func[ticks]), colorbar_title = "\n p-value", right_margin = 6mm, label = false, xlabel = "gene count", 
            markersize = 7, bottom_margin = 3mm,alpha = 0.8,color = :Blues,size = (550, 250))

savefig(p,"Julia/paper_figures/figure_4/burst_size_go.svg")
savefig(p,"Julia/paper_figures/figure_4/burst_size_go.pdf")
 

df = CSV.read("Julia/gene_enrichment/alpha_1.csv",DataFrame)

func = df[:,2]
func = [split(f,"~")[2] for f in func]
p_val = df[:,5]
all_ticks = findall(x->x<0.05,p_val)

perm = length(all_ticks) .- [0:1:length(all_ticks)-1;]
ticks = all_ticks[perm]

p = scatter(df[ticks,3], func[ticks], zcolor = p_val[ticks] .* 10^2, titlefontsize = 11, xlims = (1,6),
ylims = (0,length(ticks)), colorbar_title = "\np-value", right_margin = 4mm, label = false, xlabel = "gene count", 
markersize = 7, top_margin = 5mm,bottom_margin = 3mm,alpha = 0.8,color = :Blues,size = (500, 200))
p = annotate!([7],[3.3], "x 10⁻²",annotationfontsize = 10)

savefig(p,"Julia/paper_figures/figure_5/burst_size_go.svg")
savefig(p,"Julia/paper_figures/figure_5/burst_size_go.pdf")
 


df = CSV.read("Julia/gene_enrichment/gamma.csv",DataFrame)
col_names = names(df)
categ = unique(df[:,1])

filter_1 = findall(x->x==categ[2] || x==categ[4], df[:,1])
func = df[filter_1,2]
func = [split(f,"~")[2] for f in func]
p_val = df[filter_1,5]
all_ticks = findall(x->x<0.05,p_val)

ticks = all_ticks[[3,4,5,7,10,14,15,16,17,19,20,22,24,31]]
p = scatter(df[ticks,3], func[ticks], zcolor = p_val[ticks], titlefontsize = 11, xlims = (1,15),
            colorbar_title = "\n p-value", right_margin = 6mm, label = false, xlabel = "gene count", 
            markersize = 7, bottom_margin = 3mm,alpha = 0.8,color = :Blues,size = (550, 250))

savefig(p,"Julia/paper_figures/figure_4/decay_rate_go.svg")
savefig(p,"Julia/paper_figures/figure_4/decay_rate_go.pdf")
             
            
df = CSV.read("Julia/gene_enrichment/gamma_2.csv",DataFrame)

func = df[:,2]
func = [split(f,"~")[2] for f in func]
p_val = df[:,5]
all_ticks = findall(x->x<0.05,p_val)

perm = length(all_ticks) .- [0:1:length(all_ticks)-1;]
ticks = all_ticks[perm]

p = scatter(df[ticks,3], func[ticks], zcolor = p_val[ticks] .* 10^9, titlefontsize = 11, xlims = (1,50),
ylims = (0,length(ticks)), colorbar_title = "\np-value", right_margin = 4mm, label = false, xlabel = "gene count", 
markersize = 7, top_margin = 5mm,bottom_margin = 3mm,alpha = 0.8,color = :Blues,size = (500, 200))
p = annotate!([60],[7.7], "x 10⁻⁹",annotationfontsize = 10)

savefig(p,"Julia/paper_figures/figure_5/decay_rate_go.svg")
savefig(p,"Julia/paper_figures/figure_5/decay_rate_go.pdf")
 

