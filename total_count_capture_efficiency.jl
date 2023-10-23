using DelimitedFiles
using CSV
using DataFrames
using Statistics
using Distributions
using LsqFit
using KernelDensity
using Plots
using Plots.PlotMeasures
using Colors
using StatsPlots
using StatsBase


using("load_data.jl")


function filter_quantiles(data::Vector{Float64},lb::Float64,ub::Float64)
    lq::Float64 = quantile(data,lb)
    uq::Float64 = quantile(data,ub)
    sel = findall(x->x>=lq && x<=uq, data)
    return data[sel]
end

all_genes = Int64.(readdlm("Julia/all_data/selected_genes.txt")[:,1])
gene_id = readdlm("Julia/all_data/all_gene_ids.txt")[:,1]
all_gene_names = readdlm("Julia/all_data/all_gene_names.txt")[:,1]
genes = all_genes[Int64.(readdlm("Julia/all_data/selected_genes_main_subset.txt")[:,1])]
ids = gene_id[genes];
gene_names = replace(all_gene_names[genes], NaN => "NaN")

uu_data, us_data, lu_data, ls_data, theta, rfp, gfp, experiment, gene_id = read_all_data("Julia/all_data/",".csv")
total_data = uu_data + us_data + lu_data + ls_data
ncells,ngenes = size(us_data);

u_data = uu_data[:,genes] + us_data[:,genes];
l_data = lu_data[:,genes] + ls_data[:,genes];

t_data = u_data + l_data;

n_clusters = 5
age, age_idx, mean_age, τ_ = age_clusters(theta, n_clusters, "equidistant")
age_distribution = length.(age_idx)/ncells

pulse_idx = findall(x->x>=0 && x<=6, experiment)
chase_idx = findall(x->x>=7, experiment)

age_dist_pulse = length.([intersect(idx,pulse_idx) for idx in age_idx]) ./ length(pulse_idx)
age_dist_chase = length.([intersect(idx,chase_idx) for idx in age_idx]) ./ length(chase_idx)

cells_per_id = [findall(x->x==e,experiment) for e in sort(unique(experiment))[2:end-1]]
cells_per_age = [findall(x->x==τ,age) for τ in sort(unique(age))]
cells_age_id = [[intersect(cells_age,cells_id) for cells_age in cells_per_age] for cells_id in cells_per_id]


mean_beta_pulse = 0.1
ratios_pulse = [tc_pulse[i] / atc_pulse_age[age[pulse_idx[i]]] for i in 1:lastindex(pulse_idx)]
betas_pulse = mean_beta_pulse .* ratios_pulse

mean_beta_chase = 0.2
ratios_chase = [tc_chase[i] / atc_chase_age[age[chase_idx[i]]] for i in 1:lastindex(chase_idx)]
betas_chase = mean_beta_chase .* ratios_chase

betas = Vector{Float64}(undef,ncells)
betas[pulse_idx] = betas_pulse
betas[chase_idx] = betas_chase


#cell size 
s = sum(total_data,dims=2)[:,1];

growth(t,p) = p[1] .+ p[2] .* t
ccp = sort(theta)
fit_line = curve_fit(growth, ccp, log.(10,s[sortperm(theta)]), [1.0,1.0])
fit_param = fit_line.param

cycle_label = ["G1", "G1/S", "S", "S/G2", "G2/M"]; 
col = cgrad(:ice, 50, categorical = true)[[5,15,25,35,45]];
cells = cells_per_age[1]
p = scatter(theta[cells],log.(10,s[cells]),color = col[1], label= cycle_label[1],title = "observed scaling of total count with cell cycle", 
            ylabel="log(total count)",xlabel="cell cycle progression",legendfontsize = 8,legend = :topleft,fg_legend = :transparent, size = (600,400));
for (i,cells) in enumerate(cells_per_age[2:end])
    p = scatter!(theta[cells],log.(10,s[cells]),color = col[i+1], label= cycle_label[i+1])
end
p = plot!(theta,growth(theta,fit_param),color = :mediumorchid4,linewidth = 4,label=false)
p
savefig(p,"Julia/paper_figures/size_ccp.pdf")


growth(t,p) = p[1] .+ p[2] .* t
ccp = sort(theta)
fit_line = curve_fit(growth, ccp, s[sortperm(theta)], [1.0,1.0])
fit_param = fit_line.param




tc_pulse = sum(total_data[pulse_idx,:],dims=2)[:,1]
tc_pulse_age = [sum(total_data[intersect(cells,pulse_idx),:],dims=2)[:,1] for cells in cells_per_age]
atc_pulse_age = mean.(tc_pulse_age)

tc_chase = sum(total_data[chase_idx,:],dims=2)[:,1]
tc_chase_age = [sum(total_data[intersect(cells,chase_idx),:],dims=2)[:,1] for cells in cells_per_age]
atc_chase_age = mean.(tc_chase_age)


tc = sum(total_data,dims=2)[:,1]
tc_age = [sum(total_data[cells,:],dims=2)[:,1] for cells in cells_per_age]
atc_age = mean.(tc_age)

tc_age_id = [[sum(total_data[cells_age,:],dims=2)[:,1] for cells_age in cells_id] for cells_id in cells_age_id]
atc_age_id = [mean.(tc) for tc in tc_age_id]



p = violin([1],betas[intersect(cells_per_age[1],chase_idx)], x_ticks = ([1:5;],cycle_label), ylabel = "capture efficiency", color = col[2],
background_color_legend = nothing, fillalpha = 0.6, label = "chase cells", fg_legend = :transparent,legendfontsize = 10, xtickfontsize = 11, left_margin = 2mm, size = (450,350));
p = violin!([1],betas[intersect(cells_per_age[1],pulse_idx)], color = col[4], label = "pulse cells");
for i in 2:lastindex(cells_per_age)
    p = violin!([i],betas[intersect(cells_per_age[i],chase_idx)],fillalpha = 0.6, color = col[2],label = false)
    p = violin!([i],betas[intersect(cells_per_age[i],pulse_idx)],  color = col[4],label = false)
end
trend_pulse = mean.([betas[intersect(cells,pulse_idx)] for cells in cells_per_age])
trend_chase = mean.([betas[intersect(cells,chase_idx)] for cells in cells_per_age])
p = plot!(trend_pulse, x_ticks = ([1:5;],cycle_label),color = :lightcyan, markersize = 5, markerstrokewidth = 4,marker = :circle,label = false, linewidth = 4)
p = plot!(trend_chase, x_ticks = ([1:5;],cycle_label),color = :mediumpurple4, linealpha = 0.8, markersize = 5, markerstrokewidth = 4,marker = :circle,label = "average total count", linewidth = 4)


savefig(p,"Julia/paper_figures/capture_efficiency_ccp.png")



lb = 0.2;
ub = 0.8;

p = violin([1],filter_quantiles(tc[intersect(cells_per_age[1],chase_idx)],lb,ub), x_ticks = ([1:5;],cycle_label), ylabel = "total count", color = col[2], ylims = (0.0,30000), 
background_color_legend = nothing, fillalpha = 0.6, label = "chase cells", legend = :topleft, fg_legend = :transparent,legendfontsize = 10, xtickfontsize = 11, left_margin = 2mm, size = (450,350));
p = violin!([1],filter_quantiles(tc[intersect(cells_per_age[1],pulse_idx)],lb,ub), color = col[4], label = "pulse cells");
for i in 2:lastindex(cells_per_age)
    p = violin!([i],filter_quantiles(tc[intersect(cells_per_age[i],chase_idx)],lb,ub),fillalpha = 0.6, color = col[2],label = false)
    p = violin!([i],filter_quantiles(tc[intersect(cells_per_age[i],pulse_idx)],lb,ub),  color = col[4],label = false)
end

trend_pulse = mean.([filter_quantiles(tc[intersect(cells,pulse_idx)],lb,ub) for cells in cells_per_age])
trend_chase = mean.([filter_quantiles(tc[intersect(cells,chase_idx)],lb,ub) for cells in cells_per_age])
p = plot!(trend_pulse, x_ticks = ([1:5;],cycle_label),color = :lightcyan, markersize = 5, markerstrokewidth = 4,marker = :circle,label = "average total count", linewidth = 4)
p = plot!(trend_chase, x_ticks = ([1:5;],cycle_label),color = :mediumpurple4, linealpha = 0.8, markersize = 5, markerstrokewidth = 4,marker = :circle,label = false, linewidth = 4)

savefig(p,"Julia/paper_figures/size_ccp_distributions.pdf")







#=
#cartoons
p = scatter([1:2;],[1:2;],xlims=(3,4),ylims=(3,4),yaxis=nothing,xaxis=nothing, label = "true transcript levels", color = :steelblue, ylabel="labelled",xlabel="unlabelled", 
            title = "gene X", left_margin=5mm, size=(500,400),legend = :bottomright,legendfontsize = 11,fg_legend = :transparent)
savefig(p,"cartoon_scatter_u_vs_l_true.svg")

p = scatter!([1:2;],[1:2;],xlims=(3,4),ylims=(3,4),yaxis=nothing,xaxis=nothing, label = "observed transcript levels", color = :skyblue, ylabel="labelled",xlabel="unlabelled", 
            title = "gene X", left_margin=5mm, size=(500,400),legend = :bottomright, legendfontsize = 11,fg_legend = :transparent)
savefig(p,"cartoon_scatter_u_vs_l_observed.svg")
=#