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

function load_summary_stats(path::String,ext::String)
    pulse_data::Matrix{Float64} = readdlm(path*"pulse_data"*ext)
    pulse_se::Matrix{Float64} = readdlm(path*"pulse_se"*ext)
    chase_data::Matrix{Float64} = readdlm(path*"chase_data"*ext)
    chase_se::Matrix{Float64} = readdlm(path*"chase_se"*ext)
    ratio_data::Matrix{Float64} = readdlm(path*"ratio_data"*ext)
    ratio_se::Matrix{Float64} = readdlm(path*"ratio_se"*ext)
    mean_corr_data::Matrix{Float64} = readdlm(path*"mean_corr_data"*ext)
    mean_corr_se::Matrix{Float64} = readdlm(path*"mean_corr_se"*ext)
    corr_mean_data::Matrix{Float64} = readdlm(path*"corr_mean_data"*ext)
    corr_mean_se::Matrix{Float64} = readdlm(path*"corr_mean_se"*ext)
    return pulse_data,pulse_se,chase_data,chase_se,ratio_data,ratio_se,mean_corr_data,mean_corr_se,corr_mean_data,corr_mean_se
end

function get_mean_subset(data::Matrix{Float64})
    sub = Matrix{Float64}(undef,(Int64(size(data)[1]/2),size(data)[2]))
    for i in 1:size(sub)[1]
        sub[i,:] = data[2*i-1,:]
    end
    return sub
end

function get_ff_subset(data::Matrix{Float64})
    sub = Matrix{Float64}(undef,(Int64(size(data)[1]/2),size(data)[2]))
    for i in 1:size(sub)[1]
        sub[i,:] = data[2*i,:]
    end
    return sub
end

function get_c_nc_genes(prob::Matrix{Float64}, l_b::Matrix{Float64}, u_b::Matrix{Float64}, gene_vec::Vector{Int64}, n_particles::Vector{Int64}, min_particles::Int64, threshold::Float64)
    sel_model = Vector{Int64}(undef,length(gene_vec))
    for (i,g) in enumerate(gene_vec)
        if prob[g,:] != zeros(2) && n_particles[g] >= min_particles
            if l_b[g,1] > threshold
                sel_model[i]::Int64 = 1
            elseif u_b[g,1] <= threshold
                sel_model[i]::Int64 = 2
            else
                sel_model[i]::Int64 = 0
            end
        else
            sel_model[i]::Int64 = 3
        end
    end
    return [gene_vec[findall(x->x==1,sel_model)], gene_vec[findall(x->x==2,sel_model)], gene_vec[findall(x->x==0,sel_model)], gene_vec[findall(x->x==3,sel_model)]]
end

function cluster_c_genes(prob::Matrix{Float64}, l_b::Matrix{Float64}, u_b::Matrix{Float64}, gene_vec::Vector{Int64})
    sel_model = Int64.(zeros(length(gene_vec)))
    for (i,g) in enumerate(gene_vec)
        if prob[g,:] != zeros(size(prob)[2])
            if l_b[g,1] .>= u_b[g,2]
                sel_model[i] = 1
            elseif l_b[g,2] .>= u_b[g,1]
                sel_model[i] = 2
            end
        end
    end
    return [gene_vec[findall(x->x==1,sel_model)], gene_vec[findall(x->x==2,sel_model)], gene_vec[findall(x->x==0,sel_model)]]
end


function cluster_nc_genes(prob::Matrix{Float64}, l_b::Matrix{Float64}, u_b::Matrix{Float64}, gene_vec::Vector{Int64})
    nc_models::Vector{Vector{Int64}} = repeat([[]],length(gene_vec))
    for (i,g) in enumerate(gene_vec)
        for m1 in 1:size(prob)[2]
            for m2 in deleteat!([1:size(l_b)[2];],m1)
                if prob[g,[m1,m2]] != zeros(2) && l_b[g,m1] .>= u_b[g,m2]
                    push!(nc_models[i],m1) 
                end
            end
        end
    end
    freq::Matrix{Int64} = zeros(length(gene_vec),size(prob)[2])
    for (i,nc) in enumerate(nc_models)
        for m in unique(nc)
            freq[i,m] = length(filter(x->x==m,nc))
        end
    end
    sel::Vector{Vector{Int64}} = repeat([[0]],length(gene_vec))
    for (i,nc) in enumerate(nc_models)
        if length(nc) == size(prob)[2]-2
            sel[i] = [- nc[1]]
        elseif length(nc) == size(prob)[2]-1
            if length(unique(nc)) == 1
                sel[i] = [nc[1]]
            else
                sel[i] = nc
            end
        elseif length(nc) == size(prob)[2]
            sel[i] = findall(x->x==size(prob)[2]-1,freq[i,:])
        end
    end
    one_model_genes = [gene_vec[findall(x->x==[m],sel)] for m in 1:size(prob)[2]]
    two_model_genes = [gene_vec[findall(x->x==[m1,m2],sel)] for m1 in 1:size(prob)[2]-1 for m2 in m1+1:size(prob)[2]]
    maybe_genes = [gene_vec[findall(x->x== [-m],sel)] for m in 1:size(prob)[2]]
    und_genes = [gene_vec[findall(x->x==[0],sel)]]
    return sel, [one_model_genes, two_model_genes, maybe_genes, und_genes]
end


function cluster_nc_genes_alt(l_b::Matrix{Float64},u_b::Matrix{Float64},gene_vec::Vector{Int64})
    sel_model = Int64.(zeros(length(gene_vec)))
    for (i,g) in enumerate(gene_vec)
        for m in 1:size(l_b)[2]
            if sum(l_b[g,m] .> deleteat!(u_b[g,:],m)) == size(l_b)[2]-1
                sel_model[i] = m
            end
        end
    end
    return vcat([gene_vec[findall(x->x==m,sel_model)] for m in 1:size(l_b)[2]],[gene_vec[findall(x->x==0,sel_model)]])
end


function parallel_coordinates(c::Vector{Vector{Int64}},nc::Vector{Vector{Int64}},und::Vector{Vector{Int64}},nf::Vector{Vector{Int64}}, gene_vec::Vector{Int64})
    coord = Matrix{Float64}(undef,(length(gene_vec),length(c)))
    for i in 1:length(c)
        coord[c[i],i] = collect(range(-0.2,0.2,length(c[i])))
        coord[nc[i],i] = collect(range(0.8,1.2,length(nc[i])))
        coord[und[i],i] = collect(range(1.8,2.2,length(und[i])))
        coord[nf[i],i] = collect(range(2.8,3.2,length(nf[i])))
    end
    return coord
end

function get_n_particles(posterior::Vector{Vector{Vector{Int64}}})
    n_particles::Vector{Vector{Int64}} = [length.(p) for p in posterior]
    for (i,p_) in enumerate(posterior)
        zeros_ = findall(x->x==[0],p_)
        n_particles[i][zeros_] .= 0
    end
    return n_particles
end


function sort_by_tc(tc::Vector{Float64},gene_vec::Vector{Int64})
    return gene_vec[sortperm(tc[gene_vec],rev=true)]
end


all_genes = Int64.(readdlm("Julia/all_data/selected_genes.txt")[:,1])
gene_id = readdlm("Julia/all_data/all_gene_ids.txt")[:,1]
all_gene_names = readdlm("Julia/all_data/all_gene_names.txt")[:,1]
genes = all_genes[Int64.(readdlm("Julia/all_data/selected_genes_main_subset.txt")[:,1])]
ids = gene_id[genes];
gene_names = replace(all_gene_names[genes], NaN => "NaN")

uu_data, us_data, lu_data, ls_data, theta, rfp, gfp, experiment, gene_id = read_all_data("Julia/all_data/",".csv")
total_data = uu_data + us_data + lu_data + ls_data
ncells,ngenes = size(ls_data);
t_data = total_data[:,genes]
s = sum(t_data,dims=1)[1,:];

pulse_data,pulse_se,chase_data,chase_se,ratio_data,ratio_se,
mean_corr_data,mean_corr_se,corr_mean_data,corr_mean_se = load_summary_stats("Julia/all_data/selected_summary_stats/", ".txt");

pulse_mean = get_mean_subset(pulse_data)
pulse_mean_se = get_mean_subset(pulse_se)
pulse_ff = get_ff_subset(pulse_data)
pulse_ff_se = get_ff_subset(pulse_se)
chase_mean = get_mean_subset(chase_data)
chase_mean_se = get_mean_subset(chase_se)
chase_ff = get_ff_subset(chase_data)
chase_ff_se = get_ff_subset(chase_se)

all_ss = [log.(10,pulse_mean),log.(10,pulse_ff),log.(10,chase_mean),log.(10,chase_ff),ratio_data,mean_corr_data,corr_mean_data];
data_names = ["log(pulse_mean)","log(pulse_ff)","log(chase_mean)","log(chase_ff)","ratios","mean_corr","corr_mean"];


df = CSV.read("Julia/all_data/paper_predictions.csv", DataFrame)[:,2:end]
strategy = [findall(x->x in df[findall(x->x==s,df[:,4]),:gene_id],ids) for s in string.(unique(df[:,4]))];
genes_paper = [findall(x->x==g, ids)[1] for g in string.(df[:,1])]
#decay_estimates = Matrix(df[:,end-601:end-301]);
#synthesis_estimates = Matrix(df[:,end-300:end]);

ε = 4.8;

dir = "Julia/model_selection/";

#dir = "Julia/model_selection/alternative_error/top_particles/";
#nan_genes = Int64.(readdlm("Julia/nan_genes.txt")[:,1])

model_prob = readdlm(dir*"all/model_prob.txt")
l_b = readdlm(dir*"all/l_bound.txt") 
u_b = readdlm(dir*"all/u_bound.txt") 


c_prob = readdlm(dir*"c_/c_model_prob.txt")
c_ub = readdlm(dir*"c_/c_u_bound.txt") 
c_lb = readdlm(dir*"c_/c_l_bound.txt")

nc_prob = readdlm(dir*"nc_/nc_model_prob.txt")
nc_ub = readdlm(dir*"nc_/nc_u_bound.txt")
nc_lb = readdlm(dir*"nc_/nc_l_bound.txt")


model_name = ["const","const_const","kon","alpha","gamma"]     
posterior = Vector{Vector{Vector{Int64}}}(undef,length(model_name))
for (i,name) in enumerate(model_name)
    file = open("Julia/posteriors/$ε/particles_"*name*".txt", "r")
    posterior[i] = Vector{Vector{Int64}}(undef,length(genes_1))
    for (j,line) in enumerate(eachline(file))
        posterior[i][j] = parse.(Int64,split(line,"\t")) ## or whatever else you want to do with the line.
    end
    close(file)
end  

n_particles_model = get_n_particles(posterior)
n_particles_mat = hcat(n_particles_model...)
n_particles = sum(n_particles_mat, dims=2)[:,1]


#############################    first model selection    ############################# 
gene_vec = [1:length(genes);]
#gene_vec = [1:length(genes);]

prior_odds = 2/3;
thres = prior_odds / (1 + prior_odds)
min_particles = 5
c_, nc_, und_, nf_ = get_c_nc_genes(model_prob,l_b,u_b,gene_vec,n_particles,min_particles,thres)

freq_c = length(c_)
freq_nc = length(nc_)
freq_und = length(und_)
freq_nf = length(nf_)


#=
writedlm("Julia/model_selection/all/nc_genes_$ε.txt",nc_)
writedlm("Julia/model_selection/all/c_genes_$ε.txt",c_)

writedlm(dir*"all/c_gene_ids.txt",ids[c_])
writedlm(dir*"all/nc_gene_ids.txt",ids[nc_])
=#


model_classes = ["constant","non-constant","undetermined"];
col = [:skyblue4, :darkgoldenrod2, :azure]   #cgrad(:deep, 30, categorical = true, rev = true)[[4,15,30]]
b = bar([1,2,3],[freq_c,freq_nc,freq_und+freq_nf] ./ length(gene_vec),labels=false,xticks = ([1,2,3],model_classes), color = col[1:3], fillalpha = 0.95,
            ylabel = "relative frequency", xtickfontsize = 11, size = (450,200),leftmargin=2mm, topmargin=5mm)
b = annotate!([1,2,3],[freq_c,freq_nc,freq_und+freq_nf] ./ length(gene_vec),string.(Int64.(round.(([freq_c,freq_nc,freq_und+freq_nf] ./ length(gene_vec)) .*100,digits=0))) .* "%",:bottom)

b = bar([1,2,3],[freq_c,freq_nc,freq_und+freq_nf] ./ length(gene_vec),labels=false,xticks = ([1,2,3],model_classes), color = col[1:3], fillalpha = 0.95,
            ylabel = "relative frequency", xtickfontsize = 11, size = (450,200),leftmargin=2mm, topmargin=5mm)
b = annotate!([1,2,3],[freq_c,freq_nc,freq_und+freq_nf]./ length(gene_vec),"n = " .*string.([freq_c,freq_nc,freq_und+freq_nf]),:bottom)

savefig(b,"Julia/paper_figures/figure_2/const_vs_non.svg")



#######################    (conditional) constant model selection    #######################

c_genes = cluster_c_genes(c_prob,c_lb,c_ub,c_)

#=
writedlm("Julia/model_selection/c_/const_scaling_ids.txt",ids[c_genes[1]])
writedlm("Julia/model_selection/c_/const_non_scaling_ids.txt",ids[c_genes[2]])

writedlm("Julia/model_selection/const_genes.txt",c_genes[1])

writedlm("Julia/model_selection/const_const_genes.txt",c_genes[2])
=#
freq_c_genes = length.(c_genes)



col = [:skyblue4, :lightcyan3, :white]
model_classes = ["scaling","non-scaling","undetermined"];
b = bar([1,2,3],freq_c_genes ./ length(c_) ,labels=false,xticks = ([1,2,3],model_classes), color = col[1:3], 
            ylabel = "relative frequency", fillalpha = 1.0, xtickfontsize = 11, size = (450,200),leftmargin=2mm, topmargin=5mm)
b = annotate!([1,2,3],freq_c_genes ./ length(c_),string.(Int64.(round.((freq_c_genes ./ length(c_)) .*100,digits=0))) .* "%",:bottom)
savefig(b,"Julia/paper_figures/figure_2/const_genes_freq.svg")

b = bar([1,2,3],freq_c_genes ./ length(c_) ,labels=false,xticks = ([1,2,3],model_classes), color = col[1:3], 
            ylabel = "relative frequency", fillalpha = 1.0, xtickfontsize = 11, size = (450,200),leftmargin=2mm, topmargin=5mm)
b = annotate!([1,2,3],freq_c_genes ./ length(c_),"n = " .* string.(freq_c_genes),:bottom)
savefig(b,"Julia/paper_figures/figure_2/const_genes_freq.svg")


#######################    (conditional) non-constant model selection    #######################
nc_sel, nc_genes = cluster_nc_genes(nc_prob,nc_lb,nc_ub,nc_);

freq_nc_genes = [length.(nc) for nc in nc_genes]

#=
writedlm(dir*"nc_/kon_gene_idx.txt",ids[nc_genes[1][1]])
writedlm(dir*"nc_/alpha_gene_ids.txt",ids[nc_genes[1][2]])
writedlm(dir*"nc_/gamma_gene_ids.txt",ids[nc_genes[1][3]])
writedlm(dir*"nc_/kon_alpha_gene_ids.txt",ids[nc_genes[2][1]])
writedlm(dir*"nc_/kon_gamma_gene_ids.txt",ids[nc_genes[2][2]])
writedlm(dir*"nc_/alpha_gamma_gene_ids.txt",ids[nc_genes[2][3]])

writedlm("Julia/model_selection/kon_genes.txt",nc_genes[1][1])

writedlm("Julia/model_selection/alpha_genes.txt",nc_genes[1][2])

writedlm("Julia/model_selection/gamma_genes.txt",nc_genes[1][3])
=#

plot_freq = vcat([freq_nc_genes[1],freq_nc_genes[2][1:2],[freq_nc_genes[2][3] + sum(freq_nc_genes[3])]]...) ./ length(nc_)
nc_labels = ["burst frequency","burst size","decay rate","burst frequency or burst size","burst frequency or decay rate","undetermined"];
col = [:skyblue4, :seagreen, :gold, :azure, :azure, :azure]    #cgrad(:GnBu, length(nc_labels),rev = true, categorical = true)[[1:length(nc_labels);]];
b = bar([1:length(plot_freq);], plot_freq, label = false, xticks = ([1:length(plot_freq);],nc_labels), fillalpha = 0.8, xtickfontsize = 11, xrotation = 20, top_margin = 5mm, 
    left_margin = 3mm, bottom_margin = 15mm,color = col, fg_legend = :transparent, ylabel = "relative frequency");
b = annotate!([1:length(plot_freq);],0.015 .+ plot_freq,"n = " .* string.(Int64.(plot_freq .* length(nc_))),:bottom, size = (650,250))

savefig(b,"Julia/paper_figures/figure_2/non_const_genes_freq.svg")



##################### summarise model selection results in a DataFrame ##########################
model_sel_df = DataFrame(i = Int64[], gene_name = String[], model = Int64[], model_name = String[])

for c in c_genes[1]
    push!(model_sel_df, [c, gene_names[c], 1, "const_scaling"])
end
for c in c_genes[2]
    push!(model_sel_df, [c, gene_names[c], 2, "const_non_scaling"])
end
for c in c_genes[3]
    push!(model_sel_df, [c, gene_names[c], 0, "const_undet"])
end
for nc in nc_genes[1][1]
    push!(model_sel_df, [nc, gene_names[nc], 3, "burst_frequency"])
end
for nc in nc_genes[1][2]
    push!(model_sel_df, [nc, gene_names[nc], 4, "burst_size"])
end
for nc in nc_genes[1][3]
    push!(model_sel_df, [nc, gene_names[nc], 5, "decay_rate"])
end
for nc in nc_genes[2][1]
    push!(model_sel_df, [nc, gene_names[nc], 0, "burst_freq_or_size"])
end
for nc in nc_genes[2][2]
    push!(model_sel_df, [nc, gene_names[nc], 0, "burst_freq_or_decay"])
end
for nc in nc_genes[2][3]
    push!(model_sel_df, [nc, gene_names[nc], 0, "burst_size_or_decay"])
end
for nc in vcat(nc_genes[3]...,nc_genes[4]...)
    push!(model_sel_df, [nc, gene_names[nc], 0, "non_const_undet"])
end
for g in und_
    push!(model_sel_df, [g, gene_names[g], 0, "undet"])
end
for g in nf_
    push!(model_sel_df, [g, gene_names[g], 0, "not_fitted"])
end

CSV.write("Julia/model_selection/model_selection_results.txt", model_sel_df)

