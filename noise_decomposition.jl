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
using BenchmarkTools
using Profile
using LsqFit

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


function load_all_s_data(path::String, model_idx::Int64)
    n::Int64 = 10^6;
    model_name::String = ["const","const_const","kon","alpha","gamma"][model_idx]
    s_pulse = Matrix{Float64}(undef,(4*n,5))
    s_chase = Matrix{Float64}(undef,(4*n,5))
    s_ratios = Matrix{Float64}(undef,(2*n,11))
    s_mean_corr = Matrix{Float64}(undef,(2*n,11))
    s_corr_mean = Matrix{Float64}(undef,(2*n,11))
    if model_idx <= 2
        sets = Matrix{Float64}(undef,(2*n,5))
    else
        sets = Matrix{Float64}(undef,(2*n,9))
    end
    for i in 1:2
        s_pulse[2*n*(i-1)+1:2*n*i,:] = readdlm(path*"/simulations_$i/"*model_name*"/s_pulse_"*model_name*".txt")
        s_chase[2*n*(i-1)+1:2*n*i,:] = readdlm(path*"/simulations_$i/"*model_name*"/s_chase_"*model_name*".txt")
        s_ratios[n*(i-1)+1:n*i,:] = readdlm(path*"/simulations_$i/"*model_name*"/s_ratios_"*model_name*".txt")
        s_mean_corr[n*(i-1)+1:n*i,:] = readdlm(path*"/simulations_$i/"*model_name*"/s_mean_corr_"*model_name*".txt")
        s_corr_mean[n*(i-1)+1:n*i,:] = readdlm(path*"/simulations_$i/"*model_name*"/s_corr_mean_"*model_name*".txt")
        sets[n*(i-1)+1:n*i,:] = readdlm(path*"/simulations_$i/"*model_name*"/sets_"*model_name*".txt")
    end
    return s_pulse,s_chase,s_ratios,s_mean_corr,s_corr_mean,sets
end

function load_s_data(path::String, model_idx::Int64, data_idx::Vector{Int64})
    n::Int64 = 10^6;
    model_name::String = ["const","const_const","kon","alpha","gamma"][model_idx]
    data_name::Vector{String} = ["s_pulse_","s_chase_","s_ratios_","s_mean_corr_","s_corr_mean_"][data_idx]
    s_data = Vector{Matrix{Float64}}(undef,length(data_idx))
    for (j,idx) in enumerate(data_idx)
        if idx<=2
            s_data[j] = Matrix{Float64}(undef,(4*n,5))
            k::Int64 = 2
        else
            s_data[j] = Matrix{Float64}(undef,(2*n,11))
            k = 1
        end
        for i in 1:2
            s_data[j][k*n*(i-1)+1:k*n*i,:] = readdlm(path*"/simulations_$i/"*model_name*"/"*data_name[j]*model_name*".txt")
        end
    end
    return s_data
end

function load_param_sets(path::String, model_idx::Int64, n_groups::Int64)
    n::Int64 = 10^6;
    model_name::String = ["const","const_const","kon","alpha","gamma"][model_idx]
    if model_idx <= 2
        sets = Matrix{Float64}(undef,(n_groups*n,5))
    else
        sets = Matrix{Float64}(undef,(n_groups*n,9))
    end
    for i in 1:n_groups
        sets[n*(i-1)+1:n*i,:] = readdlm(path*"/simulations_$i/"*model_name*"/sets_"*model_name*".txt")
    end
    return sets
end

function get_n_particles(posterior::Vector{Vector{Vector{Int64}}})
    n_particles::Vector{Vector{Int64}} = [length.(p) for p in posterior]
    for (i,p_) in enumerate(posterior)
        zeros_ = findall(x->x==[0],p_)
        n_particles[i][zeros_] .= 0
    end
    return n_particles
end

function filter_quantiles(data::Matrix{Float64},q::Float64)
    sel = Vector{Vector}(undef,size(data)[2])
    for j in 1:size(data)[2]
        sel[j] = filter(x->x>=quantile(data[:,j],1-q) && x<=quantile(data[:,j],q), data[:,j])
    end
    l = maximum(length.(sel))
    for i in 1:lastindex(sel)
        if length(sel[i]) != l
            x::Int64 = l - length(sel[i])
            sel[i] = vcat(minimum(sel[i]) .* ones(floor(Int,x/2)), sel[i], maximum(sel[i]) .* ones(ceil(Int,x/2)))
        end
    end
    return hcat(sel...)
end


function weighted_cov(x::Vector{Float64},y::Vector{Float64}, w::Vector{Float64})
    return sum(w .* ((x .- sum(w .* x)) .* (y .- sum(w .* y))))
end

all_genes = Int64.(readdlm("data/selected_genes.txt")[:,1])
gene_id = readdlm("data/all_gene_ids.txt")[:,1]
all_gene_names = readdlm("data/all_gene_names.txt")[:,1]

genes = all_genes[Int64.(readdlm("data/selected_genes_main_subset.txt")[:,1])]

ids = gene_id[genes];

gene_names = replace(all_gene_names[genes], NaN => "NaN")


uu_data, us_data, lu_data, ls_data, theta, rfp, gfp, experiment, gene_id = read_all_data("data/",".csv");
total_data = uu_data + us_data + lu_data + ls_data;
ncells,ngenes = size(ls_data);
pulse_idx = findall(x->x<=6, experiment)
chase_idx = findall(x->x>=7, experiment)
n_clusters = 5
age, cluster_idx, mean_age, τ_ = age_clusters(theta, n_clusters, "equidistant")
age_distribution = length.(cluster_idx) ./ ncells


u_data = uu_data[:,genes] + us_data[:,genes];
l_data = lu_data[:,genes] + ls_data[:,genes];

t_data = u_data + l_data;

cells_per_id = [findall(x->x==e,experiment) for e in sort(unique(experiment))[2:end-1]]
cells_per_age = [findall(x->x==τ,age) for τ in sort(unique(age))]
cells_age_id = [[intersect(cells_age,cells_id) for cells_age in cells_per_age] for cells_id in cells_per_id]
cells_id_age = [[intersect(cells_age,cells_id) for cells_id in cells_per_id] for cells_age in cells_per_age]

n_cells_id = sum(length.(cells_per_id))
age_id_distribution = [length.(cells) ./ n_cells_id for cells in cells_age_id]

id_distribution = length.(cells_per_id) ./ n_cells_id


mean_data = hcat([mean(t_data[cells,:], dims = 1)[1,:] for cells in cells_per_age]...);
var_data = hcat([var(t_data[cells,:], dims = 1)[1,:] for cells in cells_per_age]...);


ε = 4.8;
model_names = ["const","const_const","kon","alpha","gamma"];     
posterior = Vector{Vector{Vector{Int64}}}(undef,length(model_names))
for (i,name) in enumerate(model_names)
    file = open("data/posteriors/$ε/particles_"*name*".txt", "r")
    posterior[i] = Vector{Vector{Int64}}(undef,length(genes))
    for (j,line) in enumerate(eachline(file))
        posterior[i][j] = parse.(Int64,split(line,"\t")) ## or whatever else you want to do with the line.
    end
    close(file)
end  

n_particles_model = get_n_particles(posterior)
n_particles_mat = hcat(n_particles_model...)
n_particles = sum(n_particles_mat, dims=2)[:,1]

n_groups_sims = 1
sets = [load_param_sets("data/large_scale_simulations/",m,n_groups_sims) for m in 1:5]


maps = Int64.(readdlm("data/posteriors/$ε/map_gene_per_model.txt")[:,1:5])


sel_genes = [Int64.(readdlm("data/model_selection/"*mn*"_genes.txt")[:,1]) for mn in model_names]

map_sets = [sets[m][maps[sel_genes[m],m],:] for m in 1:lastindex(model_names)];


col = [:skyblue4, :lightcyan3, :coral3, :mediumpurple3, :goldenrod1]

################################## Based on data ##################################

gene_vec = vcat(sel_genes...)


mn = [mean(t_data[:,sg],dims=1)[1,:] for sg in sel_genes]
v = [var(t_data[:,sg],dims=1)[1,:] for sg in sel_genes]
cv = [v[i] ./ (mn[i].^2) for i in 1:length(mn)] 

y = vcat(minimum(vcat(mn...))-0.2,sort(vcat(mn...)), maximum(vcat(mn...))+2.0)

p = scatter(log.(10,mn[1]), log.(10,cv[1]), color = :gray,  ylabel = "log₁₀(CV²)", xlims = (-1,3.5),ylims = (-2,2), xlabel = "log₁₀(mean)", label = "all genes (data)", 
    markersize = 5, markeralpha = 0.6, legendfontsize = 9,background_color_legend = :transparent, fg_legend = :transparent,xlabelfontsize = 13, ylabelfontsize = 13, size = (300,250), dpi=300)
p = scatter!(log.(10,mn[2]), log.(10,cv[2]), color = :gray, markersize = 5, markeralpha = 0.6,label = false)
p = scatter!(log.(10,mn[3]), log.(10,cv[3]), color = :gray, markersize = 5,markeralpha = 0.6,label = false)
p = scatter!(log.(10,mn[4]), log.(10,cv[4]), color = :gray, markersize = 5,markeralpha = 0.6,label = false);
p = scatter!(log.(10,mn[5]), log.(10,cv[5]), color = :gray, markersize = 5,markeralpha = 0.6,label =  false);
p = plot!(log.(10,y), log.(10, 1 ./ y), label = "1/mean scaling", linewidth = 4, linestyle = :dash, color = :lightblue3)
p = plot!(log.(10,y), ones(length(y)) .* log(10,minimum(minimum.(cv))), label = "extrinsic noise limit", linewidth = 4, linestyle = :dash, color = :coral3, linealpha = 0.75)

savefig(p, "data/paper_figures/figure_3/data_based/noise_landscape_gray.svg")

########################################################################################################################################
col = [:skyblue4, :lightcyan3, :coral3, :mediumpurple3, :goldenrod1]

id_labels = ["pulse_15", "pulse_30", "pulse_45",
            "pulse_60", "pulse_120", "pulse_180", "chase_0",
            "chase_60", "chase_120", "chase_240", "chase_360"]


u_mean = [[hcat([mean(u_data[cell,sel],dims=1)[1,:] for cell in cells]...) for cells in cells_age_id] for sel in sel_genes]
l_mean = [[hcat([mean(l_data[cell,sel],dims=1)[1,:] for cell in cells]...) for cells in cells_age_id] for sel in sel_genes]
tot_mean = u_mean .+ l_mean

u_var = [[hcat([var(u_data[cell,sel],dims=1)[1,:] for cell in cells]...) for cells in cells_age_id] for sel in sel_genes]
l_var = [[hcat([var(l_data[cell,sel],dims=1)[1,:] for cell in cells]...) for cells in cells_age_id] for sel in sel_genes]

ul_cov = Vector{Vector{Matrix{Float64}}}(undef,length(sel_genes))
for i in 1:length(sel_genes)
    ul_cov[i] = Vector{Matrix{Float64}}(undef,length(id_labels))
    for j in 1:length(id_labels)
        ul_cov[i][j] = Matrix{Float64}(undef,(length(sel_genes[i]),5))
        for k in 1:length(sel_genes[i])
            ul_cov[i][j][k,:]::Vector{Float64} = [cov(u_data[cells_age_id[j][c],k],l_data[cells_age_id[j][c],k]) for c in 1:5]
        end
    end
end

mean_ccp = [map(.+,[m .* id_distribution[j] for (j,m) in enumerate(m_)]...) for m_ in tot_mean]


#intrinsic noise components
u_m = [sum(hcat([sum(transpose(transpose(u) .* age_id),dims=2)[:,1] for (u,age_id) in zip(u_,age_id_distribution)]...),dims=2)[:,1] for u_ in u_mean]
l_m = [sum(hcat([sum(transpose(transpose(l) .* age_id),dims=2)[:,1] for (l,age_id) in zip(l_,age_id_distribution)]...),dims=2)[:,1] for l_ in l_mean]
m = u_m .+ l_m
# intrinsic noise components
mean_var_u = [sum(hcat([sum(transpose(transpose(u) .* age_id),dims=2)[:,1] for (u,age_id) in zip(u_,age_id_distribution)]...),dims=2)[:,1] for u_ in u_var]
mean_var_l = [sum(hcat([sum(transpose(transpose(l) .* age_id),dims=2)[:,1] for (l,age_id) in zip(l_,age_id_distribution)]...),dims=2)[:,1] for l_ in l_var]
mean_var = mean_var_u .+ mean_var_l

mean_cov = [sum(hcat([sum(transpose(transpose(ul) .* age_id),dims=2)[:,1] for (ul,age_id) in zip(u_l_,age_id_distribution)]...),dims=2)[:,1] for u_l_ in ul_cov]

int_v = mean_var .+ (2 .* mean_cov)

# extrinsic noise components
ext_v = Vector{Vector{Float64}}(undef,length(sel_genes));
for i in 1:length(sel_genes)
    ext_v[i] = Vector{Float64}(undef,length(sel_genes[i]))
    for j in 1:length(sel_genes[i])
        ext_v[i][j] = weighted_cov(vcat([tot_mean[i][k][j,:] for k in 1:length(id_labels)]...),vcat([tot_mean[i][k][j,:] for k in 1:length(id_labels)]...),vcat(age_id_distribution...))
    end
end


var_mean_ccp = Vector{Vector{Float64}}(undef,length(sel_genes))
for i in 1:length(sel_genes)
    var_mean_ccp[i] = Vector{Float64}(undef,length(sel_genes[i]))
    for j in 1:length(sel_genes[i])
        var_mean_ccp[i][j] = weighted_cov(mean_ccp[i][j,:],mean_ccp[i][j,:],age_dist)  #weighted_cov(s_mean_ccp[i][j,:],s_mean_ccp[i][j,:],age_distribution)
    end
end


v =  int_v .+ ext_v


tot_noise = [v_ ./ (m .^2) for (v_,m) in zip(v,m)]
kinetics_noise = [v_ ./ (m_ .^2) for (v_,m_) in zip(int_v,m)]
ext_noise = [v_ ./ (m_ .^2) for (v_,m_) in zip(ext_v,m)]

tx_deg_noise = [v_ ./ (m_ .^2) for (v_,m_) in zip(mean_var,m)]
bursty_noise = kinetics_noise .- tx_deg_noise

ccd_noise = [v_ ./ (m_ .^2) for (v_,m_) in zip(var_mean_ccp,m)]
time_noise = ext_noise .- ccd_noise



col = [:skyblue4, :lightcyan4]
p = density(filter_quantiles(reshape(ccd_noise[1],:,1),0.99), linewidth = 5, color = col[1], label = "scaling genes", title = "cell-cycle dependence of constant genes ", titlefontsize = 12,
 legendfontsize = 11, background_color_legend = :transparent,fg_legend = :transparent, size = (400,300), ylabel = "probability density", xlabel = "observed cell-cycle dependence noise")
p = density!(ccd_noise[2], linewidth = 5, linealpha = 0.5, color = col[2], label = "non-scaling genes")

savefig(p, "data/paper_figures/figure_4/observed_ccd_noise.svg")

q = 1.0
p = boxplot(filter_quantiles(reshape(ccd_noise[1],:,1),q), linewidth = 1, color = col[1], label = false, 
ylabelfontsize = 10,xtickfontsize = 11, background_color_legend = :transparent,fg_legend = :transparent, size = (400,300), xticks = ([1,2],["scaling genes","non-scaling genes"]), ylabel = "observed cell-cycle dependence noise  ", outliers = false)
p = boxplot!(filter_quantiles(reshape(ccd_noise[2],:,1),q), linewidth = 1, color = col[2], label = false, outliers = false)

savefig(p, "data/paper_figures/figure_4/observed_ccd_noise_1.svg")



#################################      based on recovered statistics      #################################
n_cells_id = length(vcat(cells_per_id...))
age_dist = [sum(length.(cells)) for cells in cells_id_age] ./ n_cells_id    #[sum(length.(cells)) for cells in cells_id_age] ./ n_cells_id

col = [:skyblue4, :lightcyan3, :coral3, :mediumpurple3, :goldenrod1]

id_labels = ["pulse_15", "pulse_30", "pulse_45",
            "pulse_60", "pulse_120", "pulse_180", "chase_0",
            "chase_60", "chase_120", "chase_240", "chase_360"]

s_u_mean = [[readdlm("data/recovered_statistics/"*mn*"/"*id*"/mean_u.txt") for id in id_labels] for mn in model_names]
s_l_mean = [[readdlm("data/recovered_statistics/"*mn*"/"*id*"/mean_l.txt") for id in id_labels] for mn in model_names]
s_mean = [s_u .+ s_l for (s_u,s_l) in zip(s_u_mean,s_l_mean)]

#s_mean_ccp = [map(.+,[transpose(transpose(s) .* age_id_distribution[j]) for (j,s) in enumerate(sm)]...) for sm in s_mean]s_mean_ccp_1 = [map(.+,[s .* id_distribution[j] for (j,s) in enumerate(sm)]...) for sm in s_mean]

s_mean_ccp = [sm[7] for sm in s_mean]

s_u_var = [[readdlm("data/recovered_statistics/"*mn*"/"*id*"/var_u.txt") for id in id_labels] for mn in model_names]
s_l_var = [[readdlm("data/recovered_statistics/"*mn*"/"*id*"/var_l.txt") for id in id_labels] for mn in model_names]

s_ul_cov = [[readdlm("data/recovered_statistics/"*mn*"/"*id*"/cov_ul.txt") for id in id_labels] for mn in model_names]

s_var_ccp = [s_uv[7] .+ s_lv[7] .+ 2 .* s_c[7] for (s_uv,s_lv,s_c) in zip(s_u_var,s_l_var,s_ul_cov)]

s_noise_ccp = [sv ./ (sm .^2) for (sv,sm) in zip(s_var_ccp,s_mean_ccp)]
#intrinsic noise components
s_u = [sum(hcat([sum(transpose(transpose(su) .* age_id),dims=2)[:,1] for (su,age_id) in zip(s_u_,age_id_distribution)]...),dims=2)[:,1] for s_u_ in s_u_mean]
s_l = [sum(hcat([sum(transpose(transpose(sl) .* age_id),dims=2)[:,1] for (sl,age_id) in zip(s_l_,age_id_distribution)]...),dims=2)[:,1] for s_l_ in s_l_mean]
s_m = s_u .+ s_l

#s_m_2 = [sum(hcat([sum(transpose(transpose(s .^2) .* age_id),dims=2)[:,1] for (s,age_id) in zip(sm,age_id_distribution)]...),dims=2)[:,1] for sm in s_mean]
#s_m2_ = [sum(transpose(transpose(sm .^2) .* age_dist),dims=2)[:,1] for sm in s_mean_ccp]
#s_var_mean_time = s_m_2 .- s_m2_


s_mean_var_u = [sum(hcat([sum(transpose(transpose(su) .* age_id),dims=2)[:,1] for (su,age_id) in zip(s_u_,age_id_distribution)]...),dims=2)[:,1] for s_u_ in s_u_var]
s_mean_var_l = [sum(hcat([sum(transpose(transpose(sl) .* age_id),dims=2)[:,1] for (sl,age_id) in zip(s_l_,age_id_distribution)]...),dims=2)[:,1] for s_l_ in s_l_var]
s_mean_var = s_mean_var_u .+ s_mean_var_l

s_mean_cov = [sum(hcat([sum(transpose(transpose(s) .* age_id),dims=2)[:,1] for (s,age_id) in zip(s_,age_id_distribution)]...),dims=2)[:,1] for s_ in s_ul_cov]


# extrinsic noise components
s_ext_v = Vector{Vector{Float64}}(undef,length(sel_genes));
for i in 1:length(sel_genes)
    s_ext_v[i] = Vector{Float64}(undef,length(sel_genes[i]))
    for j in 1:length(sel_genes[i])
        s_ext_v[i][j] = weighted_cov(vcat([s_mean[i][k][j,:] for k in 1:length(id_labels)]...),vcat([s_mean[i][k][j,:] for k in 1:length(id_labels)]...),vcat(age_id_distribution...))
    end
end

s_int_v = s_mean_var .+ (2 .* s_mean_cov)


s_var_mean_ccp = Vector{Vector{Float64}}(undef,length(sel_genes))
for i in 1:length(sel_genes)
    s_var_mean_ccp[i] = Vector{Float64}(undef,length(sel_genes[i]))
    for j in 1:length(sel_genes[i])
        s_var_mean_ccp[i][j] = weighted_cov(s_mean_ccp[i][j,:],s_mean_ccp[i][j,:],age_dist)  #weighted_cov(s_mean_ccp[i][j,:],s_mean_ccp[i][j,:],age_distribution)
    end
end


s_var_mean_time = s_ext_v .- s_var_mean_ccp
s_ext_v = s_var_mean_ccp
s_v =  s_int_v .+ s_ext_v

#CV^2
s_tot_noise = [v_ ./ (m .^2) for (v_,m) in zip(s_v,s_m)]
s_int_noise = [v_ ./ (m .^2) for (v_,m) in zip(s_int_v,s_m)]
s_ext_noise = [v_ ./ (m .^2) for (v_,m) in zip(s_ext_v,s_m)]

s_tx_deg_noise = [v_ ./ (m .^2) for (v_,m) in zip(s_mean_var,s_m)]

s_bursty_noise = s_int_noise .- s_tx_deg_noise;

s_ccd_noise = [v_ ./ (m .^2) for (v_,m) in zip(s_var_mean_ccp,s_m)]
#s_time_noise = s_ext_noise .- s_ccd_noise

s_tx_deg_prop = [tx_deg ./ tot for (tx_deg,tot) in zip(s_tx_deg_noise,s_tot_noise)]
s_bursty_prop = [bursty ./ tot for (bursty,tot) in zip(s_bursty_noise,s_tot_noise)]
s_ccd_prop = [ccd ./ tot for (ccd,tot) in zip(s_ccd_noise,s_tot_noise)]
#s_time_prop = [ccd ./ tot for (ccd,tot) in zip(s_time_noise,s_tot_noise)]


#Fano factor
s_tot_ff = [v_ ./ m for (v_,m) in zip(s_v,s_m)]
s_int_ff = [v_ ./ m for (v_,m) in zip(s_int_v,s_m)]
s_ext_ff = [v_ ./ m for (v_,m) in zip(s_ext_v,s_m)]
s_tx_deg_ff = [v_ ./ m for (v_,m) in zip(s_mean_var,s_m)]

s_bursty_ff = s_int_ff .- s_tx_deg_ff

s_ccd_ff = [v_ ./ m for (v_,m) in zip(s_var_mean_ccp,s_m)]

#s_time_ff = s_ext_ff .- s_ccd_ff


x = vcat([minimum(vcat(s_m...)) - 2.0],sort(vcat(s_m...)))
#title = "total noise vs mean",
p = scatter(log.(10,s_m[1]), log.(10,s_tot_noise[1]), color = col[1],  ylabel = "log₁₀(CV²)",xlims = (-1,3.5),ylims = (-2,2), xlabel = "log₁₀(mean)", xlabelfontsize = 13, ylabelfontsize = 13,
    label = "constant scaling genes", markersize = 5, markeralpha = 0.9, legendfontsize = 9, background_color_legend = :transparent, legend = :topright, fg_legend = :transparent, size = (450,350), dpi = 300);
p = scatter!(log.(10,s_m[2]), log.(10,s_tot_noise[2]), color = col[2], markersize = 5, markeralpha = 0.9,label = "constant non-scaling genes");
p = scatter!(log.(10,s_m[3]), log.(10,s_tot_noise[3]), color = col[3], markersize = 5,markeralpha = 0.9,label = "burst frequency genes");
p = scatter!(log.(10,s_m[4]), log.(10,s_tot_noise[4]), color = col[4], markersize = 5,markeralpha = 0.9,label = "burst size genes");
p = scatter!(log.(10,s_m[5]), log.(10,s_tot_noise[5]), color = col[5], markersize = 5,markeralpha = 0.9,label = "decay rate genes");
p = plot!(log.(10,x), log.(10, 1 ./ x), label = "1/mean scaling", linewidth = 4, linealpha = 0.8, linestyle = :dash, color = :lightblue3)

savefig(p, "data/paper_figures/figure_3/model_based/total_noise_mean.svg")



x = sort(vcat(s_m...))
#title = "total Fano factor vs mean"
p = scatter(log.(10,s_m[1]), log.(10,s_tot_ff[1]), color = col[1],  ylabel = "log₁₀(Fano factor)",xlims = (0,3.5), ylims = (-3,3), xlabel = "log₁₀(mean)", xlabelfontsize = 13,ylabelfontsize = 13,
    label = "constant scaling genes", markersize = 5, markeralpha = 0.9, legendfontsize = 9, background_color_legend = :transparent, legend = :bottomright, fg_legend = :transparent, size = (450,350), dpi = 300);
p = scatter!(log.(10,s_m[2]), log.(10,s_tot_ff[2]), color = col[2], markersize = 5, markeralpha = 0.9,label = "constant non-scaling genes");
p = scatter!(log.(10,s_m[3]), log.(10,s_tot_ff[3]), color = col[3], markersize = 5,markeralpha = 0.9,label = "burst frequency genes");
p = scatter!(log.(10,s_m[4]), log.(10,s_tot_ff[4]), color = col[4], markersize = 5,markeralpha = 0.9,label = "burst size genes");
p = scatter!(log.(10,s_m[5]), log.(10,s_tot_ff[5]), color = col[5], markersize = 5,markeralpha = 0.9,label = "decay rate genes");
p = plot!(log.(10,x), log.(10,ones(length(x))), label = "Poisson limit", linewidth = 4, linealpha = 0.8, linestyle = :dash, color = :lightblue3)

savefig(p, "data/paper_figures/figure_3/model_based/total_ff_mean.svg")


#title = "transcription and decay noise contribution       "
p = scatter(log.(10,s_m[1]), s_tx_deg_prop[1] .* 100, color = col[1], titlefontsize = 12, ylabel = "% of total noise",
    xlims = (0,3.5), ylims = (-5,100.1), xlabel = "log₁₀(mean)", label = false, markersize = 5, markeralpha = 0.9, legendfontsize = 8, legend = :topright, fg_legend = :transparent, size = (400,300), dpi = 300);
p = scatter!(log.(10,s_m[2]), s_tx_deg_prop[2] .* 100, color = col[2], markersize = 5, markeralpha = 0.9,label = false);
p = scatter!(log.(10,s_m[3]), s_tx_deg_prop[3] .* 100, color = col[3], markersize = 5,markeralpha = 0.9,label = false);
p = scatter!(log.(10,s_m[4]), s_tx_deg_prop[4] .* 100, color = col[4], markersize = 5,markeralpha = 0.9,label = false);
p = scatter!(log.(10,s_m[5]), s_tx_deg_prop[5] .* 100, color = col[5], markersize = 5,markeralpha = 0.9,label = false);
p
savefig(p, "data/paper_figures/figure_3/model_based/tx_decay_proportion.svg")

# title = "bursty promoter noise contribution   "
p = scatter(log.(10,s_m[1]), s_bursty_prop[1] .* 100, color = col[1], ylabel = "% of total noise", xlabel = "log₁₀(mean)", 
xlims = (0,3.5), ylims = (-5,100.1),label = false, markersize = 5, markeralpha = 0.9, legendfontsize = 8, legend = :topright, fg_legend = :transparent, size = (400,300), dpi = 300);
p = scatter!(log.(10,s_m[2]), s_bursty_prop[2] .* 100, color = col[2], markersize = 5, markeralpha = 0.9,label = false);
p = scatter!(log.(10,s_m[3]), s_bursty_prop[3] .* 100, color = col[3], markersize = 5,markeralpha = 0.9,label = false);
p = scatter!(log.(10,s_m[4]), s_bursty_prop[4] .* 100, color = col[4], markersize = 5,markeralpha = 0.9,label = false);
p = scatter!(log.(10,s_m[5]), s_bursty_prop[5] .* 100, color = col[5], markersize = 5,markeralpha = 0.9,label = false);
#p = plot!(log.(10,x), log.(10, 1 ./ x), label = "1/mean scaling", linewidth = 4, linestyle = :dash, color = :lightblue3)
p
savefig(p, "data/paper_figures/figure_3/model_based/bursty_promoter_proportion.svg")

#title = "cell-cycle dependence noise contribution       ",
p = scatter(log.(10,s_m[1]), s_ccd_prop[1] .* 100, color = col[1],  titlefontsize = 13, ylabel = "% of total noise", xlabel = "log₁₀(mean)", 
    xlims = (0,3.5), ylims = (-5,100.5), label = false, markersize = 5, markeralpha = 0.9, legendfontsize = 8, legend = :topright, fg_legend = :transparent, size = (400,300), dpi = 300);
p = scatter!(log.(10,s_m[2]), s_ccd_prop[2] .* 100, color = col[2], markersize = 5, markeralpha = 0.9,label = false);
p = scatter!(log.(10,s_m[3]), s_ccd_prop[3] .* 100, color = col[3], markersize = 5,markeralpha = 0.9,label = false);
p = scatter!(log.(10,s_m[4]), s_ccd_prop[4] .* 100, color = col[4], markersize = 5,markeralpha = 0.9,label = false);
p = scatter!(log.(10,s_m[5]), s_ccd_prop[5] .* 100, color = col[5], markersize = 5,markeralpha = 0.9,label = false);
#p = plot!(log.(10,x), log.(10, 1 ./ x), label = "1/mean scaling", linewidth = 4, linestyle = :dash, color = :lightblue3)
p

savefig(p, "data/paper_figures/figure_3/model_based/ccd_proportion.svg")




############################     kinetics regulation on noise-mean plots      ############################

burst_size = vcat(map_sets[1][:,3] .- map_sets[1][:,2],map_sets[2][:,3] .- map_sets[2][:,2])
x = vcat(s_m[1:2]...)[sortperm(burst_size)]
y = vcat(s_bursty_noise[1:2]...)[sortperm(burst_size)]

#x1 = vcat([minimum(vcat(s_m[1:2]...)) - 2.0],sort(vcat(s_m[1:2]...)))

p = scatter(log.(10,x), log.(10,y), zcolor = sort(burst_size), color = :viridis, title = "bursty promoter noise vs mean", ylabel = "log₁₀(noise)",
        xlims = (-1,3.5),ylims = (-10,3.5), xlabel = "log₁₀(mean)", label = nothing, colorbar_title = "log₁₀(burst size)", markersize = 5, markeralpha = 0.9,
         legendfontsize = 9, fg_legend = :transparent, size = (500,400), dpi = 300)
#p = plot!(log.(10,x1), log.(10, 1 ./ x1), label = "1/mean scaling", linewidth = 4, linestyle = :dash, color = :lightblue3)

savefig(p, "data/paper_figures/figure_3/model_based/bursty_noise_burst_size")


x = vcat(s_m[1:2]...)
y = vcat(s_bursty_noise[1:2]...)[sortperm(x)]

p = scatter(burst_size[sortperm(x)], log.(10,y),zcolor = log.(10,sort(x)), color = :viridis, title = "bursty promoter noise vs burst size", ylabel = "log₁₀(CV²)", 
        xlims = (-2.5,4.0),ylims = (-3.0,1.0),xlabel = "log₁₀(burst size)", label = nothing, colorbar_title = "\n log₁₀(mean)", markersize = 5, markeralpha = 0.9,
         legendfontsize = 9, fg_legend = :transparent, right_margin = 5mm, size = (500,400), dpi = 300)
l = [-2.5:0.1:4.0;]
p = plot!(l, l, label = false, linewidth = 4, linealpha = 0.7, linestyle = :dash, color = :lightblue3)
savefig(p, "data/paper_figures/figure_3/model_based/burst_size_bursty_noise")


burst_freq = vcat(map_sets[1][:,1],map_sets[2][:,1])
x = vcat(s_m[1:2]...)[sortperm(burst_freq)]
y = vcat(s_bursty_noise[1:2]...)[sortperm(burst_freq)]

p = scatter(log.(10,x), log.(10,y), zcolor = sort(burst_freq), color = :viridis, title = "bursty promoter noise vs mean", ylabel = "log₁₀(noise)", 
        xlims = (-1,3.5),ylims = (-10,2.0),clims = (-4,2.0), xlabel = "log₁₀(mean)", label = nothing, colorbar_title = "log₁₀(burst frequency)", markersize = 5, markeralpha = 0.9,
         legendfontsize = 9, fg_legend = :transparent, size = (500,400), dpi = 300)
#p = plot!(log.(10,x1), log.(10, 1 ./ x1), label = "1/mean scaling", linewidth = 4, linestyle = :dash, color = :lightblue3)

savefig(p, "data/paper_figures/figure_3/model_based/bursty_noise_burst_freq")


x = vcat(s_m[1:2]...)
y = vcat(s_bursty_noise[1:2]...)[sortperm(x)]

p = scatter(burst_freq[sortperm(x)], log.(10,y),zcolor = log.(10,sort(x)), color = :viridis, title = "bursty promoter noise vs burst frequency", ylabel = "log₁₀(CV²)", 
        xlims = (-3.0,3.0),ylims = (-3.0,1.0),xlabel = "log₁₀(burst frequency)", label = nothing, colorbar_title = "\n log₁₀(mean)", markersize = 5, markeralpha = 0.9,
         legendfontsize = 9, fg_legend = :transparent, right_margin = 5mm, size = (500,400), dpi = 300)
l = [-3.0:0.1:3.0;]
p = plot!(l, -l, label = false, linewidth = 4, linealpha = 0.7, linestyle = :dash, color = :lightblue3)
savefig(p, "data/paper_figures/figure_3/model_based/burst_freq_bursty_noise")

x = vcat(s_m[1:2]...)[sortperm(burst_freq)]
y = vcat(s_bursty_noise[1:2]...)[sortperm(burst_freq)]

p = scatter(log.(10,x), log.(10,y),zcolor = sort(burst_freq), color = :viridis, ylabel = "log₁₀(bursty promoter noise)", 
        xlabel = "log₁₀(mean)", label = nothing, colorbar_title = "\n log₁₀(burst frequency)", markersize = 5, markeralpha = 0.9,
         legendfontsize = 9, fg_legend = :transparent, right_margin = 5mm, size = (500,400), dpi = 300)
savefig(p, "bursty_noise_mean_burst_freq")

x = vcat(s_m[1:2]...)[sortperm(burst_size)]
y = vcat(s_bursty_noise[1:2]...)[sortperm(burst_size)]

p = scatter(log.(10,x), log.(10,y),zcolor = sort(burst_size), color = :viridis, ylabel = "log₁₀(bursty promoter noise)", 
        xlabel = "log₁₀(mean)", label = nothing, colorbar_title = "\n log₁₀(burst size)", markersize = 5, markeralpha = 0.9,
         legendfontsize = 9, fg_legend = :transparent, right_margin = 5mm, size = (500,400), dpi = 300)
savefig(p, "bursty_noise_mean_burst_size")

y = vcat(s_bursty_noise[1:2]...)
p = scatter(burst_size[sortperm(y)], burst_freq[sortperm(y)], zcolor = log.(10,sort(y)), color = :viridis, ylabel = "log₁₀(burst frequency)", 
        xlims = (-2.5,4.5),ylims = (-3.0,3.5), xlabel = "log₁₀(burst size)", label = nothing, colorbar_title = "log₁₀(bursty promoter noise)", markersize = 5, markeralpha = 0.9,
         legendfontsize = 9, fg_legend = :transparent, size = (500,400), dpi = 300)
l = [-2.5:0.1:4.0;]
p = plot!(l, -l, label = false, linewidth = 4, linealpha = 0.7, linestyle = :dash, color = :lightblue3)
savefig(p, "data/paper_figures/figure_3/model_based/burst_freq_bursty_size_1.svg")
savefig(p, "burst_freq_bursty_size_1.png")

y = vcat(s_bursty_noise[1:2]...)[sortperm(burst_freq)]
p = scatter(burst_size[sortperm(burst_freq)], log.(10,y),  zcolor = sort(burst_freq), color = :viridis, ylabel = "log₁₀(bursty promoter noise)", 
        xlims = (-3,4.5), xlabel = "log₁₀(burst size)", label = nothing, colorbar_title = "log₁₀(burst frequency)", markersize = 5, markeralpha = 0.9,
         legendfontsize = 9, fg_legend = :transparent, size = (500,400), dpi = 300)
#l = [-2.5:0.1:4.0;]
#p = plot!(l, l, label = false, linewidth = 4, linealpha = 0.7, linestyle = :dash, color = :lightblue3)
savefig(p, "data/paper_figures/figure_3/model_based/burst_freq_bursty_size_2.svg")
savefig(p, "burst_freq_bursty_size_2.png")


x = vcat(s_m[1:2]...)[sortperm(burst_size)]
y = vcat(s_tx_deg_noise[1:2]...)[sortperm(burst_size)]
p = scatter(log.(10,x), log.(10,y),  zcolor = sort(burst_size), color = :viridis, ylabel = "log₁₀(bursty promoter noise)", 
        xlims = (0,4), xlabel = "log₁₀(burst size)", label = nothing, colorbar_title = "log₁₀(burst frequency)", markersize = 5, markeralpha = 0.9,
         legendfontsize = 9, fg_legend = :transparent, size = (500,400), dpi = 300)


burst_freq = vcat(map_sets[1][:,1] .- map_sets[1][:,4],map_sets[2][:,1] .- map_sets[2][:,4])
x = vcat(s_m[1:2]...)[sortperm(burst_freq)]
y = vcat(s_tx_deg_noise[1:2]...)[sortperm(burst_freq)]

p = plot(log.(10,[minimum(x/10):1.0:maximum(x);]), -log.(10,[minimum(x/10):1.0:maximum(x);]), label = false,linewidth = 4,  color = :lightblue3, linestyle = :dash);
p = scatter!(log.(10,x), log.(10,y), zcolor = sort(burst_freq),color = :viridis, ylabel = "log₁₀(transcription & decay noise)", background_color_legend = :transparent,
        xlims = (0,3.5),ylims = (-3,1),xlabel = "log₁₀(mean)", label = nothing, colorbar_title = "log₁₀(relative burst frequency)", markersize = 5, markeralpha = 0.9,
         legendfontsize = 9, fg_legend = :transparent, size = (500,400), dpi = 300)

savefig(p, "data/paper_figures/figure_3/model_based/burst_freq_tx_deg_noise.svg")

burst_size = vcat(map_sets[1][:,3] .- map_sets[1][:,2],map_sets[2][:,3] .- map_sets[2][:,2])
x = vcat(s_m[1:2]...)[sortperm(burst_size)]
y = vcat(s_tx_deg_noise[1:2]...)[sortperm(burst_size)]

p = plot(log.(10,[minimum(x/10):1.0:maximum(x);]), -log.(10,[minimum(x/10):1.0:maximum(x);]), label = false,linewidth = 4,  color = :lightblue3, linestyle = :dash);
p = scatter!(log.(10,x), log.(10,y), zcolor = sort(burst_size), color = :viridis, ylabel = "log₁₀(transcription & decay noise)", background_color_legend = :transparent,
        xlims = (0,3.5),ylims = (-3,1),xlabel = "log₁₀(mean)", label = nothing, colorbar_title = "log₁₀(burst size)", markersize = 5, markeralpha = 0.9,
         legendfontsize = 9, fg_legend = :transparent, size = (500,400), dpi = 300)

savefig(p, "data/paper_figures/figure_3/model_based/burst_size_tx_deg_noise.svg")
    

y = vcat(s_tx_deg_noise[1:2]...)[sortperm(burst_size)]
p = scatter(burst_freq[sortperm(burst_size)], log.(10,y),  zcolor = sort(burst_size), color = :viridis, ylabel = "log₁₀(bursty promoter noise)", 
        xlims = (-3,4.5), xlabel = "log₁₀(burst size)", label = nothing, colorbar_title = "log₁₀(burst frequency)", markersize = 5, markeralpha = 0.9,
         legendfontsize = 9, fg_legend = :transparent, size = (500,400), dpi = 300)
#l = [-2.5:0.1:4.0;]
#p = plot!(l, l, label = false, linewidth = 4, linealpha = 0.7, linestyle = :dash, color = :lightblue3)
savefig(p, "data/paper_figures/figure_3/model_based/burst_freq_bursty_size_2.svg")
savefig(p, "burst_freq_bursty_size_2.png")



y = vcat(s_tx_deg_noise[1:2]...)[sortperm(burst_freq)]
p = scatter(burst_size[sortperm(burst_freq)], log.(10,y),  zcolor = sort(burst_freq), color = :viridis, ylabel = "bursty promoter noise", 
        xlims = (-3,4.5), xlabel = "log₁₀(burst size)", label = nothing, colorbar_title = "log₁₀(burst frequency)", markersize = 5, markeralpha = 0.9,
         legendfontsize = 9, fg_legend = :transparent, size = (500,400), dpi = 300)
#l = [-2.5:0.1:4.0;]
savefig(p, "burst_freq_bursty_size_3.png")

y = vcat(s_int_noise[1:2]...)[sortperm(burst_size)]
p = scatter(burst_freq[sortperm(burst_size)], log.(10,y),  zcolor = sort(burst_freq), color = :viridis, ylabel = "bursty promoter noise", 
        xlims = (-3,4.5), xlabel = "log₁₀(burst size)", label = nothing, colorbar_title = "log₁₀(burst frequency)", markersize = 5, markeralpha = 0.9,
         legendfontsize = 9, fg_legend = :transparent, size = (500,400), dpi = 300)
#l = [-2.5:0.1:4.0;]
savefig(p, "burst_freq_bursty_size_3.png")



x = vcat(s_m[1:2]...)
p = scatter(burst_size[sortperm(x)], burst_freq[sortperm(x)], zcolor = log.(10,sort(x)), color = :viridis, ylabel = "log₁₀(burst frequency)", 
        xlims = (-2.5,4.5),ylims = (-3.0,3.5), xlabel = "log₁₀(burst size)", label = nothing, colorbar_title = "\n log₁₀(mean)", markersize = 5, markeralpha = 0.9,
         legendfontsize = 9, fg_legend = :transparent, right_margin = 5mm, size = (500,400), dpi = 300)
l = [-2.5:0.1:4.0;]
p = plot!(l, -l, label = false, linewidth = 4, linealpha = 0.7, linestyle = :dash, color = :lightblue3)
savefig(p, "data/paper_figures/figure_3/model_based/burst_freq_bursty_size_mean")


decay_rate = vcat(map_sets[1][:,4],map_sets[2][:,4])
x = vcat(s_m[1:2]...)[sortperm(decay_rate)]
y = vcat(s_tx_deg_noise[1:2]...)[sortperm(decay_rate)]

x1 = vcat([minimum(vcat(s_m[1:2]...)) - 2.0],sort(vcat(s_m[1:2]...)))
p = scatter(log.(10,x), log.(10,y), zcolor = sort(10 .^decay_rate), color = :viridis, ylabel = "log₁₀(transcription & decay noise)", 
        xlims = (0.0,3.5),ylims = (-3,1),xlabel = "log₁₀(mean)", label = nothing, colorbar_title = "\ndecay rate", markersize = 5, markeralpha = 0.9,
         legendfontsize = 9, fg_legend = :transparent, right_margin = 4mm, size = (500,400), dpi = 300)
p = plot!(log.(10,x1), log.(10, 1 ./ x1), linewidth = 4, label = false,linealpha = 0.6,linestyle = :dash, color = :lightblue3)

savefig(p, "data/paper_figures/figure_3/model_based/decay_rate_tx_deg_noise_1.svg")

y = vcat(s_bursty_noise[1:2]...)[sortperm(decay_rate)]
p = scatter(log.(10,x), log.(10,y), zcolor = sort(10 .^decay_rate), color = :viridis, ylabel = "log₁₀(bursty promoter noise)", 
        xlims = (0.0,3.5),ylims = (-10,1),xlabel = "log₁₀(mean)", label = nothing, colorbar_title = "\ndecay rate", markersize = 5, markeralpha = 0.9,
         legendfontsize = 9, fg_legend = :transparent, right_margin = 4mm, size = (500,400), dpi = 300)
#p = plot!(log.(10,x1), log.(10, 1 ./ x1), linewidth = 4, label = false,linealpha = 0.6,linestyle = :dash, color = :lightblue3)

savefig(p, "data/paper_figures/figure_3/model_based/decay_rate_bursty_noise.svg")


y = vcat(s_ccd_noise[1:2]...)[sortperm(decay_rate)]
p = scatter(log.(10,x), log.(10,y), zcolor = sort(10 .^decay_rate), color = :viridis, ylabel = "log₁₀(cell cycle variation)", 
        ylims = (-3.0,-1.4), xlims = (0.0,3.5), xlabel = "log₁₀(mean)", label = nothing, colorbar_title = "\ndecay rate", markersize = 5, markeralpha = 0.9,
         legendfontsize = 9, fg_legend = :transparent, right_margin = 4mm, size = (500,400), dpi = 300)
#p = plot!(log.(10,x1), log.(10, 1 ./ x1), linewidth = 4, label = false,linealpha = 0.6,linestyle = :dash, color = :lightblue3)

savefig(p, "data/paper_figures/figure_3/model_based/decay_rate_ccd_noise.svg")

x = vcat(s_m[1:2]...)
y = vcat(s_bursty_noise[1:2]...)[sortperm(x)]
p = scatter(10 .^decay_rate[sortperm(x)], log.(10,y), zcolor = log.(10,sort(x)), color = :viridis, ylabel = "log₁₀(cell-cycle dependence noise)", 
        xlabel = "decay rate", label = nothing, colorbar_title = "log₁₀(mean)", markersize = 5, markeralpha = 0.9,
         legendfontsize = 9, fg_legend = :transparent, right_margin = 4mm, size = (500,400), dpi = 300)

savefig(p, "data/paper_figures/figure_3/model_based/decay_rate_tx_deg_noise_1.svg")

x = vcat(s_m[1:2]...)
y = vcat(s_ccd_noise[1:2]...)[sortperm(x)]
p = scatter(10 .^decay_rate[sortperm(x)], log.(10,y), zcolor = log.(10,sort(x)), color = :viridis, ylabel = "log₁₀(cell-cycle dependence noise)", 
        xlabel = "decay rate", label = nothing, colorbar_title = "log₁₀(mean)", markersize = 5, markeralpha = 0.9,
         legendfontsize = 9, fg_legend = :transparent, right_margin = 4mm, size = (500,400), dpi = 300)

savefig(p, "data/paper_figures/figure_3/model_based/decay_rate_tx_deg_noise_1.svg")




x = vcat(s_m[1:2]...)
y = vcat(s_tx_deg_noise[1:2]...)[sortperm(x)]

p = scatter(10 .^decay_rate[sortperm(x)], y,zcolor = log.(10,sort(x)), color = :viridis, ylabel = "log₁₀(transcription & decay noise)", 
   xlabel = "log₁₀(decay rate)", label = nothing, colorbar_title = "log₁₀(mean)", markersize = 5, markeralpha = 0.9,
         legendfontsize = 9, fg_legend = :transparent, size = (500,400), dpi = 300)
p = plot!([-3:0.1:0.5], [-3:0.1:0.5], label = false, linewidth = 4, linealpha = 0.7, linestyle = :dash, color = :lightblue3)

savefig(p, "data/paper_figures/figure_3/model_based/decay_rate_tx_deg_noise.svg")



col = [:skyblue4, :lightcyan4]
p = boxplot(ccd_noise[1], linewidth = 1, color = col[1], label = false, 
ylabelfontsize = 10,xtickfontsize = 11, background_color_legend = :transparent,fg_legend = :transparent, size = (500,300), bottom_margin=2mm,
xticks = ([1,2,3,4],["scaling \ngenes","non-scaling \ngenes","scaling \ngenes","non-scaling \ngenes"]), ylabel = "cell cycle variation", outliers = false)
p = boxplot!(ccd_noise[2], linewidth = 1, color = col[2], label = false, outliers = false)
p = boxplot!(s_ccd_noise[1], linewidth = 1, color = col[1], label = false, 
ylabelfontsize = 10,xtickfontsize = 11, outliers = false)
p = boxplot!(s_ccd_noise[2], linewidth = 1, color = col[2], label = false, outliers = false)

savefig(p, "data/paper_figures/figure_4/observed_ccd_noise_2.svg")





#cell cycle dependent noise
cycle_label = ["G1", "G1/S", "S", "S/G2", "G2/M"]; 
col = cgrad(:Blues,categorical=true,5);
decay_rate = vcat(map_sets[5][:,vary_maps[5][4]])
p = scatter(decay_rate[:,1],log.(10,s_noise_ccp[5][:,1]), color = col[1], ylabel = "log₁₀(noise)", 
    xlims = (-3,2),ylims = (-3,2),xlabel = "log₁₀(decay rate)", label = cycle_label[1], markersize = 5, markeralpha = 0.9,
 legendfontsize = 9, legend = :topleft,background_color_legend = :transparent, fg_legend = :transparent, size = (500,400), dpi = 300);
for j in 2:5
    p = scatter!(decay_rate[:,j],log.(10,s_noise_ccp[5][:,j]),color = col[j], label = cycle_label[j],legendfontsize = 9, markersize = 5, markeralpha = 0.9,)
end


cycle_label = ["G1", "G1/S", "S", "S/G2", "G2/M"]; 
col = cgrad(:Blues,categorical=true,5);
decay_rate = vcat(map_sets[5][:,vary_maps[5][4]])
p = scatter(log.(10,s_mean_ccp[5][:,1]),log.(10,s_noise_ccp[5][:,1]), color = col[1], ylabel = "log₁₀(noise)", 
    xlims = (-3,2),ylims = (-3,2),xlabel = "log₁₀(decay rate)", label = cycle_label[1], markerstrokewidth = 0.1, markersize = 3 .* decay_rate[:,1], markeralpha = 0.9,
 legendfontsize = 9, legend = :topleft,background_color_legend = :transparent, fg_legend = :transparent, size = (500,400), dpi = 300);
for j in 2:5
    p = scatter!(log.(10,s_mean_ccp[5][:,j]),log.(10,s_noise_ccp[5][:,j]), color = col[j], label = cycle_label[j],legendfontsize = 9,
             markersize = 3 .* decay_rate[:,j], markerstrokewidth = 0.1,markeralpha = 0.9)
end


