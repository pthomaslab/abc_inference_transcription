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

function get_subset(data::Matrix{Float64},idx::Vector{Int64})
    sub = Matrix{Float64}(undef,(length(idx)*2,size(data)[2]))
    for (i,k) in enumerate(idx)
        sub[2*i-1:2*i,:] = data[2*k-1:2*k,:]
    end
    return sub
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

function load_s_data(path::String, model_idx::Int64, n_sims::Int64, data_idx::Vector{Int64})
    n::Int64 = 10^6;
    model_name::String = ["const","const_const","kon","alpha","gamma"][model_idx]
    data_name::Vector{String} = ["s_pulse_","s_chase_","s_ratios_","s_mean_corr_","s_corr_mean_"][data_idx]
    s_data = Vector{Matrix{Float64}}(undef,length(data_idx))
    for (j,idx) in enumerate(data_idx)
        if idx<=2
            s_data[j] = Matrix{Float64}(undef,(2*n_sims*n,5))
            k::Int64 = 2
        else
            s_data[j] = Matrix{Float64}(undef,(n_sims*n,11))
            k = 1
        end
        for i in 1:n_sims
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


function filter_quantiles(data::AbstractArray{Float64},q::Float64)
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


function get_burst_kinetics(sets::Matrix{Float64}, m::Int64)
    #kinetics_names::Vector{String} = ["kon","burst size","decay rate"]
    vary_flag::Vector{Int64} = [[0,0,0,0],[0,0,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,1]][m];
    vary_map::Vector{Any} = get_vary_map(vary_flag,5)
    #kinetics = Matrix{Float64}(undef,(size(map_sets)[1],sum(length.(vary_map))))
    kinetics::Matrix{Float64} = hcat(sets[:,vary_map[1]],sets[:,vary_map[3]],sets[:,vary_map[3]] .- sets[:,vary_map[2]],sets[:,vary_map[4]])
    #return hcat(kinetics,sets[:,end])
    return kinetics
end

function plot_const_posteriors(θ::Matrix{Float64},m::Int64)
    rate_names::Vector{String} = ["burst frequency","koff","synthesis rate","burst size","decay rate","λ"]
    lims = [(-3.0,3.0),(-3.0,3.0),(-3.0,3.0),(-6.0,6.0),(-3.0,2.0),(-0.7,0.0)]
    model_label = ["constant scaling", "constant non-scaling"][m];
    col = cgrad(:acton, 30, rev=true, categorical = true)[[25,5]][m]
    p1 = [histogram(θ[:,i], normalize = :probability, tickfontsize = 5, xlims = lims[i], xticks = round.(collect(range(lims[i][1],lims[i][end],5)),digits=2), 
        label = false, title = rate_names[i], titlefontsize = 11, color = col) for i in 1:lastindex(rate_names)]
    #=p = [scatter(θ[sortperm(errors,rev=true),i], θ[sortperm(errors,rev=true),j], , color = cgrad(:viridis,rev=true), xlims = lims[i], ylims = lims[j],
        xticks = collect(range(lims[i][1],lims[i][end],5)),yticks = collect(range(lims[j][1],lims[j][end],5)),
        xlabel = rate_names[i], ylabel = rate_names[j], markersize = 3, tickfontsize = 5, xguidefontsize = 7,
        yguidefontsize = 7, label = false, legend = false) for i in 1:size(θ)[2] for j in 1:size(θ)[2]]
    =#
    p = [scatter(θ[:,i], θ[:,j], color = col, xlims = lims[i], ylims = lims[j],
        xticks = round.(collect(range(lims[i][1],lims[i][end],5)),digits=2),yticks = round.(collect(range(lims[j][1],lims[j][end],5)),digits=2),
        xlabel = rate_names[i], ylabel = rate_names[j], markersize = 3,markerstrokewidth = 0.0, tickfontsize = 5, xguidefontsize = 8,
        yguidefontsize = 8, label = false, legend = false) for i in 1:size(θ)[2] for j in 1:size(θ)[2]]
    diag_idx = [7*i+1 for i in 0:5]
    for (i,k) in enumerate(diag_idx)
        p[k] = p1[i]
    end
    pairplot = plot(p..., layout = (size(θ)[2],size(θ)[2]),plot_title = model_label*" gene MAPs", size = (900,800));
    return pairplot
    #savefig(scat, "lambda_scatter_"*rate*".pdf")
end



function plot_vary_posteriors(θ::Matrix{Float64},m::Int64)
    rate_names::Vector{String} = ["burst freq","koff","burst size","decay rate","λ"]
    model_label::String = ["burst frequency","burst size","decay rate"][m-2]
    lims = [(-3.0,3.0),(-3.0,3.0),(-6.0,6.0),(-3.0,2.0),(-0.7,0.0)]
    n_steps::Int64 = 5
    vary_flag::Vector{Int64} = [[1,0,0,0],[0,0,1,0],[0,0,0,1]][m-2]
    vary_map::Vector{Any} = get_vary_map(vary_flag,n_steps)
    vary_idx::Int64 = findall(x->length(x)>1,vary_map)[1]
    names = Vector{String}(undef,n_steps+length(vary_map)-1)
    priors = Vector{Tuple{Float64,Float64}}(undef,n_steps+length(vary_map)-1)
    for i in 1:lastindex(vary_map)
        if i != vary_idx
            names[vary_map[i]] = rate_names[i]
            priors[vary_map[i]] = round.(lims[i];digits=2)
        else
            names[vary_map[i]] = [rate_names[i]*" ($j)" for j in 1:n_steps]
            priors[vary_map[i]] = [round.(lims[i];digits=2) for _ in 1:n_steps]
        end
    end
    names = push!(names,rate_names[end])
    priors = push!(priors,lims[end])
    p1 = [density(θ[:,i], linewidth = 5, tickfontsize = 5, xlims = priors[i], xticks = round.(collect(range(priors[i][1],priors[i][end],4)),digits=2), 
    label = false, title = names[i], titlefontsize = 9, colour = :steelblue) for i in 1:size(θ)[2]]
    #=p = [scatter(θ[sortperm(errors,rev=true),i], θ[sortperm(errors,rev=true),j], zcolor = sort(errors,rev=true), color = cgrad(:viridis,rev=true), xlims = priors[i], ylims = priors[j],
        xticks = collect(range(priors[i][1],priors[i][end],5)),yticks = collect(range(priors[j][1],priors[j][end],5)),
        xlabel = names[i], ylabel = names[j], markersize = 3, tickfontsize = 5, xguidefontsize = 7,
        yguidefontsize = 7, label = false, legend = false) for i in 1:size(θ)[2] for j in 1:size(θ)[2]]
    =#
    p = [scatter(θ[:,i], θ[:,j], color = :steelblue, xlims = priors[i], ylims = priors[j],
        xticks = round.(collect(range(priors[i][1],priors[i][end],4)),digits=2),yticks = round.(collect(range(priors[j][1],priors[j][end],4)),digits=2),
        xlabel = names[i], ylabel = names[j], markersize = 3, markerstrokewidth = 0.0, tickfontsize = 5, xguidefontsize = 7,
        yguidefontsize = 7, label = false, legend = false) for i in 1:size(θ)[2] for j in 1:size(θ)[2]]
    diag_idx = [10*i+1 for i in 0:8]
    for (i,k) in enumerate(diag_idx)
        p[k] = p1[i]
    end
    pairplot = plot(p..., layout = (size(θ)[2],size(θ)[2]), plot_title = model_label*" gene MAPs", size = (1600,1400))
    return pairplot
end



function plot_fits(data::Matrix{Float64}, data_se::Matrix{Float64}, s_mean::Matrix{Float64}, s_ff::Matrix{Float64}, gene_name::String, m::Int64)
    data_col = :mediumpurple4
    col = [:royalblue4, :royalblue4, :steelblue, :seagreen4, :darkgoldenrod1][m]
    model_label = ["constant rates (scaling)", "constant rates (non-scaling)", 
                "varying burst frequency", "varying burst size", "varying decay rate"][m]
    agevec = [0.1:0.2:0.9;]
    local mean_plot::Plots.Plot{Plots.GRBackend}
    local ff_plot::Plots.Plot{Plots.GRBackend}
    mean_plot = plot(agevec,
        hcat(data[:,1] - data_se[:,1],data[:,1] + data_se[:,1]),
        fillrange = data[:,1],
        fillalpha = 0.15,
        fillcolor = [data_col data_col],
        linealpha = [0.15 0.15],
        linecolor = [data_col data_col],
        labels = [nothing nothing],
        shape = [:none :none])
    mean_plot = plot!(agevec,
            data[:,1],
            title = gene_name,
            xlabel = "cell cycle progression",
            ylabel = "mean transcript levels",
            label = "data",
            shape = :circle,
            legend = :topleft,
            legendfontsize = 11,
            fg_legend = :transparent,
            linecolor = data_col,
            linewidth = 5.0,
            markercolor = data_col)
    mean_plot = plot!(agevec,
            s_mean[1,:],
            linecolor = col,
            label = model_label,
            linewidth = 2.0)
    mean_plot = plot!(agevec,
            s_mean[2:end,:]',
            linecolor = col,
            label = false,
            linewidth = 2.0)
    ff_plot = plot(agevec,
        hcat(data[:,2] - data_se[:,2],data[:,2] + data_se[:,2]),
        fillrange = data[:,2],
        fillalpha = 0.15,
        fillcolor = [data_col data_col],
        linealpha = [0.15 0.15],
        linecolor = [data_col data_col],
        labels = [nothing nothing],
        shape = [:none :none])
    ff_plot = plot!(agevec,
            data[:,2],
            title = gene_name,
            xlabel = "cell cycle progression",
            ylabel = "Fano factor of transcript levels",
            label = "data",
            shape = :circle,
            legend = :topleft,
            legendfontsize = 11,
            fg_legend = :transparent,
            linecolor = data_col,
            linewidth = 5.0,
            markercolor = data_col)
    ff_plot = plot!(agevec,
            s_ff[1,:],
            linecolor = col,
            linalpha = 0.1,
            label = model_label,
            linewidth = 5.0)
    ff_plot = plot!(agevec,
            s_ff[2:end,:]',
            linecolor = col,
            linalpha = 0.1,
            label = false,
            linewidth = 5.0)   
    return mean_plot,ff_plot
end



#################################################################################################
#################################################################################################
#################################################################################################

all_genes = Int64.(readdlm("data/selected_genes.txt")[:,1])
gene_id = readdlm("data/all_gene_ids.txt")[:,1]
all_gene_names = readdlm("data/all_gene_names.txt")[:,1]
genes = all_genes[Int64.(readdlm("data/selected_genes_main_subset.txt")[:,1])]
ids = gene_id[genes];
gene_names = replace(all_gene_names[genes], NaN => "NaN")


pulse_data,pulse_se,chase_data,chase_se,ratio_data,ratio_se,
mean_corr_data,mean_corr_se,corr_mean_data,corr_mean_se = load_summary_stats("data/selected_summary_stats/", ".txt");

pulse_mean = get_mean_subset(pulse_data)
pulse_mean_se = get_mean_subset(pulse_se)
pulse_ff = get_ff_subset(pulse_data)
pulse_ff_se = get_ff_subset(pulse_se)
chase_mean = get_mean_subset(chase_data)
chase_mean_se = get_mean_subset(chase_se)
chase_ff = get_ff_subset(chase_data)
chase_ff_se = get_ff_subset(chase_se)



uu_data, us_data, lu_data, ls_data, theta, rfp, gfp, experiment, gene_id = read_all_data("data/",".csv");
total_data = uu_data + us_data + lu_data + ls_data;
ncells,ngenes = size(ls_data);
pulse_idx = findall(x->x<=6, experiment)

chase_idx = findall(x->x>=7, experiment)
n_clusters = 5
age, cluster_idx, mean_age, τ_ = age_clusters(theta, n_clusters, "equidistant")
age_distribution = length.(cluster_idx)/ncells;

u_data = uu_data[:,genes] + us_data[:,genes];
l_data = lu_data[:,genes] + ls_data[:,genes];

t_data = u_data + l_data;

cells_per_id = [findall(x->x==e,experiment) for e in sort(unique(experiment))[2:end-1]]
cells_per_age = [findall(x->x==τ,age) for τ in sort(unique(age))]
cells_age_id = [[intersect(cells_age,cells_id) for cells_age in cells_per_age] for cells_id in cells_per_id]

mean_data = hcat([mean(t_data[cells,:], dims = 1)[1,:] for cells in cells_per_age]...)
var_data = hcat([var(t_data[cells,:], dims = 1)[1,:] for cells in cells_per_age]...);

#################### load gene-specific kinetics posterior distributions for each gene model ####################
model_names = ["const","const_const","kon","alpha","gamma"];     
posterior = Vector{Vector{Vector{Int64}}}(undef,length(model_names))
for (i,name) in enumerate(model_names)
    file = open("data/posteriors/particles_"*name*".txt", "r")
    posterior[i] = Vector{Vector{Int64}}(undef,length(genes))
    for (j,line) in enumerate(eachline(file))
        posterior[i][j] = parse.(Int64,split(line,"\t")) ## or whatever else you want to do with the line.
    end
    close(file)
end  

n_particles_model = get_n_particles(posterior)
n_particles_mat = hcat(n_particles_model...)
n_particles = sum(n_particles_mat, dims=2)[:,1]

   
errs = Vector{Vector{Vector{Float64}}}(undef,length(model_names));
for (i,name) in enumerate(model_names)
    file = open("data/posteriors/$ε/errors_"*name*".txt", "r")
    errs[i] = Vector{Vector{Float64}}(undef,length(genes))
    for (j,line) in enumerate(eachline(file))
        errs[i][j] = parse.(Float64,split(line,"\t")) ## or whatever else you want to do with the line.
    end
    close(file)
end  

n_groups_sims = 1
sets = [load_param_sets("data/large_scale_simulations/",m,n_groups_sims) for m in 1:5]


# for each gene, get  the best performing particle per model
#=
maps = Matrix{Int64}(undef,(length(genes),length(model_names)))
for (j,post) in enumerate(posterior)
    maps[:,j] = [post[g][1] for g in 1:length(genes_1)]
end
writedlm("data/posteriors/map_gene_per_model.txt",maps)
=#

maps = Int64.(readdlm("data/posteriors/map_gene_per_model.txt")[:,1:5])


sel_genes = [Int64.(readdlm("data/model_selection/"*mn*"_genes.txt")[:,1]) for mn in model_names]

#sel_gene_ids = [ids[sel] for sel in sel_genes]

map_sets = [sets[m][maps[sel_genes[m],m],:] for m in 1:lastindex(model_names)];

ms_df = CSV.read("data/model_selection/model_selection_results.txt", DataFrame);
ms = Int64.(Matrix(ms_df[:,[1,3]]));


#################   Compare constant scaling vs non-scaling genes   ################
m = 1
map_kinetics = get_burst_kinetics(map_sets[m],m)

#p = plot_const_posteriors(map_kinetics,m)

col = [:skyblue4, :lightcyan4]


mat_1 = mean_data[sel_genes[1],:] ./ mean_data[sel_genes[1],1]
mat_2 = mean_data[sel_genes[2],:] ./ mean_data[sel_genes[2],1]


q = 0.99

h0 = density(filter_quantiles(mat_1[:,end],q), linewidth = 5, color = col[1], title = "cell growth ratio of constant genes",xlabel = "cell growth ratio", ylabel = "probability density", 
            label = "scaling genes", legendfontsize = 10, fg_legend = :transparent, size = (400,300))
h0 = density!(filter_quantiles(mat_2[:,end],q),linewidth = 5, color = col[2],label = "non-scaling genes")

savefig(h0,"data/paper_figures/figure_4/cell_growth_const_genes.pdf")


h0 = boxplot(mat_1[:,end], linewidth = 1, color = col[1], label = false, legendfontsize = 11, xtickfontsize = 11,
background_color_legend = :transparent,legend = :topright,fg_legend = :transparent, size = (400,300), xticks = ([1,2],["scaling genes","non-scaling genes"]), ylabel = "mean expression fold change", outliers = false)
h0 = boxplot!(mat_2[:,end], linewidth = 1, color = col[2], label = false, outliers = false)

savefig(h0,"data/paper_figures/figure_4/fold_change_const_genes.svg")


mat_1 = mean_data[sel_genes[1],:]
mat_2 = mean_data[sel_genes[2],:] 

der_1 = Vector{Float64}(undef,size(mat_1)[1]);
der_2 = Vector{Float64}(undef,size(mat_2)[1]);
for i in 1:size(mat_1)[1]
    der_1[i] = mean(diff(mat_1[i,:]))
end

for i in 1:size(mat_2)[1]
    der_2[i] = mean(diff(mat_2[i,:]))
end

curv_1 = Vector{Float64}(undef,size(mat_1)[1]);
curv_2 = Vector{Float64}(undef,size(mat_2)[1]);
for i in 1:size(mat_1)[1]
    curv_1[i] = mean(diff(diff(mat_1[i,:])))
end

for i in 1:size(mat_2)[1]
    curv_2[i] = mean(diff(diff(mat_2[i,:])))
end



p = boxplot(der_1, linewidth = 1, color = col[1], label = false, legendfontsize = 11,
background_color_legend = :transparent,legend = :topright,fg_legend = :transparent, size = (400,300), xtickfontsize = 11,xticks = ([1,2],["scaling genes","non-scaling genes"]), ylabel = "mean rate of change",outliers=false)
p = boxplot!(der_2, linewidth = 1, color = col[2], label = false,outliers=false)
savefig(p, "data/paper_figures/figure_4/mean_deriv.svg")


p = boxplot(curv_1, linewidth = 1, color = col[1], label = false, legendfontsize = 11,
background_color_legend = :transparent,legend = :topright,fg_legend = :transparent, size = (400,300), xtickfontsize = 11,xticks = ([1,2],["scaling genes","non-scaling genes"]), ylabel = "mean second derivative",outliers=false)
p = boxplot!(curv_2, linewidth = 1, color = col[2], label = false,outliers=false)
savefig(p, "data/paper_figures/figure_4/mean_2nd_deriv.pdf")


h0 = density(ccd_noise[1], color = col[1], linewidth = 5, xlabel = "CV²", ylabel = "probability density",label = "scaling genes",fg_legend = :transparent,
                     xlims = (0.0,0.016),legend = :topleft, legendfontsize = 10, title = "cell-cycle dependence noise of constant genes      ", size = (500,300))
h0 = density!(ccd_noise[2],color = col[2], linewidth = 5, label = "non-scaling genes")
savefig(h0, "data/paper_figures/figure_4/cc_dep_noise.pdf")



p = violin(map_sets[1][:,1], linewidth = 5, color = col[1], label = false, title = "burst frequency of constant genes    ", legendfontsize = 11,
background_color_legend = :transparent,legend = :topright,fg_legend = :transparent, size = (400,300), ylabel = "probability density", xlabel = "log₁₀(burst frequency)")
p = violin!(map_sets[2][:,1], linewidth = 5, linealpha = 0.5, color = col[2], label = false)
savefig(p, "data/paper_figures/figure_4/burst_freq.svg")



p = density(map_sets[1][:,3] .- map_sets[1][:,2], linewidth = 5, color = col[1], label = false, title = "burst size of constant genes", legendfontsize = 11,
background_color_legend = :transparent,legend = :topright,fg_legend = :transparent, size = (400,300), ylabel = "probability density", xlabel = "log₁₀(burst size)")
p = density!(map_sets[2][:,3] .- map_sets[2][:,2], linewidth = 5, linealpha = 0.5, color = col[2], label = false)
savefig(p, "data/paper_figures/figure_4/burst_size.svg")


p = density(map_sets[1][:,4], linewidth = 5, color = col[1], label = "scaling genes", title = "decay rate of constant genes", legendfontsize = 11,
background_color_legend = :transparent,legend = :topleft,fg_legend = :transparent, size = (400,300), ylabel = "probability density", xlabel = "log₁₀(decay rate)")
p = density!(map_sets[2][:,4], linewidth = 5, linealpha = 0.5, color = col[2], label = "non-scaling genes")

savefig(p, "data/paper_figures/figure_4/decay_rate.svg")

p = boxplot(10 .^map_sets[1][:,4], linewidth = 1, color = col[1], label = false, legendfontsize = 11,
background_color_legend = :transparent,legend = :topright,fg_legend = :transparent, size = (400,300), xtickfontsize = 11,xticks = ([1,2],["scaling genes","non-scaling genes"]), ylabel = "decay rate")
p = boxplot!(10 .^map_sets[2][:,4], linewidth = 1, color = col[2], label = false)

savefig(p, "data/paper_figures/figure_4/decay_rate_1.svg")

p = boxplot(10 .^map_sets[1][:,4], linewidth = 1, color = col[1], label = false, legendfontsize = 11,ylabel = "decay rate",xticks = ([1,2],["scaling genes","non-scaling genes"]),
background_color_legend = :transparent,legend = :topright,fg_legend = :transparent, size = (400,300), xtickfontsize = 11, outliers = false)
p = boxplot!(10 .^map_sets[2][:,4], linewidth = 1,  color = col[2], label = false, outliers = false)

savefig(p, "data/paper_figures/figure_4/decay_rate_2.svg")

p = boxplot(map_sets[1][:,1], linewidth = 1, color = col[1], label = false, title = "burst frequency of constant genes   ", legendfontsize = 11, background_color_legend = :transparent,
    legend = :topright,fg_legend = :transparent, size = (400,300), xticks = ([1,2],["scaling genes","non-scaling genes"]), ylabel = "log₁₀(burst frequency)", outliers = false)
p = boxplot!(map_sets[2][:,1], linewidth = 1, linealpha = 0.5, color = col[2], label = false, outliers = false)
savefig(p, "data/paper_figures/figure_4/burst_frequency.png")

p = boxplot(map_sets[1][:,3] .- map_sets[1][:,2], linewidth = 1, color = col[1], label = false, title = "burst size of constant genes", legendfontsize = 11, background_color_legend = :transparent,
legend = :topright,fg_legend = :transparent, size = (400,300), xticks = ([1,2],["scaling genes","non-scaling genes"]), ylabel = "log₁₀(burst size)", outliers = false)
p = boxplot!(map_sets[2][:,3] .- map_sets[2][:,2], linewidth = 1, linealpha = 0.5, color = col[2], label = false, outliers = false)
savefig(p, "data/paper_figures/figure_4/burst_size.png")



###########################     distributions of posterior estimates across genes     ##########################

vary_flags::Vector{Vector{Int64}} = [[0,0,0,0],[0,0,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,1]]
vary_maps::Vector{Vector{Any}} = [get_vary_map(vary_flag,5) for vary_flag in vary_flags];

col = [:skyblue4, :seagreen4, :gold]

p = violin(vcat(map_sets[1][:,vary_maps[1][4]],map_sets[2][:,vary_maps[2][4]]), xticks = ([1:4;],["constant \ngenes","burst frequency \ngenes","burst size \ngenes","decay rate \ngenes"]), 
    linewidth = 1, linealpha = 0.7,color = :grey, label = false, outliers = false);
p = violin!(map_sets[3][:,vary_maps[3][4]], linewidth = 1, color = col[1], label = false, legendfontsize = 9, 
background_color_legend = :transparent,legend = :topleft,fg_legend = :transparent, size = (400,300), ylabel = "log₁₀(decay rate)", outliers = false)
p = violin!(map_sets[4][:,vary_maps[4][4]], linewidth = 1,linealpha = 0.7, color = col[2], label = false, outliers = false)

p = violin!(vec(map_sets[5][:,vary_maps[5][4]]), linewidth = 1, linealpha = 0.7,color = col[3], label = false, outliers = false)
savefig(p, "data/paper_figures/figure_2/decay_rates.svg")
savefig(p, "data/paper_figures/figure_2/decay_rates.png")


p = violin(vcat(map_sets[1][:,vary_maps[1][1]],map_sets[2][:,vary_maps[2][1]]), xticks = ([1:4;],["constant \ngenes","burst frequency \ngenes","burst size \ngenes","decay rate \ngenes"]), 
        linewidth = 1, linealpha = 0.7,color = :grey, label = false, outliers = false);
p = violin!(vec(map_sets[3][:,vary_maps[3][1]]), linewidth = 1, color = col[1], label = false, legendfontsize = 9, 
background_color_legend = :transparent,legend = :topleft,fg_legend = :transparent, size = (400,300), ylabel = "log₁₀(burst frequency)", outliers = false)
p = violin!(map_sets[4][:,vary_maps[4][1]] , linewidth = 1,linealpha = 0.7, color = col[2], label = false, outliers = false)
p = violin!(map_sets[5][:,vary_maps[5][1]], linewidth = 1, linealpha = 0.7,color = col[3], label = false, outliers = false)

savefig(p, "data/paper_figures/figure_2/burst_freq.svg")
savefig(p, "data/paper_figures/figure_2/burst_freq.png")


p = violin(vcat(map_sets[1][:,vary_maps[1][3]] .- map_sets[1][:,vary_maps[1][2]],map_sets[2][:,vary_maps[2][3]] .- map_sets[2][:,vary_maps[2][2]]), xticks = ([1:4;],["constant \ngenes","burst frequency \ngenes","burst size \ngenes","decay rate \ngenes"]), 
        linewidth = 1, linealpha = 0.7,color = :grey, label = false, outliers = false);
p = violin!(map_sets[3][:,vary_maps[3][3]] .- map_sets[3][:,vary_maps[3][2]], linewidth = 1, color = col[1], label = false, legendfontsize = 9, 
background_color_legend = :transparent,legend = :topleft,fg_legend = :transparent, size = (400,300), ylabel = "log₁₀(burst size)", outliers = false)
p = violin!(vec(map_sets[4][:,vary_maps[4][3]] .- map_sets[4][:,vary_maps[4][2]]), linewidth = 1,linealpha = 0.7, color = col[2], label = false, outliers = false)
p = violin!(map_sets[5][:,vary_maps[5][3]] .- map_sets[5][:,vary_maps[5][2]], linewidth = 1, linealpha = 0.7,color = col[3], label = false, outliers = false)

savefig(p, "data/paper_figures/figure_2/burst_size.svg")
savefig(p, "data/paper_figures/figure_2/burst_size.png")


x = mean(t_data,dims=1)[1,:]
p = violin(vcat(x[sel_genes[1]],x[sel_genes[2]]), xticks = ([1:4;],["constant \ngenes","burst frequency \ngenes","burst size \ngenes","decay rate \ngenes"]), 
    linewidth = 1, linealpha = 0.7,color = :grey, label = false, ylims = (0.0,10), outliers = false);
p = violin!(x[sel_genes[3]], linewidth = 1, color = col[1], label = false, legendfontsize = 9, 
background_color_legend = :transparent,legend = :topleft,fg_legend = :transparent, size = (400,300), ylabel = "mean expression", outliers = false)
p = violin!(x[sel_genes[4]], linewidth = 1,linealpha = 0.7, color = col[2], label = false, outliers = false)
p = violin!(x[sel_genes[5]], linewidth = 1, linealpha = 0.7,color = col[3], label = false, outliers = false)
savefig(p, "data/paper_figures/figure_2/mean_expression.svg")
savefig(p, "data/paper_figures/figure_2/mean_expression.png")




###############################################################################################
############################### mRNA half-life analysis #######################################
###############################################################################################

vary_flags::Vector{Vector{Int64}} = [[0,0,0,0],[0,0,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,1]]
vary_maps::Vector{Vector{Any}} = [get_vary_map(vary_flag,5) for vary_flag in vary_flags];

cycle = 20.0
dilution_rate = 0.0; # log(2) / cycle;

half_life = [log(2) ./ (10 .^m_s[:,vary_maps[i][4]]) for (i,m_s) in enumerate(map_sets[1:4])];
half_life = push!(half_life, log(2) ./ ((10 .^mean(map_sets[5][:,vary_maps[5][4]],dims=2)[:,1]) .+ dilution_rate))

col = [:skyblue4, :lightcyan4]
p = boxplot(half_life[1], linewidth = 1, color = col[1], label = false, legendfontsize = 11, ylims = (0.0,20.0),
background_color_legend = :transparent,legend = :topright,fg_legend = :transparent, size = (400,300), xtickfontsize = 11,xticks = ([1,2],["scaling genes","non-scaling genes"]), ylabel = "mRNA half-life (h)", outliers = false)
p = boxplot!(half_life[2], linewidth = 1, color = col[2], label = false, dpi = 300, outliers = false)
savefig(p, "data/paper_figures/figure_4/half_life.svg")
savefig(p, "data/paper_figures/figure_4/half_life.png")

savefig(p, "data/paper_figures/figure_4/half_life_no_deg.svg")
savefig(p, "data/paper_figures/figure_4/half_life_deg.png")


col = [:skyblue4, :seagreen4, :gold]
p = boxplot(half_life[3], linewidth = 1, color = col[1], label = false, legendfontsize = 11, ylims = (0.0,40),
background_color_legend = :transparent,legend = :topright,fg_legend = :transparent, size = (400,300), xtickfontsize = 11, bottom_margin = 2mm,
xticks = ([1,2,3],["burst frequency \ngenes","burst size \ngenes", "decay rate \ngenes"]), ylabel = "mRNA half-life (h)", dpi = 300, outliers = false)
p = boxplot!(half_life[4], linewidth = 1, color = col[2], label = false, outliers = false)
p = boxplot!(half_life[5], linewidth = 1, color = col[3], label = false, outliers = false)
savefig(p, "data/paper_figures/figure_5/half_life.svg")
savefig(p, "data/paper_figures/figure_5/half_life.png")

#compare with scSLAM-seq
scslam = CSV.read("scSLAMseq_hf.csv",DataFrame)
gene_vec = vcat(sel_genes[1:2]...)
idx = [];
idx_1 = [];
for (i,g) in enumerate(uppercase.(gene_names[gene_vec]))
    x = findall(x->x==g,uppercase.(scslam[:,1]))
    if x != []
        push!(idx,i)
        push!(idx_1,x[1])
    end
end


hl = vcat(half_life[1:2]...)[idx]
hl_1 = scslam[idx_1,2]
    
p = scatter(hl_1, hl, size = (400,300), label = false,ylabel = "mRNA half-life (scEU-seq)",xlabel = "mRNA half-life (scSLAM-seq)",
                ylims = (-3,100), markeralpha = 0.6, color = :skyblue4);
rho = StatsBase.corspearman(hl,hl_1)
p = plot!([0:maximum(hl_1);],[0:maximum(hl_1);],label = false,linewidth = 4, linestyle = :dash, top_margin = 2mm, color = :darkorange2, linealpha=0.75, dpi = 300)
p = annotate!(35,60,text("Spearman's ρ = $(round(rho,digits = 2))", 11, color = :darkorange2))
savefig(p, "data/paper_figures/figure_4/half_life_scSLAM.svg")
savefig(p, "data/paper_figures/figure_4/half_life_scSLAM.png")


gene_vec = vcat(sel_genes[3:5]...)
idx = [];
idx_1 = [];
for (i,g) in enumerate(uppercase.(gene_names[gene_vec]))
    x = findall(x->x==g,uppercase.(scslam[:,1]))
    if x != []
        push!(idx,i)
        push!(idx_1,x[1])
    end
end

hl = vcat(half_life[3:5]...)[idx]
hl_1 = scslam[idx_1,2]
  
p = scatter(hl_1, hl, size = (400,300), label = false,ylabel = "mRNA half-life (scEU-seq)",xlabel = "mRNA half-life (scSLAM-seq)",
                 ylims = (-3,40), markeralpha = 0.6, color = :skyblue4);
rho = StatsBase.corspearman(hl,hl_1)
p = plot!([0:maximum(hl_1);],[0:maximum(hl_1);],label = false,linewidth = 4, linestyle = :dash, top_margin = 2mm, color = :darkorange2, linealpha=0.75, dpi = 300)
p = annotate!(7,25,text("Spearman's ρ = $(round(rho,digits = 2))", 11, color = :darkorange2))
savefig(p, "data/paper_figures/figure_4/half_life_nc_scSLAM.svg")
savefig(p, "data/paper_figures/figure_4/half_life_nc_scSLAM.png")



#compare with sci-fate
scifate = CSV.read("scifate_hf.csv",DataFrame)
gene_vec = vcat(sel_genes[1:2]...)
idx = [];
idx_1 = [];
ids_1 = [split(id,".")[1] for id in scifate[:,1]]
for (i,g) in enumerate(ids[gene_vec])
    x = findall(x->x==g,ids_1)
    if x != []
        push!(idx,i)
        push!(idx_1,x[1])
    end
end



hl = vcat(half_life[1:2]...)[idx]
hl_1 = scifate[idx_1,3]
   
p = scatter(hl_1, hl, size = (400,300), label = false,ylabel = "mRNA half-life (scEU-seq)",xlabel = "mRNA half-life (sci-fate)",
               ylims = (-3,100), markeralpha = 0.6, color = :skyblue4);
rho = StatsBase.corspearman(hl,hl_1)
p = plot!([0:maximum(hl_1);],[0:maximum(hl_1);],label = false,linewidth = 4, linestyle = :dash, top_margin = 2mm, color = :darkorange2, linealpha=0.75, dpi = 300)
p = annotate!(88,40,text("Spearman's ρ = $(round(rho,digits = 2))", 11, color = :darkorange2))
#p = annotate!(88,47,text("$(length(idx)) common genes", 11))

savefig(p, "data/paper_figures/figure_4/half_life_scifate.svg")
savefig(p, "data/paper_figures/figure_4/half_life_scifate.png")

gene_vec = vcat(sel_genes[3:5]...)
idx = [];
idx_1 = [];
ids_1 = [split(id,".")[1] for id in scifate[:,1]]
for (i,g) in enumerate(ids[gene_vec])
    x = findall(x->x==g,ids_1)
    if x != []
        push!(idx,i)
        push!(idx_1,x[1])
    end
end

hl = vcat(half_life[3:5]...)[idx]
hl_1 = scifate[idx_1,3]
p = scatter(hl_1, hl, size = (400,300), label = false,ylabel = "mRNA half-life (scEU-seq)",xlabel = "mRNA half-life (sci-fate)",
                ylims = (-3,40), markeralpha = 0.6, color = :skyblue4);
rho = StatsBase.corspearman(hl,hl_1)
p = plot!([0:maximum(hl_1);],[0:maximum(hl_1);],label = false,linewidth = 4, linestyle = :dash, top_margin = 2mm, color = :darkorange2, linealpha=0.75, dpi = 300)
p = annotate!(15,23,text("Spearman's ρ = $(round(rho,digits = 2))", 11, color = :darkorange2))
#p = annotate!(15,52,text("$(length(idx)) common genes", 11))

savefig(p, "data/paper_figures/figure_4/half_life_nc_scifate.svg")
savefig(p, "data/paper_figures/figure_4/half_life_nc_scifate.png")



########################################################################################################################
#####################################    Plot posterior model fits    ##################################################
########################################################################################################################
function plot_fit(data::Matrix{Float64}, data_se::Matrix{Float64}, s_mean::Matrix{Float64}, s_ff::Matrix{Float64}, gene_name::String, m::Int64)
    data_col = :mediumpurple4
    col = [:skyblue4, :lightcyan4, :skyblue4, :seagreen4, :gold][m]
    model_label = ["constant scaling", "constant non-scaling", 
                "varying burst frequency", "varying burst size", "varying decay rate"][m]
    agevec = [0.1:0.2:0.9;]
    s_mean_min = s_mean[argmin(sum(s_mean, dims=2)[:,1]),:]    #minimum(s_mean,dims=1)[1,:]
    s_mean_mean = mean(s_mean,dims=1)[1,:]
    s_mean_max = s_mean[argmax(sum(s_mean, dims=2)[:,1]),:]    #maximum(s_mean,dims=1)[1,:]
    s_ff_min = s_ff[argmin(sum(s_ff, dims=2)[:,1]),:]         #minimum(s_ff,dims=1)[1,:]
    s_ff_mean = mean(s_ff,dims=1)[1,:]
    s_ff_max = s_ff[argmax(sum(s_ff, dims=2)[:,1]),:]          #maximum(s_ff,dims=1)[1,:]
    lims = [(0.0,maximum(vcat(data[:,1] + data_se[:,1], s_mean_max)) + 2.0), (0.0,maximum(vcat(data[:,2] + data_se[:,2], s_ff_max,10.0)))];
    #lims = [(0.0,maximum(vcat(data[:,1] + data_se[:,1], s_mean_max)) + 1.0), (0.0,5.0)];
    local mean_plot::Plots.Plot{Plots.GRBackend}
    local ff_plot::Plots.Plot{Plots.GRBackend}
    mean_plot = plot(agevec,
        hcat(data[:,1] - data_se[:,1],data[:,1] + data_se[:,1]),
        fillrange = data[:,1],
        fillalpha = 0.25,
        fillcolor = [data_col data_col],
        linealpha = [0.0 0.0],
        linecolor = [data_col data_col],
        labels = [nothing nothing],
        shape = [:none :none])
    mean_plot = plot!(agevec,
            data[:,1],
            title = gene_name,
            xlabel = "cell cycle progression",
            ylabel = "mean transcript levels",
            label = "data",
            shape = :circle,
            legend = :topleft,
            legendfontsize = 11,
            fg_legend = :transparent,
            linecolor = data_col,
            linewidth = 5.0,
            markercolor = data_col)
    mean_plot = plot!(agevec,
            hcat(s_mean_min,s_mean_max),
            fillrange = s_mean_mean,
            fillalpha = 0.25,
            fillcolor = [col col],
            linecolor = [col col],
            labels = [nothing nothing],
            shape = [:none :none],
            linewidth = 0.1)
    mean_plot = plot!(agevec,
            s_mean_mean,
            ylims = lims[1],
            linecolor = col,
            linestyle = :dash,
            label = model_label*" model fit",
            linewidth = 5.0,
            size = (400,300),
            background_color_legend = nothing)
    ff_plot = plot(agevec,
        hcat(data[:,2] - data_se[:,2],data[:,2] + data_se[:,2]),
        fillrange = data[:,2],
        fillalpha = 0.25,
        fillcolor = [data_col data_col],
        linealpha = [0.0 0.0],
        linecolor = [data_col data_col],
        labels = [nothing nothing],
        shape = [:none :none])
    ff_plot = plot!(agevec,
            data[:,2],
            title = gene_name,
            xlabel = "cell cycle progression",
            ylabel = "Fano factor of transcript levels",
            label = "data",
            shape = :circle,
            legend = :topleft,
            legendfontsize = 11,
            fg_legend = :transparent,
            linecolor = data_col,
            linewidth = 5.0,
            markercolor = data_col)
    ff_plot = plot!(agevec,
            hcat(s_ff_min,s_ff_max),
            fillrange = s_ff_mean,
            fillalpha = 0.25,
            fillcolor = [col col],
            linecolor = [col col],
            labels = [nothing nothing],
            shape = [:none :none],
            linewidth = 0.1)
    ff_plot = plot!(agevec,
            s_ff_mean,
            ylims = lims[2],
            linecolor = col,
            linestyle = :dash,
            label = model_label*" model fit",
            linewidth = 5.0,
            size = (400,300),
            background_color_legend = nothing)   
    return mean_plot,ff_plot
end

function plot_id_specific_fits(data::Vector{Float64},se::Vector{Float64}, s_data::Matrix{Float64}, gene_name::String, m::Int64, which_data::Int64)
    data_label = ["mean labelled vs total ratio", "mean correlation", "correlation of means"][which_data]
    col = [:skyblue4, :lightcyan4, :skyblue4, :seagreen4, :gold][m]
    model_label = ["constant scaling", "constant non-scaling", 
                "varying burst frequency", "varying burst size", "varying decay rate"][m]
    s_min = s_data[argmin(sum(s_data, dims=2)[:,1]),:]    #minimum(s_mean,dims=1)[1,:]
    s_mean = mean(s_data,dims=1)[1,:]
    s_max = s_data[argmax(sum(s_data, dims=2)[:,1]),:] 
    condition_id = hcat([1/4,1/2,3/4,1,2,3,22,22,22,22,22], [0,0,0,0,0,0,0,1,2,4,6])
    ticks = vcat(condition_id[1:6,1], condition_id[7:end,1] .+ condition_id[7:end,2] .- 17)
    x_tick_labels = vcat(["$j" for j in condition_id[1:3,1]], ["$(round(Int64,j))" for j in condition_id[4:6,1]], ["$(round(Int64,i))"* " + " *"$(round(Int64,j))" for i in condition_id[7:end,1] for j in condition_id[7:end,2]])
    idx::Vector{Int64} = [1,4,5,6,7,8,9,10,11]
    data_col = :mediumpurple4
    p = vspan([0,3], color = [:lightcyan3], alpha = 0.35, label = "pulse conditions")
    p = vspan!([5,11], color = [:mediumpurple4], alpha = 0.15, label = "chase conditions")
    p = plot!(ticks, hcat(data .+ se, data .- se), fillrange = data, fillalpha = 0.25, xrotation=30, linealpha = [0.0 0.0],
            xticks = (ticks[idx], x_tick_labels[idx]), legend = :topleft,fg_legend = :transparent,background_color_legend = :transparent, legendfontsize = 9,
            fillcolor = [data_col data_col], linecolor = [data_col data_col], labels = [nothing nothing], ylims = (-0.05,1.0), size = (400,300), dpi = 300)
    p = plot!(ticks, data, xlabel = "labelling time (h)", ylabel = data_label, title = gene_name,
             linewidth = 5, markersize = 5, label = "data", linecolor = data_col, margin = 2mm)
    p = plot!(ticks, hcat(s_min,s_max), fillrange = s_mean,fillalpha = 0.25, fillcolor = [col col],linecolor = [col col],
        labels = [nothing nothing], shape = [:none :none], linewidth = 0.1)
    p = plot!(ticks, s_mean,linecolor = col,label = model_label*" model fit", linestyle = :dash, linewidth = 5)
    return p
end

model_names = ["const","const_const","kon","alpha","gamma"];

n_sim_groups = 1
s_pulse = Vector{Matrix{Float64}}(undef,length(model_names))
s_chase = Vector{Matrix{Float64}}(undef,length(model_names))
for m in 1:5 
    s_pulse[m], s_chase[m] = load_s_data("data/large_scale_simulations/",m,n_sim_groups,[1,2])
end

s_ratios = Vector{Matrix{Float64}}(undef,length(model_names))
s_mean_corr = Vector{Matrix{Float64}}(undef,length(model_names))
s_corr_mean = Vector{Matrix{Float64}}(undef,length(model_names))
for m in 1:5 
    s_ratios[m], s_mean_corr[m], s_corr_mean[m] = load_s_data("data/large_scale_simulations/",m,n_sim_groups,[3,4,5])
end

ms_df = CSV.read("data/model_selection/model_selection_results.txt", DataFrame);
ms = Int64.(Matrix(ms_df[:,[1,3]]));

#map_sets = [get_posterior_estimate(sets[m],posterior[m],nc_genes[1][m-2],"map") for m in 3:5]
#post_mean_sets = [get_posterior_estimate(sets[m],posterior[m],nc_genes[1][m-2],"mean") for m in 3:5]

#pick model and gene
m = 2
model_name = model_names[m]

g = 30              #const_scaling: 51, 102    #non_scaling:     #burst_freq: 44    #burst_size: 10  # decay rate:  58, 8
g_name = String(ms_df.gene_name[ms_df.model .== m][g])


n_particles_mat[ms[ms[:,2] .== m,1][g],m]
#ms_df.gene_name[ms_df.model .== m]

#get data
p_mean = pulse_mean[ms[ms[:,2] .== m,1][g],:]
p_ff = pulse_ff[ms[ms[:,2] .== m,1][g],:]

c_mean = chase_mean[ms[ms[:,2] .== m,1][g],:]
c_ff = chase_ff[ms[ms[:,2] .== m,1][g],:]

p_mean_se = pulse_mean_se[ms[ms[:,2] .== m,1][g],:]
p_ff_se = pulse_ff_se[ms[ms[:,2] .== m,1][g],:]
c_mean_se = chase_mean_se[ms[ms[:,2] .== m,1][g],:]
c_ff_se = chase_ff_se[ms[ms[:,2] .== m,1][g],:]


#get posterior fits
prop = 1.0;
particles = posterior[m][ms[ms[:,2] .== m,1][g]][1:Int64(floor(n_particles_mat[ms[ms[:,2] .== m,1][g],m] * prop))]

q = 1.0;

s_p_mean = filter_quantiles(get_mean_subset(s_pulse[m])[particles,:],q)
s_p_ff = filter_quantiles(get_ff_subset(s_pulse[m])[particles,:],q)
s_c_mean = filter_quantiles(get_mean_subset(s_chase[m])[particles,:],q)
s_c_ff = filter_quantiles(get_ff_subset(s_chase[m])[particles,:],q)


p1, p2 = plot_fit(hcat(p_mean,p_ff),hcat(p_mean_se,p_ff_se),s_p_mean,s_p_ff,g_name,m);
p1
p2

p1, p2 = plot_fit(hcat(c_mean,c_ff),hcat(c_mean_se,c_ff_se),s_c_mean,s_c_ff,g_name,m);
p1
p2

savefig(p1,"data/model_fits/"*model_name*"/"*lowercase(g_name)*"_mean_fit.svg")
savefig(p2,"data/model_fits/"*model_name*"/"*lowercase(g_name)*"_ff_fit.svg")


# other summary statistics=
r_data = ratio_data[ms[ms[:,2] .== m,1][g],:]
r_se = ratio_se[ms[ms[:,2] .== m,1][g],:]

m_c_data = mean_corr_data[ms[ms[:,2] .== m,1][g],:]
m_c_se = mean_corr_se[ms[ms[:,2] .== m,1][g],:]

c_m_data = corr_mean_data[ms[ms[:,2] .== m,1][g],:]
c_m_se = corr_mean_se[ms[ms[:,2] .== m,1][g],:]

q = 1.0;

s_r = filter_quantiles(s_ratios[m][particles,:],q)
s_m_c = filter_quantiles(s_mean_corr[m][particles,:],q)
s_c_m = filter_quantiles(s_corr_mean[m][particles,:],q)

p1 = plot_id_specific_fits(r_data,r_se,s_r,g_name, m, 1)

p2 = plot_id_specific_fits(m_c_data,m_c_se,s_m_c,g_name, m, 2)

p3 = plot_id_specific_fits(c_m_data,c_m_se,s_c_m,g_name,m, 3)

savefig(p1,"data/model_fits/"*model_name*"/"*lowercase(g_name)*"_ratio_fit.svg")
savefig(p2,"data/model_fits/"*model_name*"/"*lowercase(g_name)*"_mean_corr_fit.svg")
savefig(p3,"data/model_fits/"*model_name*"/"*lowercase(g_name)*"_corr_mean_fit.svg")



#plot posterior kinetic rates

vary_flag = [[0,0,0,0],[0,0,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,1]][m]
scaling = 1 * (!=(m,2))
n_steps = 5
cycle = 20.0

θ = sets[m][particles,:]
q = 0.75;          #constant genes: q = 0.75 # burst size gene: 0.7

p = plot_posterior_rate(filter_quantiles(θ,q),vary_flag,n_steps,scaling,cycle,g_name,m);
p

savefig(p,"data/model_fits/"*model_name*"/"*lowercase(g_name)*"_rate.svg")



p = plot_rates(θ,vary_flag,n_steps,scaling,cycle,g_name,m);
p

savefig(p,"data/model_fits/"*model_name*"/"*lowercase(g_name)*"_all_rates.svg")



function plot_posterior_rate(θ::AbstractArray{Float64}, vary_flag::Vector{Int64}, n_steps::Int64, scaling::Int64, cycle::Float64, gene_name::String, m::Int64)
    col = [:skyblue4, :lightcyan4, :skyblue4, :seagreen4, :gold][m]
    rate_label = ["synthesis rate", "synthesis rate", 
                "burst frequency", "burst size", "decay rate"][m]          
    tspan::Vector{Float64} = [0.0:0.1:cycle-0.1;]
    t_steps::Vector{Float64} = [0.0:cycle/n_steps:cycle;]
    vary_map = get_vary_map(vary_flag,n_steps)
    if m >= 3
        vary_idx::Vector{Int64} = findall(x->x>0,vary_flag)
    else
        vary_idx = [3]
    end
    θ_lims = [minimum(θ,dims=1)[1,:],maximum(θ,dims=1)[1,:]]
    local denom::Vector{Float64}
    local p::Plots.Plot{Plots.GRBackend}
    y = Matrix{Float64}(undef,(length(vary_flag),length(tspan)))
    for (k,t) in enumerate(tspan)
        y[:,k] = 10 .^(get_rate(θ[1,1:end-1], vary_map, scaling, t_steps, cycle, t))
    end
    p = plot(tspan ./ cycle, y[vary_idx[1],:],
                    linewidth = 5.0,
                    xlabel = "cell cycle progression",
                    ylabel = rate_label,
                    label = false,
                    title = gene_name,
                    size = (400, 300),
                    linecolor = col)
    y1 = Matrix{Float64}(undef,(length(vary_flag),length(tspan)))
    y2 = Matrix{Float64}(undef,(length(vary_flag),length(tspan)))
    for (k,t) in enumerate(tspan)
        y1[:,k] = 10 .^(get_rate(θ_lims[1][1:end-1], vary_map, scaling, t_steps, cycle, t))
        y2[:,k] = 10 .^(get_rate(θ_lims[2][1:end-1], vary_map, scaling, t_steps, cycle, t))
    end
    p = plot!(tspan ./ cycle, 
        hcat(y1[vary_idx[1],:],y2[vary_idx[1],:]),
        fillrange = y[vary_idx[1],:],
        ylims = (0.0, maximum(y2[vary_idx[1],:]) + 5.0),
        fillalpha = 0.25,
        fillcolor = [col col],
        linealpha = [0.0 0.0],
        linecolor = [col col],
        labels = [nothing nothing],
        shape = [:none :none])
    return p
end

function plot_rates(θ::AbstractArray{Float64}, vary_flag::Vector{Int64}, n_steps::Int64, scaling::Int64, cycle::Float64, gene_name::String, m::Int64)
    u = median(θ[:,1:end-1],dims=1)[1,:]
    col = reshape([:skyblue4, :seagreen4, :gold],1,:)         
    tspan::Vector{Float64} = [0.0:0.1:cycle-0.1;]
    t_steps::Vector{Float64} = [0.0:cycle/n_steps:cycle;]
    vary_map = get_vary_map(vary_flag,n_steps)
    local p::Plots.Plot{Plots.GRBackend}
    y = Matrix{Float64}(undef,(length(vary_flag),length(tspan)))
    for (k,t) in enumerate(tspan)
        y[:,k] = get_rate(u, vary_map, scaling, t_steps, cycle, t)
    end
    #y[vary_map[3],:] = y[vary_map[3],:] .- y[vary_map[2],:]
    p = plot(tspan ./ cycle, y[[1,3,4],:]',
                    ylims = (-3,2),
                    linewidth = 5.0,
                    xlabel = "cell cycle progression",
                    ylabel = "log₁₀(kinetic rates)",
                    label = ["burst frequency" "synthesis rate" "decay rate"],
                    title = gene_name,
                    size = (400, 300),
                    #left_margin = 4Plots.mm,
                    #bottom_margin = 5Plots.mm,
                    legendfontsize = 10,
                    legend = :right,
                    background_color_legend = :transparent,
                    fg_legend = :transparent,
                    linecolor = col)
    return p
end


################################# Posterior distribution of a single gene #######################################


function plot_const_posterior(θ::Matrix{Float64},err::Vector{Float64},m::Int64)
    rate_names::Vector{String} = ["burst frequency","synthesis rate","burst size","decay rate"]
    lims = [(-3,3),(-3,3),(-6,6),(-3,2)]
    col = [:skyblue4,:lightcyan4][m]
    p1 = [density(θ[:,i], linewidth = 4, xlims = lims[i], 
        label = false, ylabel = "probability density", title = rate_names[i], titlefontsize = 11, color = col, size = (350,250), dpi=300) for i in 1:lastindex(rate_names)]
    # xticks = round.(collect(range(lims[i][1],lims[i][end],5))
    #yticks = round.(collect(range(lims[j][1],lims[j][end],5))
    p2 = [scatter(θ[sortperm(err,rev=true),i], θ[sortperm(err,rev=true),j], zcolor = sort(err,rev=true), color = :viridis, xlims = lims[i], ylims = lims[j],colorbar_ticks = nothing,
        xlabel = rate_names[i], ylabel = rate_names[j], markersize = 5,markeralpha = 0.8,markerstrokewidth = 0.0, label = false, size = (350,250),right_margin = 2mm, dpi=300) for i in 1:size(θ)[2] for j in 1:i-1]
    #pairplot = plot(p..., layout = (size(θ)[2],size(θ)[2]),plot_title = model_label*" gene MAPs", size = (900,800));
    #return pairplot
    return p1,p2
end

function plot_const_posterior_density(θ::Matrix{Float64},m::Int64)
    rate_names::Vector{String} = ["burst frequency","burst size","synthesis rate","decay rate"]
    lims = [(-3,3),(-6,6),(-3,3),(-3,2)]
    col = [:skyblue4,:lightcyan4][m]
    p1 = [density(θ[:,i], linewidth = 4, xlims = lims[i], 
        label = false, ylabel = "probability density", title = rate_names[i], titlefontsize = 11, color = col, size = (350,250), dpi=300) for i in 1:lastindex(rate_names)]
    # xticks = round.(collect(range(lims[i][1],lims[i][end],5))
    #yticks = round.(collect(range(lims[j][1],lims[j][end],5))
    p2 = [plot(KernelDensity.kde((θ[:,i], θ[:,j])), linewidth = 3, color = :viridis, xlims = lims[i], ylims = lims[j],colorbar_ticks = nothing,
        xlabel = rate_names[i], ylabel = rate_names[j], markersize = 5,markeralpha = 0.8,markerstrokewidth = 0.0, label = false, colorbar_title = "\nprobability density",
        size = (350,250),right_margin = 3mm, dpi=300) for i in 1:size(θ)[2] for j in 1:i-1]
    #pairplot = plot(p..., layout = (size(θ)[2],size(θ)[2]),plot_title = model_label*" gene MAPs", size = (900,800));
    #return pairplot
    return p1,p2
end


function get_burst_kinetics(sets::Matrix{Float64}, m::Int64)
    #kinetics_names::Vector{String} = ["kon","burst size","decay rate"]
    vary_flag::Vector{Int64} = [[0,0,0,0],[0,0,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,1]][m];
    vary_map::Vector{Any} = get_vary_map(vary_flag,5)
    #kinetics = Matrix{Float64}(undef,(size(map_sets)[1],sum(length.(vary_map))))
    kinetics::Matrix{Float64} = hcat(sets[:,vary_map[1]],sets[:,vary_map[3]] .- sets[:,vary_map[2]],sets[:,vary_map[3]],sets[:,vary_map[4]])
    #return hcat(kinetics,sets[:,end])
    return kinetics
end

m = 5
model_name = model_names[m]

g = 16        #const_scaling: 51, 102    #non_scaling: 146, 316, 80, 77   #burst_freq: 44    #burst_size: 10  # decay rate:  58, 8
g_name = String(ms_df.gene_name[ms_df.model .== m][g])


n_particles_mat[ms[ms[:,2] .== m,1][g],m]

prop = 1.0;
particles = posterior[m][ms[ms[:,2] .== m,1][g]][1:Int64(floor(n_particles_mat[ms[ms[:,2] .== m,1][g],m] * prop))]


θ = sets[m][particles,:]
kinetics = get_burst_kinetics(θ,m)
lambdas = θ[:,end]

p = density(10 .^lambdas, xlabel = "labelling rate", ylabel = "probabiliy density", title = g_name, xlims = (10^(-0.8), 1.1), linewidth = 5, color = :skyblue4, size = (400,300), dpi = 100,label = false)
savefig(p,"gene_lambda_$m.png")


p1,p2 = plot_const_posterior_density(kinetics,m);
[savefig(p1[i],"data/paper_figures/supplement/fig_3/posterior_densities_$i.svg") for i in 1:lastindex(p1)];
[savefig(p2[i],"data/paper_figures/supplement/fig_3/posterior_scatters_$i.svg") for i in 1:lastindex(p2)];



#################################  laballing rate and gene length  #######################################


col = [:skyblue4, :lightcyan3, :coral3, :mediumpurple3, :goldenrod1]


genome_length_df = CSV.read("data/genes_loc_H.csv",DataFrame)


all_sel_genes = vcat(sel_genes...)
sub_genes = all_sel_genes[sort(findall(x->x in genome_df[:,5], gene_names[all_sel_genes]))]

lambdas = vcat([m_s[:,end] for m_s in map_sets]...)[sort(findall(x->x in genome_df[:,5], gene_names[all_sel_genes]))]


loc_idx = [];
for g in gene_names[sub_genes]
    push!(loc_idx,findall(x->x==g,genome_df[:,5])[1])
end

gene_length = genome_df[loc_idx,4] .- genome_df[loc_idx,3]

rho = corspearman(10 .^lambdas,log.(10,gene_length))


# group by model
sub_genes = [sel[sort(findall(x->x in genome_df[:,5], gene_names[sel]))] for sel in sel_genes]

lambdas = [map_sets[m][sort(findall(x->x in genome_df[:,5], gene_names[sel_genes[m]])),end] for m in 1:5]


loc_idx = [[],[],[],[],[]];
for m in 1:5
    for g in gene_names[sub_genes[m]]
        push!(loc_idx[m],findall(x->x==g,genome_df[:,5])[1])
    end
end


gene_length = [genome_df[loc_idx[m],4] .- genome_df[loc_idx[m],3] for m in 1:5]

labs = ["constant scaling genes","constant non-scaling genes","burst frequency genes","burst size genes","decay rate genes"];
p = scatter(log.(10,gene_length[1]), 10 .^lambdas[1],  y_lims = (0.1,1.1), x_lims = (2.0,7.5), ylabel = "labelling efficiency", xlabel = "log₁₀(gene length)", legendfontsize = 7, xlabelfontsize = 13, ylabelfontsize = 13,
        color = col[1], label = false, fg_legend = :transparent, markeralpha = 0.9, markersize = 5, background_color_legend = nothing, size = (450,350), dpi = 300);
for i in 2:5
    p = scatter!(log.(10,gene_length[i]), 10 .^lambdas[i], color = col[i], markeralpha = 0.9, markersize = 5, label = false)
end
p = annotate!(4.5,0.15,text("overall Spearman's ρ = $(round(rho,digits = 2))", 12, color = :lightcyan4))

savefig(p,"data/paper_figures/figure_3/gene_length_lambdas_per_m.svg")
