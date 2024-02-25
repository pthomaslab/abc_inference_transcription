###########################     distributions of MAP estimates across genes     ##########################
vary_flags::Vector{Vector{Int64}} = [[0,0,0,0],[0,0,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,1]]
vary_maps::Vector{Vector{Any}} = [get_vary_map(vf,5) for vf in vary_flags];

col = [:skyblue4, :seagreen4, :gold]
x_labs = ["constant \ngenes","burst frequency \ngenes","burst size \ngenes","decay rate \ngenes"];

p = violin(vcat(map_sets[1][:,vary_maps[1][4]],map_sets[2][:,vary_maps[2][4]]), xticks = ([1:4;],x_labs), 
    linewidth = 1, linealpha = 0.7,color = :grey, label = false, outliers = false);
p = violin!(map_sets[3][:,vary_maps[3][4]], linewidth = 1, color = col[1], label = false, legendfontsize = 9, 
background_color_legend = :transparent,legend = :topleft,fg_legend = :transparent, size = (400,300), ylabel = "log₁₀(decay rate)", outliers = false)
p = violin!(map_sets[4][:,vary_maps[4][4]], linewidth = 1,linealpha = 0.7, color = col[2], label = false, outliers = false)
p = violin!(vec(map_sets[5][:,vary_maps[5][4]]), linewidth = 1, linealpha = 0.7,color = col[3], label = false, outliers = false)
savefig(p, "data/paper_figures/figure_2/decay_rates.svg")

p = violin(vcat(map_sets[1][:,vary_maps[1][1]],map_sets[2][:,vary_maps[2][1]]), xticks = ([1:4;],x_labs), 
        linewidth = 1, linealpha = 0.7,color = :grey, label = false, outliers = false);
p = violin!(vec(map_sets[3][:,vary_maps[3][1]]), linewidth = 1, color = col[1], label = false, legendfontsize = 9, 
background_color_legend = :transparent,legend = :topleft,fg_legend = :transparent, size = (400,300), ylabel = "log₁₀(burst frequency)", outliers = false)
p = violin!(map_sets[4][:,vary_maps[4][1]] , linewidth = 1,linealpha = 0.7, color = col[2], label = false, outliers = false)
p = violin!(map_sets[5][:,vary_maps[5][1]], linewidth = 1, linealpha = 0.7,color = col[3], label = false, outliers = false)

savefig(p, "data/paper_figures/figure_2/burst_freq.svg")


p = violin(vcat(map_sets[1][:,vary_maps[1][3]] .- map_sets[1][:,vary_maps[1][2]],map_sets[2][:,vary_maps[2][3]] .- map_sets[2][:,vary_maps[2][2]]), xticks = ([1:4;],x_labs), 
        linewidth = 1, linealpha = 0.7,color = :grey, label = false, outliers = false);
p = violin!(map_sets[3][:,vary_maps[3][3]] .- map_sets[3][:,vary_maps[3][2]], linewidth = 1, color = col[1], label = false, legendfontsize = 9, 
background_color_legend = :transparent,legend = :topleft,fg_legend = :transparent, size = (400,300), ylabel = "log₁₀(burst size)", outliers = false)
p = violin!(vec(map_sets[4][:,vary_maps[4][3]] .- map_sets[4][:,vary_maps[4][2]]), linewidth = 1,linealpha = 0.7, color = col[2], label = false, outliers = false)
p = violin!(map_sets[5][:,vary_maps[5][3]] .- map_sets[5][:,vary_maps[5][2]], linewidth = 1, linealpha = 0.7,color = col[3], label = false, outliers = false)

savefig(p, "data/paper_figures/figure_2/burst_size.svg")


x = mean(t_data,dims=1)[1,:]
p = violin(vcat(x[sel_genes[1]],x[sel_genes[2]]), xticks = ([1:4;],x_labs), 
    linewidth = 1, linealpha = 0.7,color = :grey, label = false, ylims = (0.0,10), outliers = false);
p = violin!(x[sel_genes[3]], linewidth = 1, color = col[1], label = false, legendfontsize = 9, 
background_color_legend = :transparent,legend = :topleft,fg_legend = :transparent, size = (400,300), ylabel = "mean expression", outliers = false)
p = violin!(x[sel_genes[4]], linewidth = 1,linealpha = 0.7, color = col[2], label = false, outliers = false)
p = violin!(x[sel_genes[5]], linewidth = 1, linealpha = 0.7,color = col[3], label = false, outliers = false)
savefig(p, "data/paper_figures/figure_2/mean_expression.svg")


############################### mRNA half-lives #######################################
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
