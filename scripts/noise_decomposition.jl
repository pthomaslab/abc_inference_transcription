function weighted_cov(x::Vector{Float64},y::Vector{Float64}, w::Vector{Float64})
    return sum(w .* ((x .- sum(w .* x)) .* (y .- sum(w .* y))))
end

cells_per_id = [findall(x->x==e,experiment) for e in sort(unique(experiment))[2:end-1]];
cells_per_age = [findall(x->x==τ,age) for τ in sort(unique(age))];
cells_age_id = [[intersect(cells_age,cells_id) for cells_age in cells_per_age] for cells_id in cells_per_id];
cells_id_age = [[intersect(cells_age,cells_id) for cells_id in cells_per_id] for cells_age in cells_per_age];

n_cells_id = sum(length.(cells_per_id));
age_id_distribution = [length.(cells) ./ n_cells_id for cells in cells_age_id];
age_dist = [sum(length.(cells)) for cells in cells_id_age] ./ n_cells_id;  
id_distribution = length.(cells_per_id) ./ n_cells_id;

# noise landscape based on raw data 
mn = [mean(t_data[:,sg],dims=1)[1,:] for sg in sel_genes];
v = [var(t_data[:,sg],dims=1)[1,:] for sg in sel_genes];
cv = [v[i] ./ (mn[i].^2) for i in 1:length(mn)];

y = vcat(minimum(vcat(mn...))-0.2,sort(vcat(mn...)), maximum(vcat(mn...))+2.0);
p = scatter(log.(10,mn[1]), log.(10,cv[1]), color = :gray,  ylabel = "log₁₀(CV²)", xlims = (-1,3.5),ylims = (-2,2), xlabel = "log₁₀(mean)", label = "all genes (data)", 
    markersize = 5, markeralpha = 0.6, legendfontsize = 9,background_color_legend = :transparent, fg_legend = :transparent,xlabelfontsize = 13, ylabelfontsize = 13, size = (300,250), dpi=300);
p = scatter!(log.(10,mn[2]), log.(10,cv[2]), color = :gray, markersize = 5, markeralpha = 0.6,label = false);
p = scatter!(log.(10,mn[3]), log.(10,cv[3]), color = :gray, markersize = 5,markeralpha = 0.6,label = false);
p = scatter!(log.(10,mn[4]), log.(10,cv[4]), color = :gray, markersize = 5,markeralpha = 0.6,label = false);
p = scatter!(log.(10,mn[5]), log.(10,cv[5]), color = :gray, markersize = 5,markeralpha = 0.6,label =  false);
p = plot!(log.(10,y), log.(10, 1 ./ y), label = "1/mean scaling", linewidth = 4, linestyle = :dash, color = :lightblue3);
p = plot!(log.(10,y), ones(length(y)) .* log(10,minimum(minimum.(cv))), label = "extrinsic noise limit", linewidth = 4, linestyle = :dash, color = :coral3, linealpha = 0.75);
savefig(p, "data/paper_figures/figure_3/data_based/data_noise_mean.svg");


id_labels = ["pulse_15", "pulse_30", "pulse_45",
            "pulse_60", "pulse_120", "pulse_180", "chase_0",
            "chase_60", "chase_120", "chase_240", "chase_360"];


u_mean = [[hcat([mean(u_data[cell,sel],dims=1)[1,:] for cell in cells]...) for cells in cells_age_id] for sel in sel_genes];
l_mean = [[hcat([mean(l_data[cell,sel],dims=1)[1,:] for cell in cells]...) for cells in cells_age_id] for sel in sel_genes];
tot_mean = u_mean .+ l_mean;

u_var = [[hcat([var(u_data[cell,sel],dims=1)[1,:] for cell in cells]...) for cells in cells_age_id] for sel in sel_genes];
l_var = [[hcat([var(l_data[cell,sel],dims=1)[1,:] for cell in cells]...) for cells in cells_age_id] for sel in sel_genes];

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

mean_ccp = [map(.+,[m .* id_distribution[j] for (j,m) in enumerate(m_)]...) for m_ in tot_mean];
u_m = [sum(hcat([sum(transpose(transpose(u) .* age_id),dims=2)[:,1] for (u,age_id) in zip(u_,age_id_distribution)]...),dims=2)[:,1] for u_ in u_mean];
l_m = [sum(hcat([sum(transpose(transpose(l) .* age_id),dims=2)[:,1] for (l,age_id) in zip(l_,age_id_distribution)]...),dims=2)[:,1] for l_ in l_mean];
tot_m = u_m .+ l_m;

var_mean_ccp = Vector{Vector{Float64}}(undef,length(sel_genes));
for i in 1:length(sel_genes)
    var_mean_ccp[i] = Vector{Float64}(undef,length(sel_genes[i]))
    for j in 1:length(sel_genes[i])
        var_mean_ccp[i][j] = weighted_cov(mean_ccp[i][j,:],mean_ccp[i][j,:],age_dist) 
    end
end
ccd_noise = [v_ ./ (m_ .^2) for (v_,m_) in zip(var_mean_ccp,tot_m)];



# noise landscape based on recovered statistics 

col = [:skyblue4, :lightcyan3, :coral3, :mediumpurple3, :goldenrod1]

s_u_mean = [[readdlm("Julia/recovered_statistics/"*mn*"/"*id*"/mean_u.txt") for id in id_labels] for mn in model_names];
s_l_mean = [[readdlm("Julia/recovered_statistics/"*mn*"/"*id*"/mean_l.txt") for id in id_labels] for mn in model_names];
s_mean = [s_u .+ s_l for (s_u,s_l) in zip(s_u_mean,s_l_mean)];

s_mean_ccp = [sm[7] for sm in s_mean];

s_u_var = [[readdlm("Julia/recovered_statistics/"*mn*"/"*id*"/var_u.txt") for id in id_labels] for mn in model_names];
s_l_var = [[readdlm("Julia/recovered_statistics/"*mn*"/"*id*"/var_l.txt") for id in id_labels] for mn in model_names];
s_ul_cov = [[readdlm("Julia/recovered_statistics/"*mn*"/"*id*"/cov_ul.txt") for id in id_labels] for mn in model_names];

#intrinsic noise components
s_u = [sum(hcat([sum(transpose(transpose(su) .* age_id),dims=2)[:,1] for (su,age_id) in zip(s_u_,age_id_distribution)]...),dims=2)[:,1] for s_u_ in s_u_mean];
s_l = [sum(hcat([sum(transpose(transpose(sl) .* age_id),dims=2)[:,1] for (sl,age_id) in zip(s_l_,age_id_distribution)]...),dims=2)[:,1] for s_l_ in s_l_mean];
s_m = s_u .+ s_l;

s_mean_var_u = [sum(hcat([sum(transpose(transpose(su) .* age_id),dims=2)[:,1] for (su,age_id) in zip(s_u_,age_id_distribution)]...),dims=2)[:,1] for s_u_ in s_u_var];
s_mean_var_l = [sum(hcat([sum(transpose(transpose(sl) .* age_id),dims=2)[:,1] for (sl,age_id) in zip(s_l_,age_id_distribution)]...),dims=2)[:,1] for s_l_ in s_l_var];
s_mean_var = s_mean_var_u .+ s_mean_var_l;

s_mean_cov = [sum(hcat([sum(transpose(transpose(s) .* age_id),dims=2)[:,1] for (s,age_id) in zip(s_,age_id_distribution)]...),dims=2)[:,1] for s_ in s_ul_cov];

s_int_v = s_mean_var .+ (2 .* s_mean_cov);

# extrinsic noise components
s_ext_v = Vector{Vector{Float64}}(undef,length(sel_genes));
for i in 1:length(sel_genes)
    s_ext_v[i] = Vector{Float64}(undef,length(sel_genes[i]))
    for j in 1:length(sel_genes[i])
        s_ext_v[i][j] = weighted_cov(vcat([s_mean[i][k][j,:] for k in 1:length(id_labels)]...),vcat([s_mean[i][k][j,:] for k in 1:length(id_labels)]...),vcat(age_id_distribution...))
    end
end

s_var_mean_ccp = Vector{Vector{Float64}}(undef,length(sel_genes));
for i in 1:length(sel_genes)
    s_var_mean_ccp[i] = Vector{Float64}(undef,length(sel_genes[i]))
    for j in 1:length(sel_genes[i])
        s_var_mean_ccp[i][j] = weighted_cov(s_mean_ccp[i][j,:],s_mean_ccp[i][j,:],age_dist)  
    end
end

s_var_mean_time = s_ext_v .- s_var_mean_ccp;  # --> ~zero contribution

s_ext_v = s_var_mean_ccp;
s_v =  s_int_v .+ s_ext_v;

s_tot_noise = [v_ ./ (m .^2) for (v_,m) in zip(s_v,s_m)];
s_int_noise = [v_ ./ (m .^2) for (v_,m) in zip(s_int_v,s_m)];
s_ext_noise = [v_ ./ (m .^2) for (v_,m) in zip(s_ext_v,s_m)];

s_tx_deg_noise = [v_ ./ (m .^2) for (v_,m) in zip(s_mean_var,s_m)];
s_bursty_noise = s_int_noise .- s_tx_deg_noise;
s_ccd_noise = [v_ ./ (m .^2) for (v_,m) in zip(s_var_mean_ccp,s_m)];

s_tx_deg_prop = [tx_deg ./ tot for (tx_deg,tot) in zip(s_tx_deg_noise,s_tot_noise)];
s_bursty_prop = [bursty ./ tot for (bursty,tot) in zip(s_bursty_noise,s_tot_noise)];
s_ccd_prop = [ccd ./ tot for (ccd,tot) in zip(s_ccd_noise,s_tot_noise)];


# total noise vs mean
x = vcat([minimum(vcat(s_m...)) - 2.0],sort(vcat(s_m...)));
p = scatter(log.(10,s_m[1]), log.(10,s_tot_noise[1]), color = col[1],  ylabel = "log₁₀(CV²)",xlims = (-1,3.5),ylims = (-2,2), xlabel = "log₁₀(mean)", xlabelfontsize = 13, ylabelfontsize = 13,
    label = "constant scaling genes", markersize = 5, markeralpha = 0.9, legendfontsize = 9, legend = :topright, size = (450,350), dpi = 200);
p = scatter!(log.(10,s_m[2]), log.(10,s_tot_noise[2]), color = col[2], markersize = 5, markeralpha = 0.9,label = "constant non-scaling genes");
p = scatter!(log.(10,s_m[3]), log.(10,s_tot_noise[3]), color = col[3], markersize = 5,markeralpha = 0.9,label = "burst frequency genes");
p = scatter!(log.(10,s_m[4]), log.(10,s_tot_noise[4]), color = col[4], markersize = 5,markeralpha = 0.9,label = "burst size genes");
p = scatter!(log.(10,s_m[5]), log.(10,s_tot_noise[5]), color = col[5], markersize = 5,markeralpha = 0.9,label = "decay rate genes");
p = plot!(log.(10,x), log.(10, 1 ./ x), label = "1/mean scaling", linewidth = 4, linealpha = 0.8, linestyle = :dash, color = :lightblue3);
savefig(p, "data/paper_figures/figure_3/model_based/total_noise_mean.svg")

# transcription and degradation contribution
p = scatter(log.(10,s_m[1]), s_tx_deg_prop[1] .* 100, color = col[1], titlefontsize = 12, ylabel = "% of total noise",
    xlims = (0,3.5), ylims = (-5,100.1), xlabel = "log₁₀(mean)", label = false, markersize = 5, markeralpha = 0.9, size = (400,300), dpi = 200);
p = scatter!(log.(10,s_m[2]), s_tx_deg_prop[2] .* 100, color = col[2], markersize = 5, markeralpha = 0.9,label = false);
p = scatter!(log.(10,s_m[3]), s_tx_deg_prop[3] .* 100, color = col[3], markersize = 5,markeralpha = 0.9,label = false);
p = scatter!(log.(10,s_m[4]), s_tx_deg_prop[4] .* 100, color = col[4], markersize = 5,markeralpha = 0.9,label = false);
p = scatter!(log.(10,s_m[5]), s_tx_deg_prop[5] .* 100, color = col[5], markersize = 5,markeralpha = 0.9,label = false);
savefig(p, "data/paper_figures/figure_3/model_based/tx_decay_proportion.svg")

# bursty promoter noise contribution 
p = scatter(log.(10,s_m[1]), s_bursty_prop[1] .* 100, color = col[1], ylabel = "% of total noise", xlabel = "log₁₀(mean)", 
xlims = (0,3.5), ylims = (-5,100.1),label = false, markersize = 5, markeralpha = 0.9, size = (400,300), dpi = 200);
p = scatter!(log.(10,s_m[2]), s_bursty_prop[2] .* 100, color = col[2], markersize = 5, markeralpha = 0.9,label = false);
p = scatter!(log.(10,s_m[3]), s_bursty_prop[3] .* 100, color = col[3], markersize = 5,markeralpha = 0.9,label = false);
p = scatter!(log.(10,s_m[4]), s_bursty_prop[4] .* 100, color = col[4], markersize = 5,markeralpha = 0.9,label = false);
p = scatter!(log.(10,s_m[5]), s_bursty_prop[5] .* 100, color = col[5], markersize = 5,markeralpha = 0.9,label = false);
savefig(p, "data/paper_figures/figure_3/model_based/bursty_promoter_proportion.svg")

# cell-cycle dependence noise contribution 
p = scatter(log.(10,s_m[1]), s_ccd_prop[1] .* 100, color = col[1],  titlefontsize = 13, ylabel = "% of total noise", xlabel = "log₁₀(mean)", 
    xlims = (0,3.5), ylims = (-5,100.5), label = false, markersize = 5, markeralpha = 0.9, size = (400,300), dpi = 200);
p = scatter!(log.(10,s_m[2]), s_ccd_prop[2] .* 100, color = col[2], markersize = 5, markeralpha = 0.9,label = false);
p = scatter!(log.(10,s_m[3]), s_ccd_prop[3] .* 100, color = col[3], markersize = 5,markeralpha = 0.9,label = false);
p = scatter!(log.(10,s_m[4]), s_ccd_prop[4] .* 100, color = col[4], markersize = 5,markeralpha = 0.9,label = false);
p = scatter!(log.(10,s_m[5]), s_ccd_prop[5] .* 100, color = col[5], markersize = 5,markeralpha = 0.9,label = false);
savefig(p, "data/paper_figures/figure_3/model_based/ccd_proportion.svg")

# compare observed and recovered cell cycle-dependent variation
col = [:skyblue4, :lightcyan4];
p = boxplot(ccd_noise[1], linewidth = 1, color = col[1], label = false, ylabelfontsize = 10,xtickfontsize = 11, background_color_legend = :transparent,fg_legend = :transparent, size = (500,300), bottom_margin=2mm,
    xticks = ([1,2,3,4],["scaling \ngenes","non-scaling \ngenes","scaling \ngenes","non-scaling \ngenes"]), ylabel = "cell cycle variation", outliers = false);
p = boxplot!(ccd_noise[2], linewidth = 1, color = col[2], label = false, outliers = false);
p = boxplot!(s_ccd_noise[1], linewidth = 1, color = col[1], label = false, ylabelfontsize = 10,xtickfontsize = 11, outliers = false);
p = boxplot!(s_ccd_noise[2], linewidth = 1, color = col[2], label = false, outliers = false);
savefig(p, "data/paper_figures/figure_4/observed_recovered_ccd_noise.svg")


############################     kinetics regulation on noise-mean plots      ############################
burst_size = vcat(map_sets[1][:,3] .- map_sets[1][:,2],map_sets[2][:,3] .- map_sets[2][:,2])
x = vcat(s_m[1:2]...)[sortperm(burst_size)]
y = vcat(s_bursty_noise[1:2]...)[sortperm(burst_size)]



burst_freq = vcat(map_sets[1][:,1],map_sets[2][:,1])
x = vcat(s_m[1:2]...)[sortperm(burst_freq)]
y = vcat(s_bursty_noise[1:2]...)[sortperm(burst_freq)]


y = vcat(s_bursty_noise[1:2]...)[sortperm(burst_freq)]
p = scatter(burst_size[sortperm(burst_freq)], log.(10,y),  zcolor = sort(burst_freq), color = :viridis, ylabel = "log₁₀(bursty promoter noise)", 
        xlims = (-3,4.5), xlabel = "log₁₀(burst size)", label = nothing, colorbar_title = "log₁₀(burst frequency)", markersize = 5, markeralpha = 0.9,
         legendfontsize = 9, fg_legend = :transparent, size = (500,400), dpi = 300)
#l = [-2.5:0.1:4.0;]
#p = plot!(l, l, label = false, linewidth = 4, linealpha = 0.7, linestyle = :dash, color = :lightblue3)
savefig(p, "data/paper_figures/figure_3/model_based/burst_freq_bursty_size.svg")


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
    


decay_rate = vcat(map_sets[1][:,4],map_sets[2][:,4])
x = vcat(s_m[1:2]...)[sortperm(decay_rate)]
y = vcat(s_tx_deg_noise[1:2]...)[sortperm(decay_rate)]

x1 = vcat([minimum(vcat(s_m[1:2]...)) - 2.0],sort(vcat(s_m[1:2]...)))
p = scatter(log.(10,x), log.(10,y), zcolor = sort(10 .^decay_rate), color = :viridis, ylabel = "log₁₀(transcription & decay noise)", 
        xlims = (0.0,3.5),ylims = (-3,1),xlabel = "log₁₀(mean)", label = nothing, colorbar_title = "\ndecay rate", markersize = 5, markeralpha = 0.9,
         legendfontsize = 9, fg_legend = :transparent, right_margin = 4mm, size = (500,400), dpi = 300)
p = plot!(log.(10,x1), log.(10, 1 ./ x1), linewidth = 4, label = false,linealpha = 0.6,linestyle = :dash, color = :lightblue3)

savefig(p, "data/paper_figures/figure_3/model_based/decay_rate_tx_deg_noise.svg")

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






