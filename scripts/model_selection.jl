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

function get_n_particles(posterior::Vector{Vector{Vector{Int64}}})
    n_particles::Vector{Vector{Int64}} = [length.(p) for p in posterior]
    for (i,p_) in enumerate(posterior)
        zeros_ = findall(x->x==[0],p_)
        n_particles[i][zeros_] .= 0
    end
    return n_particles
end

ε = 4.8;

dir = "data/model_selection/";

model_prob = readdlm(dir*"all/model_prob.txt")
l_b = readdlm(dir*"all/l_bound.txt") 
u_b = readdlm(dir*"all/u_bound.txt") 

c_prob = readdlm(dir*"c_/c_model_prob.txt")
c_ub = readdlm(dir*"c_/c_u_bound.txt") 
c_lb = readdlm(dir*"c_/c_l_bound.txt")

nc_prob = readdlm(dir*"nc_/nc_model_prob.txt")
nc_ub = readdlm(dir*"nc_/nc_u_bound.txt")
nc_lb = readdlm(dir*"nc_/nc_l_bound.txt")

 
posterior = Vector{Vector{Vector{Int64}}}(undef,length(model_names));
for (i,name) in enumerate(model_names)
    file = open("data/posteriors/particles_"*name*".txt", "r")
    posterior[i] = Vector{Vector{Int64}}(undef,length(genes))
    for (j,line) in enumerate(eachline(file))
        posterior[i][j] = parse.(Int64,split(line,"\t")) 
    end
    close(file)
end  

n_particles_model = get_n_particles(posterior);
n_particles_mat = hcat(n_particles_model...);
n_particles = sum(n_particles_mat, dims=2)[:,1];

#############################    first model selection    ############################# 
gene_vec = [1:length(genes);];
prior_odds = 2/3;
thres = prior_odds / (1 + prior_odds);
min_particles = 5;
c_, nc_, und_, nf_ = get_c_nc_genes(model_prob,l_b,u_b,gene_vec,n_particles,min_particles,thres);


writedlm(dir*"all/c_gene_ids.txt",gene_ids[c_]);
writedlm(dir*"all/nc_gene_ids.txt",gene_ids[nc_]);


freq_c = length(c_);
freq_nc = length(nc_);
freq_und = length(und_);
freq_nf = length(nf_);
model_classes = ["constant","non-constant","undetermined"];
col = [:skyblue4, :darkgoldenrod2, :azure]; 
b = bar([1,2,3],[freq_c,freq_nc,freq_und+freq_nf] ./ length(gene_vec),labels=false,xticks = ([1,2,3],model_classes), color = col[1:3], fillalpha = 0.95,
            ylabel = "relative frequency", xtickfontsize = 11, size = (450,200),leftmargin=2mm, topmargin=5mm);
b = annotate!([1,2,3],[freq_c,freq_nc,freq_und+freq_nf]./ length(gene_vec),"n = " .*string.([freq_c,freq_nc,freq_und+freq_nf]),:bottom);
savefig(b,"data/paper_figures/figure_2/const_vs_non_genes_freq.svg")


#######################    (conditional) constant model selection    #######################

c_genes = cluster_c_genes(c_prob,c_lb,c_ub,c_);

writedlm(dir*"c_/const_scaling_ids.txt",gene_ids[c_genes[1]]);
writedlm(dir*"c_/const_non_scaling_ids.txt",gene_ids[c_genes[2]]);

writedlm("const_genes.txt",c_genes[1]);
writedlm("const_const_genes.txt",c_genes[2]);


freq_c_genes = length.(c_genes);
col = [:skyblue4, :lightcyan3, :white];
model_classes = ["scaling","non-scaling","undetermined"];
b = bar([1,2,3],freq_c_genes ./ length(c_) ,labels=false,xticks = ([1,2,3],model_classes), color = col[1:3], 
            ylabel = "relative frequency", fillalpha = 1.0, xtickfontsize = 11, size = (450,200),leftmargin=2mm, topmargin=5mm);
b = annotate!([1,2,3],freq_c_genes ./ length(c_),"n = " .* string.(freq_c_genes),:bottom);
savefig(b,"data/paper_figures/figure_2/const_genes_freq.svg");


#######################    (conditional) non-constant model selection    #######################
nc_sel, nc_genes = cluster_nc_genes(nc_prob,nc_lb,nc_ub,nc_);

writedlm(dir*"nc_/kon_gene_idx.txt",gene_ids[nc_genes[1][1]]);
writedlm(dir*"nc_/alpha_gene_ids.txt",gene_ids[nc_genes[1][2]]);
writedlm(dir*"nc_/gamma_gene_ids.txt",gene_ids[nc_genes[1][3]]);

writedlm("kon_genes.txt",nc_genes[1][1]);
writedlm("alpha_genes.txt",nc_genes[1][2]);
writedlm("gamma_genes.txt",nc_genes[1][3]);


freq_nc_genes = [length.(nc) for nc in nc_genes]
plot_freq = vcat([freq_nc_genes[1],freq_nc_genes[2][1:2],[freq_nc_genes[2][3] + sum(freq_nc_genes[3]) + sum(freq_nc_genes[4])]]...) ./ length(nc_)
nc_labels = ["burst frequency","burst size","decay rate","burst frequency or burst size","burst frequency or decay rate","undetermined"];
col = [:skyblue4, :seagreen, :gold, :azure, :azure, :azure];  
b = bar([1:length(plot_freq);], plot_freq, label = false, xticks = ([1:length(plot_freq);],nc_labels), fillalpha = 0.8, xtickfontsize = 11, xrotation = 20, 
    top_margin = 5mm, left_margin = 3mm, bottom_margin = 15mm,color = col, ylabel = "relative frequency");
b = annotate!([1:length(plot_freq);],0.015 .+ plot_freq,"n = " .* string.(Int64.(plot_freq .* length(nc_))),:bottom, size = (650,250));
savefig(b,"data/paper_figures/figure_2/non_const_genes_freq.svg")



##################### summarise model selection results in a DataFrame ##########################
model_sel_df = DataFrame(i = Int64[], gene_name = String[], model = Int64[], model_name = String[]);

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
for nc in vcat(nc_genes[2]...,nc_genes[3]...,nc_genes[4]...)
    push!(model_sel_df, [nc, gene_names[nc], 0, "non_const_undet"])
end
for g in und_
    push!(model_sel_df, [g, gene_names[g], 0, "undet"])
end
for g in nf_
    push!(model_sel_df, [g, gene_names[g], 0, "not_fitted"])
end

CSV.write("data/model_selection/model_selection_results.txt", model_sel_df);

