function filter_genes(lu_data::Matrix{Float64},ls_data::Matrix{Float64},total_data::Matrix{Float64},cells_per_age::Vector{Vector{Int64}})
    #1st filtering step
    #select genes with at least 500 labelled counts
    tc_l = sum(lu_data .+ ls_data,dims=1)[1,:]
    all_genes = findall(x->x>=500.0, tc_l)
    #2nd filtering step
    s = sum(total_data[:,all_genes],dims=1)[1,:]
    min_tc = size(total_data)[1];
    age_freq = length.(cells_per_age);
    gene_tc_per_age = [sum(total_data[cells,all_genes],dims=1)[1,:] for cells in cells_per_age];
    #select genes that have on average at least *min_count* counts in at least one cell-cycle phase
    sel_1 = [];
    min_count = 1;
    for g in 1:lastindex(all_genes)
        if sum([gene_tc_per_age[i][g] >= min_count*age_freq[i] for i in 1:n_clusters]) >= 1
            push!(sel_1,g)
        end
    end
    sel_2 = [];
    for g in 1:lastindex(all_genes)
        if s[g] >= min_tc
            push!(sel_2,g)
        end
    end
    return all_genes[Int64.(union(sel_1,sel_2))]
end
    
sel = filter_genes(lu_data,ls_data,total_data,cells_per_age);

writedlm("data/selected_genes.txt",sel);

