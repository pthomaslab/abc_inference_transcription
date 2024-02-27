function read_all_data(dir::String,ext::String)
    uu_df = CSV.read(dir*"scEU_unlabeled_unspliced"*ext, DataFrame)
    us_df = CSV.read(dir*"scEU_unlabeled_spliced"*ext, DataFrame)
    lu_df = CSV.read(dir*"scEU_labeled_unspliced"*ext, DataFrame)
    ls_df = CSV.read(dir*"scEU_labeled_spliced"*ext, DataFrame)
    cells_df = CSV.read(dir*"scEU_cell_cycle"*ext, DataFrame)
    uu_data = transpose(Matrix(select(uu_df, Not(:Gene_Id))))
    us_data = transpose(Matrix(select(us_df, Not(:Gene_Id))))
    lu_data = transpose(Matrix(select(lu_df, Not(:Gene_Id))))
    ls_data = transpose(Matrix(select(ls_df, Not(:Gene_Id))))
    cells = Matrix(select(cells_df, Not(:Cell_Id)))
    cell_cycle = parse.(Float64, cells[end,:])
    rfp = parse.(Float64, cells[5,:])
    gfp = parse.(Float64, cells[4,:])
    experiment = Vector{Int64}(undef, length(cell_cycle))
    conditions = ["Pulse_dmso", "Pulse_15", "Pulse_30", "Pulse_45",
                  "Pulse_60", "Pulse_120", "Pulse_180", "Chase_0",
                  "Chase_60", "Chase_120", "Chase_240", "Chase_360", "Chase_dmso"]
    for i in 0:length(conditions)-1
        idx = findall(x->x==conditions[i+1], cells[2,:])
        experiment[idx] .= i
    end
    return copy(uu_data), copy(us_data), copy(lu_data), copy(ls_data), cell_cycle, rfp, gfp, experiment
end

function get_gene_names(dir::String, id::Vector{String})
    local gene_name::String
    local vec::Vector{FASTX.FASTA.Record}
    for ensembl in id
        vec = fetchseq(ensembl)
        if length(vec)>0
            descr = description(vec[1])
            pos_1 = findfirst("GN=",descr)
            pos_2 = findfirst("PE=",descr)
            if (!=(pos_1,nothing)) && (!=(pos_2,nothing))
                gene_name = descr[pos_1[end]+1:pos_2[1]-2]
            else
                gene_name = "nan"
            end
        else
            gene_name = "nan"
        end
        open("dir/gene_names.txt", "a") do io
           println(io, gene_name)
        end
    end
end


function age_clusters(theta::Vector{Float64}, n_clusters::Int64)
    clusters = Vector{Int64}(undef, length(theta))
    τ_::Vector{Float64} = [(2k+1)/(2*n_clusters) for k in 0:n_clusters-1]
    for k in 0:n_clusters-1
        clusters[findall(x-> x > k/n_clusters && x <= (k+1)/n_clusters, theta)] .= k+1
    end
    return clusters, τ_
end

uu_data, us_data, lu_data, ls_data, theta, rfp, gfp, experiment = read_all_data(dir,".csv");

#total counts matrix:
total_data = uu_data + us_data + lu_data + ls_data;

########## Pre-process data ##########

#cluster cells with respect to cell cycle position
n_age_clusters = 5;
age, τ_ = age_clusters(theta, n_age_clusters);

#distinguishing between pulse-treated and chase-treated cells
pulse_idx = findall(x->x>=0 && x<=6, experiment);
chase_idx = findall(x->x>=7, experiment);

#empirical distributions of cells across the 5 cell cycle stages and 11 labelling conditions
cells_per_id = [findall(x->x==e,experiment) for e in sort(unique(experiment))[2:end-1]];
cells_per_age = [findall(x->x==τ,age) for τ in sort(unique(age))];
cells_age_id = [[intersect(cells_age,cells_id) for cells_age in cells_per_age] for cells_id in cells_per_id];

n_cells_id = sum(length.(cells_per_id));
age_id_distribution = [length.(cells) ./ n_cells_id for cells in cells_age_id]
age_dist_pulse = length.([intersect(cells,pulse_idx) for cells in cells_per_age]) ./ length(pulse_idx);
age_dist_chase = length.([intersect(cells,chase_idx) for cells in cells_per_age]) ./ length(chase_idx);

