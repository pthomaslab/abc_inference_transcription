include("load_data.jl")
include("model.jl")


function run_part_sim(θ::Vector{Float64}, iv::Vector{Float64}, agevec::Vector{Float64}, cycle::Float64, pulsevec::Vector{Float64}, chasevec::Vector{Float64}, t0::Float64, vary_map::Vector{Any},
    scaling::Int64, n_steps::Int64, downsampling::Bool, betas::Vector{Float64}, age::Vector{Int64})
    s = Vector{Matrix{Float64}}(undef,length(pulsevec))
    ss_iv::Vector{Float64} = get_steady_state_iv(θ,vary_map,scaling,n_steps,iv,cycle) 
    for i in 1:length(pulsevec) 
        s[i] = Matrix{Float64}(undef,(length(agevec),5)) 
        s[i] = get_synthetic_data(θ,vary_map,scaling,n_steps,ss_iv,agevec,cycle,pulsevec[i],chasevec[i],t0,downsampling,betas,age)
    end
    return s
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


################################### Load data ##########################################
uu_data, us_data, lu_data, ls_data, theta, rfp, gfp, experiment, gene_id = read_all_data("data/",".csv")
total_data = uu_data + us_data + lu_data + ls_data

ncells, ngenes = size(us_data)

n_clusters = 5
age, age_idx, mean_age, τ_ = age_clusters(theta, n_clusters, "equidistant")
age_distribution = length.(age_idx)/ncells

cells_per_id = [findall(x->x==e,experiment) for e in sort(unique(experiment))[2:end-1]]
cells_per_age = [findall(x->x==τ,age) for τ in sort(unique(age))]
cells_age_id = [[intersect(cells_age,cells_id) for cells_age in cells_per_age] for cells_id in cells_per_id]

pulse_idx = findall(x->x>=0 && x<=6, experiment)
tc_pulse = sum(total_data[pulse_idx,:],dims=2)[:,1]
tc_pulse_age = [sum(total_data[intersect(cells,pulse_idx),:],dims=2)[:,1] for cells in cells_per_age]
atc_pulse_age = mean.(tc_pulse_age)

chase_idx = findall(x->x>=7, experiment)
tc_chase = sum(total_data[chase_idx,:],dims=2)[:,1]
tc_chase_age = [sum(total_data[intersect(cells,chase_idx),:],dims=2)[:,1] for cells in cells_per_age]
atc_chase_age = mean.(tc_chase_age)

mean_beta_pulse = 0.1
ratios_pulse = [tc_pulse[i] / atc_pulse_age[age[pulse_idx[i]]] for i in 1:lastindex(pulse_idx)]
betas_pulse = mean_beta_pulse .* ratios_pulse

mean_beta_chase = 0.2
ratios_chase = [tc_chase[i] / atc_chase_age[age[chase_idx[i]]] for i in 1:lastindex(chase_idx)]
betas_chase = mean_beta_chase .* ratios_chase

betas = Vector{Float64}(undef,ncells)
betas[pulse_idx] = betas_pulse
betas[chase_idx] = betas_chase


m = 5
model_name = ["const","const_const","kon","alpha","gamma"][m]

n_groups_sims = 1
sets = load_param_sets("data/large_scale_simulations/",m,n_groups_sims)

maps = Int64.(readdlm("data/posteriors/map_gene_per_model.txt")[:,m])

sel_genes = Int64.(readdlm("data/model_selection/"*model_name*"_genes.txt")[:,1])

map_sets = sets[maps[non_sel_genes],:]

##################################################################################################
################################# Define model framework #########################################
##################################################################################################
#matrix of experimental conditions: 1st col. - pulse, 2nd col. - chase
condition_id = hcat([1/4,1/2,3/4,1,2,3,22,22,22,22,22],[0,0,0,0,0,0,0,1,2,4,6])
cond_idx = [1:11;]

n_steps = length(τ_);
iv = zeros(9)   #zeros(14)   # zeros(27)
iv[1] = 1/2
cycle = 20.0
t0 = -3 * cycle
agevec = τ_ .* cycle                   #[log(s_c[c]) for c in 1:10] ./ growth_rate
pulsevec = condition_id[cond_idx,1]
chasevec = condition_id[cond_idx,2]      #[0.05:0.1:1;] .* cycle
downsampling = false

vary_flag = [[0,0,0,0],[0,0,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,1]][m]
vary_map = get_vary_map(vary_flag,n_steps)
scaling = 1 * !=(m,2)
###################################################################################################


id_labels = ["pulse_15", "pulse_30", "pulse_45",
                "pulse_60", "pulse_120", "pulse_180", "chase_0",
                         "chase_60", "chase_120", "chase_240", "chase_360"]


@time for i in 1:length(sel_genes)
    s::Vector{Matrix{Float64}} = run_part_sim(map_sets[i,:],iv,agevec,cycle,pulsevec,chasevec,t0,vary_map,scaling,n_steps,downsampling,betas,age)
    for (k,id) in enumerate(id_labels)
        open("data/recovered_statistics/"*model_name*"/"*id*"/mean_u.txt", "a") do io
            writedlm(io, reshape(s[k][:,1],(1,:)))
        end 
        open("data/recovered_statistics/"*model_name*"/"*id*"/mean_l.txt", "a") do io
            writedlm(io, reshape(s[k][:,2],(1,:)))
        end 
        open("data/recovered_statistics/"*model_name*"/"*id*"/var_u.txt", "a") do io
            writedlm(io, reshape(s[k][:,3],(1,:)))
        end 
        open("data/recovered_statistics/"*model_name*"/"*id*"/cov_ul.txt", "a") do io
            writedlm(io, reshape(s[k][:,4],(1,:)))
        end 
        open("data/recovered_statistics/"*model_name*"/"*id*"/var_l.txt", "a") do io
            writedlm(io, reshape(s[k][:,5],(1,:)))
        end 
    end
end


 