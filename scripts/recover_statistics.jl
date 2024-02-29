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

model_name = ["const","const_const","kon","alpha","gamma"][m];
maps = readdlm("data/posterior_estimates/map_sets_"*model_name*".txt");

condition_id = hcat([1/4,1/2,3/4,1,2,3,22,22,22,22,22],[0,0,0,0,0,0,0,1,2,4,6]);
cond_idx = [1:11;];
n_steps = length(τ_);
iv = zeros(9); 
iv[1] = 1/2;
cycle = 20.0;
t0 = -3 * cycle;
agevec = τ_ .* cycle;                   
pulsevec = condition_id[cond_idx,1];
chasevec = condition_id[cond_idx,2];    
downsampling = false;

vary_flag = [[0,0,0,0],[0,0,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,1]][m];
vary_map = get_vary_map(vary_flag,n_steps);
scaling = 1 * !=(m,2);
###################################################################################################

id_labels = ["pulse_15", "pulse_30", "pulse_45", "pulse_60", "pulse_120", "pulse_180", "chase_0", "chase_60", "chase_120", "chase_240", "chase_360"];

@time for i in 1:size(maps)[1]
    s::Vector{Matrix{Float64}} = run_part_sim(maps[i,:],iv,agevec,cycle,pulsevec,chasevec,t0,vary_map,scaling,n_steps,downsampling,betas,age)
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


 
