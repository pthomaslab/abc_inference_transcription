using DifferentialEquations
using Sundials
using DelimitedFiles
using CSV
using DataFrames
using Statistics
using Distributions


include("load_data.jl")
include("model.jl")

function fix_params(vary_map::Vector{Any}, N::Int64)
    kon::Matrix{Float64} = rand(Uniform(-3.0, 3.0), (N,length(vary_map[1])))
    #burstsize::Matrix{Float64} = rand(Uniform(-1.0, 2.7), (N,length(vary_map[2])))
    koff::Matrix{Float64} = rand(Uniform(-3.0, 3.0), (N,length(vary_map[2])))
    alpha::Matrix{Float64} = rand(Uniform(-3.0, 3.0), (N,length(vary_map[3])))
    gamma::Matrix{Float64} = rand(Uniform(-3.0, 2.0), (N,length(vary_map[4])))
    lambda::Matrix{Float64} = rand(Uniform(-0.7, 0.0), (N,1))
    return hcat(kon, koff, alpha, gamma, lambda)
end

function abc_sim(θ::Vector{Float64}, iv::Vector{Float64}, agevec::Vector{Float64}, cycle::Float64, pulsevec::Vector{Float64}, chasevec::Vector{Float64}, t0::Float64, vary_map::Vector{Any}, scaling::Int64, n_steps::Int64,
    downsampling::Bool, betas::Vector{Float64}, age::Vector{Int64}, pulse_idx::Vector{Int64}, chase_idx::Vector{Int64}, age_dist::Matrix{Float64}, rate_name::String, n::Int64)
    local io::IOStream
    s = Matrix{Float64}(undef,(length(agevec),5))
    s_pulse = Matrix{Float64}(undef,(length(agevec),2))
    s_chase = Matrix{Float64}(undef,(length(agevec),2))
    s_ratios = Vector{Float64}(undef,length(pulsevec))
    s_mean_corr = Vector{Float64}(undef,length(pulsevec))    #mean covariance over age distribution
    s_corr_mean = Vector{Float64}(undef,length(pulsevec))
    ss_iv::Vector{Float64} = get_steady_state_iv(θ,vary_map,scaling,n_steps,iv,cycle)
    for j in 1:6
        s = get_synthetic_data(θ,vary_map,scaling,n_steps,ss_iv,agevec,cycle,pulsevec[j],chasevec[j],t0,downsampling,betas[pulse_idx],age[pulse_idx])
        s_ratios[j] = sum(age_dist[:,j] .* s[:,2]) ./ (sum(age_dist[:,j] .* s[:,1]) + sum(age_dist[:,j] .* s[:,2]))
        stds::Float64 = sqrt(abs((sum(age_dist[:,j] .* s[:,3]) + weighted_cov(s[:,1],s[:,1],age_dist[:,j])) * (sum(age_dist[:,j] .* s[:,5]) + weighted_cov(s[:,2],s[:,2],age_dist[:,j]))))
        s_mean_corr[j] = sum(age_dist[:,j] .* s[:,4]) / stds
        s_corr_mean[j] = weighted_cov(s[:,1],s[:,2],age_dist[:,j]) / stds
        if j == 6
            ε::Vector{Float64} = [0.0 * (>(x,0.0)) + 0.0001 *(==(x,0.0)) for x in s[:,1] + s[:,2]]
            s_pulse[:,1] = s[:,1] + s[:,2]
            s_pulse[:,2] = (s[:,3] + 2*s[:,4] + s[:,5]) ./ (s[:,1] + s[:,2] + ε)
        end
    end
    for j in 7:length(pulsevec)
        s = get_synthetic_data(θ,vary_map,scaling,n_steps,ss_iv,agevec,cycle,pulsevec[j],chasevec[j],t0,downsampling,betas[chase_idx],age[chase_idx])
        s_ratios[j] = sum(age_dist[:,j] .* s[:,2]) ./ (sum(age_dist[:,j] .* s[:,1]) + sum(age_dist[:,j] .* s[:,2]))
        stds::Float64 = sqrt(abs((sum(age_dist[:,j] .* s[:,3]) + weighted_cov(s[:,1],s[:,1],age_dist[:,j])) * (sum(age_dist[:,j] .* s[:,5]) + weighted_cov(s[:,2],s[:,2],age_dist[:,j]))))
        s_mean_corr[j] = sum(age_dist[:,j] .* s[:,4]) / stds
        s_corr_mean[j] = weighted_cov(s[:,1],s[:,2],age_dist[:,j]) / stds
        if j == 7
            ε::Vector{Float64} = [0.0 * (>(x,0.0)) + 0.0001 *(==(x,0.0)) for x in s[:,1] + s[:,2]]
            s_chase[:,1] = s[:,1] + s[:,2]
            s_chase[:,2] = (s[:,3] + 2*s[:,4] + s[:,5]) ./ (s[:,1] + s[:,2] + ε)
        end
    end
    open("Julia/new_simulations/"*rate_name*"/s_pulse_"*rate_name*"_$n.txt", "a") do io
        writedlm(io, transpose(s_pulse))
    end
    open("Julia/new_simulations/"*rate_name*"/s_chase_"*rate_name*"_$n.txt", "a") do io
        writedlm(io, transpose(s_chase))
    end
    open("Julia/new_simulations/"*rate_name*"/s_ratios_"*rate_name*"_$n.txt", "a") do io
        writedlm(io, transpose(s_ratios))
    end
    open("Julia/new_simulations/"*rate_name*"/s_mean_corr_"*rate_name*"_$n.txt", "a") do io
        writedlm(io, transpose(s_mean_corr))
    end
    open("Julia/new_simulations/"*rate_name*"/s_corr_mean_"*rate_name*"_$n.txt", "a") do io
        writedlm(io, transpose(s_corr_mean))
    end
end


uu_data, us_data, lu_data, ls_data, theta, rfp, gfp, experiment, gene_id = read_all_data("data/",".csv")
total_data = uu_data + us_data + lu_data + ls_data

n_clusters = 5
age, age_idx, mean_age, τ_ = age_clusters(theta, n_clusters, "equidistant")
#age_distribution = length.(age_idx)/ncells

cells_per_age = [findall(x->x==τ,age) for τ in sort(unique(age))]
cells_per_id = [findall(x->x==e,experiment) for e in sort(unique(experiment))[2:end-1]]
cells_age_id = [[intersect(cells_age,cells_id) for cells_age in cells_per_age] for cells_id in cells_per_id]

age_id_distribution = hcat([length.(cells) / sum(length.(cells)) for cells in cells_age_id]...)


pulse_idx = findall(x->x<=6, experiment)
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


#matrix of experimental conditions: 1st col. - pulse, 2nd col. - chase
condition_id = hcat([1/4,1/2,3/4,1,2,3,22,22,22,22,22], [0,0,0,0,0,0,0,1,2,4,6])
cond_idx = [1:11;]

##################################################################################################
#####################################  model configuration #######################################
##################################################################################################
iv = zeros(9)   
iv[2] = 1/2
cycle = 20.0
t0 = -3 * cycle
agevec = τ_ .* cycle                  
pulsevec = condition_id[cond_idx,1]
chasevec = condition_id[cond_idx,2]      

downsampling = true

###################################################################################################
# m: modelling hypothesis index  --> {1: constant rates - scaling with size
#                                     2: constant rates - non-scaling
#                                     3: varying burst frequency - scaling with size
#                                     4: varying burst size - scaling with size
#                                     5: varying decay rate - scaling with size }
#                                       
 

m = 1   #2, 3, 4, 5
model_name = ["const","const_const","kon","alpha","gamma"][m]
vary_flag = [[0,0,0,0],[0,0,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,1]][m]
vary_map = get_vary_map(vary_flag,n_clusters)
scaling = 1 * (!=(m,2))

#number of simulations in current run
n_trials = 250000

#submit file on high-performance computing cluster
submit = 1

@time for i in 1:n_trials
    open("Julia/simulations/"*model_name*"/progress_"*model_name*"_$submit.txt", "a") do io
        writedlm(io, i)
    end
    θ = fix_params(vary_map,1)
    open("Julia/simulations/"*model_name*"/sets_"*model_name*"_$submit.txt", "a") do io
        writedlm(io, θ)
    end
    abc_sim(θ[1,:],iv,agevec,cycle,pulsevec,chasevec,t0,vary_map,scaling,n_clusters,downsampling,betas,age,pulse_idx,chase_idx,age_id_distribution,model_name,submit)
end
