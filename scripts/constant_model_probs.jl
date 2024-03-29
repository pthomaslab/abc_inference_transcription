function get_model_probs(err::Matrix{Float64},ε::Float64)
    stats = Vector{Float64}(undef, 100)
    alpha::Float64 = 0.95
    n_bootstraps::Int64 = 100
    #l = number of particles below ε for each model separately
    l::Vector{Int64} = [length(filter(x->x<=ε, err[:,i])) for i in 1:size(err)[2]]
    if sum(l) == 0.0
        model_prob::Vector{Float64} = zeros(size(err)[2])
        l_b::Vector{Float64} = zeros(size(err)[2])
        u_b::Vector{Float64} = zeros(size(err)[2])
    else
        model_prob = l ./ sum(l)
        model_vec = Vector{Int64}(undef, Int64(sum(l)))
        model_vec[1:l[1]] .= 1
        model_vec[l[1]+1:l[1]+l[2]] .= 2
        stats = Matrix{Float64}(undef, (n_bootstraps, size(err)[2]))
        for b in 1:n_bootstraps
            sample_idx::Vector{Int64} = Distributions.sample([1:1:sum(l);], sum(l); replace = true)
            sample_model_vec::Vector{Int64} = model_vec[sample_idx]
            sample_l::Vector{Int64} = [length(findall(x->x==i, sample_model_vec)) for i in 1:lastindex(l)]
            stats[b,:] = sample_l ./ sum(sample_l)
        end
        sorted_stats::Matrix{Float64} = sort(stats,dims=1)
        l_b = [quantile(sorted_stats[:,j], 1.0 - alpha) for j in 1:size(stats)[2]]
        u_b = [quantile(sorted_stats[:,j], alpha) for j in 1:size(stats)[2]]
    end
    return model_prob, l_b, u_b
end

f_1 = JDFFile("data/errors/error_const.jdf");
f_2 = JDFFile("data/errors/error_const_const.jdf");

#set acceptance threshold ε
ε = 4.8


for g in 1:length(genes)
    err_mat = Matrix{Float64}(undef,(1000000,2))
    err_mat[:,1] = f_1["x"*string(g)]  
    err_mat[:,2] = f_2["x"*string(g)] 
    err_mat = sort(err_mat,dims=1)
    which = [1,2][[x<=ε for x in err_mat[1,:]]]
    if which == []
        model_prob::Vector{Float64} = zeros(2)
        l_bound::Vector{Float64} = zeros(2)
        u_bound::Vector{Float64} = zeros(2)
    elseif length(which) == 1
        model_prob = zeros(2)
        l_bound = zeros(2)
        u_bound = zeros(2)
        model_prob[which[1]],l_bound[which[1]],u_bound[which[1]] = ones(3)
    else
        model_prob,l_bound,u_bound = get_model_probs(err_mat,ε)
    end
    open("data/model_selection/c_/c_model_prob.txt", "a") do io
        writedlm(io, reshape(model_prob,(1,:)))
    end
    open("data/model_selection/c_/c_l_bound.txt", "a") do io
        writedlm(io, reshape(l_bound,(1,:)))
    end
    open("data/model_selection/c_/c_u_bound.txt", "a") do io
        writedlm(io, reshape(u_bound,(1,:)))
    end
end
