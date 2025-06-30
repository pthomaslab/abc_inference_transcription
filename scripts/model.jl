function get_rate(θ::Vector{Float64}, vary_map::Vector{Any}, scaling::Int64, t_steps::Vector{Float64}, cycle::Float64, t::Float64)::Vector{Float64}
    p = Vector{Float64}(undef,length(vary_map))
    for i in 1:length(vary_map)
        if i==3 
            txn_size = scaling * log10(size_scaling(cycle,t))
        else
            txn_size = 0.0
        end
        if length(vary_map[i]) == 1
            p[i] = θ[vary_map[i]] + txn_size
        else
            for j in 1:length(t_steps)-1
                if mod(t,cycle) > t_steps[j] && mod(t,cycle) <= t_steps[j+1]
                    p[i] =  θ[vary_map[i][j]] + txn_size
                elseif mod(t,cycle) == 0.0
                    p[i] =  θ[vary_map[i][1]] + txn_size
                end
            end
        end
    end
    return p
end


function size_scaling(cycle::Float64,t::Float64)::Float64
    return 1.0 + (mod(t,cycle)/cycle) * (<(t,cycle)) + (==(t,cycle))
end


function get_vary_map(vary_flag::Vector{Int64},n_steps::Int64)
    local keys = []
    local k = 1
    for i in 1:length(vary_flag)
        if vary_flag[i] == 0
            push!(keys,k)
            k += 1
        else
            push!(keys,[k:1:k+n_steps-1;])
            k = k + n_steps
        end
    end
    return keys
end

#=
function get_vary_matrix(vary_map::Vector{Any},n_steps::Int64)
    mat = Matrix{Int64}(undef,(length(vary_map),n_steps))
    for (i,v) in enumerate(vary_map) 
        if length(v) > 1
            mat[i,:] = v
        else
            mat[i,:] = ones(n_steps) .* v
        end
    end
    return mat
end
=#
function labelling(lambda::Float64, texp::Float64, pulse::Float64, t::Float64)::Float64
    if t >= texp && t <= texp + pulse
        return 10^lambda
    else
        return 0.0
    end
end


function momentodes(dy, y, θ, t, cycle, texp, pulse, vary_map, t_steps, scaling)
    p = 10 .^(get_rate(θ[1:end-1], vary_map, scaling, t_steps, cycle, t))
    λ = labelling(θ[end], texp, pulse, t)
    dy .= f(y,p,λ)
end


function f(y::Vector{Float64},p::Vector{Float64},λ::Float64)
    dy = Vector{Float64}(undef,9)
    dy[1] = p[1]*(1-y[1]) - p[2]*y[1]
    dy[2] = (1-λ)*p[3]*y[1] - p[4]*y[2]
    dy[3] = λ*p[3]*y[1] - p[4]*y[3]
    dy[4] = p[2]*y[1] + p[1]*(1-y[1]) - 2*(p[1]+p[2])*y[4]
    dy[5] = p[3]*(1-λ)*y[4] - (p[1]+p[2]+p[4])*y[5]
    dy[6] = p[3]*λ*y[4] - (p[1]+p[2]+p[4])*y[6]
    dy[7] = p[3]*(1-λ)*y[1] + p[4]*y[2] + 2*p[3]*(1-λ)*y[5] - 2*p[4]*y[7]
    dy[8] = p[3]*λ*y[5] + p[3]*(1-λ)*y[6] - 2*p[4]*y[8]
    dy[9] = p[3]*λ*y[1] + p[4]*y[3] + 2*p[3]*λ*y[6] - 2*p[4]*y[9]
    return dy
end

#ode model definition conditioning on labelling periods
function model(paramvals::Vector{Float64}, vary_map::Vector{Any}, scaling::Int64, t_steps::Vector{Float64}, iv::Vector{Float64}, tmin::Float64,
    tmax::Float64, cycle::Float64, texp::Float64, pulse::Float64)
    tspan::Tuple{Float64, Float64} = (tmin, tmax)
    _moment_odes(dμ, μ, θ, t) = momentodes(dμ, μ, θ, t, cycle, texp, pulse, vary_map, t_steps, scaling)
    odeprob = ODEProblem(_moment_odes, iv, tspan, paramvals)
    sol = solve(odeprob, CVODE_BDF())
    return sol.u[end]
end

function periodic_boundary(endpoint::Vector{Float64})
    v::Vector{Float64} = copy(endpoint)
    #mRNA - gene covariances
    gene_covidx::Vector{Int64} = [5,6]
    #mRNA means and covariances
    meanidx::Vector{Int64} = [2,3]
    varidx::Vector{Int64} = [7,9]
    covidx::Vector{Int64} = [8] 
    v[gene_covidx] = endpoint[gene_covidx]/2
    v[varidx] = endpoint[varidx]/4 + endpoint[meanidx]/4
    v[covidx] = endpoint[covidx]/4
    v[meanidx] = endpoint[meanidx]/2
    return v
end


function transient_phase(paramvals::Vector{Float64}, vary_map::Vector{Any}, scaling::Int64, t_steps::Vector{Float64}, iv::Vector{Float64}, cycle::Float64)
    #convergence threshold
    ε::Float64 = 0.01
    #simply ensuring that there is no labelling in transient phase
    texp::Float64 = -1.0
    pulse::Float64 = 0.1
    #select indices of mean, variance & covariance for species of interest - gene & unlabelled mRNA
    check_idx::Vector{Int64} = [1,2,4,5,7]
    τ0::Float64 = 0.0
    τf::Float64 = cycle
    endpoint_2 = Vector{Float64}(undef,length(iv))
    k::Int64 = 0
    convergence::Bool = false
    max_iter::Int64 = 100
    endpoint_1::Vector{Float64} = model(paramvals, vary_map, scaling, t_steps, iv, τ0, τf, cycle, texp, pulse)
    while convergence == false && k <= max_iter
        k += 1
        #set boundary conditions at cell division according to binomial partitioning
        iv = periodic_boundary(endpoint_1)
        endpoint_2 = model(paramvals, vary_map, scaling, t_steps, iv, τ0, τf, cycle, texp, pulse)
        #check if mean and variance of mRNA have reached steady state
        nz_check_idx = check_idx[endpoint_1[check_idx] .> 0.0]
        if sum(abs.((endpoint_1[nz_check_idx] .- endpoint_2[nz_check_idx]) ./ endpoint_1[nz_check_idx]) .<= ε) == length(nz_check_idx)
            convergence = true
        end
        endpoint_1 = endpoint_2
    end
    return periodic_boundary(endpoint_2)
end
 


function trajectories(paramvals::Vector{Float64}, vary_map::Vector{Any}, scaling::Int64, t_steps::Vector{Float64}, iv::Vector{Float64}, age::Float64,
    cycle::Float64, pulse::Float64, t0::Float64, texp::Float64)
    local endpoint::Vector{Float64}
    local nsols::Int64
    #mRNA means and covariances
    meanidx::Vector{Int64} = [2,3]
    allcovidx::Vector{Int64} = [7:9;] 
    if age > 0 && age < cycle
        nsols = floor((age-t0)/cycle) + 1
    else
        nsols = floor((age-t0)/cycle)
    end
    τ = t0
    k = 1
    for k in 1:nsols
        if τ < 0
            tf = τ + cycle
        else
            tf = age
        end
        if τ == t0
            endpoint = model(paramvals, vary_map, scaling, t_steps, iv, τ, tf, cycle, texp, pulse)
        else
            endpoint = model(paramvals, vary_map, scaling, t_steps, periodic_boundary(endpoint), τ, tf, cycle, texp, pulse)
        end
        τ = tf
    end
    return endpoint[meanidx], endpoint[allcovidx]
end

function syntheticdata(paramvals::Vector{Float64}, vary_map::Vector{Any}, scaling::Int64, n_steps::Int64, iv::Vector{Float64}, agevec::Vector{Float64},
    cycle::Float64, pulse::Float64, chase::Float64, t0::Float64)
    local texp::Float64
    mean_data = Matrix{Float64}(undef,length(agevec),2)
    cov_data = Matrix{Float64}(undef,length(agevec),3)
    t_steps::Vector{Float64} = [0.0:cycle/n_steps:cycle;]
    for (i,age) in enumerate(agevec)
        texp = age - pulse - chase
        mean_data[i,:], cov_data[i,:] = trajectories(paramvals,vary_map,scaling,t_steps,iv,age,cycle,pulse,t0,texp)
    end
    return hcat(mean_data, cov_data)
end



function get_synthetic_data(θ::Vector{Float64}, vary_map::Vector{Any}, scaling::Int64, n_steps::Int64, iv::Vector{Float64}, agevec::Vector{Float64},
    cycle::Float64, pulse::Float64, chase::Float64, t0::Float64, downsampling::Bool, betas::Vector{Float64}, age_clusters::Vector{Int64})
    s_data::Matrix{Float64} = syntheticdata(θ, vary_map, scaling, n_steps, iv, agevec, cycle, pulse, chase, t0)
    if downsampling == true
        return downsample(s_data, betas, age_clusters)
    else
        return s_data
    end
end

function get_steady_state_iv(paramvals::Vector{Float64}, vary_map::Vector{Any}, scaling::Int64, n_steps::Int64, iv::Vector{Float64}, cycle::Float64)
    t_steps::Vector{Float64} = [0.0:cycle/n_steps:cycle;]
    #run transient dynamics to determine steady state
    ss_iv::Vector{Float64} = transient_phase(paramvals,vary_map,scaling,t_steps,iv,cycle)
    return ss_iv
end 


function get_total_synthetic_data(θ::Vector{Float64}, vary_map::Vector{Any}, scaling::Int64, n_steps::Int64, iv::Vector{Float64}, agevec::Vector{Float64},
    cycle::Float64, pulse::Float64, chase::Float64, t0::Float64, downsampling::Bool, betas::Vector{Float64}, clusters::Vector{Int64})
    s_data::Matrix{Float64} = syntheticdata(θ, vary_map, scaling, n_steps, iv, agevec, cycle, pulse, chase, t0)
    ε::Vector{Float64} = [0.0 * (>(x,0.0)) + 0.0001 *(==(x,0.0)) for x in s_data[:,1] + s_data[:,2]]
    if downsampling == true
        return downsample_total(hcat(s_data[:,1] + s_data[:,2] + ε, s_data[:,3] + s_data[:,4] + s_data[:,5]), betas, clusters)
    else
        return hcat(s_data[:,1] + s_data[:,2], (s_data[:,3] + s_data[:,4] + s_data[:,5]) ./ (s_data[:,1] + s_data[:,2] + ε))
    end
end


function downsample(data::Matrix{Float64}, betas::Vector{Float64}, clusters::Vector{Int64})
    meanidx::Vector{Int64} = [1,2]
    varidx::Vector{Int64} = [3,5]
    covidx::Vector{Int64} = [4]
    new_data = Matrix{Float64}(undef,size(data))
    mean_betas = Vector{Float64}(undef, length(unique(clusters)))
    m2_betas = Vector{Float64}(undef, length(unique(clusters)))  #second moments of betas
    var_betas = Vector{Float64}(undef, length(unique(clusters)))
    for c in sort(unique(clusters))
        idx = findall(x->x==c, clusters)
        mean_betas[c] = mean(betas[idx])
        m2_betas[c] = mean(betas[idx].^2)
        var_betas[c] = var(betas[idx])
    end
    new_data[:,meanidx] = data[:,meanidx] .* mean_betas
    new_data[:,varidx] = ((mean_betas .- m2_betas) .* data[:,meanidx]) .+ (var_betas .* (data[:,meanidx].^2 .+ data[:,varidx])) .+ (mean_betas.^2 .* data[:,varidx])
    new_data[:,covidx] = (var_betas .* (data[:,meanidx[1]] .* data[:,meanidx[2]] .+ data[:,covidx])) .+ (mean_betas.^2 .* data[:,covidx])
    return new_data
end

