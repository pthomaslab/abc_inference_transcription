function get_mean_subset(data::Matrix{Float64})
    sub = Matrix{Float64}(undef,(Int64(size(data)[1]/2),size(data)[2]))
    for i in 1:size(sub)[1]
        sub[i,:] = data[2*i-1,:]
    end
    return sub
end

function get_ff_subset(data::Matrix{Float64})
    sub = Matrix{Float64}(undef,(Int64(size(data)[1]/2),size(data)[2]))
    for i in 1:size(sub)[1]
        sub[i,:] = data[2*i,:]
    end
    return sub
end
    
function load_s_data(path::String, model_name::String, ext::String)
    s_pulse::Matrix{Float64} = readdlm(path*model_name*"/s_pulse_"*model_name*ext)
    s_chase::Matrix{Float64} = readdlm(path*model_name*"/s_chase_"*model_name*ext)
    s_pulse_mean = get_mean_subset(s_pulse)
    s_pulse_ff = get_ff_subset(s_pulse)
    s_chase_mean = get_mean_subset(s_chase)
    s_chase_ff = get_ff_subset(s_chase)
    s_ratios::Matrix{Float64} = readdlm(path*model_name*"/s_ratios_"*model_name*ext)
    s_mean_corr::Matrix{Float64} = readdlm(path*model_name*"/s_mean_corr_"*model_name*ext)
    s_corr_mean::Matrix{Float64} = readdlm(path*model_name*"/s_corr_mean_"*model_name*ext)
    return s_pulse_mean,s_pulse_ff,s_chase_mean,s_chase_ff,s_ratios,s_mean_corr,s_corr_mean
end

function nlsqerror_part(data::Vector{Float64},se::Vector{Float64},s_data::Vector{Float64},n_summary_stats::Int64)::Float64
    local ε::Float64
    σ::Float64 = 0.1
    err::Float64 = 0.0
    for i in 1:length(data)
        if se[i] + data[i] != 0.0
            ε = 0.0
        else
            ε = 0.0001
        end
        err += (data[i] - s_data[i])^2 / (se[i]^2 + σ^2 * data[i]^2 + ε)
    end
    return err / n_summary_stats
end

function compute_trunc_errors(pulse_mean::Matrix{Float64},pulse_mean_se::Matrix{Float64},pulse_ff::Matrix{Float64},pulse_ff_se::Matrix{Float64},
    chase_mean::Matrix{Float64},chase_mean_se::Matrix{Float64},chase_ff::Matrix{Float64},chase_ff_se::Matrix{Float64},
    ratio_data::Matrix{Float64},ratio_se::Matrix{Float64},mean_corr_data::Matrix{Float64},mean_corr_se::Matrix{Float64}, 
    corr_mean_data::Matrix{Float64},corr_mean_se::Matrix{Float64},s_pulse_mean::Matrix{Float64},s_pulse_ff::Matrix{Float64},
    s_chase_mean::Matrix{Float64},s_chase_ff::Matrix{Float64},s_ratios::Matrix{Float64},s_mean_corr::Matrix{Float64},s_corr_mean::Matrix{Float64},model_name::String)
    local io::IOStream
    n_summary_stats::Int64 = 4*size(pulse_data)[2] + 3*size(ratio_data)[2]
    for i in 1:size(s_ratios)[1]
        err = Vector{Float64}(undef,size(ratio_data)[1])
        for j in 1:size(ratio_data)[1]
            data_::Vector{Array{Float64}} = [pulse_mean[j,:],pulse_ff[j,:],chase_mean[j,:],chase_ff[j,:],ratio_data[j,:],mean_corr_data[j,:],corr_mean_data[j,:]]
            se_::Vector{Array{Float64}} = [pulse_mean_se[j,:],pulse_ff_se[j,:],chase_mean_se[j,:],chase_ff_se[j,:],ratio_se[j,:],mean_corr_se[j,:],corr_mean_se[j,:]]
            s_data_::Vector{Array{Float64}} = [s_pulse_mean[i,:],s_pulse_ff[i,:],s_chase_mean[i,:],s_chase_ff[i,:],s_ratios[i,:],s_mean_corr[i,:],s_corr_mean[i,:]]
            err[j] = 0.0
            for l in 1:7
                err[j] += nlsqerror_part(data_[l],se_[l],s_data_[l],n_summary_stats)
            end
            if err[j] > 10.0
                err[j] = 10.0
            end
        end
        open("/rds/general/user/dv19/ephemeral/errors/error_"*model_name*".txt", "a") do io
            DelimitedFiles.writedlm(io, reshape(err,1,:))
        end
    end
end

################################################################################################################
#load model-generated summary statistics
s_pulse_mean,s_pulse_ff,s_chase_mean,s_chase_ff,s_ratios,s_mean_corr,s_corr_mean = load_s_data("data/simulations/",model_name,".txt");

#compute errors
@time (compute_trunc_errors(pulse_mean,pulse_mean_se,pulse_ff,pulse_ff_se,chase_mean,chase_mean_se,chase_ff,chase_ff_se,ratio_data,ratio_se,
    mean_corr_data,mean_corr_se,corr_mean_data,corr_mean_se,s_pulse_mean,s_pulse_ff,s_chase_mean,s_chase_ff,s_ratios,s_mean_corr,s_corr_mean,model_name));



