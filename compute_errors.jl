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


function load_summary_stats(path::String,ext::String)
    pulse_data::Matrix{Float64} = readdlm(path*"pulse_data"*ext)
    pulse_se::Matrix{Float64} = readdlm(path*"pulse_se"*ext)
    chase_data::Matrix{Float64} = readdlm(path*"chase_data"*ext)
    chase_se::Matrix{Float64} = readdlm(path*"chase_se"*ext)
    ratio_data::Matrix{Float64} = readdlm(path*"ratio_data"*ext)
    ratio_se::Matrix{Float64} = readdlm(path*"ratio_se"*ext)
    mean_corr_data::Matrix{Float64} = readdlm(path*"mean_corr_data"*ext)
    mean_corr_se::Matrix{Float64} = readdlm(path*"mean_corr_se"*ext)
    corr_mean_data::Matrix{Float64} = readdlm(path*"corr_mean_data"*ext)
    corr_mean_se::Matrix{Float64} = readdlm(path*"corr_mean_se"*ext)
    return pulse_data,pulse_se,chase_data,chase_se,ratio_data,ratio_se,mean_corr_data,mean_corr_se,corr_mean_data,corr_mean_se
end

function load_s_data(path::String, model_name::String, ext::String)
    s_pulse::Matrix{Float64} = readdlm(path*model_name*"/s_pulse_"*model_name*ext)
    s_chase::Matrix{Float64} = readdlm(path*model_name*"/s_chase_"*model_name*ext)
    s_ratios::Matrix{Float64} = readdlm(path*model_name*"/s_ratios_"*model_name*ext)
    s_mean_corr::Matrix{Float64} = readdlm(path*model_name*"/s_mean_corr_"*model_name*ext)
    s_corr_mean::Matrix{Float64} = readdlm(path*model_name*"/s_corr_mean_"*model_name*ext)
    return s_pulse,s_chase,s_ratios,s_mean_corr,s_corr_mean
end

function compute_trunc_errors(pulse_data::Matrix{Float64}, pulse_se::Matrix{Float64}, chase_data::Matrix{Float64}, chase_se::Matrix{Float64},
    ratio_data::Matrix{Float64}, ratio_se::Matrix{Float64}, mean_corr_data::Matrix{Float64}, mean_corr_se::Matrix{Float64}, 
    corr_mean_data::Matrix{Float64}, corr_mean_se::Matrix{Float64}, s_pulse::Matrix{Float64}, s_chase::Matrix{Float64}, 
    s_ratios::Matrix{Float64}, s_mean_corr::Matrix{Float64}, s_corr_mean::Matrix{Float64}, model_name::String)
    local io::IOStream
    n_summary_stats::Int64 = 4*size(pulse_data)[2] + 3*size(ratio_data)[2]
    for i in 1:size(s_ratios)[1]
        err = Vector{Float64}(undef,size(ratio_data)[1])
        for j in 1:size(ratio_data)[1]
            data_::Vector{Array{Float64}} = [pulse_data[j*2-1,:],pulse_data[j*2,:],chase_data[j*2-1,:],chase_data[j*2,:],ratio_data[j,:],mean_corr_data[j,:],corr_mean_data[j,:]]
            se_::Vector{Array{Float64}} = [pulse_se[j*2-1,:],pulse_se[j*2,:],chase_se[j*2-1,:],chase_se[j*2,:],ratio_se[j,:],mean_corr_se[j,:],corr_mean_se[j,:]]
            s_data_::Vector{Array{Float64}} = [s_pulse[i*2-1,:],s_pulse[i*2,:],s_chase[i*2-1,:],s_chase[i*2,:],s_ratios[i,:],s_mean_corr[i,:],s_corr_mean[i,:]]
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

# m: modelling hypothesis index  --> {1: constant rates - scaling with size
#                                     2: constant rates - non-scaling
#                                     3: varying burst frequency - scaling with size
#                                     4: varying burst size - scaling with size
#                                     5: varying decay rate - scaling with size }
#                                       
 

m = 1   #2, 3, 4, 5

model_name = ["const","const_const","kon","alpha","gamma"][m]

#load model outputs
s_pulse,s_chase,s_ratios,s_mean_corr,s_corr_mean = load_s_data("data/large_scale_simulations/simulations_1/",model_name,".txt")

#load data summary statistics
pulse_data,pulse_se,chase_data,chase_se,ratio_data,ratio_se,
mean_corr_data,mean_corr_se,corr_mean_data,corr_mean_se = load_summary_stats("data/summary_stats/", ".txt");

@time (compute_trunc_errors(pulse_data,pulse_se,chase_data,chase_se,ratio_data,ratio_se,
        mean_corr_data,mean_corr_se,corr_mean_data,corr_mean_se,s_pulse,s_chase,s_ratios,
        s_mean_corr,s_corr_mean,model_name));



