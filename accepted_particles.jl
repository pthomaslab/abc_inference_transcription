f_1 = JDFFile("data/large_scale_simulations/errors/error_const.jdf");
f_2 = JDFFile("data/large_scale_simulations/errors/error_const_const.jdf");
f_3 = JDFFile("data/large_scale_simulations/errors/error_kon.jdf");
f_4 = JDFFile("data/large_scale_simulations/errors/error_alpha.jdf");
f_5 = JDFFile("data/large_scale_simulations/errors/error_gamma.jdf");

model_names = ["const","const_const","kon","alpha","gamma"];

#set acceptance threshold ε
ε = 4.8;

for g in 1:length(genes)
    err_mat = Matrix{Float64}(undef,(1000000,length(model_names)))
    err_mat[:,1] = f_1["x"*string(g)]    
    err_mat[:,2] = f_2["x"*string(g)]   
    err_mat[:,3] = f_3["x"*string(g)]    
    err_mat[:,4] = f_4["x"*string(g)]    
    err_mat[:,5] = f_5["x"*string(g)]    
    for j in 1:lastindex(model_name)
        v::Vector{Int64} = findall(x->x<=ε,err_mat[:,j])
        if length(v) > 0
            sort_idx::Vector{Int64} = sortperm(err_mat[v,j])
            open("data/posteriors/particles_"*model_names[j]*".txt", "a") do io
                DelimitedFiles.writedlm(io, reshape(v[sort_idx],1,:))
            end
        else
            open("data/posteriors/particles_"*model_names[j]*".txt", "a") do io
                DelimitedFiles.writedlm(io, [0])
            end
        end
    end
end
