using JDF

for model_name in ["const","const_const","kon","alpha","gamma"]
    error_data = readdlm("/rds/general/user/dv19/ephemeral/errors/error_"*model_name*".txt")
    df = DataFrame(data,:auto)
    JDF.save("data/large_scale_simulations/errors/error_"*model_name*".jdf", df)
end
