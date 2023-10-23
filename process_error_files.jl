using JDF
using DelimitedFiles
using DataFrames


# m: modelling hypothesis index  --> {1: constant rates - scaling with size
#                                     2: constant rates - non-scaling
#                                     3: varying burst frequency - scaling with size
#                                     4: varying burst size - scaling with size
#                                     5: varying decay rate - scaling with size }
#                                       
 

m = 1   #2, 3, 4, 5

model_name = ["const","const_const","kon","alpha","gamma"][m]


error_data = readdlm("/rds/general/user/dv19/ephemeral/errors/error_"*model_name*".txt")
println("Read.")

df = DataFrame(data,:auto)
println("Dataframe.")

JDF.save("Julia/large_scale_simulations/errors/error_const.jdf", df)
println("JDF done.")
