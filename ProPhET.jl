### -*- Mode: Julia -*-

### ProPhET.jl
using RData
using NamedArrays
include("libs/homeoplasia.jl")
include("libs/entropy.jl")
include("libs/bioModel.jl")

data = load(ARGS[1], convert = true)
data_names = load(ARGS[2], convert = true)["names"]

# P = NamedArray(load(ARGS[3], convert = true)["processed_variants"],
#                (data_names["P_rows"], data_names["P_cols"]))

B = NamedArray(data["inference"]["B"],
               (data_names["B_rows"], data_names["B_cols"]))

C = NamedArray(data["inference"]["C"]["Experiment_1"][:,1],
               data_names["C_rows"])

G = NamedArray(data["inference"]["corrected_genotypes"],
               (data_names["G_rows"], data_names["G_cols"]))

α = data["inference"]["error_rates"]["alpha"]
β = data["inference"]["error_rates"]["beta"]
# D = buildD(B, C, P)
D = NamedArray(load(ARGS[3], convert=true)["clonal_variants_1"],
               (data_names["C_rows"], data_names["C_cols"]))

## print(findSNV(B, C))
## print(allnames(buildD(B, C, P)))
# print(compare(D, G))
Main.ErrorStats.check_HP_violation(Main.ErrorStats.compare(D,G), α, β)

clone_error = Main.ErrorStats.error_distribution(G, C, D)
println(clone_error)
E = Main.Entropy.prob_distribution(clone_error, C, α, β)
println(E)
entropyV = Main.Entropy.entropy_array(E)
println(entropyV)

Main.BioModel.inference(entropyV, 10000, 1)

# for (k,v) in clone_error
#     println("Clone $k")
#     println(v)
#     println()
# end
### end of file -- ProPhET.jl
