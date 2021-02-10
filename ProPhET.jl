### -*- Mode: Julia -*-

### ProPhET.jl
using RData
using Plots
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
# D = NamedArray(load(ARGS[3], convert=true)["clonal_variants_1"],
#                (data_names["C_rows"], data_names["C_cols"]))
clonal_variants = NamedArray(load(ARGS[3], convert=true)["clonal_variants"],
                             (data_names["CL_rows"], data_names["CL_cols"]))
D = Main.ErrorStats.buildD_clonal(clonal_variants, 0.030)
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

# inf = Main.BioModel.inference(entropyV, 10000, 1)
# Main.BioModel.save_plot(inf, "./sample.png")

# println(Main.BioModel.inference(entropyV, 1000, 1))
a, b, t = Main.BioModel.optimal(entropyV, 10000, 1, 0.01)
println(a)
println(b)
savefig(plot(a), "./a.png")
savefig(plot(b), "./b.png")
savefig(plot(t), "./t.png")
### end of file -- ProPhET.jl
