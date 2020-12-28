using RData
using NamedArrays

function findSNV(B, C)
    # C = samples
    for i in 1:size(B,1)
        row = B[i,:]
        println(row)

    end
end



data = load(ARGS[1], convert=true)
names = load(ARGS[2], convert=true)["names"]
B = NamedArray(data["inference"]["B"], (names["B_rows"], names["B_cols"]))
# C = Dict(zip(data["inference"]["C"]["Experiment_1"], names))
C = NamedArray(data["inference"]["C"]["Experiment_1"][:,1], names["C_rows"])
α = data["inference"]["error_rates"]["alpha"]
β = data["inference"]["error_rates"]["beta"]

findSNV(B,C)
