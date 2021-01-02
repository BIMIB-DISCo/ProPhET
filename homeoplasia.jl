using RData
using NamedArrays

function findSNV(B, C)
    genSample = Dict()
    for i in 2:size(B,1)
        push!(genSample, i-1 => [k for (k,v) in enamerate(C) if v==i])
    end
    return genSample
end

function buildD(B, C, P)
    # Get mutation and sample names
    mutation_name = filter!(e->e!="Root", allnames(B)[2])
    sample_name = allnames(C)[1]

    # Build D getting values from P
    # result = reshape([P[sample,mutation] for sample in allnames(C) for mutation in allnames(B)[2]], (size(C)[1], size(B)[2]))
    result = reshape([P[sample,mutation] for sample in sample_name for mutation in mutation_name], (length(sample_name), length(mutation_name)))

    # Set D colnames and rownames
    result = NamedArray(result, (sample_name, mutation_name))

    # element >= 0.9 -> 1
    # element <  0.9 -> 0
    function repl(x)
        if ismissing(x)
            return missing
        elseif x>=0.9
            return 1
        elseif x<0.9
            return 0
        end
    end


    return map(x->repl(x), result)
end


data = load(ARGS[1], convert=true)
names = load(ARGS[2], convert=true)["names"]
P = NamedArray(load(ARGS[3], convert=true)["processed_variants"], (names["P_rows"], names["P_cols"]))
B = NamedArray(data["inference"]["B"], (names["B_rows"], names["B_cols"]))
# C = Dict(zip(data["inference"]["C"]["Experiment_1"], names))
C = NamedArray(data["inference"]["C"]["Experiment_1"][:,1], names["C_rows"])
α = data["inference"]["error_rates"]["alpha"]
β = data["inference"]["error_rates"]["beta"]
D = buildD(B, C, P)
print(findSNV(B,C))
print(buildD(B, C, P))
