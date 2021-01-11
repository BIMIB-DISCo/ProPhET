### -*- Mode: Julia -*-

### homeoplasia.jl

using RData
using NamedArrays


function findSNV(B, C)
    genSample = Dict()
    for i in 2:size(B,1)
        push!(genSample, i-1 => [k for (k, v) in enumerate(C) if v == i])
    end
    return genSample
end


function buildD(B, C, P)
    ## Get mutation and sample names
    mutation_name = filter!(e->e != "Root", allnames(B)[2])
    sample_name = allnames(C)[1]

    ## Build D getting values from P
    ## result = reshape([P[sample,mutation] for sample in allnames(C)
    ## for mutation in allnames(B)[2]], (size(C)[1], size(B)[2]))
    result = reshape([P[sample,mutation]
                      for sample in sample_name
                      for mutation in mutation_name],
                     (length(sample_name), length(mutation_name))) 

    ## Set D colnames and rownames
    result = NamedArray(result, (sample_name, mutation_name))

    ## element >= 0.9 -> 1
    ## element <  0.9 -> 0
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


function compare(D, G)
    ## compare D (original data) with G (corrected genotypes)
    total = Dict()
    for j in allnames(D)[2]
        fp = 0
        fn = 0
        for i in allnames(D)[1]
            if !ismissing(D[i,j]) && !ismissing(G[i,j])
                if D[i,j] != G[i,j]
                    if G[i,j] == 0
                        fp += 1 #/size(G)[1]
                    elseif G[i,j] == 1
                        fn += 1 #/size(G)[1]
                    end
                end
            end
        end
        push!(total, j => (fp, fn))
    end

    return total
end


function check_HP_violation(total, α, β)
    ## find suspicious mutation given alpha and beta
    ## cand = {}
    for (k,v) in total
        if v[1] >= α
            println("FP: $k")
            # println(α)
            println(v[1])
        end
        if v[2] >= β
            println("FN: $k")
            # println(β)
            println(v[1])
        end
    end
end


data = load(ARGS[1], convert=true)
names = load(ARGS[2], convert=true)["names"]

P = NamedArray(load(ARGS[3], convert = true)["processed_variants"],
               (names["P_rows"], names["P_cols"]))

B = NamedArray(data["inference"]["B"],
               (names["B_rows"], names["B_cols"]))

C = NamedArray(data["inference"]["C"]["Experiment_1"][:,1],
               names["C_rows"])

G = NamedArray(data["inference"]["corrected_genotypes"],
               (names["G_rows"], names["G_cols"]))

α = data["inference"]["error_rates"]["alpha"]
β = data["inference"]["error_rates"]["beta"]
D = buildD(B, C, P)

## print(findSNV(B, C))
## print(allnames(buildD(B, C, P)))
print(compare(D, G))
## check_HP_violation(compare(D,G), α, β)


### end of file -- homeoplasia.jl

