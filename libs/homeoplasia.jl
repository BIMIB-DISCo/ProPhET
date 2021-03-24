### -*- Mode: Julia -*-

### homeoplasia.jl

module ErrorStats

using NamedArrays
using Distributions

function findSNV(B::NamedArray, C::NamedArray)
    genSample = Dict()
    for i in 2 : size(B, 1)
        push!(genSample, i - 1 => [k for (k, v) in enumerate(C) if v == i])
    end
    return genSample
end


function buildD(B, C, P)
    ## Get mutation and sample names
    mutation_name = filter!(e -> e != "Root", allnames(B)[2])
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
        elseif x >= 0.9
            return 1
        elseif x < 0.9
            return 0
        end
    end

    return map(repl, result)
end

function buildD_clonal(cl::NamedArray, p::Float64)
    # manual skip missing
    s = sum(map(x -> ismissing(x) ? 0 : x, cl), dims=1)[1,:] / size(cl)[1]

    # filter! is not supported by NamedArray
    # s = filter!(x -> x> 0.030, e)
    #
    # workaround:
    # get mutations names
    mut = [allnames(s)[1][i] for i in 1 : length(s) if s[i] > p]

    # foldl cl columns for the mutations names in mut
    return NamedArray(foldl(hcat, [cl[:,i] for i in mut]),
                      (allnames(cl)[1], mut))
end


function compare(D::NamedArray, G::NamedArray)
    ## compare D (original data) with G (corrected genotypes)
    total = Dict()
    for j in allnames(D)[2]
        fp = 0
        fn = 0
        for i in allnames(D)[1]
            if !ismissing(D[i,j]) && !ismissing(G[i,j])
                if D[i,j] != G[i,j]
                    if G[i,j] == 0
                        fp += 1/size(G)[1]
                    elseif G[i,j] == 1
                        fn += 1/size(G)[1]
                    end
                end
            end
        end
        push!(total, j => (fp, fn))
    end

    return total
end


function check_HP_violation(total::Dict{}, α::Float64, β::Float64)
    ## find suspicious mutation given alpha and beta
    ## cand = {}
    println("Beta = $β")
    println("Alpha = $α")
    println()
    for (k,v) in total
        if v[1] >= α
            println(k)
            println("FP: $v\n")
        end
        if v[2] >= β
            println(k)
            println("FN: $v\n")
        end
    end
end


function error_distribution(G::NamedArray, C::NamedArray, D::NamedArray)
    ## compare G with C/B
    ## total = Dict()
    rows = sort(map(x -> string(x), unique(C)))
    cols = allnames(D)[2]
    tot_fp = NamedArray(zeros(length(rows), length(cols)),
                        (rows, cols))

    tot_fn = NamedArray(zeros(length(rows), length(cols)),
                        (rows, cols))

    for i in rows
        for j in cols
            fp = 0
            fn = 0
            cluster = allnames(filter(x -> x == parse(Int32, i), C))[1]
            for c in cluster
                if !(ismissing(G[c,j]) || ismissing(D[c,j]))
                    if G[c,j] != D[c,j]
                        if G[c,j] == 0
                            fp += 1
                        end
                        if G[c,j] == 1
                            fn += 1
                        end
                    end
                end
                tot_fn[i,j] = Int64(fn) # / length(cluster)
                tot_fp[i,j] = Int64(fp) # / length(cluster)
            end
        end
    end

    return (tot_fp, tot_fn)
end

end

### end of file -- homeoplasia.jl
