### -*- Mode: Julia -*-

### entropy.jl

module Entropy
using NamedArrays

### Define single node distribution
function node_prob(n::Integer, fp::Integer, fn::Integer, α::Float64, β::Float64)
    bin = ((x, fx, ξ) -> binomial(big(n), fx) * ξ ^ fx * (1 - ξ) ^ (n - fx))

    return bin(n, fp, α) * bin(n, fn, β)
end

### Build error matrix distribution given C, α, β
function prob_distribution(error_m::Tuple{NamedArray, NamedArray},
                           C::NamedArray, α::Float64, β::Float64)
    s = size(error_m[1])
    rows = allnames(error_m[1])[1]
    result = NamedArray(zeros(s[1], s[2]),
                        (allnames(error_m[1])[1], allnames(error_m[1])[2]))

    # number of samples
    n = Dict([(i, count(x -> x == i, C)) for i in unique(C)])
    for i in 1 : s[1]
        for j in 1 : s[2]
            current = parse(Int, rows[i])
            result[i,j] = node_prob(n[current], Int(error_m[1][i,j]),
                                     Int(error_m[2][i,j]), α, β)
        end
    end
    # for i in allnames(error_m[1])[1]
    #     for j in allnames(error_m[1])[2]
    #         current = parse(Int, i)
    #         result[i,j] = node_prob(n[current], Int(error_m[1][i,j]),
    #                                 Int(error_m[2][i,j]), α, β)
    #     end
    # end
    return result
end


### Generic single node entropy function
### error_p = probability matrix
### node = rowname
function entropy(error_p::NamedArray, node)
    # e = ((pᵢ, m) -> pᵢ <= pᵢ / m ? pᵢ :  1 / (1 - m) * pᵢ - 1 / (1 - m))
    e = ((pᵢ, m) -> pᵢ > 1 / m ? (1 / (1 - m)) * ( pᵢ - 1) : pᵢ)
    m = size(error_p)[2]
    row = error_p[node, :]

    return sum([e(i, m) for i in row])
end


function entropy_array(error_p::NamedArray)
    result = NamedArray([entropy(error_p, i) for i in allnames(error_p)[1]],
                        (allnames(error_p))[1])
    return result
end

end

### end of file -- entropy.jl
