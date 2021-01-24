### -*- Mode: Julia -*-

### entropy.jl

module Entropy
using NamedArrays

# define single node distribution
function node_prob(n, fp, fn, α, β)
    bin = ((x, fx, ξ) -> binomial(n, fx) * ξ ^ fx * (1 - ξ) ^ (n - fx))

    return bin(n, fp, α) * bin(n, fn, β)
end

function prob_distribution(error_m, C, α, β)
    s = size(error_m[1])
    rows = allnames(error_m[1])[1]
    result = NamedArray(zeros(s[1], s[2]),
                        (allnames(error_m[1])[1], allnames(error_m[1])[2]))

    # number of samples
    n = Dict([(i, count(x -> x == i, C)) for i in unique(C)])
    for i in 1 : s[1]
        for j in 1 : s[2]
            current = parse(Int, rows[i])
            # println(i," ",j)
            result[i,j] = node_prob(n[current], Int(error_m[1][i,j]),
                                    Int(error_m[2][i,j]), α, β)
            println(result[i,j])
        end
    end
    return result
end

end
### end of file -- entropy.jl
