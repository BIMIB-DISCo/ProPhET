### -*- Mode: Julia -*-

### bioModel.jl

module BioModel

using NamedArrays
using Turing
using Distributions
ENV["GKS_ENCODING"]="utf-8"
using Plots
using StatsPlots


@model function mutation(y, λ)
    α ~ Exponential(λ)
    β ~ Exponential(λ)
    for i in 1 : length(y)
        y[i] ~ Beta(α, β)
    end
end

function inference(data, n, λ)
    return sample(mutation(data, λ), Gibbs(MH(:α), MH(:β)), n)
end

function save_plot(inf, filename)
    gr()
    savefig(plot(inf, fontfamily = "Symbol"), "./sample.png")
end
end
### end of file -- bioModel.jl
