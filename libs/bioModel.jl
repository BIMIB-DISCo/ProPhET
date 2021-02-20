### -*- Mode: Julia -*-

### bioModel.jl

module BioModel

using NamedArrays
using Turing
using Distributions
ENV["GKS_ENCODING"]="utf-8"
using Plots
using StatsPlots
using SpecialFunctions
using Random

@model function mutation(y, λ)
    α ~ Exponential(λ)
    β ~ Exponential(λ)
    for i in 1 : length(y)
        y[i] ~ Beta(α, β)
    end
end


function opt(data, n, λ, sampler, figure)
    if sampler == "NUTS"
        sa = sample(mutation(data, λ), NUTS(), n)
        sb = sample(mutation(data, λ), NUTS(), n)
    elseif sampler == "MH"
        sa = sample(mutation(data, λ), Gibbs(MH(:α), MH(:β)), n)
        sb = sample(mutation(data, λ), Gibbs(MH(:α), MH(:β)), n)
    end
    save_plot(sa, string(figure, "a.png"))
    save_plot(sb, string(figure, "b.png"))
    return [sa, sb]
end

function get_δ(a, b)
    sample_α = a[:α]
    sample_β = b[:β]
    println(mean(sample_α))
    println(mean(sample_β))
    return mean(map((x, y) -> y / x, sample_α, sample_β))
end


function save_plot(inf, filename)
    gr()
    savefig(plot(inf, fontfamily = "Symbol"), filename)
end
end

### end of file -- bioModel.jl
