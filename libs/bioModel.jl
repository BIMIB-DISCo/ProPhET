### -*- Mode: Julia -*-

### bioModel.jl

module BioModel

using Turing
using Distributions
using SpecialFunctions
ENV["GKS_ENCODING"]="utf-8"
using Plots
using StatsPlots
using Statistics

@model function mutation(y::Array{Float64,1}, λ::Float64)
    α ~ Exponential(λ)
    β ~ Exponential(λ)
    for i in 1 : length(y)
        y[i] ~ Beta(α, β)
    end
end


function opt(data::Array{Float64,1}, n::Integer, λ::Float64, sampler,
    figure::String)
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

function save_plot(inf, filename::String)
    gr()
    savefig(plot(inf, fontfamily = "Symbol"), filename)
end

function print_stats(values::Array{Any,1})
    return [Statistics.mean(values), Statistics.std(values)]
end

function variation(data::Array{Float64,1}, λ::Float64, n::Integer,
                   sampler::String, figure::String, exc::Integer)
    va = []
    vb = []
    vδ = []
    for i in 1:exc
        println(sampler)
        sa, sb = opt(data, n, λ, sampler, figure)
        push!(va, mean(sa[:α]))
        push!(vb, mean(sb[:β]))
        push!(vδ, get_δ(sa,sb))
    end
    savefig(plot([va, vb, vδ], label = ["α" "β" "δ"]), "./var_a.png")
    return map(x -> print_stats(x), [va, vb, vδ])
end

end

### end of file -- bioModel.jl
