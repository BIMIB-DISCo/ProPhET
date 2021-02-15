### -*- Mode: Julia -*-

### bioModel.jl

module BioModel

using NamedArrays
using Turing
using Distributions
ENV["GKS_ENCODING"]="utf-8"
using Plots
using StatsPlots

@model function opt_α(y, λ, β)
    α ~ Exponential(λ)
    dist ~ Beta(α, β)
    s = sum([log(i) for i in y])
    w = (dist ^ - length(y)) * exp((- α * (λ - s)))
    Turing.@addlogprob! w
end

@model function opt_β(y, λ, α)
    β ~ Exponential(λ)
    dist ~ Beta(α, β)
    s = sum([log(1 - i) for i in y])
    w = (dist ^ - length(y)) * exp((- β * (λ - s)))
    Turing.@addlogprob! w
end

@model function opt(y, λ)
    α ~ Exponential(λ)
    β ~ Exponential(λ)
    dist ~ Beta(α, β)
    sa = sum([log(i) for i in y])
    sb = sum([log(1-i) for i in y])
    w = (dist ^ - length(y)) * exp(- α * (λ - sa)) * exp(- β * (λ - sb))
    Turing.@addlogprob! w
end

function sampling(var, model, par, data, n, λ)
    #s = ((samplers, vars) -> map(x -> (foldl(∘, samplers))(x), vars))
    return sample(model(data, λ, par), Gibbs(MH(:dist), MH(var)), n)
end

function sampling_NUTS(var, model, par, data, n, λ)
    return sample(model(data, λ, par), NUTS(), n)
end

function opt(data, n, λ, figure)
    # sa = sample(opt(data, λ), Gibbs(MH(:α), MH(:β), MH(:dist)), n)
    # sb = sample(opt(data, λ), Gibbs(MH(:α), MH(:β), MH(:dist)), n)
    sa = sample(opt(data, λ), NUTS(), n)
    sb = sample(opt(data, λ), NUTS(), n)
    println(mean(sa[:α]))
    println(mean(sa[:β]))
    println(mean(sb[:α]))
    println(mean(sb[:β]))
    save_plot(sa, string(figure, "a.png"))
    save_plot(sb, string(figure, "b.png"))
    return [sa, sb]
end

function generic_optimal(model, var, par, data, n, λ, error, figure)
    ξ = 10
    tmp_ξ = 1
    steps = []
    actual_sample = 0
    while abs(ξ - tmp_ξ) > error * tmp_ξ
        tmp_ξ = ξ
        actual_sample = sampling(var, model, tmp_ξ, data, n, λ)
        println(ξ)
        # actual_sample = sampling_NUTS(var, model, tmp_ξ, data, n, λ)
        # valore posterior maggiore
        
        ξ = mean(actual_sample[var])
        push!(steps, ξ)
    end
    save_plot(actual_sample, figure)
    return steps
end

function optimal(data, n, λ, error)
    # tmp = rand(uniform)
    tmp_β = mean(sample(mutation(data, λ), Gibbs(MH(:α), MH(:β)), n)[:β])
    tmp_α = mean(sample(opt_α(data, λ, tmp_β), Gibbs(MH(:α), MH(:dist)), n)[:α])
    β_steps = generic_optimal(opt_β, :β, tmp_α, data, n, λ, error,
    "./betaDist.png") 
    α_steps = generic_optimal(opt_α, :α, last(β_steps), data, n, λ, error,
    "./alphaDist.png") 
    return [α_steps, β_steps]
end

function save_plot(inf, filename)
    gr()
    savefig(plot(inf, fontfamily = "Symbol"), filename)
end
end
### end of file -- bioModel.jl
