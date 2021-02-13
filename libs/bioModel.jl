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
    # β ~ Exponential(λ)
    dist ~ Beta(α, β)
    s = sum([log(i) for i in y])
    w = (dist ^ - length(y)) ^ (- α * (λ - s))
    Turing.@addlogprob! w
end

@model function opt_β(y, λ, α)
    # α ~ Exponential(λ)
    β ~ Exponential(λ)
    dist ~ Beta(α, β)
    s = sum([log(1 - i) for i in y])
    w = (dist ^ - length(y)) ^ (- β * (λ - s))
    Turing.@addlogprob! w
end

@model function mutation(y, λ)
    α ~ Exponential(λ)
    β ~ Exponential(λ)
    for i in 1 : length(y)
        y[i] ~ Beta(α, β)
    end
end

function inference(data, n, λ)
    # β is the mean of the sample, it could be a random numeric value
    temp_β = mean(sample(mutation(data, λ), Gibbs(MH(:α), MH(:β)), n)[:β])
    # evaluate α (where β is fixed)
    α = mean(sample(opt_α(data, λ, temp_β), Gibbs(MH(:α)), n)[:α])
    # evaluate β (where α is fixed)
    β = mean(sample(opt_β(data, λ, α), Gibbs(MH(:β)), n)[:β])

    return [α, β]
end

function sampling(var, model, par, data, n, λ)
    #s = ((samplers, vars) -> map(x -> (foldl(∘, samplers))(x), vars))
    return sample(model(data, λ, par), Gibbs(MH(:dist), MH(var)), n)
end

function sampling_NUTS(var, model, par, data, n, λ)
    return sample(model(data, λ, par), NUTS(), n)
end

function generic_optimal(model, var, par, data, n, λ, error, figure)
    ξ = 10
    tmp_ξ = 1
    steps = []
    actual_sample = 0
    while abs(ξ - tmp_ξ) > error * tmp_ξ
        tmp_ξ = ξ
        # actual_sample = sampling(var, model, tmp_ξ, data, n, λ)
        actual_sample = sampling_NUTS(var, model, tmp_ξ, data, n, λ)
        ξ = mean(actual_sample[var])
        push!(steps, ξ)
    end
    save_plot(actual_sample, figure)
    return steps
end

function optimal(data, n, λ, error)
    tmp_β = mean(sample(mutation(data, λ), Gibbs(MH(:α), MH(:β)), n)[:β])
    tmp_α = mean(sample(opt_α(data, λ, tmp_β), Gibbs(MH(:α), MH(:dist)), n)[:α])
    β_steps = generic_optimal(opt_β, :β, tmp_α, data, n, λ, error, "./betaDist.png")
    α_steps = generic_optimal(opt_α, :α, last(β_steps), data, n, λ, error, "./alphaDist.png")
    return [α_steps, β_steps]
end

function save_plot(inf, filename)
    gr()
    savefig(plot(inf, fontfamily = "Symbol"), filename)
end
end
### end of file -- bioModel.jl
