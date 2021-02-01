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
    dist ~ Beta(α, mean(β))
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

function optimal(data, n, λ, error)
    α = β = 10
    temp_β = mean(sample(mutation(data, λ), Gibbs(MH(:α), MH(:β)), n)[:β])
    temp_α = mean(sample(opt_α(data, λ, temp_β), Gibbs(MH(:α)), n)[:α])
    # tmp = 0
    # while tmp <= 2
    lα = []
    lβ = []
    lt = []
    # end
    while abs(α - temp_α) > error * temp_α && abs(β - temp_β) > error * temp_β
        temp_α = α
        temp_β = β
        # evaluate α (where β is fixed)
        α = mean(sample(opt_α(data, λ, temp_β), Gibbs(MH(:α)), n)[:α])
        # evaluate β (where α is fixed)
        β = mean(sample(opt_β(data, λ, temp_α), Gibbs(MH(:β)), n)[:β])
        println("α: ", α, " β:", β)
        append!(lα, α)
        append!(lβ, β)
        append!(lt, abs(α - temp_α))
    end
    return [lα, lβ, lt]

end


function tmp_α(beta)
end

function save_plot(inf, filename)
    gr()
    savefig(plot(inf, fontfamily = "Symbol"), "./sample.png")
end
end
### end of file -- bioModel.jl
