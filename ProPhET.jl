### -*- Mode: Julia -*-

### ProPhET.jl
using RData
using NamedArrays
using Fire
include("libs/homeoplasia.jl")
include("libs/entropy.jl")
include("libs/bioModel.jl")

@main function inference(inference::AbstractString,
                         clonalVariants::AbstractString,
                         dataNames::AbstractString;
                         sampler::AbstractString="NUTS",
                         n::Integer=40000,
                         lambda::Float64=0.1,
                         error::Float64=0.030,
                         pltName::AbstractString="plot_")
    data = load(inference, convert = true)
    data_names = load(dataNames, convert = true)["names"]

    # P = NamedArray(load(ARGS[3], convert = true)["processed_variants"],
    #                (data_names["P_rows"], data_names["P_cols"]))

    B = NamedArray(data["inference"]["B"],
                   (data_names["B_rows"], data_names["B_cols"]))

    C = NamedArray(data["inference"]["C"]["Experiment_1"][:,1],
                   data_names["C_rows"])

    G = NamedArray(data["inference"]["corrected_genotypes"],
                   (data_names["G_rows"], data_names["G_cols"]))

    α = data["inference"]["error_rates"]["alpha"]
    β = data["inference"]["error_rates"]["beta"]
    # D = buildD(B, C, P)
    # D = NamedArray(load(ARGS[3], convert=true)["clonal_variants_1"],
    #                (data_names["C_rows"], data_names["C_cols"]))
    clonal_variants = NamedArray(load(clonalVariants,
                                      convert=true)["clonal_variants"],
                                 (data_names["CL_rows"], data_names["CL_cols"]))
    D = Main.ErrorStats.buildD_clonal(clonal_variants, error)
    ## print(findSNV(B, C))
    ## print(allnames(buildD(B, C, P)))
    # print(compare(D, G))
    Main.ErrorStats.check_HP_violation(Main.ErrorStats.compare(D,G), α, β)

    clone_error = Main.ErrorStats.error_distribution(G, C, D)
    println(clone_error)
    E = Main.Entropy.prob_distribution(clone_error, C, α, β)
    println(E)
    entropyV = Main.Entropy.entropy_array(E)
    println(entropyV)

    sa, sb = Main.BioModel.opt(Array(entropyV), n, lambda, sampler, pltName)
    println(string("δ = ", Main.BioModel.get_δ(sa, sb)))
end
### end of file -- ProPhET.jl
