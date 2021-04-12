# ProPhET

ProPhET is a probabilistic programming tool written in [Julia](https://julialang.org) and [Turing](https://turing.ml) that can be used to estimate the measure of
a *phylogenetic signal*.  The ispiration for the tool came form the work our group did on the inference of SARS-CoV-2 sublineage variants identification, [*VERSO: A comprehensive framework for the inference of robust phylogenies and the quantification of intra-host genomic diversity of viral samples*](https://doi.org/10.1016/j.patter.2021.100212) (cfr., the [VERSO](https://github.com/BIMIB-DISCo/VERSO) repository). 
The original Bayesian model used in ProPhET is described in [*Measuring phylogenetic signal between categorical traits and phylogenies*](https://doi.org/10.1093/bioinformatics/bty800).


## Input

ProPhET takes in input the following files:
- **inference**: the VERSO step one output
- **clonal_variants**: the VERSO input
- **names**: exported names of rows and columns of the following objects:
    - B matrix
    - C matrix
    - corrected_genotypes
    - clonal_variants


### Input preparation

To build **names** it is possible to export the data from R:
```R
names = list("B_rows" = rownames(inference$B),
             "B_cols" = colnames(inference$B),
             "C_rows" = rownames(inference$C$Experiment_1),
             "G_rows" = rownames(inference$corrected_genotypes),
             "G_cols" = colnames(inference$corrected_genotypes),
             "CL_rows" = rownames(clonal_variants),
             "CL_cols" = colnames(clonal_variants))
save(names, file="~/names.RData")
```

## Running ProPhET

**ProPheT** can be run form the CLI as follows:

```bash
$ julia ProPhET.jl inference.RData clonal_variants.RData names.RData
```
Is also possible to specify the following optional parameters:

| Parameter | Effect                                     | Type                | Default |
|-----------|--------------------------------------------|---------------------|---------|
| `--sampler` | Set sampler                                | String (NUTS or MH) |    NUTS |
| `--n`       | Set the number of generated samples        | Integer             |   40000 |
| `--lambda`  | Set the lambda hyperparameter              | Float64             |     0.1 |
| `--error`   | Set the error rate for mutations selection | Float64             |   0.030 |
| `--pltName` | Set the plots filename prefix              | String              |   plot_ |

## Inference

ProPhET computes the entropy value of each node and uses the entropy vector to
perform a Bayesian inference. The sampling process can be built with the
following methods:

| --sampler value | Sampler                            |
|-----------------|------------------------------------|
| MH              | ```Gibbs(MH(:alpha), MH(:beta))``` |
| NUTS            | ```NUTS()```                       |

## Output

The output of ProPhET is a real number $\delta$. An high value of $\delta$ means
that the phylogenetic signals are more informative. A low value of $\delta$ 
could suggest that a case of homoplasia could be present in the phylogenetic 
tree.
