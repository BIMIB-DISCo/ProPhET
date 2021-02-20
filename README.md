# ProPhET
ProPhET is a probabilistic programming tool able to estimate the measure of
phylogenetic signal generated by VERSO. The original Bayesian model used in
ProPhET can be found at:
```10.1093/bioinformatics/bty800```


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
names = list("B_rows" = rownames(inference$B), "B_cols"=colnames(inference$B),
"C_rows"=rownames(inference$C$Experiment_1),
"G_rows"=rownames(inference$corrected_genotypes),
"G_cols"=colnames(inference$corrected_genotypes),
"CL_rows"=rownames(clonal_variants), "CL_cols"=colnames(clonal_variants))
save(names, file="~/names.RData")
```

## Running ProPhET
```bash
$ julia ProPhET.jl inference.RData clonal_variants.RData names.RData
```
Is also possible to specify the following optional parameters:
| Parameter | Effect                                     |                Type | Default |
|-----------|--------------------------------------------|---------------------|---------|
| --sampler | Set samper                                | String (NUTS or MH) |    NUTS |
| --n       | Set the number of generated samples        |               40000 |   40000 |
| --lambda  | Set the lambda hyperparameter              |                 0.1 |     0.1 |
| --error   | Set the error rate for mutations selection |               0.030 |   0.030 |
| --pltName | Set the plots filename prefix              |               plot_ |   plot_ |

## Inference
ProPhET computes the entropy value of each node and uses the entropy vector to
perform a Bayesian inference. The sampling process is built with
```julia
Gibbs(MH(:alpha), MH(:beta))
```

## Output
The output of ProPhET is a real number $\delta$. An high value of $\delta$ means
that the phylogenetic signals are more informative. A low value of $\delta$ 
could suggest that a case of homoplasia could be present in the phylogenetic 
tree.
