# Calculate top ligands using a PCA approach from receptor genes given a database

Calculate top ligands using a PCA approach from receptor genes given a
database

## Usage

``` r
pca_limma(eset, design, db, obs_id = NULL, correlation = NULL)
```

## Arguments

- eset:

  Expression Set object containing gene expression data.

- design:

  Design matrix generated from create_design()

- db:

  Ligand-receptor database

- obs_id:

  Optional: provide a vector of sample IDs making sure the order matches
  with the eset

- correlation:

  Optional: input the correlation consensus between the samples to
  evaluate if it is paired data

## Value

a data frame of differentially expressed ligands ordered by p-values

## Examples

``` r
if (FALSE) { # \dontrun{
# pca_limma_res <- pca_limma(eset, design, db)
} # }
```
