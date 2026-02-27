# Perform PCA and fit a multivariate regression to measure ligand activity

Perform PCA and fit a multivariate regression to measure ligand activity

## Usage

``` r
pca_plsda(eset, treatment, db)
```

## Arguments

- eset:

  Expression Set object containing gene expression data.

- treatment:

  Treatment response variable

- db:

  ligand-receptor database

## Value

Table of ranked ligands by coef order (largest to smallest)

## Examples

``` r
if (FALSE) { # \dontrun{
# Load the database
# data(dbs_all)
# Take one database and run the method on the data and condition
# pca_plsda(eset, treatment, dbs_all$baderlab)
} # }
```
