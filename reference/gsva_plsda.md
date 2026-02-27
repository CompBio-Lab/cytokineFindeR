# Calculate coefficients from variable selection using GSVA and PLS-DA

Uses a classification approach when computing receptor weights with
Partial Least Squares - Discriminant Analysis and pools these receptors
together to compute a value for the ligand or cytokine given a gene set.
As limma treats all variables as equal importance, PLSDA provides a
multiple regression approach

## Usage

``` r
gsva_plsda(eset, treatment, db)
```

## Arguments

- eset:

  Expression Set object containing gene expression data.

- treatment:

  Treatment response variable

- db:

  ligand-receptor database

## Value

A named vector of ligands indicating importance

## Examples

``` r
if (FALSE) { # \dontrun{
# Load the database
# data(dbs_all)
# Take one database and run the method on the data and condition for an unpaired dataset
gsva_plsda(eset, treatment, dbs_all$baderlab)
} # }
```
