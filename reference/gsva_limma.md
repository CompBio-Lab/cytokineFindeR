# Gene Set Variation Analysis for cytokines (cGSVA)

Gene Set Variation Analysis for cytokines (cGSVA)

## Usage

``` r
gsva_limma(eset, design, db, obs_id = NULL, correlation = NULL)
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

top table of all ligands

## Examples

``` r
if (FALSE) { # \dontrun{
# Load the database
# data(dbs_all)
# Take one database and run the method on the data and condition for an unpaired dataset
# gsva_limma(eset, treatment, dbs_all$baderlab)
} # }
```
