# Gene Set Enrichment Analysis

Gene Set Enrichment Analysis

## Usage

``` r
cfgsea(eset, design, db, obs_id = NULL, correlation = NULL)
```

## Arguments

- eset:

  Expression Set object containing gene expression data.

- design:

  Design matrix generated from create_design()

- db:

  Ligand-receptor database

- obs_id:

  A vector of observation IDs to indicate if it's paired data

- correlation:

  Add a correlation block based on the dupcor package for paired
  analysis

## Value

a table with GSEA results. Each row corresponds to a ligand

## Examples

``` r
# This is part of a series of enrichment analysis methods
# Basic usage:
if (FALSE) { # \dontrun{
#gsea_res <- cfgsea(eset, design, dbs)
} # }
```
