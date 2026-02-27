# Given an expression set, make sure that data is cleaned by gene name using a provided gene list data frame

Given an expression set, make sure that data is cleaned by gene name
using a provided gene list data frame

## Usage

``` r
clean_eset(eset, gene_list_df)
```

## Arguments

- eset:

  Expression Set as a numeric matrix

- gene_list_df:

## Value

cleaned up expression matrix

## Examples

``` r
if (FALSE) { # \dontrun{
# create a small gene_symbol and ensembl_id df from the GEO dataset example (golimumab TNF-targeted treatment) if annotation df exists
# gensym <- sapply(strsplit(golimumab$GSE92415_series_matrix.txt.gz$annotations$`Gene Symbol`, "///"), trimws) 
# probe2gene_df <- tibble(probeids = rep(rownames(golimumab$GSE92415_series_matrix.txt.gz$annotations), sapply(gensym, length)), gensym = unlist(gensym)) 
# clean_eset(eset, probe2gene_df)
} # }
```
