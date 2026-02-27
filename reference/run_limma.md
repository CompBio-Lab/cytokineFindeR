# Run limma in unpaired or paired mode

Run limma in unpaired or paired mode

## Usage

``` r
run_limma(eset, design, obs_id = NULL, correlation = NULL)
```

## Arguments

- eset:

  expression matrix

- design:

  design matrix using the create_design function which indicates blocks
  for paired samples

- obs_id:

  required to indicate the blocks that map biological replicates to the
  same sample

- correlation:

  the average estimated inter-duplicate correlation. The average is the
  trimmed mean of the individual correlations on the atanh-transformed
  scale.

## Value

top table that runs differential gene expression analysis. Main purpose
for this is to get logFC of all genes for CytoSig input.

## Examples

``` r
if (FALSE) { # \dontrun{
# Load up the datasets
# pert_dats <- readRDS("pert_dats.rds")
# runs all pert data using limma for CytoSig
# cytosig_pert_input <- lapply(pert_dats, run_limma)
} # }
```
