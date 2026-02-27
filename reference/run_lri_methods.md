# The core function to run benchmarking of several methods

The core function to run benchmarking of several methods

## Usage

``` r
run_lri_methods(
  eset,
  design,
  dbs,
  methods,
  treatment = NULL,
  obs_id = NULL,
  correlation = NULL
)
```

## Arguments

- eset:

  An expression set (numeric matrix) of genes x samples

- design:

  The design matrix for the data set used to generate the model

- dbs:

  The databases in a list of list format (for package default, use
  dbs_all)

- methods:

  A vector of methods contained in this package

- treatment:

  A vector containing the treatment (specific to the demo data set, this
  is to analyze differentially expressed genes between week 0 and
  week 6) with the drug gollimumab on Ulcerative colitis patients)

- obs_id:

  A vector of sample IDs

- correlation:

  the average estimated inter-duplicate correlation for

## Value

A large BenchmarkResults object containing a nested list of methods and
the results

## Examples

``` r
# This is the core function for running benchmarks
# Basic usage:
if (FALSE) { # \dontrun{
result <- run_lri_methods(eset, design, dbs, methods = c("fgsea", "pca_limma"))
} # }
```
