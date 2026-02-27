# Wrapper function to run complete cytokinefinder benchmarking workflow with preprocessing

Wrapper function to run complete cytokinefinder benchmarking workflow
with preprocessing

## Usage

``` r
run_cytokinefinder(study_data, databases, methods)
```

## Arguments

- study_data:

  Named list containing study data with required elements:

  - qc_eset: Quality-controlled expression set matrix

  - cond: Treatment/condition vector

  - obs_id: Sample IDs for paired experiments (optional)

- databases:

  List of LRI databases for benchmarking

- methods:

  Vector of method names to benchmark

## Value

Original study_data with added 'benchmarks', 'design', and
'preprocessing' elements

## Examples

``` r
# Single study example
study_data <- list(
  qc_eset = matrix(rnorm(1000), nrow = 100, ncol = 10),
  cond = rep(c("control", "treatment"), each = 5)
)

# With paired samples
study_data_paired <- list(
  qc_eset = matrix(rnorm(1000), nrow = 100, ncol = 10),
  cond = rep(c("control", "treatment"), each = 5),
  obs_id = rep(1:5, 2)
)

if (FALSE) { # \dontrun{
# Load databases and methods
data(dbs_all)  # or however you load your databases
methods <- c("gsea", "ssgsea", "cytosig_custom_ridge")

# Run workflow
results <- run_cytokinefinder_workflow(study_data, dbs_all, methods)

# Access results
benchmark_results <- results$benchmarks
design_matrix <- results$design
} # }
```
