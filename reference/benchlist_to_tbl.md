# Convert nested benchmark results to tibble format

Convert nested benchmark results to tibble format

## Usage

``` r
benchlist_to_tbl(results_list, study_type, has_benchmarks_layer = TRUE)
```

## Arguments

- results_list:

  Nested list of benchmark results (cytokine -\> benchmarks -\> method
  -\> database -\> data) or CytoSig Results (cytokine -\> method -\>
  data)

- study_type:

  Character string identifying the study type

- has_benchmarks_layer:

  boolean to handle benchmark and cytosig structures

## Value

A tibble with flattened benchmark results Convert nested results to
tibble format (handles both benchmark and cytosig structures)
