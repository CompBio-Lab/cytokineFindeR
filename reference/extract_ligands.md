# Extract specific ligands and merge results based on a chosen metric

This function extracts specified ligands from a `BenchmarkResults`
object and merges the results into a single data frame using the
specified metric.

## Usage

``` r
extract_ligands(benchmark_results, ligands, metrics = c("padj", "coef"))

process_method_db(df, method_name, db_name, metrics)

summarize_df(results_df, metric)

reshape_metric(df, metric, type)
```

## Arguments

- benchmark_results:

  A BenchmarkResults object containing nested results

- ligands:

  A vector of ligands to filter BenchmarkResults against

- df:

  Input data frame

- method_name:

  Name of the method

- db_name:

  Name of the database

- results_df:

  Combined results data frame

- metric:

  Column name to reshape

- type:

  Type label for the metric

## Value

A data frame containing extracted results for specified ligands

Processed data frame with ranking

Summary data frame in long format

Long format data frame with value and type columns

## Examples

``` r
if (FALSE) { # \dontrun{
results_df <- extract_ligands(
    benchmark_results = results,
    ligands = c("LigandA", "LigandB"),
    metric = "coef")
} # }
```
