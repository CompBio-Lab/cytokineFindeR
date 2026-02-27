# Create ensemble results combining LRI-based and CytoSig methods

This function combines results from LRI-based methods with CytoSig
results using a rank-based ensemble approach. It performs an inner join
between LRI and CytoSig results and computes ensemble rankings based on
p-values.

## Usage

``` r
create_ensemble_results(
  master_tbl,
  ensemble_method = "mean_rank",
  pval_col_lri = "pval",
  pval_col_cytosig = "pval"
)
```

## Arguments

- master_tbl:

  A tibble containing benchmark results from benchlist_to_tbl(), with
  columns: study_type, cytokine, method, database, class, ligand_tables

- ensemble_method:

  Character string specifying ensemble method. Currently supports
  "mean_rank" (default)

- pval_col_lri:

  Character string specifying the p-value column name in LRI method
  results. Default is "pval"

- pval_col_cytosig:

  Character string specifying the p-value column name in CytoSig method
  results. Default is "pval"

## Value

A tibble with additional ensemble columns:

- overlap_count: Number of overlapping ligands between LRI and CytoSig

- ensemble_table: Combined LRI and CytoSig results with rankings

- ensemble_rank: Mean rank score combining LRI and CytoSig rankings

- lri_rank: Percentile rank from LRI method (100 = best, 0 = worst)

- cytosig_rank: Percentile rank from CytoSig method (100 = best, 0 =
  worst)

## Details

The function:

1.  Filters for LRI-based methods (class == "LRI")

2.  Joins with CytoSig results (class == "CytoSig_Web") by study_type
    and cytokine

3.  Computes overlap between LRI and CytoSig ligand sets

4.  Creates ensemble rankings using percentile ranks (lower p-value =
    higher rank)

5.  Falls back to CytoSig-only results when no LRI data is available

Ranking system: Converts p-values to percentile ranks where 100
represents the best (lowest p-value) and 0 represents the worst (highest
p-value).

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming you have benchmark results
results <- run_cytokinefinder_workflow(study_data, databases, methods)

# Convert to tibble format  
master_tbl <- benchlist_to_tbl(results$benchmarks, "my_study", FALSE)

# Create ensemble results
ensemble_results <- create_ensemble_results(
  master_tbl, 
  ensemble_method = "mean_rank",
  pval_col_lri = "pval",
  pval_col_cytosig = "pval"
)

# View ensemble rankings
ensemble_results %>% 
  select(cytokine, method, database, ensemble_rank, lri_rank, cytosig_rank)
} # }
```
