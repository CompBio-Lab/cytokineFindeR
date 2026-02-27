# Helper fun to create the design matrix

Helper fun to create the design matrix

## Usage

``` r
create_design(treatment, obs_id = NULL, eset = NULL)
```

## Arguments

- treatment:

  Treatment response variable

- obs_id:

  Observation ID: some samples may have unique IDs but come from the
  same tissue of origin, if that exists, provide a vector of this to
  make sure the Expression matrix accounts for this to avoid incorrect
  DEA input.

## Value

Either a design matrix or a design list with components

## Examples

``` r
if (FALSE) { # \dontrun{
# data(golimumab)
# For paired datasets (biological replicates), we can use dupcor blocks
# design <- create_design(treatment = golimumab$cond, obs_id = golimumab$obs_id, eset = golimumab$qc_eset)

# If unpaired datasets, we can ignore the obs_id argument
# design <- create_design(treatment = golimumab$cond, eset = golimumab$qc_eset)
} # }
```
