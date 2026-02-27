# Run CytoSig Custom Ridge Regression function

Run CytoSig Custom Ridge Regression function

## Usage

``` r
cytosig_custom_ridge(
  eset,
  design,
  obs_id = NULL,
  correlation = NULL,
  beta_coef = cytosig_beta
)
```

## Arguments

- eset:

  Expression matrix of microarray or RNAseq experiments

- design:

  design matrix list that indicates conditions to compare and if paired
  design

- obs_id:

  if paired design, indicate observation ID to map biological replicates
  to sample of origin if paired

- correlation:

  if paired design, indicate correlation blocks for paired data

- beta_coef:

  beta matrix for ridge regression. Defaults to cytosig_beta included in
  package.

## Value

CytoSig ridge regression results table

## Details

Cytokine Signaling Analyzer or CytoSig (Jiang et al, 2021) is a ridge
regression predictive model of cytokine signaling cascades trained on
20,591 transcriptomic profiles of perturbations relevant to cytokine
activity. The following is our adoption of CytoSig using the beta
coefficients published using logFC results generated differential
expression through our run_limma function based on the standard of
practice, which can be done through cytokineFindeR's functions.

## Examples

``` r
if (FALSE) { # \dontrun{
# Ridge regression with default CytoSig beta coefficients
result <- cytosig_custom_ridge(eset, design)

# Ridge regression with custom beta coefficients
result <- cytosig_custom_ridge(eset, design, beta_coef = custom_beta)
} # }
```
