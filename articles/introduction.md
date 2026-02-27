# cytokineFindeR-introduction

The following is a brief introduction on the usage of `CytokineFindeR`.

This vignette will outline the following procedures: 1) Setting up the
design matrix for the linear models 2) Loading in the demo data
`golimumab` from `GSE92415` on GEO. The data set is subset specifically
for ulcerative colitis cell isolates looking at the gene expression
difference for golimumab, which is an anti-TNF drug, between week 0 and
week 6. 3) Running the benchmark analysis in parallel using `future` 4)
Data visualization which shows p-values

## Load the R packages

``` r
library(cytokineFindeR)
#> Setting options('download.file.method.GEOquery'='auto')
#> Setting options('GEOquery.inmemory.gpl'=FALSE)
```

Load demo data (GSE92415) which holds a cleaned up expression set and
meta data subset for the treatment (golimumab) on ulcerative colitus
patients

``` r
data(golimumab)

obs_id <- golimumab$obs_id
treatment <- golimumab$cond

eset <- golimumab$qc_eset
```

## Pre-filtered databases against CytoSig cytokines

The full archive of LRI reference databases is `data("dbs_all")`

For this demo, we will be using a subset containing 3 databases specific
to cytosig cytokines.

``` r
# load the default databases list of lists
data("dbs_cytosig")

# 
dbs <- dbs_cytosig[names(dbs_cytosig) %in% c("baderlab", "lianaplus", "cellchat")]
```

## Set up the design matrix for demo data

``` r
# Input treatment: difference between week 0 and week 6
# input obs_id if paired data (in this case we do have paired data), GSM IDs different for same sample IDs due to different time points

design <- create_design(treatment, obs_id, eset)
#> Within-block correlation: 0.209764916604386
design_mat <- design$design
dupcor <- design$dupcor
```

## run preprocessing step

``` r
preprocess_list <- preprocess_eset(eset, dbs)
eset_f <- preprocess_list$eset_f
dbs_f <- preprocess_list$dbs_f
```

## Call the databases, list out the methods to benchmark, run `cytokinefinder()`

The result should return a nested list of benchmark results

``` r
methods <- c("gsva_limma",
             "pca_limma")

results <- run_lri_methods(eset_f,
                           design_mat,
                           dbs_f,
                           methods,
                           treatment = treatment,
                           obs_id = obs_id,
                           correlation = dupcor$consensus
                          )
```

## To reconstruct the list of benchmarks into a tibble:

`benchlist_to_tbl` will rank the cytokine of interest based on the name
assignment provided in `golimumab_treatment`, in this case, as the
benchmark results are named by `TNF`, the function will output the
ranking of TNF across all method+db combinations.

``` r
golimumab_treatment <- list(TNF = list(benchmarks = results))
golimumab_tibble <- benchlist_to_tbl(results_list = golimumab_treatment,
                                     study_type = "treatment")
golimumab_tibble
```

## CytoSig Custom ridge can be run through the following function:

The `cytosig_custom_ridge` function will run differential expression
using limma and parse the logFC of all genes in the top table.

The ridge regression implementation takes the common genes based on the
beta coefficients from CytoSig and computes the cytokine signatures
accordingly.

``` r
data("cytosig_beta")
cytosig <- cytosig_custom_ridge(eset, 
                     design_mat, 
                     obs_id = obs_id,
                     correlation = dupcor$consensus, 
                     beta_coef = cytosig_beta)
#> fitting model with paired samples.
cytosig
#>                 ligand          coef
#> 1            Activin.A -1.001954e-03
#> 2                 BDNF -6.778455e-03
#> 3                 BMP2  3.905990e-03
#> 4                 BMP4  6.320596e-05
#> 5                 BMP6 -1.412075e-03
#> 6                 BMP7  2.462472e-03
#> 7                 CCL2  6.566563e-03
#> 8                CD40L  7.630706e-03
#> 9               CXCL12 -2.300357e-05
#> 10 Dihydrotestosterone -2.017470e-03
#> 11                 EGF  2.875290e-02
#> 12                 EPO  1.643687e-03
#> 13           Estradiol -8.203645e-03
#> 14                FGF2  2.694235e-03
#> 15                GCSF -5.566465e-03
#> 16               GDF11  3.248006e-03
#> 17               GMCSF -2.503910e-03
#> 18                 HGF -6.671887e-03
#> 19               HMGB1 -4.528049e-04
#> 20                IFN1  3.414549e-03
#> 21                IFNG -4.049134e-04
#> 22                IFNL -3.581396e-03
#> 23                IGF1 -8.871237e-03
#> 24                IL10 -6.302162e-03
#> 25                IL12  1.673317e-03
#> 26                IL13  6.652357e-03
#> 27                IL15 -2.406481e-03
#> 28               IL17A -6.323613e-03
#> 29                IL18  8.895798e-04
#> 30                IL1A -2.298248e-03
#> 31                IL1B -3.505999e-03
#> 32                 IL2  1.805892e-03
#> 33                IL21 -2.170944e-03
#> 34                IL22 -6.243744e-03
#> 35                IL27 -1.114670e-02
#> 36                 IL3 -4.140105e-03
#> 37                IL33  3.598250e-03
#> 38                IL36  1.038565e-02
#> 39                 IL4 -1.658351e-03
#> 40                 IL6  1.476195e-03
#> 41                 INS -1.268250e-03
#> 42                 LIF -4.977668e-03
#> 43                 LTA -1.388151e-03
#> 44                MCSF -3.296464e-03
#> 45            N.3.PUFA -5.857839e-03
#> 46                  NO -4.416151e-03
#> 47                NRG1 -3.086280e-03
#> 48                 OSM -5.394475e-03
#> 49               PDGFB  1.717544e-03
#> 50               PDGFD  2.178155e-03
#> 51                PGE2  1.441978e-03
#> 52       Palmitic.acid  4.216780e-03
#> 53        Progesterone -1.415484e-03
#> 54                TGFA -3.016584e-03
#> 55               TGFB1 -6.223328e-03
#> 56               TGFB3 -8.011084e-04
#> 57                TNFA -1.679128e-02
#> 58               TRAIL  1.786749e-03
#> 59                TSLP -1.957181e-03
#> 60               TWEAK -2.537485e-03
#> 61        Testosterone  3.860417e-03
#> 62               VEGFA  8.014912e-03
#> 63               WNT3A -1.439883e-03
#> 64               WNT5A -3.481603e-03
```
