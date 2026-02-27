# Golimumab Treatment Dataset (GSE92415)

A pre-processed subset of GEO accession GSE92415 containing bulk gene
expression data from ulcerative colitis patients treated with golimumab,
an anti-TNF biologic. Samples represent paired measurements at Week 0
(pre-treatment) and Week 6 (post-treatment) from colonic tissue
isolates. This dataset is intended for use in package examples and
vignettes.

## Usage

``` r
data(golimumab)
```

## Format

A named list with three elements:

- qc_eset:

  A numeric matrix of normalized gene expression values (genes x
  samples), where rows are gene symbols and columns are sample IDs.

- cond:

  A character vector indicating the treatment time point for each
  sample. Values are `"week0"` (pre-treatment) or `"week6"`
  (post-treatment).

- obs_id:

  A character vector of subject IDs used to match paired samples across
  time points for use with
  [`limma::duplicateCorrelation`](https://rdrr.io/pkg/limma/man/dupcor.html).

## Source

GEO accession GSE92415:
<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92415>
