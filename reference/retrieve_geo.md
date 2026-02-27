# Retrieve GEO data using the Bioconductor package GEOquery, clean it up, and store in a df list

Retrieve GEO data using the Bioconductor package GEOquery, clean it up,
and store in a df list

## Usage

``` r
retrieve_geo(geo_id)
```

## Arguments

- geo_ids:

  A character string representing the GEO ID

## Value

A combined list of data frames saved by GEO_ID series matrix containing
metadata, the eset, and annotations (if available)

## Examples

``` r
if (FALSE) { # \dontrun{
geo_datasets <- retrieve_geo("GSE179478")
GSE179478$GSE179478_series_matrix.txt.gz$metadata
} # }
```
