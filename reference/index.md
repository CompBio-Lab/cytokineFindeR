# Package index

## All functions

- [`benchlist_to_tbl()`](https://compbio-lab.github.io/cytokineFinder/reference/benchlist_to_tbl.md)
  : Convert nested benchmark results to tibble format
- [`cfgsea()`](https://compbio-lab.github.io/cytokineFinder/reference/cfgsea.md)
  : Gene Set Enrichment Analysis
- [`clean_eset()`](https://compbio-lab.github.io/cytokineFinder/reference/clean_eset.md)
  : Given an expression set, make sure that data is cleaned by gene name
  using a provided gene list data frame
- [`create_db_space()`](https://compbio-lab.github.io/cytokineFinder/reference/create_db_space.md)
  : Create an LRI database space
- [`create_design()`](https://compbio-lab.github.io/cytokineFinder/reference/create_design.md)
  : Helper fun to create the design matrix
- [`create_ensemble_results()`](https://compbio-lab.github.io/cytokineFinder/reference/create_ensemble_results.md)
  : Create ensemble results combining LRI-based and CytoSig methods
- [`cytosig_beta`](https://compbio-lab.github.io/cytokineFinder/reference/cytosig_beta.md)
  : CytoSig Betas
- [`cytosig_custom_ridge()`](https://compbio-lab.github.io/cytokineFinder/reference/cytosig_custom_ridge.md)
  : Run CytoSig Custom Ridge Regression function
- [`dbs_all`](https://compbio-lab.github.io/cytokineFinder/reference/dbs_all.md)
  : Database List unfiltered
- [`dbs_cytosig`](https://compbio-lab.github.io/cytokineFinder/reference/dbs_cytosig.md)
  : Database List to make comparable benchmarks against CytoSig
- [`dbs_subset`](https://compbio-lab.github.io/cytokineFinder/reference/dbs_subset.md)
  : Small Database example
- [`extract_ligands()`](https://compbio-lab.github.io/cytokineFinder/reference/extract_ligands.md)
  [`process_method_db()`](https://compbio-lab.github.io/cytokineFinder/reference/extract_ligands.md)
  [`summarize_df()`](https://compbio-lab.github.io/cytokineFinder/reference/extract_ligands.md)
  [`reshape_metric()`](https://compbio-lab.github.io/cytokineFinder/reference/extract_ligands.md)
  : Extract specific ligands and merge results based on a chosen metric
- [`golimumab`](https://compbio-lab.github.io/cytokineFinder/reference/golimumab.md)
  : Golimumab Dataset
- [`gsva_limma()`](https://compbio-lab.github.io/cytokineFinder/reference/gsva_limma.md)
  : Gene Set Variation Analysis for cytokines (cGSVA)
- [`gsva_plsda()`](https://compbio-lab.github.io/cytokineFinder/reference/gsva_plsda.md)
  : Calculate coefficients from variable selection using GSVA and PLS-DA
- [`pca_limma()`](https://compbio-lab.github.io/cytokineFinder/reference/pca_limma.md)
  : Calculate top ligands using a PCA approach from receptor genes given
  a database
- [`pca_plsda()`](https://compbio-lab.github.io/cytokineFinder/reference/pca_plsda.md)
  : Perform PCA and fit a multivariate regression to measure ligand
  activity
- [`preprocess_eset()`](https://compbio-lab.github.io/cytokineFinder/reference/preprocess_eset.md)
  : Preprocess the eset
- [`retrieve_geo()`](https://compbio-lab.github.io/cytokineFinder/reference/retrieve_geo.md)
  : Retrieve GEO data using the Bioconductor package GEOquery, clean it
  up, and store in a df list
- [`run_cytokinefinder()`](https://compbio-lab.github.io/cytokineFinder/reference/run_cytokinefinder.md)
  : Wrapper function to run complete cytokinefinder benchmarking
  workflow with preprocessing
- [`run_limma()`](https://compbio-lab.github.io/cytokineFinder/reference/run_limma.md)
  : Run limma in unpaired or paired mode
- [`run_lri_methods()`](https://compbio-lab.github.io/cytokineFinder/reference/run_lri_methods.md)
  : The core function to run benchmarking of several methods
