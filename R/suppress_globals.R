# Suppress R CMD check notes for dplyr column name variables and package data
# objects used as default arguments or within non-standard evaluation contexts.
utils::globalVariables(c(
  # dplyr / tidyr column names used in NSE (non-standard evaluation)
  ".", "ligand", "value", "method", "database", "rank", "metric_type",
  "type", "genesym",
  # limma topTable column names
  "P.Value", "adj.P.Val", "genes", "logFC",
  # fgsea column name
  "pathway",
  # create_ensemble_results column names
  "study_type", "cytokine", "ligand_tables", "cytosig_table",
  "lri_rank", "cytosig_rank", "ensemble_data",
  # package data objects used as default argument values
  "cytosig_beta"
))
