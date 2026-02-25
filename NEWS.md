# cytokineFindeR 0.99.0

* Initial Bioconductor submission.
* Implements seven statistical methods for ligand activity inference from bulk transcriptomics: `gsva_limma`, `pca_limma`, `cfgsea`, `run_limma`, `gsva_plsda`, `pca_plsda`, and `cytosig_custom_ridge`.
* Includes curated ligand-receptor interaction (LRI) databases: `dbs_all`, `dbs_subset`, `dbs_cytosig`.
* Includes CytoSig beta coefficient matrix (`cytosig_beta`) for ridge regression.
* Includes golimumab demo dataset (`golimumab`) derived from GEO accession GSE92415.
* Supports paired and unpaired experimental designs via `limma::duplicateCorrelation`.
* Parallel execution via `future` and `future.apply`.
* Rank-based ensemble scoring combining LRI-based methods with CytoSig via `create_ensemble_results()`.
