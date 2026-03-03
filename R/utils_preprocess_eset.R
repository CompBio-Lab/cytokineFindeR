#' Preprocess the eset
#'
#' @param eset the input eset typically the raw gene expression matrix
#' @param dbs the list of databases in the list of lists that can be called from the package's internal data
#'
#' @return The filtered expression set that contains 0 variability ligands
#' @export
#'
#' @examples
#' set.seed(42)
#' genes <- paste0("GENE", 1:20)
#' eset  <- matrix(rnorm(160), nrow = 20, ncol = 8,
#'                 dimnames = list(genes, paste0("S", 1:8)))
#' dbs   <- list(db1 = list(LigandA = genes[1:5], LigandB = genes[6:10]))
#' result <- preprocess_eset(eset, dbs)
#' names(result)

preprocess_eset <- function(eset, dbs) {
  # filter the dbs against the eset 
  filtered_dbs <- filter_db_against_eset(eset, dbs)
  
  # check 0 variable ligands from dbs and intersect with eset 
  # Run PCA to get the first PC
  filtered_dbs_zeroVariance <- remove_zero_variance_ligands(eset, filtered_dbs)
  
  unique_receptors <- unique(unlist(filtered_dbs_zeroVariance))
  filtered_eset <- eset[intersect(rownames(eset), unique_receptors),]

  return(list(eset_f = filtered_eset, 
              dbs_f = filtered_dbs_zeroVariance))
}

#' Filter database and extract ligands' receptors list that match the eset data (internal)
#'
#' @description
#' Provide a list of databases and match against receptors of interest.
#' Remove ligands that do not have matching receptors in the eset.
#'
#' @details
#' Iterate through the expression matrix input and match a vector of genes
#' of interest.
#'
#' @param eset A numeric matrix that represents the expression set
#' @param dbs A nested list of lists containing the databases, the ligand genes,
#' and the receptor genes for each ligand gene
#' @keywords internal
#' @noRd
#'
#' @return A list of lists showing a set of ligand genes where each ligand
#' contains a list of genes that were subset defined as "receptor" genes

filter_db_against_eset <- function(eset, dbs) {
  # For all databases, search for ligand that matched cytokine and eset
  genes <- rownames(eset)
  
  # Iterate through each database and subset based on receptor genes available
  # in the eset
  output_db <- lapply(dbs, function(db) {
    # Update list of receptors against the genes
    db <- lapply(db, function(receptors) { 
      intersect(receptors, genes) 
    })
    db[sapply(db, length) > 0]
  })
  return(output_db)
}


#' Remove zero variance ligands from the database list of list (internal)
#'
#' @description
#' Filters out ligands that have zero variance across samples using PCA.
#' Ligands with insufficient variance cannot be meaningfully analyzed.
#'
#' @param eset A numeric matrix representing the expression set
#' @param dbs A nested list of lists containing the databases, the ligand genes,
#' and the receptor genes for each ligand gene
#' @keywords internal
#' @noRd
#'
#' @return A list of databases with ligands removed if there is 0 variance based on the eset

remove_zero_variance_ligands <- function(eset, dbs) {
  dbs_list_removedZeroVariance <- 
    lapply(dbs, function(db) {
      pc <- lapply(db, function(ligand){
        tryCatch({
          genexp <- t(eset[intersect(rownames(eset), ligand), , drop=FALSE])
          prcomp(genexp, center = TRUE, scale. = TRUE, rank. = 1)$x[, "PC1"]
        }, error = function(e) NA)
      })
      db[sapply(pc, length) > 1]
    })
  # remove 0 variability ligands
  return(dbs_list_removedZeroVariance)
}
