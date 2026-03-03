#' Gene Set Variation Analysis for cytokines (cGSVA)
#'
#' @param eset Expression Set object containing gene expression data.
#' @param design Design matrix generated from create_design()
#' @param db Ligand-receptor database
#' @param obs_id Optional: provide a vector of sample IDs making sure the order matches with the eset
#' @param correlation Optional: input the correlation consensus between the samples to evaluate if it is paired data
#'
#' @return top table of all ligands
#' @export
#'
#' @importFrom limma lmFit
#' @importFrom limma eBayes
#' @importFrom limma topTable
#' @importFrom GSVA gsvaParam
#' @importFrom GSVA gsva
#' @importFrom tibble enframe
#' @importFrom tibble rownames_to_column
#' 
#' @examples
#' \donttest{
#' set.seed(42)
#' genes   <- paste0("GENE", 1:50)
#' samples <- paste0("S", 1:8)
#' eset    <- matrix(rnorm(400), nrow = 50, ncol = 8,
#'                  dimnames = list(genes, samples))
#' treatment <- rep(c("ctrl", "trt"), each = 4)
#' design  <- model.matrix(~ treatment)
#' rownames(design) <- samples
#' db      <- list(LigandA = genes[1:10], LigandB = genes[11:20])
#' result  <- gsva_limma(eset, design, db)
#' }

gsva_limma <- function(eset, design, 
                       db, obs_id = NULL, 
                       correlation = NULL) {
  
  length_receptors <- sapply(db, length)
  
  gsvapar <- gsvaParam(eset, 
                       db, 
                       minSize = min(length_receptors), 
                       maxDiff = TRUE)
  gsva_eset <- gsva(gsvapar)
  
  # Run DEA
  # First check if experiment samples are paired
  if(!is.null(obs_id)) {
    fit <- eBayes(lmFit(gsva_eset, 
                        design, 
                        block = obs_id, 
                        correlation = correlation))
    message("fitting model with paired samples.")
  } else {
    fit <- eBayes(lmFit(gsva_eset, design))
    message("fitting model without paired sample consideration.")
  }
  # Get top table
  top <- topTable(fit, coef = 2, number = nrow(fit)) %>%
    rownames_to_column("ligand") %>%
    dplyr::rename(pval = P.Value, padj = adj.P.Val)
  
  return(top)
}