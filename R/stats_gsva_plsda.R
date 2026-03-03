#' Calculate coefficients from variable selection using GSVA and PLS-DA
#' 
#' Uses a classification approach when computing receptor weights with
#' Partial Least Squares - Discriminant Analysis and pools these receptors
#' together to compute a value for the ligand or cytokine given a gene set.
#' As limma treats all variables as equal importance, PLSDA provides a multiple
#' regression approach
#'
#' @param eset Expression Set object containing gene expression data.
#' @param treatment Treatment response variable
#' @param db ligand-receptor database
#'
#' @return A named vector of ligands indicating importance 
#' @export
#' 
#' @importFrom GSVA gsvaParam
#' @importFrom GSVA gsva
#' @importFrom mixOmics plsda
#' @importFrom mixOmics selectVar
#' 
#' @examples
#' \donttest{
#' set.seed(42)
#' genes     <- paste0("GENE", 1:50)
#' eset      <- matrix(rnorm(400), nrow = 50, ncol = 8,
#'                     dimnames = list(genes, paste0("S", 1:8)))
#' treatment <- rep(c("ctrl", "trt"), each = 4)
#' db        <- list(LigandA = genes[1:10], LigandB = genes[11:20])
#' result    <- gsva_plsda(eset, treatment, db)
#' }
 
gsva_plsda <- function(eset, treatment, db){
  length_receptors <- sapply(db, length)
  
  gsvapar <- gsvaParam(eset, 
                       db, 
                       minSize = min(length_receptors), 
                       maxDiff = TRUE)
  gsva_eset <- gsva(gsvapar)
  
  # Use mixomics to fit regression
  # make sure to transpose matrix from pxn to nxp
  fit <- plsda(t(gsva_eset), treatment)
  coef <- abs(selectVar(fit, comp=1)$value$value.var)
  names(coef) <- rownames(selectVar(fit, comp=1)$value)
  return(enframe(coef[order(coef, decreasing = TRUE)], 
                 name = "ligand", 
                 value = "coef"))
  }