#' Perform PCA and fit a multivariate regression to measure ligand activity
#'
#' @param eset Expression Set object containing gene expression data.
#' @param treatment Treatment response variable
#' @param db ligand-receptor database
#'
#' @return Table of ranked ligands by coef order (largest to smallest)
#' @export
#'
#' @examples
#' set.seed(42)
#' genes     <- paste0("GENE", 1:20)
#' eset      <- matrix(rnorm(160), nrow = 20, ncol = 8,
#'                     dimnames = list(genes, paste0("S", 1:8)))
#' treatment <- rep(c("ctrl", "trt"), each = 4)
#' db        <- list(LigandA = genes[1:5], LigandB = genes[6:10])
#' result    <- pca_plsda(eset, treatment, db)

pca_plsda <- function(eset, treatment, db){
  pc <- lapply(db, function(ligand){
    tryCatch({
      genexp <- t(eset[intersect(rownames(eset), ligand), , drop=FALSE])
      prcomp(genexp, center = TRUE, scale. = TRUE, rank. = 1)$x[, "PC1"]  
    }, error = function(e) NA)
  })
  # remove 0 variability ligands
  pc <- pc[sapply(pc, length) > 1] %>% 
    do.call(rbind, .)
  # fit to plsda
  fit <- mixOmics::plsda(t(pc), treatment)
  coef <- abs(mixOmics::selectVar(fit, comp=1)$value$value.var)
  names(coef) <- rownames(mixOmics::selectVar(fit, comp=1)$value)
  return(enframe(coef[order(coef, decreasing = TRUE)],
                 name = "ligand",
                 value = "coef")
         )
}
