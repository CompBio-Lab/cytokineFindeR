#' Gene Set Enrichment Analysis
#'
#' @param eset Expression Set object containing gene expression data.
#' @param design Design matrix generated from create_design()
#' @param db Ligand-receptor database
#' @param obs_id A vector of observation IDs to indicate if it's paired data
#' @param correlation Add a correlation block based on the dupcor package for paired analysis
#'
#' @return a table with GSEA results. Each row corresponds to a ligand
#' @export
#'
#' @importFrom fgsea fgsea
#' 
#' @examples
#' set.seed(42)
#' genes   <- paste0("GENE", 1:50)
#' samples <- paste0("S", 1:8)
#' eset    <- matrix(rnorm(400), nrow = 50, ncol = 8,
#'                  dimnames = list(genes, samples))
#' treatment <- rep(c("ctrl", "trt"), each = 4)
#' design  <- model.matrix(~ treatment)
#' rownames(design) <- samples
#' db <- list(LigandA = genes[1:10], LigandB = genes[11:20])
#' result <- cfgsea(eset, design, db)

cfgsea <- function(eset, design, db, 
                   obs_id = NULL, correlation = NULL) {
  # generate linear model from limma for DEA
  # First check if experiment samples are paired
  top <- run_limma(eset, design, obs_id, correlation)
  
  # create named vector of t-stats for each gene
  stats <- top$t
  names(stats) <- top$genes
  
  # Get the list of lengths of receptors (gene sets) 
  # representing each ligand
  length_receptors <- sapply(db, length)
  
  # Run fgsea
  fgsea_results <- fgsea(pathways = db,
                         stats = stats,
                         minSize = min(length_receptors),
                         maxSize = max(length_receptors)
                         )
  fgsea_results <- dplyr::rename(fgsea_results, ligand = pathway)
  return(fgsea_results)
}
