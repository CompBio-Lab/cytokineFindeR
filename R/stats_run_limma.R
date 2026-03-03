#' Run limma in unpaired or paired mode
#'
#' @param eset expression matrix
#' @param design design matrix using the create_design function which indicates blocks for paired samples
#' @param obs_id required to indicate the blocks that map biological replicates to the same sample
#' @param correlation the average estimated inter-duplicate correlation. The average is the trimmed mean of the individual correlations on the atanh-transformed scale. 
#'
#' @returns top table that runs differential gene expression analysis. Main purpose for this is to get logFC of all genes for CytoSig input.
#' @export
#'
#' @importFrom limma lmFit
#' @importFrom limma eBayes
#' @importFrom limma topTable
#' @importFrom tibble rownames_to_column
#'
#' @examples
#' set.seed(42)
#' genes   <- paste0("GENE", 1:20)
#' samples <- paste0("S", 1:8)
#' eset    <- matrix(rnorm(160), nrow = 20, ncol = 8,
#'                  dimnames = list(genes, samples))
#' treatment <- rep(c("ctrl", "trt"), each = 4)
#' design  <- model.matrix(~ treatment)
#' rownames(design) <- samples
#' result  <- run_limma(eset, design)

run_limma <- function (eset, design, obs_id = NULL, correlation = NULL) 
{
  if (!is.null(obs_id)) {
    fit <- limma::lmFit(eset, design, block = obs_id, correlation = correlation)
    message("fitting model with paired samples.")
  }
  else {
    fit <- limma::lmFit(eset, design)
    message("fitting model without paired sample consideration.")
  }
  efit <- limma::eBayes(fit)
  top <- limma::topTable(efit, coef = 2, number = nrow(efit)) %>%
    rownames_to_column(var = "genes")
  return(top)
}
