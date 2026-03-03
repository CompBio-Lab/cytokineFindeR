#' Helper fun to create the design matrix
#'
#' @param treatment Treatment response variable
#' @param obs_id Observation ID: some samples may have unique IDs but come from
#' the same tissue of origin, if that exists, provide a vector of this to make
#' sure the Expression matrix accounts for this to avoid incorrect DEA input.
#' @param eset Expression matrix (genes x samples); required when \code{obs_id}
#' is provided in order to compute duplicate correlations via
#' \code{limma::duplicateCorrelation}.
#'
#' @return Either a design matrix or a design list with components
#' 
#' @import limma
#' @importFrom stats model.matrix
#' @export
#'
#' @examples
#' # Unpaired design
#' treatment <- rep(c("ctrl", "trt"), each = 4)
#' result <- create_design(treatment)
#' result$design
#'
#' # Paired design with duplicate correlation
#' set.seed(42)
#' genes   <- paste0("GENE", 1:20)
#' samples <- paste0("S", 1:8)
#' eset    <- matrix(rnorm(160), nrow = 20, ncol = 8,
#'                  dimnames = list(genes, samples))
#' obs_id  <- rep(1:4, 2)
#' result_paired <- create_design(treatment, obs_id = obs_id, eset = eset)
#' result_paired$design

create_design <- function(treatment, obs_id = NULL, eset = NULL){
  if(is.null(obs_id)) {
    #for unpaired datasets
    design <- model.matrix(~treatment)
    
    } else{
      design <- model.matrix(~treatment)
      dupcor <- duplicateCorrelation(eset, 
                                   design = model.matrix(~treatment), 
                                   block = obs_id)
      message("Within-block correlation: ", dupcor$consensus.correlation)
    }
  
  return(list(design = design, 
              dupcor = if (!is.null(eset) && !is.null(obs_id)) dupcor else NULL))
  }
