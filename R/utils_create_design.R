#' Helper fun to create the design matrix
#'
#' @param treatment Treatment response variable
#' @param obs_id Observation ID: some samples may have unique IDs but come from
#' the same tissue of origin, if that exists, provide a vector of this to make
#' sure the Expression matrix accounts for this to avoid incorrect DEA input.
#' 
#' @return Either a design matrix or a design list with components
#' 
#' @importFrom limma duplicateCorrelation
#' @export
#'
#' @examples
#' \dontrun{
#' # data(golimumab)
#' # For paired datasets (biological replicates), we can use dupcor blocks
#' # design <- create_design(treatment = golimumab$cond, obs_id = golimumab$obs_id, eset = golimumab$qc_eset)
#' 
#' # If unpaired datasets, we can ignore the obs_id argument
#' # design <- create_design(treatment = golimumab$cond, eset = golimumab$qc_eset)
#' }

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
