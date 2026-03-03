#' Given an expression set, make sure that data is 
#' cleaned by gene name using a provided gene list data frame
#'
#' @param eset Expression Set as a numeric matrix (genes x samples).
#' @param gene_list_df A data frame with columns \code{probeids} (probe
#'   identifiers matching rownames of \code{eset}) and \code{gensym}
#'   (gene symbols to map to). Probes mapping to the same gene are averaged.
#'
#' @return cleaned up expression matrix
#' @export
#' @name clean_eset
#'
#' @examples
#' # Create a small expression matrix with probe IDs as rownames
#' set.seed(42)
#' eset <- matrix(rnorm(60), nrow = 6, ncol = 10,
#'                dimnames = list(paste0("probe", 1:6), paste0("S", 1:10)))
#'
#' # Map probes to gene symbols (probes 1-2 both map to GeneA; 5-6 to GeneD)
#' gene_list_df <- data.frame(
#'     probeids = paste0("probe", 1:6),
#'     gensym   = c("GeneA", "GeneA", "GeneB", "GeneC", "GeneD", "GeneD")
#' )
#' result <- clean_eset(eset, gene_list_df)
#' dim(result)

clean_eset <- function(eset, gene_list_df){
  # clean eset against list of probe genes
  # combine probes that bind to multiple genes
  
  X <- eset[gene_list_df$probeids, ] %>% 
    # convert to a data frame to transform probe IDs to genes
    # Take a mean of all probes that match to the Gene ID
    as.data.frame() %>% 
    mutate(genesym = gene_list_df$gensym) %>% 
    group_by(genesym) %>% 
    summarise(across(everything(), ~ mean(.x, na.rm=TRUE)))
  
  eset <- as.matrix(X[,-1])
  rownames(eset) <- X$genesym
  
  return(eset)
  
}
