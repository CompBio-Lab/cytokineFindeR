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
#' \dontrun{
#' # create a small gene_symbol and ensembl_id df from the GEO dataset example (golimumab TNF-targeted treatment) if annotation df exists
#' # gensym <- sapply(strsplit(golimumab$GSE92415_series_matrix.txt.gz$annotations$`Gene Symbol`, "///"), trimws) 
#' # probe2gene_df <- tibble(probeids = rep(rownames(golimumab$GSE92415_series_matrix.txt.gz$annotations), sapply(gensym, length)), gensym = unlist(gensym)) 
#' # clean_eset(eset, probe2gene_df)
#' }

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
