#' Retrieve GEO data using the Bioconductor package GEOquery, clean it up, and store in a df list
#'
#' @param geo_ids A character string representing the GEO ID
#'
#' @return A combined list of data frames saved by GEO_ID series matrix containing metadata, the eset, and annotations (if available)
#' @export
#' @name retrieve_geo
#'
#' @examples
#' \dontrun{
#' geo_datasets <- retrieve_geo("GSE179478")
#' GSE179478$GSE179478_series_matrix.txt.gz$metadata
#' }
#'
#' @importFrom GEOquery getGEO
#' @importFrom Biobase pData
#' @importFrom Biobase exprs

retrieve_geo <- function(geo_id){
  geo_data <- tryCatch({
    getGEO(geo_id, GSEMatrix = TRUE)
  }, error = function(e) {
    message("Error in fetching GEO data for ID: ", geo_id, 
            " - ", e$message)
    return(NULL)
  })
  e1 <- geo_data
  # handle retrieval if e1 returns a list
  combined_data <- lapply(e1, function(i){
    metadata <- pData(i)
    eset <- exprs(i)
    annotations <- i@featureData@data
    annotations <- annotations[annotations$`Gene Symbol` != "", ]
    list(metadata = metadata, eset = eset, annotations = annotations)
  })
  if (any(sum(rapply(combined_data, nrow)) == 0)) {
    message("Warning: One or more data frames           
            in the combined data list are empty.")
  }
  return(combined_data)
}