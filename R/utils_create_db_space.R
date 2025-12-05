#' Create an LRI database space
#' @details This is a function that saves the database list as an RDA in case future contributors wanted to create other variations of the LRI database compendium 
#'
#' @param ... 
#' @param filePath 
#' @param saveToFile 
#'
#' @return Saves the database list
#' @export
#'

create_db_space <- function(...,
                            filePath = "data/ligand_receptor_db.rda",
                            saveToFile = TRUE) {
  # Check if the data directory exists might need to come back to this
  data_dir <- dirname(filePath)
  if (!dir.exists(data_dir)) {
    stop(paste("The directory",
               data_dir,
               "does not exist. Please create it first."))
  }

  # Save the database list to the file if the save_to_file parameter is TRUE
  if (saveToFile) {
    # Capture the objects passed via ...
    objs <- list(...)
    # Save the list to the specified file path
    save(list = names(objs), file = filePath, envir = list2env(objs))
    message("Database list saved to ", filePath)
  } else {
    message("Database list not saved.")
  }
}