#' The core function to run benchmarking of several methods
#'
#' @param eset An expression set (numeric matrix) of genes x samples
#' @param design The design matrix for the data set used to generate the model
#' @param dbs The databases in a list of list format (for package default, use dbs_all)
#' @param methods A character vector of method names to benchmark. Must be one
#'   or more of: `"gsva_limma"`, `"pca_limma"`, `"cfgsea"`,
#'   `"gsva_plsda"`, `"pca_plsda"`, `"cytosig_custom_ridge"`.
#' @param treatment A vector containing the treatment (specific to the demo
#' data set, this is to analyze differentially expressed genes between week 0
#' and week 6) with the drug gollimumab on Ulcerative colitis patients)
#' @param obs_id A vector of sample IDs
#' @param correlation the average estimated inter-duplicate correlation for
#' @param verbose Logical; if `TRUE`, progress messages are printed for each
#'   method-database combination. Default is `FALSE`.
#'
#' @return A large BenchmarkResults object containing a nested list of methods
#' and the results
#' @export
#'
#' @importFrom future plan multicore multisession
#' @importFrom future.apply future_lapply future_sapply
#' @importFrom stats setNames
#' @examples
#' \donttest{
#' set.seed(42)
#' genes   <- paste0("GENE", 1:50)
#' samples <- paste0("S", 1:8)
#' eset    <- matrix(rnorm(400), nrow = 50, ncol = 8,
#'                  dimnames = list(genes, samples))
#' treatment <- rep(c("ctrl", "trt"), each = 4)
#' design  <- model.matrix(~ treatment)
#' rownames(design) <- samples
#' dbs     <- list(db1 = list(LigandA = genes[1:10], LigandB = genes[11:20]))
#' result  <- run_lri_methods(eset, design, dbs, methods = "cfgsea",
#'                            treatment = treatment)
#' }

run_lri_methods <- function(eset, design, dbs, methods,
                            treatment = NULL, obs_id = NULL,
                            correlation = NULL, verbose = FALSE) {
  method_registry <- list(
    gsva_limma           = gsva_limma,
    pca_limma            = pca_limma,
    cfgsea               = cfgsea,
    gsva_plsda           = gsva_plsda,
    pca_plsda            = pca_plsda,
    cytosig_custom_ridge = cytosig_custom_ridge
  )

  methods <- match.arg(methods, choices = names(method_registry),
                       several.ok = TRUE)

  # Set up the future plan (multicore on Unix, multisession on Windows)
  # Falls back to sequential if parallel plan cannot be initialised (e.g. in
  # test environments or nested subprocess contexts)
  tryCatch(
    future::plan(if (.Platform$OS.type == "unix") future::multicore else future::multisession),
    error = function(e) {
      message("Parallel execution unavailable, using sequential plan: ", e$message)
      future::plan(future::sequential)
    }
  )

  # Initialize an empty list to store results
  results <- list()

  # Iterate over each method in the methods list
  for (method_name in methods) {
    method <- method_registry[[method_name]]

    # Parallel execution for the current method
    method_results <- future_lapply(names(dbs), function(database) {
      if (verbose)
        message("Processing method: ", method_name, " with database: ", database)

      # Check if the method is PLSDA-based:
      if (grepl("plsda", method_name)) {
        result <- method(eset, treatment, dbs[[database]])
      } else {
        if (!is.null(obs_id)) {
          result <- method(eset, design,
                           dbs[[database]],
                           obs_id = obs_id,
                           correlation = correlation)
        } else {
          result <- method(eset, design, dbs[[database]])
        }
      }

      if (verbose)
        message("Finished processing method: ", method_name, " with database: ", database)
      return(list(database = database, result = result))
    }, future.seed = TRUE)

    if (verbose) message("Combining into one BenchmarkResults Object")
    method_results_named <- setNames(
      lapply(method_results, `[[`, "result"),
      sapply(method_results, `[[`, "database")
    )
    results[[method_name]] <- method_results_named
  }

  # Create and return results object
  return(structure(results, class = "BenchmarkResults"))
}