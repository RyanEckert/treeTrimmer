#' Prune Clones from Dataset
#'
#' Takes the output from \code{\link{tree_trimmer}} and a metadata table,
#' cuts the tree at the calibrated height, and selects one representative
#' sample per clonal cluster. Optionally selects the highest quality sample
#' per cluster if a quality column is provided.
#'
#' @param tree_result A \code{tree_trimmer} object (output from \code{\link{tree_trimmer}})
#' @param data A data frame containing sample metadata. Must include a column
#'   with sample names that match the labels in the original tree.
#' @param sample Character string specifying the column name in \code{data}
#'   that contains sample identifiers (must match tree labels).
#' @param quality Character string specifying the column name in \code{data}
#'   that contains quality scores, or \code{NULL} (default). If provided,
#'   the highest quality sample is selected from each cluster. If \code{NULL},
#'   the first sample (alphabetically) is selected.
#' @param suffix Character string to append to each sample name (default = "").
#'   Useful for creating file lists (e.g., \code{suffix = ".bam"} to generate
#'   BAM file names).
#'
#' @return A character vector of retained sample names (with suffix if provided).
#'   Length equals the number of unique clonal groups detected. Can be written
#'   directly to a file with \code{writeLines()}.
#'
#' @details
#' This function implements the second step of clone removal after calibrating
#' the tree cut height with \code{\link{tree_trimmer}}. It:
#' \enumerate{
#'   \item Cuts the tree at the calibrated height from \code{tree_result$height}
#'   \item Groups samples into clonal clusters
#'   \item Selects one representative per cluster based on quality (if provided)
#'      or alphabetically (if no quality column)
#'   \item Returns sample names with optional suffix
#' }
#'
#' The function preserves the order of samples as they appear in the tree,
#' making it easier to match with downstream analyses.
#'
#' @section Quality Selection:
#' When \code{quality} is specified:
#' \itemize{
#'   \item Higher values are assumed to be better quality
#'   \item The sample with maximum quality in each cluster is retained
#'   \item If multiple samples have identical max quality, the first
#'         (alphabetically) is selected
#' }
#'
#' When \code{quality = NULL}:
#' \itemize{
#'   \item The first sample (alphabetically) in each cluster is retained
#'   \item This provides deterministic behavior even without quality scores
#' }
#'
#' @examples
#' # Load example data
#' data_file <- system.file("extdata", "coral_ibs.csv", package = "treeTrimmer")
#' coral_ibs <- as.matrix(read.csv(data_file, row.names = 1))
#' 
#' # Load metadata
#' meta_file <- system.file("extdata", "coral_metadata.csv", package = "treeTrimmer")
#' coral_metadata <- read.csv(meta_file)
#'
#' # Build tree and calibrate threshold
#' dist_mat <- as.dist(1 - coral_ibs)
#' hc <- hclust(dist_mat, method = "average")
#' 
#' tech_reps <- list(
#'   OGM024 = c("OGM024.1", "OGM024.2", "OGM024.3"),
#'   OGM066 = c("OGM066.1", "OGM066.2", "OGM066.3"),
#'   OGM108 = c("OGM108.1", "OGM108.2", "OGM108.3")
#' )
#' 
#' result <- tree_trimmer(hc, tech_reps, verbose = FALSE)
#'
#' # Example 1: Simple pruning (no quality, just sample names)
#' retained <- prune_clones(result, 
#'                         data = coral_metadata,
#'                         sample = "sample")
#' head(retained)
#'
#' # Example 2: Quality-based selection
#' retained <- prune_clones(result,
#'                         data = coral_metadata,
#'                         sample = "sample",
#'                         quality = "quality_score")
#' head(retained)
#'
#' # Example 3: Generate BAM file list
#' bam_list <- prune_clones(result,
#'                         data = coral_metadata,
#'                         sample = "sample",
#'                         quality = "quality_score",
#'                         suffix = ".trim.un.bt2.bam")
#' head(bam_list)
#'
#' # Write to file
#' # writeLines(bam_list, "bamsNoClones.txt")
#'
#' @seealso
#' \code{\link{tree_trimmer}} for calibrating the tree cut height,
#' \code{\link[stats]{cutree}} for the underlying tree cutting
#'
#' @importFrom dplyr group_by arrange desc slice ungroup
#' @export
prune_clones <- function(tree_result, data, sample, quality = NULL, suffix = "") {
  
  # ---- Input validation ----
  if (!inherits(tree_result, "tree_trimmer")) {
    stop("tree_result must be a tree_trimmer object from tree_trimmer()",
         call. = FALSE)
  }
  
  if (!is.data.frame(data)) {
    stop("data must be a data frame",
         call. = FALSE)
  }
  
  if (!is.character(sample) || length(sample) != 1) {
    stop("sample must be a single character string specifying the column name",
         call. = FALSE)
  }
  
  if (!sample %in% names(data)) {
    stop("Column '", sample, "' not found in data",
         call. = FALSE)
  }
  
  if (!is.null(quality)) {
    if (!is.character(quality) || length(quality) != 1) {
      stop("quality must be a single character string or NULL",
           call. = FALSE)
    }
    if (!quality %in% names(data)) {
      stop("Column '", quality, "' not found in data",
           call. = FALSE)
    }
  }
  
  if (!is.character(suffix) || length(suffix) != 1) {
    stop("suffix must be a single character string",
         call. = FALSE)
  }
  
  if (is.na(tree_result$height)) {
    stop("tree_result$height is NA - tree_trimmer() did not find a valid threshold",
         call. = FALSE)
  }
  
  # ---- Extract cluster assignments ----
  # Get the original hclust object
  if (!"tree" %in% names(tree_result)) {
    stop("tree_result does not contain the tree object. ",
         "Make sure you're using tree_trimmer() from treeTrimmer v1.0.0 or later.",
         call. = FALSE)
  }
  
  tree_obj <- tree_result$tree
  
  if (!inherits(tree_obj, "hclust")) {
    stop("tree_result$tree is not an hclust object",
         call. = FALSE)
  }
  
  # Cut the tree at the calibrated height
  clusters <- cutree(tree_obj, h = tree_result$height)
  
  # Create a data frame with clusters
  cluster_df <- data.frame(
    sample_name = names(clusters),
    cluster_id = clusters,
    stringsAsFactors = FALSE
  )
  
  # Merge with user data
  merged_data <- merge(
    cluster_df,
    data,
    by.x = "sample_name",
    by.y = sample,
    all.x = TRUE
  )
  
  # Check for samples not in data
  missing_in_data <- merged_data$sample_name[is.na(merged_data$cluster_id)]
  if (length(missing_in_data) > 0) {
    warning("Some samples from the tree are not in data: ",
            paste(head(missing_in_data, 3), collapse = ", "),
            if (length(missing_in_data) > 3) " ...",
            call. = FALSE)
  }
  
  # ---- Select best sample per cluster ----
  if (!is.null(quality)) {
    # Quality-based selection: highest quality per cluster
    retained_df <- merged_data %>%
      dplyr::group_by(cluster_id) %>%
      dplyr::arrange(dplyr::desc(.data[[quality]]), sample_name) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(sample_name)
  } else {
    # No quality: just take first sample (alphabetically) per cluster
    retained_df <- merged_data %>%
      dplyr::group_by(cluster_id) %>%
      dplyr::arrange(sample_name) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(sample_name)
  }
  
  # ---- Add suffix and return ----
  retained_samples <- paste0(retained_df$sample_name, suffix)
  
  return(retained_samples)
}
