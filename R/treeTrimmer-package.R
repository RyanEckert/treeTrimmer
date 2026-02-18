#' treeTrimmer: Calibrate Hierarchical Clustering Thresholds
#'
#' Tools to calibrate tree cut heights for clone detection using technical
#' replicates. Particularly useful for population genetic studies with clonal
#' organisms where identifying genetic clones is critical.
#'
#' @section Main Function:
#' \code{\link{tree_trimmer}}: Determines the minimum hierarchical clustering
#' cut height needed to group known technical replicates (clones) into single
#' clusters.
#'
#' @section Typical Workflow:
#' \enumerate{
#'   \item Calculate pairwise genetic distances (e.g., 1 - IBS)
#'   \item Build hierarchical clustering tree with \code{hclust()}
#'   \item Use \code{tree_trimmer()} with known technical replicates to
#'         calibrate the cut height
#'   \item Apply the calibrated height with \code{cutree()} for clone detection
#'   \item Remove clones by selecting one representative per cluster
#' }
#'
#' @keywords internal
"_PACKAGE"

