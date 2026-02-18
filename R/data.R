#' Example Coral IBS Matrix
#'
#' A simulated identity-by-state (IBS) matrix from a coral population genetics
#' study. Contains 24 samples including 3 colonies sequenced in triplicate
#' (technical replicates) plus 15 unique genets.
#'
#' @format A 24x24 symmetric matrix with sample names as row/column names.
#'   Values represent genetic similarity (IBS), ranging from 0 to 1.
#'   \describe{
#'     \item{Technical replicates}{OGM024.1-3, OGM066.1-3, OGM108.1-3 have
#'           very high IBS (0.95-0.99) within each group}
#'     \item{Unique genets}{OGM001-005, OGM010-015, OGM020-023 have moderate
#'           IBS (0.65-0.85) to each other}
#'   }
#'
#' @details
#' This dataset simulates RADseq data from *Orbicella* corals sampled across
#' shallow and mesophotic depth zones. The technical replicates represent DNA
#' extracted and sequenced multiple times from the same colony, providing known
#' clones for calibrating the clone detection threshold.
#'
#' The IBS values were generated to mimic realistic patterns:
#' \itemize{
#'   \item Technical replicates: High IBS (0.95-0.99) with some technical variation
#'   \item Within-depth zone: Moderate-high IBS (0.75-0.85)
#'   \item Between-depth zones: Lower IBS (0.65-0.75)
#' }
#'
#' @source Simulated data based on empirical coral population genetics studies
#'
#' @name coral_ibs
#' @examples
#' # Load data from CSV (included in package)
#' data_file <- system.file("extdata", "coral_ibs.csv", package = "treeTrimmer")
#' coral_ibs <- read.csv(data_file, row.names = 1)
#' coral_ibs <- as.matrix(coral_ibs)
#'
#' # Convert to distance and build tree
#' dist_mat <- as.dist(1 - coral_ibs)
#' hc <- hclust(dist_mat, method = "average")
#'
#' # Calibrate using technical replicates
#' tech_reps <- list(
#'   OGM024 = c("OGM024.1", "OGM024.2", "OGM024.3"),
#'   OGM066 = c("OGM066.1", "OGM066.2", "OGM066.3"),
#'   OGM108 = c("OGM108.1", "OGM108.2", "OGM108.3")
#' )
#' result <- tree_trimmer(hc, tech_reps)
#' result$height
NULL

#' Example Coral Sample Metadata
#'
#' Metadata for the 24 samples in the \code{\link{coral_ibs}} dataset,
#' including depth zone, sequencing quality, and replicate information.
#'
#' @format A data frame with 24 rows and 5 columns:
#' \describe{
#'   \item{sample}{Character; sample identifier matching row names in coral_ibs}
#'   \item{depth_zone}{Character; either "Shallow" or "Mesophotic"}
#'   \item{quality_score}{Numeric; sequencing quality metric (0-1)}
#'   \item{is_tech_rep}{Logical; TRUE if sample is a technical replicate}
#'   \item{colony_id}{Character; biological colony identifier (replicates
#'         share the same colony_id)}
#' }
#'
#' @source Simulated data
#'
#' @name coral_metadata
#' @examples
#' # Load metadata from CSV
#' data_file <- system.file("extdata", "coral_metadata.csv", package = "treeTrimmer")
#' coral_metadata <- read.csv(data_file)
#'
#' # Identify technical replicates
#' tech_reps_samples <- coral_metadata$sample[coral_metadata$is_tech_rep]
#' print(tech_reps_samples)
#'
#' # After clone detection, select best quality rep per cluster
#' # (see vignette for full workflow)
NULL

#' Technical Replicate Groups
#'
#' A named list defining which samples are technical replicates of each other
#' in the \code{\link{coral_ibs}} dataset. Each element contains sample names
#' from one biological colony that was sequenced multiple times.
#'
#' @format A named list with 3 elements:
#' \describe{
#'   \item{OGM024}{Character vector of 3 replicate names}
#'   \item{OGM066}{Character vector of 3 replicate names}
#'   \item{OGM108}{Character vector of 3 replicate names}
#' }
#'
#' @source Simulated data
#'
#' @name tech_rep_groups
#' @examples
#' # Tech rep groups are provided as package data
#' tech_rep_groups <- list(
#'   OGM024 = c("OGM024.1", "OGM024.2", "OGM024.3"),
#'   OGM066 = c("OGM066.1", "OGM066.2", "OGM066.3"),
#'   OGM108 = c("OGM108.1", "OGM108.2", "OGM108.3")
#' )
#'
#' # Use directly with tree_trimmer
#' data_file <- system.file("extdata", "coral_ibs.csv", package = "treeTrimmer")
#' coral_ibs <- as.matrix(read.csv(data_file, row.names = 1))
#' dist_mat <- as.dist(1 - coral_ibs)
#' hc <- hclust(dist_mat, method = "average")
#'
#' result <- tree_trimmer(hc, tech_rep_groups)
NULL
