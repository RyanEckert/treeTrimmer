#' Calibrate Tree Cut Height Using Technical Replicates
#'
#' Determines the minimum hierarchical clustering cut height needed to group
#' technical replicates (known clones) into single clusters. This calibrated
#' threshold can then be used for clone detection across the full dataset.
#'
#' @param tree_obj An \code{hclust} object produced by \code{\link[stats]{hclust}}
#' @param rep_groups Replicate group specification, format depends on \code{method}:
#'   \describe{
#'     \item{\strong{method="list"}}{A character vector (single group) or named list
#'           of character vectors (multiple groups), where each element contains
#'           sample names belonging to that replicate group}
#'     \item{\strong{method="prefix"}}{A character vector of sample name prefixes
#'           used to auto-detect replicate groups from tree labels}
#'   }
#' @param method Character string specifying how to interpret \code{rep_groups}.
#'   Either \code{"list"} (default, explicit groups) or \code{"prefix"}
#'   (auto-detect from naming pattern).
#' @param separator Character string; required when \code{method="prefix"}.
#'   The delimiter between prefix and suffix in sample names (e.g., \code{"."}
#'   for samples like \code{"OGM024.1"}). Must be a single non-empty character.
#' @param buffer Numeric value >= 1.0 specifying the multiplicative safety margin
#'   to add above the empirical threshold (default = 1.1, i.e., 10\% buffer).
#'   The final cut height is calculated as \code{max_height * buffer}.
#'   Use \code{buffer = 1.0} for no margin, or increase for more conservative
#'   clone detection.
#' @param verbose Logical; if \code{TRUE} (default), prints diagnostic information
#'   showing sample names, unification heights, and summary statistics for each
#'   replicate group.
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{\code{height}}{Final cut height (= max_raw_height * buffer).
#'           Use this value with \code{\link[stats]{cutree}} for clone detection.}
#'     \item{\code{max_raw_height}}{Maximum height across all replicate groups
#'           before applying the buffer.}
#'     \item{\code{buffer}}{The buffer value used in the calculation.}
#'     \item{\code{group_heights}}{Named numeric vector of the minimum unification
#'           height for each replicate group.}
#'     \item{\code{group_sizes}}{Named integer vector of sample counts per group.}
#'     \item{\code{groups}}{List of character vectors showing which samples were
#'           assigned to each group (useful for verification).}
#'   }
#'
#' @details
#' \strong{Background}
#'
#' In population genetic studies of clonal organisms (corals, plants, asexual fungi),
#' identifying and removing genetic clones is essential to avoid pseudoreplication
#' in demographic and evolutionary analyses. Hierarchical clustering based on genetic
#' distance (e.g., 1 - IBS) can detect clones, but choosing an appropriate tree cut
#' height is non-trivial. This function solves that problem by using technical
#' replicates—samples that are known a priori to be clones (e.g., DNA extracted
#' multiple times from the same individual)—to empirically calibrate the threshold.
#'
#' \strong{Method: "list" (Explicit Groups)}
#'
#' Provide replicate groups as a named list where each element is a character
#' vector of sample names. This gives you full control over group membership:
#' \preformatted{
#' tech_reps <- list(
#'   colony_A = c("Sample_A_rep1", "Sample_A_rep2", "Sample_A_rep3"),
#'   colony_B = c("Sample_B_rep1", "Sample_B_rep2")
#' )
#' result <- tree_trimmer(hc, tech_reps, method = "list")
#' }
#'
#' \strong{Method: "prefix" (Auto-Detection)}
#'
#' If your samples follow a consistent naming convention (e.g., \code{PREFIX.NUMBER}),
#' the function can automatically detect replicate groups. For each prefix, it finds
#' all tree labels matching the pattern \code{prefix<separator>*}. The separator
#' must immediately follow the prefix to avoid false matches (e.g., \code{"OGM11"}
#' won't match \code{"OGM111"} if separator is \code{"."}):
#' \preformatted{
#' prefixes <- c("OGM024", "OGM066", "OGM108")
#' result <- tree_trimmer(hc, prefixes, method = "prefix", separator = ".")
#' }
#'
#' \strong{Buffer Selection}
#'
#' The buffer parameter adds a safety margin above the empirical threshold:
#' \itemize{
#'   \item \code{buffer = 1.0}: No margin (uses exact empirical threshold; risky
#'         if replicates have any technical variation)
#'   \item \code{buffer = 1.05}: Conservative 5\% margin
#'   \item \code{buffer = 1.1}: Default 10\% margin (recommended for most cases)
#'   \item \code{buffer = 1.2}: Liberal 20\% margin (may inadvertently merge
#'         distinct genotypes)
#' }
#'
#' \strong{Algorithm}
#'
#' For each replicate group:
#' \enumerate{
#'   \item Extract all unique heights from the dendrogram
#'   \item Test each height in ascending order using \code{\link[stats]{cutree}}
#'   \item Find the minimum height where all replicates in the group belong to
#'         the same cluster
#' }
#' The final threshold is the maximum height across all groups (ensuring ALL
#' groups unify at the returned threshold), multiplied by the buffer.
#'
#' @section Warning:
#' Groups with fewer than 2 samples are automatically skipped with a warning.
#' If any replicate group never unifies (highly unusual), that group returns
#' \code{NA} and a warning is issued, but the function continues processing
#' other groups.
#'
#' @examples
#' # Simulate genetic distance data
#' set.seed(123)
#' dist_mat <- dist(matrix(rnorm(100), nrow = 20))
#' hc <- hclust(dist_mat, method = "average")
#'
#' # Example 1: Explicit list method
#' tech_reps <- list(
#'   group1 = c("1", "2", "3"),
#'   group2 = c("4", "5", "6")
#' )
#' result <- tree_trimmer(hc, tech_reps)
#' print(result)
#'
#' # Use the calibrated height for clone detection
#' clone_groups <- cutree(hc, h = result$height)
#' table(clone_groups)
#'
#' # Example 2: Custom buffer
#' result2 <- tree_trimmer(hc, tech_reps, buffer = 1.05, verbose = FALSE)
#' result2$height
#'
#' # Example 3: Access individual components
#' result$group_heights
#' result$group_sizes
#'
#'
#' @seealso
#' \code{\link[stats]{hclust}} for hierarchical clustering,
#' \code{\link[stats]{cutree}} for cutting dendrograms,
#' \code{\link[stats]{dist}} for distance matrix computation
#'
#' @export
#' @importFrom stats cutree
tree_trimmer <- function(tree_obj, rep_groups,
                        method = c("list", "prefix"),
                        separator = NULL,
                        buffer = 1.1,
                        verbose = TRUE) {

  # Match method argument
  method <- match.arg(method)

  # ---- Input validation ----
  if (!inherits(tree_obj, "hclust")) {
    stop("tree_obj must be an hclust object from stats::hclust()",
         call. = FALSE)
  }

  if (!is.numeric(buffer) || length(buffer) != 1 || buffer < 1.0) {
    stop("buffer must be a single numeric value >= 1.0",
         call. = FALSE)
  }

  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("verbose must be TRUE or FALSE",
         call. = FALSE)
  }

  # ---- Method-specific processing ----
  if (method == "prefix") {
    # Prefix mode: auto-detect groups from tree labels

    if (is.null(separator)) {
      stop("separator argument is required when method='prefix'",
           call. = FALSE)
    }

    if (!is.character(separator) || length(separator) != 1 || nchar(separator) == 0) {
      stop("separator must be a single non-empty character string",
           call. = FALSE)
    }

    if (!is.character(rep_groups)) {
      stop("When method='prefix', rep_groups must be a character vector of prefixes",
           call. = FALSE)
    }

    if (verbose) {
      message("Processing ", length(rep_groups),
              " replicate group", if(length(rep_groups) > 1) "s" else "",
              " using prefix method (separator = '", separator, "')\n")
    }

    # Build groups by finding matches
    processed_groups <- list()
    for (prefix in rep_groups) {
      # Escape special regex characters in separator
      sep_escaped <- gsub("([.|()\\^{}+$*?\\[\\]])", "\\\\\\1", separator)
      # Build pattern: prefix must be followed by separator
      pattern <- paste0("^", prefix, sep_escaped)
      matches <- grep(pattern, tree_obj$labels, value = TRUE)

      if (length(matches) == 0) {
        warning("Prefix '", prefix, "' matched zero samples in tree labels",
                call. = FALSE)
      }

      processed_groups[[prefix]] <- matches
    }

  } else {
    # List mode: use groups as provided

    if (verbose) {
      n_groups <- if (is.list(rep_groups)) length(rep_groups) else 1
      message("Processing ", n_groups,
              " replicate group", if(n_groups > 1) "s" else "",
              " using list method\n")
    }

    # Convert single vector to list
    if (is.character(rep_groups)) {
      processed_groups <- list(group1 = rep_groups)
    } else if (is.list(rep_groups)) {
      processed_groups <- rep_groups
    } else {
      stop("When method='list', rep_groups must be a character vector or list of character vectors",
           call. = FALSE)
    }
  }

  # Remove empty groups
  processed_groups <- processed_groups[lengths(processed_groups) > 0]

  if (length(processed_groups) == 0) {
    stop("No valid replicate groups found",
         call. = FALSE)
  }

  # ---- Process each group ----
  group_heights <- numeric(length(processed_groups))
  group_sizes <- integer(length(processed_groups))
  names(group_heights) <- names(processed_groups)
  names(group_sizes) <- names(processed_groups)

  for (i in seq_along(processed_groups)) {
    group_name <- if (is.null(names(processed_groups)[i]) || names(processed_groups)[i] == "") {
      paste0("Group_", i)
    } else {
      names(processed_groups)[i]
    }

    rep_names <- processed_groups[[i]]
    group_sizes[i] <- length(rep_names)

    # Validate group
    if (length(rep_names) < 2) {
      if (verbose) {
        message("Group '", group_name, "' (n=", length(rep_names),
                " sample", if(length(rep_names) > 1) "s" else "",
                "): SKIPPED - need at least 2 replicates\n")
      }
      warning("Group '", group_name, "' has fewer than 2 replicates - skipping",
              call. = FALSE)
      group_heights[i] <- NA_real_
      next
    }

    if (!all(rep_names %in% tree_obj$labels)) {
      missing <- rep_names[!rep_names %in% tree_obj$labels]
      stop("Group '", group_name, "': the following samples were not found in tree labels: ",
           paste(missing, collapse = ", "),
           call. = FALSE)
    }

    # Print group members
    if (verbose) {
      message("Group '", group_name, "' (n=", length(rep_names), " samples): ",
              paste(rep_names, collapse = ", "))
    }

    # Find minimum unifying height
    heights <- sort(unique(tree_obj$height))
    found <- FALSE

    for (h in heights) {
      groups <- cutree(tree_obj, h = h)
      rep_group_ids <- groups[rep_names]

      if (length(unique(rep_group_ids)) == 1) {
        group_heights[i] <- h
        found <- TRUE

        if (verbose) {
          message("  Unified at height = ", round(h, 6), "\n")
        }
        break
      }
    }

    if (!found) {
      if (verbose) {
        message("  WARNING: Never unified\n")
      }
      warning("Group '", group_name, "' never unified across all tree heights. ",
              "This may indicate very high divergence or data quality issues.",
              call. = FALSE)
      group_heights[i] <- NA_real_
    }
  }

  # ---- Calculate final threshold ----
  if (all(is.na(group_heights))) {
    warning("No groups successfully unified - returning NA values",
            call. = FALSE)
    return(list(
      height = NA_real_,
      max_raw_height = NA_real_,
      buffer = buffer,
      group_heights = group_heights,
      group_sizes = group_sizes,
      groups = processed_groups
    ))
  }

  max_raw_height <- max(group_heights, na.rm = TRUE)
  final_height <- max_raw_height * buffer

  # ---- Print summary ----
  if (verbose) {
    message("=== SUMMARY ===")
    message("Max raw height: ", round(max_raw_height, 6))
    message("Buffer applied: ", buffer, " (", round((buffer - 1) * 100, 1), "%)")
    message("Final cut height: ", round(final_height, 6))

    if (any(is.na(group_heights))) {
      failed <- names(group_heights)[is.na(group_heights)]
      message("\nWARNING: ", length(failed), " group(s) failed to unify:")
      message("  ", paste(failed, collapse = ", "))
    }
  }

  # ---- Return structured object ----
  structure(
    list(
      height = final_height,
      max_raw_height = max_raw_height,
      buffer = buffer,
      group_heights = group_heights,
      group_sizes = group_sizes,
      groups = processed_groups
    ),
    class = c("tree_trimmer", "list")
  )
}
