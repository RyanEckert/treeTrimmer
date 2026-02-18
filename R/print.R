#' Print Method for tree_trimmer Objects
#'
#' @param x A tree_trimmer object returned by \code{\link{tree_trimmer}}
#' @param ... Additional arguments (not used)
#'
#' @return Invisibly returns the input object
#' @export
print.tree_trimmer <- function(x, ...) {
  cat("Tree Trimmer Result\n")
  cat("===================\n\n")
  cat("Final cut height:", round(x$height, 6), "\n")
  cat("  (Max raw height:", round(x$max_raw_height, 6), "× buffer:", x$buffer, ")\n\n")

  cat("Groups processed:", length(x$groups), "\n")
  cat("  Successfully unified:", sum(!is.na(x$group_heights)), "\n")
  if (any(is.na(x$group_heights))) {
    cat("  Failed to unify:", sum(is.na(x$group_heights)), "\n")
  }

  cat("\nGroup details:\n")
  for (i in seq_along(x$groups)) {
    name <- names(x$groups)[i]
    if (is.null(name) || name == "") name <- paste0("Group_", i)

    cat("  ", name, ": n=", x$group_sizes[i], sep = "")
    if (!is.na(x$group_heights[i])) {
      cat(", height=", round(x$group_heights[i], 6), sep = "")
    } else {
      cat(", FAILED")
    }
    cat("\n")
  }

  invisible(x)
}
