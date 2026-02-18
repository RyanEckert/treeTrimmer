test_that("tree_trimmer works with basic example", {
  # Create simple test data
  set.seed(123)
  dist_mat <- dist(matrix(rnorm(50), nrow = 10))
  hc <- hclust(dist_mat, method = "average")
  
  # Create test groups
  tech_reps <- list(
    group1 = c("1", "2"),
    group2 = c("3", "4")
  )
  
  # Run function
  result <- tree_trimmer(hc, tech_reps, verbose = FALSE)
  
  # Check return structure
  expect_type(result, "list")
  expect_named(result, c("height", "max_raw_height", "buffer", 
                         "group_heights", "group_sizes", "groups"))
  expect_true(is.numeric(result$height))
  expect_true(result$height > 0)
  expect_equal(result$buffer, 1.1)
})

test_that("prefix method works", {
  set.seed(123)
  # Create labels with prefix pattern
  labels <- c("OGM024.1", "OGM024.2", "OGM066.1", "OGM066.2", 
              "ILL123.1", "ILL123.2")
  dist_mat <- dist(matrix(rnorm(length(labels)^2), nrow = length(labels)))
  hc <- hclust(dist_mat, method = "average")
  hc$labels <- labels
  
  prefixes <- c("OGM024", "OGM066")
  
  result <- tree_trimmer(hc, prefixes, 
                        method = "prefix", 
                        separator = ".",
                        verbose = FALSE)
  
  expect_equal(length(result$groups), 2)
  expect_equal(names(result$groups), prefixes)
})

test_that("buffer parameter works", {
  set.seed(123)
  dist_mat <- dist(matrix(rnorm(50), nrow = 10))
  hc <- hclust(dist_mat, method = "average")
  
  tech_reps <- list(group1 = c("1", "2"))
  
  result1 <- tree_trimmer(hc, tech_reps, buffer = 1.0, verbose = FALSE)
  result2 <- tree_trimmer(hc, tech_reps, buffer = 1.2, verbose = FALSE)
  
  expect_equal(result1$max_raw_height, result2$max_raw_height)
  expect_true(result2$height > result1$height)
  expect_equal(result2$height / result1$height, 1.2)
})

test_that("errors are thrown for invalid inputs", {
  set.seed(123)
  dist_mat <- dist(matrix(rnorm(50), nrow = 10))
  hc <- hclust(dist_mat, method = "average")
  
  # Test invalid tree object
  expect_error(tree_trimmer("not a tree", list(a = c("1", "2"))))
  
  # Test invalid buffer
  expect_error(tree_trimmer(hc, list(a = c("1", "2")), buffer = 0.5))
  
  # Test missing separator with prefix method
  expect_error(tree_trimmer(hc, c("OGM024"), method = "prefix"))
  
  # Test samples not in tree
  expect_error(tree_trimmer(hc, list(a = c("999", "888"))))
})
