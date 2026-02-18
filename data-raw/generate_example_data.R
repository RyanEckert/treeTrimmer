## Code to generate example datasets for treeTrimmer package
## This simulates realistic genetic distance data from a coral population genetics study

library(MASS)  # For mvrnorm

set.seed(2024)

# ============================================================================
# Dataset 1: coral_ibs - Simulated IBS matrix from coral samples
# ============================================================================

# Sample structure:
# - 3 colonies sequenced 3 times each (technical replicates)
# - 15 additional unique genets
# Total: 24 samples

n_colonies <- 3
reps_per_colony <- 3
n_unique <- 15
n_total <- n_colonies * reps_per_colony + n_unique

# Sample names following typical coral study naming
tech_rep_names <- c(
  paste0("OGM024.", 1:3),  # Colony 1
  paste0("OGM066.", 1:3),  # Colony 2
  paste0("OGM108.", 1:3)   # Colony 3
)

unique_names <- paste0("OGM", sprintf("%03d", c(1:5, 10:15, 20:23)))

sample_names <- c(tech_rep_names, unique_names)

# Generate realistic IBS matrix
# Technical replicates should have very high IBS (0.95-0.99)
# Different genets should have lower IBS (0.60-0.90)

coral_ibs <- matrix(0, nrow = n_total, ncol = n_total)
rownames(coral_ibs) <- colnames(coral_ibs) <- sample_names

# Set diagonal to 1 (identical to self)
diag(coral_ibs) <- 1

# Technical replicates: very high IBS
for (i in 1:n_colonies) {
  idx <- ((i-1)*3 + 1):(i*3)
  coral_ibs[idx, idx] <- runif(9, 0.95, 0.99)
  diag(coral_ibs[idx, idx]) <- 1
}

# Different genets: lower IBS with population structure
# Add some population structure (shallow vs mesophotic)
shallow_idx <- c(1:9, 10:14)  # Tech reps + some unique
meso_idx <- c(15:24)

# Within shallow group
for (i in shallow_idx) {
  for (j in shallow_idx) {
    if (i < j) {
      coral_ibs[i, j] <- coral_ibs[j, i] <- runif(1, 0.70, 0.85)
    }
  }
}

# Within mesophotic group  
for (i in meso_idx) {
  for (j in meso_idx) {
    if (i < j) {
      coral_ibs[i, j] <- coral_ibs[j, i] <- runif(1, 0.70, 0.85)
    }
  }
}

# Between groups: lower IBS
for (i in shallow_idx) {
  for (j in meso_idx) {
    coral_ibs[i, j] <- coral_ibs[j, i] <- runif(1, 0.60, 0.75)
  }
}

# Restore tech rep high IBS (may have been overwritten)
for (i in 1:n_colonies) {
  idx <- ((i-1)*3 + 1):(i*3)
  for (ii in idx) {
    for (jj in idx) {
      if (ii != jj) {
        coral_ibs[ii, jj] <- runif(1, 0.95, 0.99)
      }
    }
  }
  diag(coral_ibs[idx, idx]) <- 1
}

# ============================================================================
# Dataset 2: coral_metadata - Sample metadata
# ============================================================================

coral_metadata <- data.frame(
  sample = sample_names,
  depth_zone = c(
    rep("Shallow", 9),      # All tech reps are shallow
    rep("Shallow", 5),      # First 5 unique
    rep("Mesophotic", 10)   # Last 10 unique
  ),
  quality_score = runif(n_total, 0.75, 1.0),  # Sequencing quality
  is_tech_rep = c(rep(TRUE, 9), rep(FALSE, 15)),
  colony_id = c(
    rep("OGM024", 3),
    rep("OGM066", 3),
    rep("OGM108", 3),
    unique_names
  ),
  stringsAsFactors = FALSE
)

# ============================================================================
# Dataset 3: tech_rep_groups - Known technical replicate groupings
# ============================================================================

tech_rep_groups <- list(
  OGM024 = c("OGM024.1", "OGM024.2", "OGM024.3"),
  OGM066 = c("OGM066.1", "OGM066.2", "OGM066.3"),
  OGM108 = c("OGM108.1", "OGM108.2", "OGM108.3")
)

# ============================================================================
# Save datasets
# ============================================================================

usethis::use_data(coral_ibs, overwrite = TRUE)
usethis::use_data(coral_metadata, overwrite = TRUE)
usethis::use_data(tech_rep_groups, overwrite = TRUE)

# ============================================================================
# Document datasets (documentation will go in R/data.R)
# ============================================================================

cat("Datasets created:
  - coral_ibs: 24x24 IBS matrix
  - coral_metadata: Sample information
  - tech_rep_groups: Known technical replicates
")
