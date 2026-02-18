# treeTrimmer

<!-- badges: start -->
[![R-CMD-check](https://github.com/ryanEckert/treeTrimmer/workflows/R-CMD-check/badge.svg)](https://github.com/ryanEckert/treeTrimmer/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**Calibrate Hierarchical Clustering Thresholds Using Technical Replicates**

## Overview

`treeTrimmer` provides a principled way to calibrate tree cut heights for clone detection in hierarchical clustering analyses. It uses known technical replicates (samples sequenced multiple times) to empirically determine the minimum genetic distance threshold at which true clones cluster together.

This is particularly useful for population genetic studies with clonal organisms (corals, plants, asexual fungi) where identifying and removing clonal replicates is critical for accurate demographic and evolutionary inference.

## Installation

Install from GitHub using `devtools`:

```r
# install.packages("devtools")
devtools::install_github("ryanEckert/treeTrimmer")
```

## Quick Start

```r
library(treeTrimmer)

# Build hierarchical clustering tree from genetic distances
hc <- hclust(distance_matrix, method = "average")

# Define your technical replicates (samples you KNOW are clones)
tech_reps <- list(
  colony_A = c("OGM024.1", "OGM024.2", "OGM024.3"),
  colony_B = c("OGM066.1", "OGM066.2", "OGM066.3"),
  colony_C = c("OGM108.1", "OGM108.2", "OGM108.3")
)

# Calibrate the cut height
result <- tree_trimmer(hc, tech_reps)

# Use the calibrated height for clone detection
clone_groups <- cutree(hc, h = result$height)
```

## Key Features

- **Empirical calibration** - Uses your data to determine the threshold, not arbitrary cutoffs
- **Prefix auto-detection** - Automatically finds replicate groups from sample naming patterns
- **Customizable buffer** - Add safety margins above the empirical threshold
- **Diagnostic output** - Verbose mode shows exactly what's happening
- **Structured return** - Access individual components (heights, group info, etc.)

## Example: Prefix Method

If your samples follow a consistent naming pattern:

```r
# Instead of manually listing all replicates...
prefixes <- c("OGM024", "OGM066", "OGM108")

result <- tree_trimmer(hc, prefixes,
                       method = "prefix",
                       separator = ".")

# This auto-detects:
# - OGM024.1, OGM024.2, OGM024.3
# - OGM066.1, OGM066.2, OGM066.3
# - OGM108.1, OGM108.2, OGM108.3
```

## Output

```
Processing 3 replicate groups using prefix method (separator = '.')

Group 'OGM024' (n=3 samples): OGM024.1, OGM024.2, OGM024.3
  Unified at height = 0.012345

Group 'OGM066' (n=3 samples): OGM066.1, OGM066.2, OGM066.3
  Unified at height = 0.015678

Group 'OGM108' (n=3 samples): OGM108.1, OGM108.2, OGM108.3
  Unified at height = 0.014321

=== SUMMARY ===
Max raw height: 0.015678
Buffer applied: 1.1 (10.0%)
Final cut height: 0.017246
```

## Complete Workflow

```r
# 1. Calculate genetic distances (e.g., from IBS matrix)
ibs_matrix <- read.table("samples.ibsMat")
dist_matrix <- as.dist(1 - ibs_matrix)

# 2. Build tree
hc <- hclust(dist_matrix, method = "average")

# 3. Calibrate threshold
result <- tree_trimmer(hc, tech_reps)

# 4. Detect clones
clone_groups <- cutree(hc, h = result$height)

# 5. Remove clones (keep one per group)
sample_metadata <- data.frame(
  sample = rownames(ibs_matrix),
  quality = quality_scores,
  cluster = clone_groups
)

unique_genets <- sample_metadata %>%
  group_by(cluster) %>%
  arrange(desc(quality)) %>%
  slice(1)

cat("Started with:", nrow(sample_metadata), "samples\n")
cat("Final dataset:", nrow(unique_genets), "unique genets\n")
```

## Integration with ANGSD

This package was designed for ANGSD (Analysis of Next Generation Sequencing Data) workflows:

```bash
# Generate IBS matrix
angsd -b samples.bams -GL 1 \
      -doIBS 1 -doCounts 1 -doMajorMinor 1 -makeMatrix 1 \
      -out myproject
```

```r
# Detect and remove clones in R
library(treeTrimmer)
ibs_mat <- as.matrix(read.table("myproject.ibsMat"))
dist_mat <- as.dist(1 - ibs_mat)
hc <- hclust(dist_mat, method = "average")

result <- tree_trimmer(hc, tech_reps)
clone_groups <- cutree(hc, h = result$height)

# Write clone-free BAM list
unique_samples <- # ... (see full workflow above)
writeLines(unique_samples$bam, "samples_noclones.bams")
```

```bash
# Continue with downstream ANGSD analyses
angsd -b samples_noclones.bams -GL 1 \
      -doMaf 1 -doGlf 2 # ...
```

## Documentation

- [Getting Started Vignette](vignettes/introduction.Rmd) - Basic usage and examples
- [Advanced Usage](vignettes/advanced.Rmd) - Troubleshooting and edge cases
- Function help: `?tree_trimmer`

## Citation

If you use `treeTrimmer` in your research, please cite:

```
Eckert, R. (2025). treeTrimmer: Calibrate Hierarchical Clustering
Thresholds Using Technical Replicates. R package version 1.0.0.
https://github.com/ryanEckert/treeTrimmer
```

## Issues and Contributing

- Report bugs or request features: [GitHub Issues](https://github.com/ryanEckert/treeTrimmer/issues)
- Contribute code: Submit a pull request

## License

MIT License - see [LICENSE](LICENSE) file for details

## Related Work

- [ANGSD](http://www.popgen.dk/angsd/) - Analysis of Next Generation Sequencing Data
- [pcangsd](https://github.com/Rosemeis/pcangsd) - PCA from genotype likelihoods
- [NGSadmix](https://github.com/ANGSD/NgsAdmix) - Admixture analysis from genotype likelihoods

## Authors

- R. J. Eckert (Florida Atlantic University)
