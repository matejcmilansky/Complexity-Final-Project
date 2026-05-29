library(bootnet)
library(qgraph)
library(psych)
library(networktools)
library(NetworkComparisonTest)

# Read data
no_mistreat <- read.csv("no_trauma_final.csv")
high_mistreat <- read.csv("high_trauma_final.csv")

# Draw a size-matched random sample from the no trauma group
n_high <- nrow(high_mistreat)
set.seed(45)
no_mistreat_matched <- no_mistreat[sample(nrow(no_mistreat), size = n_high, replace = FALSE), ]

# Verify sizes match
cat(sprintf("\nHigh mistreat n = %d\n", n_high))
cat(sprintf("No mistreat (matched) n = %d\n", nrow(no_mistreat_matched)))

# Estimate networks
cat("\nEstimating networks...\n")

net_no <- estimateNetwork(dat_no,   default = "EBICglasso", tuning=0.5)
net_high <- estimateNetwork(dat_high, default = "EBICglasso", tuning=0.5)

# Plot the networks; use the same layout and edge weight scale for comparability
layout_fixed <- averageLayout(net_no$graph, net_high$graph)
max_val <- max(abs(c(net_no$graph, net_high$graph)))

par(mfrow = c(1, 2))

qgraph(net_no$graph,
       layout     = layout_fixed,
       maximum    = max_val,
       title      = sprintf("No Trauma (n=%d)", nrow(dat_no)),
       labels     = vars,
       color      = "#AED6F1",
       vsize      = 6,
       esize      = 10,
       posCol     = "#2874A6",
       negCol     = "#C0392B")

qgraph(net_high$graph,
       layout     = layout_fixed,
       maximum    = max_val,
       title      = sprintf("High Trauma (n=%d)", nrow(dat_high)),
       labels     = vars,
       color      = "#F1948A",
       vsize      = 6,
       esize      = 10,
       posCol     = "#2874A6",
       negCol     = "#C0392B")

par(mfrow = c(1, 1))

# Run NCT
cat("\nRunning NCT (this may take several minutes)...\n")

# Set the No. of iterations to 2000 and test all edges
nct_result <- NCT(
  dat_no,
  dat_high,
  it            = 2000,     
  test.edges    = TRUE,       
  edges         = "all",
  progressbar   = TRUE
)

# Print NCT results
cat("\n=== NCT RESULTS ===\n")
print(summary(nct_result))

cat("\n--- Global Strength Invariance Test ---\n")
cat(sprintf("  Observed S:  %.4f\n", nct_result$glstrinv.real))
cat(sprintf("  p-value:     %.4f\n", nct_result$glstrinv.pval))

cat("\n--- Network Structure Invariance Test ---\n")
cat(sprintf("  Observed M:  %.4f\n", nct_result$nwinv.real))
cat(sprintf("  p-value:     %.4f\n", nct_result$nwinv.pval))

# Individual edges involve many comparisons and should be interpreted
# cautiously
cat("\n--- Individual Edge Tests (p < .05) ---\n")
if (!is.null(nct_result$einv.pvals)) {
  edge_pvals <- nct_result$einv.pvals
  sig_edges  <- edge_pvals[edge_pvals[, "p-value"] < 0.05, ]
  if (nrow(sig_edges) > 0) {
    print(sig_edges)
  } else {
    cat("  No individual edges significantly differ between groups.\n")
  }
}

