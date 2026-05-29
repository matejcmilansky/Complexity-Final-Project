library(bootnet)
library(qgraph)
library(psych)
library(networktools)
library(NetworkComparisonTest)

# Load raw depression scores per group
no_mistreat <- read.csv("no_trauma_final.csv")
high_mistreat <- read.csv("high_trauma_final.csv")

set.seed(129)

# Number of subsampling/bootstrap iterations
n_subsamples <- 2000

# Keep default tuning parameter (0.5) for network estimation
net_tuning <- 0.5

# Sample size of the trauma group (the smaller group)
n_trauma <- nrow(high_mistreat)  

# Number of nodes in network (No. of depression symptoms measured)
n_nodes <- ncol(no_mistreat)

# Calculate possible number of node connections for later density calculation
n_possible <- n_nodes * (n_nodes - 1) / 2

# Empty arrays for storing
densities_no <- numeric(n_subsamples)
strengths_no <- numeric(n_subsamples)
densities_trauma <- numeric(n_subsamples)
strengths_trauma <- numeric(n_subsamples)
diff_density <- numeric(n_subsamples)
diff_strength <- numeric(n_subsamples)


# Unequal group size prevents direct comparison; run a bootstrapping and
# subsampling procedure instead
cat("Running paired bootstrap (n =", n_trauma, "per group per iteration...)\n")
for (i in 1:n_subsamples) {
  
  # Randomly subsample no trauma group to trauma group size
  idx_no <- sample(nrow(no_mistreat), n_trauma, replace = FALSE)
  dat_no_i <- no_mistreat[idx_no, , drop = FALSE]
  
  # Bootstrap sample of high trauma group (with replacement)
  idx_trauma <- sample(nrow(high_mistreat), n_trauma, replace = TRUE)
  dat_trauma_i <- high_mistreat[idx_trauma, , drop = FALSE]
  
  # Estimate networks
  net_no_i <- estimateNetwork(dat_no_i, default = "EBICglasso", tuning=net_tuning)
  net_trauma_i <- estimateNetwork(dat_trauma_i, default = "EBICglasso", tuning=net_tuning)
  
  # Calculate densities
  dens_no <- sum(net_no_i$graph[lower.tri(net_no_i$graph)] != 0) / n_possible
  dens_trauma <- sum(net_trauma_i$graph[lower.tri(net_trauma_i$graph)] != 0) / n_possible
  
  # Calculate strengths
  str_no <- sum(abs(net_no_i$graph[lower.tri(net_no_i$graph)]))
  str_trauma <- sum(abs(net_trauma_i$graph[lower.tri(net_trauma_i$graph)]))
  
  # Save densities, strengths, and the difference in them between the groups
  densities_no[i] <- dens_no
  densities_trauma[i] <- dens_trauma
  strengths_no[i] <- str_no
  strengths_trauma[i] <- str_trauma
  diff_density[i] <- dens_trauma - dens_no
  diff_strength[i] <- str_trauma - str_no
  
  if (i %% 100 == 0) cat("Completed iteration", i, "\n")
}

# Print results
cat("\n=== BOOTSTRAP RESULTS (n_iter =", n_subsamples, ") ===\n")

ci_dens <- quantile(diff_density, c(0.025, 0.975))
cat(sprintf("\nDensity (Trauma - No Trauma):\n"))
cat(sprintf("  Trauma mean:    %.3f (SD = %.3f)\n", mean(densities_trauma), sd(densities_trauma)))
cat(sprintf("  No Trauma mean: %.3f (SD = %.3f)\n", mean(densities_no), sd(densities_no)))
cat(sprintf("  Mean diff:    %.4f\n", mean(diff_density)))
cat(sprintf("  95%% CI:       [%.4f, %.4f]\n", ci_dens[1], ci_dens[2]))


ci_str <- quantile(diff_strength, c(0.025, 0.975))
cat(sprintf("\nStrength (Trauma - No Trauma):\n"))
cat(sprintf("  Trauma mean:     %.3f (SD = %.3f)\n", mean(strengths_trauma), sd(strengths_trauma)))
cat(sprintf("  No Trauma mean: %.3f (SD = %.3f)\n", mean(strengths_no), sd(strengths_no)))
cat(sprintf("  Mean diff:     %.4f\n", mean(diff_strength)))
cat(sprintf("  95%% CI:        [%.4f, %.4f]\n", ci_str[1], ci_str[2]))

# Plot results
par(mfrow = c(1, 2))

hist(densities_no, col = rgb(0, 0, 1, 0.4), border = NA,
     main = "Density Distribution", xlab = "Density",
     xlim = range(c(densities_no, densities_trauma)), breaks = 20)
hist(densities_trauma, col = rgb(1, 0, 0, 0.4), border = NA, add = TRUE, breaks = 20)
legend("topright", legend = c("No Trauma", "High Trauma"),
       fill = c(rgb(0, 0, 1, 0.4), rgb(1, 0, 0, 0.4)), border = NA)
abline(v = mean(densities_no), col = "blue", lwd = 2, lty = 2)
abline(v = mean(densities_trauma), col = "red", lwd = 2, lty = 2)

hist(strengths_no, col = rgb(0, 0, 1, 0.4), border = NA,
     main = "Strength Distribution", xlab = "Global Strength",
     xlim = range(c(strengths_no, strengths_trauma)), breaks = 30)
hist(strengths_trauma, col = rgb(1, 0, 0, 0.4), border = NA, add = TRUE, breaks = 30)
legend("topright", legend = c("No Trauma", "High Trauma"),
       fill = c(rgb(0, 0, 1, 0.4), rgb(1, 0, 0, 0.4)), border = NA)
abline(v = mean(strengths_no), col = "blue", lwd = 2, lty = 2)
abline(v = mean(strengths_trauma), col = "red", lwd = 2, lty = 2)

par(mfrow = c(1, 1))

