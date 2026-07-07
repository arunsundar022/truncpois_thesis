devtools::install_github("arunsundar022/truncpois")
library(truncpois)

dir.create("Figures", showWarnings = FALSE)
dir.create("RData", showWarnings = FALSE)

png("Figures/01_pmf_visualization.png", width=900, height=700)
plottruncpois(lambda = 5, b = 10)
dev.off()

# How mean and variance converge as bounds increase
lambda <- 5
bounds <- seq(5, 50, by = 1)
means <- sapply(bounds, function(b) extruncpois(lambda = lambda, a = 0, b = b))
variances <- sapply(bounds, function(b) vartruncpois(lambda = lambda, a = 0, b = b))

png("Figures/02_moment_convergence.png", width=900, height=700)
par(mfrow = c(1, 2))

plot(bounds, means, type = "l", lwd = 2.5, col = "darkblue", 
     ylab = "Expected Value", xlab = "Upper Bound (b)",
     main = "Mean Convergence for Truncated Poisson(λ=5)")
abline(h = lambda, lty = 2, col = "red", lwd = 2)
legend("bottomright", c("Truncated Mean", "Untruncated Mean"), 
       col = c("darkblue", "red"), lty = c(1, 2), lwd = c(2.5, 2))

plot(bounds, variances, type = "l", lwd = 2.5, col = "darkgreen",
     ylab = "Variance", xlab = "Upper Bound (b)",
     main = "Variance Convergence for Truncated Poisson(λ=5)")
abline(h = lambda, lty = 2, col = "red", lwd = 2)
legend("bottomright", c("Truncated Variance", "Untruncated Variance"), 
       col = c("darkgreen", "red"), lty = c(1, 2), lwd = c(2.5, 2))

par(mfrow = c(1, 1))
dev.off()

# Compare sampling methods via microbenchmark::microbenchmark()
set.seed(42)
n_samples <- 10000

scenarios <- list(
  list(lambda = 5, a = 0, b = 10, name = "Small\n(λ=5, [0,10])"),
  list(lambda = 50, a = 40, b = 60, name = "Bounded\n(λ=50, [40,60])"),
  list(lambda = 100, a = 90, b = 110, name = "Heavy\n(λ=100, [90,110])")
)

methods <- c("direct", "inversion", "bounded")
timing_matrix <- matrix(NA, nrow = 3, ncol = 3)
colnames(timing_matrix) <- sapply(scenarios, function(s) s$name)
rownames(timing_matrix) <- methods

for (i in seq_along(scenarios)) {
  scenario <- scenarios[[i]]
  mb <- microbenchmark::microbenchmark(
    direct    = rtruncpois(n_samples, lambda = scenario$lambda, a = scenario$a, b = scenario$b, method = "direct"),
    inversion = rtruncpois(n_samples, lambda = scenario$lambda, a = scenario$a, b = scenario$b, method = "inversion"),
    bounded   = rtruncpois(n_samples, lambda = scenario$lambda, a = scenario$a, b = scenario$b, method = "bounded"),
    times = 20
  )
  med_ms <- tapply(mb$time, mb$expr, median) / 1e6  # nanoseconds -> milliseconds
  timing_matrix[, i] <- med_ms[methods]
}

png("Figures/03_sampling_efficiency.png", width=900, height=700)
barplot(timing_matrix, beside = TRUE, col = c("lightblue", "lightcoral", "lightgreen"),
        main = "Sampling Method Efficiency (10,000 samples, median of 20 microbenchmark reps)",
        ylab = "Time (milliseconds)", xlab = "Parameter Scenario",
        legend.text = methods, args.legend = list(x = "topleft"))
dev.off()

# Numerical stability with large lambdas, extended with stress-test extremes
lambdas <- c(10, 50, 100, 200, 500, 1e5)
a_vals  <- c(lambdas[1:5] - 10, 1e5 - 500)
b_vals  <- c(lambdas[1:5] + 10, 1e5 + 500)

results <- data.frame(lambda = lambdas, log_pmf = NA)

for (i in seq_along(lambdas)) {
  lambda_i <- lambdas[i]
  a_i <- a_vals[i]
  b_i <- b_vals[i]
  log_pmf <- dtruncpois(floor(lambda_i), lambda = lambda_i, a = a_i, b = b_i, log = TRUE)
  results$log_pmf[i] <- log_pmf
}

# Narrowest possible truncation window (b = a + 1)
narrow_lambda <- 5
narrow_a <- 3
narrow_b <- narrow_a + 1
narrow_log_pmf <- dtruncpois(narrow_a, lambda = narrow_lambda, a = narrow_a, b = narrow_b, log = TRUE)

png("Figures/04_numerical_stability.png", width=900, height=700)
plot(results$lambda, results$log_pmf, type = "b", pch = 19, col = "darkblue", lwd = 2, log = "x",
     main = "Log-Scale Stability for Extreme Parameters",
     xlab = "Lambda (λ, log scale)", ylab = "Log-PMF at floor(λ)",
     xlim = c(3, 1.5e5), ylim = c(-15, 0))
grid(nx = NA, ny = NULL, col = "gray", lty = "dotted")
points(narrow_lambda, narrow_log_pmf, pch = 17, col = "firebrick", cex = 1.4)
text(narrow_lambda, narrow_log_pmf, labels = "narrowest window (b = a+1)",
     pos = 4, offset = 0.7, col = "firebrick", cex = 0.8)
legend("bottomleft", c("Narrow window (±10 to ±500)", "Narrowest possible window"),
       col = c("darkblue", "firebrick"), pch = c(19, 17), lty = c(1, NA))
dev.off()

# Zero-truncated Poisson example
set.seed(123)
lambda <- 2.5
samples <- rtruncpois(5000, lambda = lambda, a = 1)

theo_mean <- extruncpois(lambda = lambda, a = 1)
theo_var <- vartruncpois(lambda = lambda, a = 1)
sample_mean <- mean(samples)
sample_var <- var(samples)

png("Figures/05_zero_truncation_application.png", width=900, height=700)
par(mfrow = c(2, 2))

hist(samples, breaks = 0:max(samples) + 0.5, freq = FALSE, 
     main = "Zero-Truncated Poisson Distribution",
     xlab = "Count", ylab = "Probability", col = "lightblue", border = "black")
x_vals <- 1:max(samples)
pmf_vals <- dtruncpois(x_vals, lambda = lambda, a = 1)
points(x_vals, pmf_vals, type = "o", col = "red", lwd = 2, pch = 16)
legend("topright", c("Sample", "Theoretical"), col = c("lightblue", "red"), pch = c(15, 16))

theoretical_quantiles <- qtruncpois((1:length(samples))/(length(samples)+1), lambda = lambda, a = 1)
sample_quantiles <- sort(samples)
plot(theoretical_quantiles, sample_quantiles,
     main = "Q-Q Plot", xlab = "Theoretical", ylab = "Sample", pch = 16, col = "darkblue")
abline(0, 1, col = "red", lwd = 2, lty = 2)

moments_data <- data.frame(Stat = c("Mean", "Variance"),
                           Theoretical = c(theo_mean, theo_var),
                           Sample = c(sample_mean, sample_var))
barplot(t(moments_data[, -1]), beside = TRUE, col = c("steelblue", "coral"),
        main = "Theoretical vs Sample Moments", ylab = "Value",
        names.arg = moments_data$Stat, legend.text = c("Theoretical", "Sample"),
        args.legend = list(x = "topleft"))

bootstrap_means <- replicate(1000, mean(sample(samples, replace = TRUE)))
hist(bootstrap_means, breaks = 30, main = "Bootstrap Distribution of Mean",
     xlab = "Sample Mean", col = "lightgreen", border = "black")
abline(v = theo_mean, col = "red", lwd = 2, lty = 2)

par(mfrow = c(1, 1))
dev.off()

# MLE parameter recovery
set.seed(2026)
lambda_true <- 6
a6 <- 2
b6 <- 15
x6 <- rtruncpois(5000, lambda = lambda_true, a = a6, b = b6)
fit6 <- mletruncpois(x6, a = a6, b = b6)

support6   <- a6:b6
tab6       <- table(factor(x6, levels = support6))
emp_prop6  <- as.numeric(tab6) / length(x6)
theo_prop6 <- dtruncpois(support6, lambda = fit6$lambda, a = a6, b = b6)

png("Figures/06_mle_parameter_recovery.png", width=900, height=700)
barplot(rbind(Empirical = emp_prop6, Theoretical = theo_prop6),
        beside = TRUE, names.arg = support6,
        col = c("grey70", "steelblue"),
        ylab = "Probability", xlab = "x",
        main = paste0("MLE Recovery: True λ=", lambda_true,
                       ", Fitted λ=", round(fit6$lambda, 3)),
        legend.text = c("Empirical", "Theoretical (fitted λ)"),
        args.legend = list(x = "topright"))
dev.off()

cat("Fitted lambda:", fit6$lambda, " | True lambda:", lambda_true, "\n")

# dtruncpois(): PMF sums to exactly 1 over the support
upper <- qtruncpois(1 - 1e-12, a = 1, lambda = 6)
pmf_sum <- sum(dtruncpois(seq(1, upper, by = 1), a = 1, lambda = 6))
cat("dtruncpois() validity check: PMF sums to", pmf_sum, "over the support\n")

# ptruncpois(): lower.tail and log.p combinations
p_lower     <- ptruncpois(7, lambda = 5, a = 2, b = 10)
p_upper     <- ptruncpois(7, lambda = 5, a = 2, b = 10, lower.tail = FALSE)
p_upper_log <- ptruncpois(7, lambda = 5, a = 2, b = 10, lower.tail = FALSE, log.p = TRUE)
p_lower_log <- ptruncpois(10, lambda = 5, a = 2, b = 10, log.p = TRUE)
cat("ptruncpois() tail/log demonstrations:\n")
print(c(lower_tail = p_lower, upper_tail = p_upper,
        upper_tail_log = p_upper_log, lower_tail_log = p_lower_log))

# qtruncpois(): round-trip check against ptruncpois(), plus log.p and upper-tail variants
p_seq        <- seq(0.1, 0.9, by = 0.2)
q_vals       <- qtruncpois(p_seq, lambda = 10, a = 5, b = 15)
cdf_at_q     <- ptruncpois(q_vals, lambda = 10, a = 5, b = 15)
q_from_logp  <- qtruncpois(log(p_seq), lambda = 10, a = 5, b = 15, log.p = TRUE)
q_upper_tail <- qtruncpois(p_seq, lambda = 10, a = 5, b = 15, lower.tail = FALSE)
cat("qtruncpois() round-trip check (cdf_at_q should be >= target p):\n")
print(rbind(target_p = p_seq, quantile = q_vals, cdf_at_q = round(cdf_at_q, 4)))
cat("qtruncpois() log.p input matches linear-scale input:", all(q_from_logp == q_vals), "\n")
cat("qtruncpois() upper-tail quantiles:", q_upper_tail, "\n")

# modtruncpois(): non-integer lambda (unique mode) vs integer lambda (tied modes, warns)
cat("modtruncpois() non-integer lambda (unique mode):\n")
print(modtruncpois(lambda = 2.5, a = 0, b = 10))
cat("modtruncpois() integer lambda (tied modes, warning expected):\n")
print(suppressWarnings(modtruncpois(lambda = 5, a = 0, b = 10)))

# plottruncpois(): cdf and quantile plot types (pmf type already shown above)
png("Figures/plottruncpois_cdf_quantile.png", width=1400, height=700)
par(mfrow = c(1, 2))
plottruncpois(lambda = 5, a = 2, b = 10, type = "cdf")
plottruncpois(lambda = 5, a = 2, b = 10, type = "quantile")
par(mfrow = c(1, 1))
dev.off()

moment_convergence <- data.frame(bound = bounds, mean = means, variance = variances)
sampling_timing_ms <- timing_matrix
mle_recovery <- list(true_lambda = lambda_true, fitted = fit6,
                      sample = x6, a = a6, b = b6)

save(moment_convergence, sampling_timing_ms, mle_recovery,
     file = "RData/results.RData")

cat("All experiments done. Figures saved to Figures/ folder, results saved to RData/results.RData.\n")