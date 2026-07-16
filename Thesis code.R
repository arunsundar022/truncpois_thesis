devtools::install_github("arunsundar022/truncpois")
library(truncpois)
library(LaplacesDemon)

dir.create("Figures", showWarnings = FALSE)

# Exp 1: PMF plot, no title, a/b shown in legend instead
png("Figures/01_pmf_visualization.png", width = 900, height = 700)
plottruncpois(lambda = 5, b = 10, main = "")
legend("topleft", legend = "a = 0, b = 10", bty = "n")
dev.off()

# Exp 2: moment convergence. No LaplacesDemon comparison - no closed-form moments there
lambda <- 5
bounds <- seq(5, 50, by = 1)
means <- sapply(bounds, function(b) extruncpois(lambda = lambda, a = 0, b = b))
variances <- sapply(bounds, function(b) vartruncpois(lambda = lambda, a = 0, b = b))

png("Figures/02_moment_convergence.png", width = 900, height = 700)
par(mfrow = c(1, 2))

plot(bounds, means,
  type = "l", lwd = 2.5, col = "darkblue",
  ylab = "Expected Value", xlab = "Upper Bound (b)", main = ""
)
abline(h = lambda, lty = 2, col = "red", lwd = 2)
legend("bottomright", c("Truncated Mean", "Untruncated Mean"),
  col = c("darkblue", "red"), lty = c(1, 2), lwd = c(2.5, 2)
)

plot(bounds, variances,
  type = "l", lwd = 2.5, col = "darkgreen",
  ylab = "Variance", xlab = "Upper Bound (b)", main = ""
)
abline(h = lambda, lty = 2, col = "red", lwd = 2)
legend("bottomright", c("Truncated Variance", "Untruncated Variance"),
  col = c("darkgreen", "red"), lty = c(1, 2), lwd = c(2.5, 2)
)

par(mfrow = c(1, 1))
dev.off()

# Exp 3: sampling efficiency, 5 diverse scenarios + LaplacesDemon (a corrected to a-1)
set.seed(42)
n_samples <- 10000

scenarios <- list(
  list(lambda = 3, a = 1, b = Inf, name = "Left-trunc (ZTP)\n(λ=3, a=1)"),
  list(lambda = 5, a = 0, b = 50, name = "Right-trunc, wide\n(λ=5, [0,50])"),
  list(lambda = 50, a = 40, b = 60, name = "Double, centered\n(λ=50, [40,60])"),
  list(lambda = 50, a = 48, b = 52, name = "Double, near edge\n(λ=50, [48,52])"),
  list(lambda = 20, a = 40, b = 60, name = "Double, λ outside\n(λ=20, [40,60])")
)

methods <- c("direct", "inversion", "bounded", "LaplacesDemon")
timing_matrix <- matrix(NA, nrow = length(methods), ncol = length(scenarios))
colnames(timing_matrix) <- sapply(scenarios, function(s) s$name)
rownames(timing_matrix) <- methods

for (i in seq_along(scenarios)) {
  scenario <- scenarios[[i]]
  mb <- microbenchmark::microbenchmark(
    direct        = rtruncpois(n_samples, lambda = scenario$lambda, a = scenario$a, b = scenario$b, method = "direct"),
    inversion     = rtruncpois(n_samples, lambda = scenario$lambda, a = scenario$a, b = scenario$b, method = "inversion"),
    bounded       = rtruncpois(n_samples, lambda = scenario$lambda, a = scenario$a, b = scenario$b, method = "bounded"),
    LaplacesDemon = rtrunc(n_samples, spec = "pois", a = scenario$a - 1, b = scenario$b, lambda = scenario$lambda),
    times = 20
  )
  med_ms <- tapply(mb$time, mb$expr, median) / 1e6 # nanoseconds -> milliseconds
  timing_matrix[, i] <- med_ms[c("direct", "inversion", "bounded", "LaplacesDemon")]
}

png("Figures/03_sampling_efficiency.png", width = 1100, height = 700)
barplot(timing_matrix,
  beside = TRUE, col = c("lightblue", "lightcoral", "lightgreen", "grey60"),
  main = "", ylab = "Time (milliseconds)", xlab = "Parameter Scenario",
  legend.text = methods, args.legend = list(x = "topleft", inset = 0.025)
)
dev.off()

sampling_timing_ms <- as.data.frame(timing_matrix)
print(round(sampling_timing_ms, 2))

# Exp 4: numerical stability, extreme lambda / narrow windows, vs LaplacesDemon
lambdas <- c(10, 50, 100, 200, 500, 1e5)
a_vals <- c(lambdas[1:5] - 10, 1e5 - 500)
b_vals <- c(lambdas[1:5] + 10, 1e5 + 500)

stability <- data.frame(
  lambda = lambdas, a = a_vals, b = b_vals,
  log_pmf_truncpois = NA, log_pmf_laplacesdemon_naive_a = NA
)

for (i in seq_along(lambdas)) {
  lambda_i <- lambdas[i]
  a_i <- a_vals[i]
  b_i <- b_vals[i]
  stability$log_pmf_truncpois[i] <- dtruncpois(floor(lambda_i), lambda = lambda_i, a = a_i, b = b_i, log = TRUE)
  # a passed naively, no -1 correction
  stability$log_pmf_laplacesdemon_naive_a[i] <- dtrunc(floor(lambda_i), spec = "pois", a = a_i, b = b_i, lambda = lambda_i, log = TRUE)
}

# Narrowest possible truncation window (b = a + 1)
narrow_lambda <- 5
narrow_a <- 3
narrow_b <- narrow_a + 1
narrow_log_pmf <- dtruncpois(narrow_a, lambda = narrow_lambda, a = narrow_a, b = narrow_b, log = TRUE)

png("Figures/04_numerical_stability.png", width = 900, height = 700)
plot(stability$lambda, stability$log_pmf_truncpois,
  type = "b", pch = 19, col = "darkblue", lwd = 2, log = "x",
  main = "",
  xlab = "Lambda (λ, log scale)", ylab = "Log-PMF at floor(λ)",
  xlim = c(3, 1.5e5), ylim = c(-15, 0)
)
grid(nx = NA, ny = NULL, col = "gray", lty = "dotted")
points(narrow_lambda, narrow_log_pmf, pch = 17, col = "firebrick", cex = 1.4)
text(narrow_lambda, narrow_log_pmf,
  labels = "narrowest window (b = a+1)",
  pos = 4, offset = 0.7, col = "firebrick", cex = 0.8
)
legend("bottomleft", c("Narrow window (±10 to ±500)", "Narrowest possible window"),
  col = c("darkblue", "firebrick"), pch = c(19, 17), lty = c(1, NA),
  inset = 0.025
)
dev.off()

print(stability)

# lower.tail = FALSE gives the same result as the default - ignored
tail_check <- data.frame(
  lower_tail_arg = c(TRUE, FALSE),
  ptruncpois = c(
    ptruncpois(7, lambda = 5, a = 2, b = 10, lower.tail = TRUE),
    ptruncpois(7, lambda = 5, a = 2, b = 10, lower.tail = FALSE)
  ),
  laplacesdemon_ptrunc = c(
    ptrunc(7, spec = "pois", a = 1, b = 10, lambda = 5, lower.tail = TRUE),
    ptrunc(7, spec = "pois", a = 1, b = 10, lambda = 5, lower.tail = FALSE)
  )
)
print(tail_check)

# same bug in qtrunc, plus no log.p support at all
qtail_check <- data.frame(
  lower_tail_arg = c(TRUE, FALSE),
  qtruncpois = c(
    qtruncpois(0.3, lambda = 5, a = 2, b = 10, lower.tail = TRUE),
    qtruncpois(0.3, lambda = 5, a = 2, b = 10, lower.tail = FALSE)
  ),
  laplacesdemon_qtrunc = c(
    qtrunc(0.3, spec = "pois", a = 1, b = 10, lambda = 5, lower.tail = TRUE),
    qtrunc(0.3, spec = "pois", a = 1, b = 10, lambda = 5, lower.tail = FALSE)
  )
)
print(qtail_check)
qtrunc_logp_error <- tryCatch(
  {
    qtrunc(log(0.3), spec = "pois", a = 1, b = 10, lambda = 5, log.p = TRUE)
    "no error"
  },
  error = function(e) conditionMessage(e)
)
print(c(qtrunc_log.p_TRUE_result = qtrunc_logp_error))

# Exp 5: zero-truncated simulation
set.seed(123)
lambda <- 2.5
samples <- rtruncpois(5000, lambda = lambda, a = 1)

theo_mean <- extruncpois(lambda = lambda, a = 1)
theo_var <- vartruncpois(lambda = lambda, a = 1)
sample_mean <- mean(samples)
sample_var <- var(samples)

png("Figures/05_zero_truncation_application.png", width = 900, height = 700)
par(mfrow = c(2, 2))

hist(samples,
  breaks = 0:max(samples) + 0.5, freq = FALSE,
  main = "", xlab = "Count", ylab = "Probability", col = "lightblue", border = "black"
)
x_vals <- 1:max(samples)
pmf_vals <- dtruncpois(x_vals, lambda = lambda, a = 1)
points(x_vals, pmf_vals, type = "o", col = "red", lwd = 2, pch = 16)
legend("topright", c("Sample", "Theoretical"), col = c("lightblue", "red"), pch = c(15, 16), inset = 0.025)

theoretical_quantiles <- qtruncpois((1:length(samples)) / (length(samples) + 1), lambda = lambda, a = 1)
sample_quantiles <- sort(samples)
plot(theoretical_quantiles, sample_quantiles,
  main = "", xlab = "Theoretical", ylab = "Sample", pch = 16, col = "darkblue"
)
abline(0, 1, col = "red", lwd = 2, lty = 2)

moments_data <- data.frame(
  Stat = c("Mean", "Variance"),
  Theoretical = c(theo_mean, theo_var),
  Sample = c(sample_mean, sample_var)
)
barplot(t(moments_data[, -1]),
  beside = TRUE, col = c("steelblue", "coral"),
  main = "", ylab = "Value",
  names.arg = moments_data$Stat, legend.text = c("Theoretical", "Sample"),
  args.legend = list(x = "topleft", inset = 0.025)
)

bootstrap_means <- replicate(1000, mean(sample(samples, replace = TRUE)))
hist(bootstrap_means,
  breaks = 30, main = "",
  xlab = "Sample Mean", col = "lightgreen", border = "black"
)
abline(v = theo_mean, col = "red", lwd = 2, lty = 2)

par(mfrow = c(1, 1))
dev.off()

# Exp 6: MLE recovery
set.seed(2026)
lambda_true <- 6
a6 <- 2
b6 <- 15
x6 <- rtruncpois(5000, lambda = lambda_true, a = a6, b = b6)
fit6 <- mletruncpois(x6, a = a6, b = b6)

support6 <- a6:b6
tab6 <- table(factor(x6, levels = support6))
emp_prop6 <- as.numeric(tab6) / length(x6)
theo_prop6 <- dtruncpois(support6, lambda = fit6$lambda, a = a6, b = b6)

png("Figures/06_mle_parameter_recovery.png", width = 900, height = 700)
barplot(rbind(Empirical = emp_prop6, Theoretical = theo_prop6),
  beside = TRUE, names.arg = support6,
  col = c("grey70", "steelblue"),
  ylab = "Probability", xlab = "x", main = "",
  legend.text = c("Empirical", "Theoretical (fitted λ)"),
  args.legend = list(x = "topright", inset = 0.025)
)
dev.off()

print(fit6)

# function checks

# dtruncpois: PMF sums to 1
upper <- qtruncpois(1 - 1e-12, a = 1, lambda = 6)
pmf_sum <- sum(dtruncpois(seq(1, upper, by = 1), a = 1, lambda = 6))
print(c(pmf_sum_over_support = pmf_sum))

# ptruncpois: lower.tail / log.p combos
ptruncpois_demo <- c(
  lower_tail     = ptruncpois(7, lambda = 5, a = 2, b = 10),
  upper_tail     = ptruncpois(7, lambda = 5, a = 2, b = 10, lower.tail = FALSE),
  upper_tail_log = ptruncpois(7, lambda = 5, a = 2, b = 10, lower.tail = FALSE, log.p = TRUE),
  lower_tail_log = ptruncpois(10, lambda = 5, a = 2, b = 10, log.p = TRUE)
)
print(ptruncpois_demo)

# qtruncpois: round-trip check against ptruncpois
p_seq <- seq(0.1, 0.9, by = 0.2)
q_vals <- qtruncpois(p_seq, lambda = 10, a = 5, b = 15)
cdf_at_q <- ptruncpois(q_vals, lambda = 10, a = 5, b = 15)
q_from_logp <- qtruncpois(log(p_seq), lambda = 10, a = 5, b = 15, log.p = TRUE)
q_upper_tail <- qtruncpois(p_seq, lambda = 10, a = 5, b = 15, lower.tail = FALSE)
qtruncpois_roundtrip <- rbind(
  target_p = p_seq, quantile = q_vals, cdf_at_q = round(cdf_at_q, 4),
  quantile_from_logp = q_from_logp, quantile_upper_tail = q_upper_tail
)
print(qtruncpois_roundtrip)
print(c(logp_matches_linear_input = all(q_from_logp == q_vals)))

# modtruncpois: non-integer vs integer lambda (tied modes)
print(modtruncpois(lambda = 2.5, a = 0, b = 10))
print(suppressWarnings(modtruncpois(lambda = 5, a = 0, b = 10)))

# plottruncpois: cdf and quantile plot types
png("Figures/plottruncpois_cdf_quantile.png", width = 1400, height = 700)
par(mfrow = c(1, 2))
plottruncpois(lambda = 5, a = 2, b = 10, type = "cdf")
plottruncpois(lambda = 5, a = 2, b = 10, type = "quantile")
par(mfrow = c(1, 1))
dev.off()

message("All experiments done. Figures saved to Figures/ folder.")
