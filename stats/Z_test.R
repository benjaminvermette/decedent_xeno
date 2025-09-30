n1 <- 1743 # Total Number of XDRTCCs
p1 <- 0.1927711 # Value of Nominator

n2 <- 25604 # Total Number of Non-XDRTCCs
p2 <- 0.3371739 # Value of Denominator

# Pooled proportion
p_pooled <- (n1 * p1 + n2 * p2) / (n1 + n2)

# Calculate Z-statistic
z_stat <- (p1 - p2) / sqrt(p_pooled * (1 - p_pooled) * (1/n1 + 1/n2))


# Two-tailed test, find critical Z-value for alpha = 0.05
critical_value <- qnorm(1 - 0.05/2)

# Compare Z-statistic to critical value
if (abs(z_stat) > critical_value) {
  cat("Reject the null hypothesis\n")
} else {
  cat("Fail to reject the null hypothesis\n")
}

# Alternatively, calculate the p-value
p_value <- 2 * (1 - pnorm(abs(z_stat)))

# Compare p-value to alpha
if (p_value < 0.05) {
  cat("Reject the null hypothesis\n")
} else {
  cat("Fail to reject the null hypothesis\n")
}
