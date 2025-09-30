# Define the observed counts
apples_time_A <- 68
apples_time_B <- 232
oranges_time_A <- 14819
oranges_time_B <- 1658

# Define the contingency table
contingency_table <- matrix(c(apples_time_A, apples_time_B, oranges_time_A, oranges_time_B), nrow = 2, byrow = TRUE)

# Perform Fisher's exact test
fisher_result <- fisher.test(contingency_table)

# Print the results
print(fisher_result)
