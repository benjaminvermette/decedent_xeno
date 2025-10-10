# Example totals
N_pre <- 19783   # total templates at PreTx
N_post <- 60142  # total CD8 templates at POD33 or 49

# Hyperexpanded counts
hyper_pre <- round(0.00 * N_pre) # % should be cumulative frequency of hyperexpanded clones
hyper_post <- round(0.2793686 * N_post) # % should be cumulative frequency of hyperexpanded clones

# Not hyperexpanded counts
nonhyper_pre <- N_pre - hyper_pre
nonhyper_post <- N_post - hyper_post

# Contingency table
mat <- matrix(c(hyper_pre, nonhyper_pre,
                hyper_post, nonhyper_post),
              nrow = 2, byrow = TRUE,
              dimnames = list(Time = c("PreTx","POD33"),
                              Status = c("Hyper","NotHyper")))

print(mat)

# Fisher’s exact test
fisher.test(mat)
