# Inputs
N1 <- 30       # XDRTCC aa clones
r1 <- 1.1     # nt/aa ratio for XDRTCC
N2 <- 3635    # NonXDRTCC aa clones
r2 <- 1.01     # nt/aa ratio for NonXDRTCC

# Compute converged counts
conv1 <- round((r1 - 1) * N1)
conv2 <- round((r2 - 1) * N2)

notconv1 <- N1 - conv1
notconv2 <- N2 - conv2

mat <- matrix(c(conv1, notconv1,
                conv2, notconv2),
              nrow = 2, byrow = TRUE,
              dimnames = list(
                Group = c("XDRTCC","NonXDRTCC"),
                Status = c("Converged","NotConverged")
              ))

print(mat)
fisher.test(mat)
