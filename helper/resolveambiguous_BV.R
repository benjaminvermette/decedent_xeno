# resolveambiguous function adapted from Obradovic et. al 2022 (doi: 10.1016/j.simpa.2021.100142)
# that keeps adequately promiscuous clones present > cutoff (1e-5 default) in both CD4s & CD8s

normalize <- function(data) {
  nc <- ncol(data)
  for (i in seq_len(nc)) {
    s <- sum(data[, i], na.rm = TRUE)
    if (s > 0) data[, i] <- data[, i] / s
  }
  data
}

resolveambiguous <- function(data, c1, c2, ratio = 5, cutoff = 1e-5, do_normalize = TRUE) {
  # Normalize the two columns IN PLACE (before resolving -- optional)
  if (do_normalize) {
    data[, c(c1, c2)] <- normalize(as.matrix(data[, c(c1, c2)]))
  }
  
  both_pos <- which(data[, c1] > 0 & data[, c2] > 0)
  
  # Dominance: zero the weaker column when one is >= ratio * the other
  c1_dom <- both_pos[data[both_pos, c1] > ratio * data[both_pos, c2]]
  c2_dom <- both_pos[data[both_pos, c2] > ratio * data[both_pos, c1]]
  if (length(c1_dom)) data[c1_dom, c2] <- 0
  if (length(c2_dom)) data[c2_dom, c1] <- 0
  
  # Ambiguous: both > 0 but neither dominates
  ambi <- setdiff(both_pos, c(c1_dom, c2_dom))
  if (length(ambi)) {
    # If both > cutoff → keep both; if either < cutoff → set both to 0
    low <- ambi[data[ambi, c1] < cutoff | data[ambi, c2] < cutoff]
    if (length(low)) data[low, c(c1, c2)] <- 0
  }
  
  # Final normalization of the two columns
  data[, c(c1, c2)] <- normalize(as.matrix(data[, c(c1, c2)]))
  
  data
}

