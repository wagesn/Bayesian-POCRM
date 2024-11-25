library(stats)

bpocrm <- function(p0, p_skel, ttr, cohort_size, n_cohort, n_stop, start_comb, cs) {
  # Ensure p_skel is a matrix
  if (is.vector(p_skel)) {
    p_skel <- matrix(p_skel, nrow = 1)
  }
  
  nord_tox <- nrow(p_skel)
  mprior_tox <- rep(1 / nord_tox, nord_tox)  # Prior for each toxicity ordering
  asd <- 1.34
  
  bcrmh <- function(a, p, y, n) {
    lik <- exp(-0.5 * a^2 / asd)
    for (j in seq_along(p)) {
      pj <- p[j]^exp(a)
      lik <- lik * pj^y[j] * (1 - pj)^(n[j] - y[j])
    }
    return(lik)
  }
  
  bcrmht <- function(a, p, y, n) {
    lik <- a * exp(-0.5 * a^2 / asd)
    for (j in seq_along(p)) {
      pj <- p[j]^exp(a)
      lik <- lik * pj^y[j] * (1 - pj)^(n[j] - y[j])
    }
    return(lik)
  }
  
  bcrmht2 <- function(a, p, y, n) {
    lik <- a^2 * exp(-0.5 * a^2 / asd)
    for (j in seq_along(p)) {
      pj <- p[j]^exp(a)
      lik <- lik * pj^y[j] * (1 - pj)^(n[j] - y[j])
    }
    return(lik)
  }
  
  ncomb <- ncol(p_skel)
  y <- rep(0, ncomb)  # Toxicities at each dose level
  n <- rep(0, ncomb)  # Patients treated at each dose level
  comb_curr <- start_comb  # Current dose level
  comb_select <- rep(0, ncomb)  # Indicators for dose selection
  stop <- 0
  
  for (i in seq_len(n_cohort)) {
    # Generate data for a new cohort
    y[comb_curr] <- y[comb_curr] + rbinom(1, cohort_size, p0[comb_curr])
    n[comb_curr] <- n[comb_curr] + cohort_size
    
    if (any(n > n_stop)) {
      stop <- 1
      break
    }
    
    marginal_tox <- rep(0, nord_tox)
    est_tox <- rep(0, nord_tox)
    e2_tox <- rep(0, nord_tox)
    
    for (k in seq_len(nord_tox)) {
      marginal_tox[k] <- integrate(bcrmh, -Inf, Inf, p = p_skel[k, ], y = y, n = n)$value
      est_tox[k] <- integrate(bcrmht, -10, 10, p = p_skel[k, ], y = y, n = n)$value / marginal_tox[k]
      e2_tox[k] <- integrate(bcrmht2, -10, 10, p = p_skel[k, ], y = y, n = n)$value / marginal_tox[k]
    }
    
    postprob_tox <- (marginal_tox * mprior_tox) / sum(marginal_tox * mprior_tox)
    mtox_sel <- if (nord_tox > 1) which.max(postprob_tox) else 1
    ptox_hat <- p_skel[mtox_sel, ]^exp(est_tox[mtox_sel])
    post_var_tox <- e2_tox[mtox_sel] - est_tox[mtox_sel]^2
    crit_tox <- qnorm(0.5 + cs / 2)
    lb_tox <- est_tox[mtox_sel] - crit_tox * sqrt(post_var_tox)
    ub_tox <- est_tox[mtox_sel] + crit_tox * sqrt(post_var_tox)
    ptoxL <- p_skel[mtox_sel, ]^exp(ub_tox)
    ptoxU <- p_skel[mtox_sel, ]^exp(lb_tox)
    
    if (ptoxL[1] > ttr) {
      stop <- 1
      break
    }
    
    comb_curr <- which.min(abs(ptox_hat - ttr))
  }
  
  if (stop == 0) {
    comb_select[comb_curr] <- comb_select[comb_curr] + 1
  }
  
  return(list(comb_select = comb_select, tox_data = y, pt_allocation = n, stop = stop))
}

bpocrm_sim <- function(p0, p_skel, ttr, cohort_size, n_cohort, n_stop, n_trial, start_comb, cs) {
  ncomb <- length(p0)
  comb_select <- matrix(0, nrow = n_trial, ncol = ncomb)
  y <- matrix(0, nrow = n_trial, ncol = ncomb)
  n <- matrix(0, nrow = n_trial, ncol = ncomb)
  nstop <- 0
  
  for (i in seq_len(n_trial)) {
    result <- bpocrm(p0, p_skel, ttr, cohort_size, n_cohort, n_stop, start_comb, cs)
    comb_select[i, ] <- result$comb_select
    y[i, ] <- result$tox_data
    n[i, ] <- result$pt_allocation
    nstop <- nstop + result$stop
  }
  
  cat("True tox probability: ", round(p0, 3), "\n")
  cat("Selection percentage: ", formatC(colMeans(comb_select) * 100, digits = 1, format = "f"), "\n")
  cat("Number of toxicities: ", formatC(colMeans(y), digits = 1, format = "f"), "\n")
  cat("Number of patients treated: ", formatC(colMeans(n), digits = 1, format = "f"), "\n")
  cat("Percentage of stop: ", nstop / n_trial * 100, "\n")
}

# Example parameters
d <- 6
s <- 6

orders <- matrix(c(
  1, 2, 3, 4, 5, 6,
  1, 2, 3, 5, 4, 6,
  1, 2, 3, 5, 6, 4,
  1, 2, 5, 3, 4, 6,
  1, 2, 5, 3, 6, 4,
  1, 2, 5, 6, 3, 4
), nrow = s, byrow = TRUE)

skeleton <- c(0.07, 0.18, 0.33, 0.49, 0.63, 0.75)
p_skel <- t(apply(orders, 1, function(order) skeleton[order]))

p0 <- c(0.12, 0.21, 0.34, 0.50, 0.66, 0.79)
ttr <- 0.33
cohort_size <- 1
n_cohort <- 22
start_comb <- 3
n_stop <- 100
n_trial <- 100
cs <- 0.90

set.seed(580)
bpocrm_sim(p0, p_skel, ttr, cohort_size, n_cohort, n_stop, n_trial, start_comb, cs)
