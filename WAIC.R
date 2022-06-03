# Here is the code to estimate elpd via WAIC in section 5

# WAIC
waic <- numeric(length(eta.cand))
lppd <- numeric(length(eta.cand))
elpd_waic <- numeric(length(eta.cand))
for(i in 1:length(eta.cand)){
  up <- covid_sampling(N1=100000, M, U, N, u, n, eta=eta.cand[i], alpha, beta, alpha_q, beta_q)
  pi <- up$pi[10*(1:10000)]
  q <- up$q[10*(1:10000)]
  waic[i] <- var(dbinom(n, ceiling(pi*M), q, log = TRUE)+dbinom(N-n, M-ceiling(pi*M), (N-n)/M, log = TRUE))
  lppd[i] <- log(mean(dbinom(n, ceiling(pi*M), q)*dbinom(N-n, M-ceiling(pi*M), (N-n)/M)))
  elpd_waic[i] <- lppd[i]-waic[i]
}

# Fig.6
plot(eta.cand, elpd_waic, type="l", xlab="eta", ylab="elpd")
