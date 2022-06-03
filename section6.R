# Here is the code to produce Fig.7 in section 6.

# week number where U is not 0
valid_index <- (southeast_df$Nr>0)*(1:length(southeast_df$Nr))

# save the optimal eta value
opt_eta <- numeric(length(valid_index))

for(g in 1:length(valid_index)){
  
  data_num <- valid_index[g]
  M <- southeast_df$M[data_num]
  N <- southeast_df$Nt[data_num]
  n <- southeast_df$nt[data_num]
  U <- southeast_df$Nr[data_num]
  u <- southeast_df$nr[data_num]
  p.hat <- (N-n)/M
  
  L <- n/M; K <- 1000
  x <- L + delta^(-(1:K)); max_x_index <- 0
  for(k in 1:K){
    if(x[k] <= (M-N+n)/M) max_x_index <- max_x_index+1
  }
  x <- x[K:(K-max_x_index+1)]; K <- max_x_index
  
  inte_cut <- 0
  for(i in 1:(K-1)){inte_cut <- inte_cut + p_cut(x[i])*(x[i+1]-x[i])}
  cdf_cut <- cumsum((p_cut(x[1:(K-1)])/inte_cut)*(x[2:K]-x[1:(K-1)]))
  
  d_nobs_eta <- numeric(n.eta)
  for(i in 1:n.eta){
    eta <- eta.cand[i]
    inte <- 0; for(z in 1:(K-1)){inte <- inte + f(x[z])*(x[z+1]-x[z])}
    cdf_em <- cumsum((f(x[1:(K-1)])/inte)*(x[2:K]-x[1:(K-1)]))
    pi_smi_sample <- pi_inver(N=1000) 
    d_nobs_eta[i] <- mean(f_eta_pi_n(pi_smi_sample, inte, inte_cut))
  }
  
  pi_cut_sample <- pi_inver_cut(100) # sample pi from the cut posterior
  q_sample <- rbeta(100, alpha_q, beta_q) # sample q from prior
  n_sample <- numeric(100)
  for(i in 1:100){ # we use S = 1
    n_sample[i] <- n_mcmc(N1=100000, pi_cut_sample[i], q_sample[i])[100000]
  }
  
  eta <- 1; d_n <- numeric(100)
  for(i in 1:100){
    n <- n_sample[i]
    L <- n/M; K <- 1000
    x <- L + delta^(-(1:K)); max_x_index <- 0
    for(k in 1:K){
      if(x[k] <= (M-N+n)/M) max_x_index <- max_x_index+1
    }
    x <- x[K:(K-max_x_index+1)]; K <- max_x_index
    inte_full <- 0
    for(z in 1:(K-1)){inte_full <- inte_full + f(x[z])*(x[z+1]-x[z])}
    cdf_em <- cumsum((f(x[1:(K-1)])/inte_full)*(x[2:K]-x[1:(K-1)]))
    pi_full_sample <- pi_inver(1000)
    inte_cut_new_n <- 0
    for(z in 1:(K-1)){inte_cut_new_n <- inte_cut_new_n + p_cut(x[z])*(x[z+1]-x[z])}
    d_n[i] <- mean(f_eta_pi_n(pi_full_sample, inte_full, inte_cut_new_n))
  }
  
  pr <- numeric(n.eta)
  for(i in 1:n.eta){
    pr[i] <- mean(d_n > d_nobs_eta[i])
  }
  
  count <- 0
  for(i in 1:n.eta){
    if(pr[i] >= 0.05) count <- count+1
  }
  
  opt_eta[g] <- eta.cand[count]
}


# repeat but under q uniform prior
opt_eta_unifq <- numeric(length(valid_index))

f_unifq <- function(pi){ 
  exp(lbeta(n*eta+1, ceiling(pi*M)*eta-n*eta+1) +
        dhyper(u, ceiling(pi*M), M-ceiling(pi*M), U, log = TRUE) + 
        eta*dbinom(N-n, M-ceiling(pi*M), (N-n)/M, log = TRUE) + eta*lchoose(ceiling(pi*M),n) 
      + (alpha-1)*log(pi) + (beta-1)*log(1-pi))
}

f_eta_pi_n_unifq <- function(pi, inte, inte_cut){
  lbeta(n*eta+1, ceiling(pi*M)*eta-n*eta+1) +
    dhyper(u, ceiling(pi*M), M-ceiling(pi*M), U, log = TRUE) + 
    eta*dbinom(N-n, M-ceiling(pi*M), (N-n)/M, log = TRUE) + eta*lchoose(ceiling(pi*M),n) 
  + (alpha-1)*log(pi) + (beta-1)*log(1-pi) - dhyper(u, ceiling(pi*M), M-ceiling(pi*M), U, log = TRUE) 
  - dbeta(pi, alpha, beta, log=TRUE) - log(inte) + log(inte_cut)
}

for(g in 23:length(valid_index)){
  
  data_num <- valid_index[g]
  M <- southeast_df$M[data_num]
  N <- southeast_df$Nt[data_num]
  n <- southeast_df$nt[data_num]
  U <- southeast_df$Nr[data_num]
  u <- southeast_df$nr[data_num]
  p.hat <- (N-n)/M
  
  L <- n/M; K <- 1000
  x <- L + delta^(-(1:K)); max_x_index <- 0
  for(k in 1:K){
    if(x[k] <= (M-N+n)/M) max_x_index <- max_x_index+1
  }
  x <- x[K:(K-max_x_index+1)]; K <- max_x_index
  
  inte_cut <- 0
  for(i in 1:(K-1)){inte_cut <- inte_cut + p_cut(x[i])*(x[i+1]-x[i])}
  cdf_cut <- cumsum((p_cut(x[1:(K-1)])/inte_cut)*(x[2:K]-x[1:(K-1)]))
  
  d_nobs_eta <- numeric(n.eta)
  for(i in 1:n.eta){
    eta <- eta.cand[i]
    inte <- 0; for(z in 1:(K-1)){inte <- inte + f_unifq(x[z])*(x[z+1]-x[z])}
    cdf_em <- cumsum((f_unifq(x[1:(K-1)])/inte)*(x[2:K]-x[1:(K-1)]))
    pi_smi_sample <- pi_inver(N=1000) 
    d_nobs_eta[i] <- mean(f_eta_pi_n_unifq(pi_smi_sample, inte, inte_cut))
  }
  
  pi_cut_sample <- pi_inver_cut(100) # sample pi from the cut posterior
  q_sample <- runif(100) # sample q from prior
  n_sample <- numeric(100)
  for(i in 1:100){ # we use S = 1
    n_sample[i] <- n_mcmc(N1=100000, pi_cut_sample[i], q_sample[i])[100000]
  }
  
  eta <- 1; d_n <- numeric(100)
  for(i in 1:100){
    n <- n_sample[i]
    L <- n/M
    x <- L + delta^(-(1:K)); x <- x[K:1]
    inte_full <- 0
    for(z in 1:(K-1)){inte_full <- inte_full + f_unifq(x[z])*(x[z+1]-x[z])}
    cdf_em <- cumsum((f_unifq(x[1:(K-1)])/inte_full)*(x[2:K]-x[1:(K-1)]))
    pi_full_sample <- pi_inver(1000)
    inte_cut_new_n <- 0
    for(z in 1:(K-1)){inte_cut_new_n <- inte_cut_new_n + p_cut(x[z])*(x[z+1]-x[z])}
    d_n[i] <- mean(f_eta_pi_n_unifq(pi_full_sample, inte_full, inte_cut_new_n))
  }
  
  pr <- numeric(n.eta)
  for(i in 1:n.eta){
    pr[i] <- mean(d_n > d_nobs_eta[i])
  }
  
  count <- 0
  for(i in 1:n.eta){
    if(pr[i] >= 0.05) count <- count+1
  }
  
  opt_eta_unifq[g] <- eta.cand[count]
}

# u>0 among U>0 (u in the denominator)
valid_index1 <- (southeast_df$nr[valid_index]>0)*(1:length(valid_index))
estimated.q1 <- southeast_df$nt[valid_index1]*southeast_df$Nr[valid_index1]/
  (southeast_df$M[valid_index1]*southeast_df$nr[valid_index1])
estimated.q <- c(estimated.q1[1:19], 0, 0, estimated.q1[22:40]) # set estimated q=0 when u=0

## Fig.7(a)
plot(1:40, opt_eta, type="l", xlab="week", ylab="optimal eta")
lines(1:40, opt_eta_unifq, lty=2, col="grey")
lines(1:40, estimated.q, col="type")

## Fig.7(b)
plot(estimated.q, opt_eta, xlab="estimated q", ylab="optimal eta")