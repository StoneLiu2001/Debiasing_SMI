# synthetic data to test elpd code (Appendix C)

# synthetic data
M_s <- 9180135
pi_s <- rbeta(1, alpha, beta)
q_s <- rbeta(1, alpha_q, beta_q)
U_s <- southeast_df$Nr[c(1,4,5,6,9,10,11,13,14,15)] 
N_s <- southeast_df$Nt[c(1,4,5,6,9,10,11,13,14,15)]
u_s <- rhyper(10, ceiling(pi_s*M_s), M_s-ceiling(pi_s*M_s), U_s)
n_s <- n_mcmc(N1=200000, pi_s, q_s)[10000*(11:20)]
p.hat <- (N_s-n_s)/M_s

update_pow <- function(pi, q.aux, M, U, N, u, n, eta, alpha, beta, alpha_q, beta_q){
  
  # update pi (MH)
  logdensity0 <- pi_loglik(pi, M, U, N, u, n, eta, alpha, beta, alpha_q, beta_q)
  pi.new <- rbeta(1, alpha, beta)
  if((pi.new >= max(n_s/M_s)) & (pi.new <= min((M_s-N_s+n_s)/M_s))){
    logdensity1 <- pi_loglik(pi.new, M, U, N, u, n, eta, alpha, beta, alpha_q, beta_q)
    r <- logdensity1 - logdensity0
    if(log(runif(1)) < r) pi <- pi.new
  }
  
  # update q
  q.aux <- rbeta(1, sum(n)*eta+alpha_q, sum(pi*M-n)*eta+beta_q)
  
  # return new samples
  c(pi, q.aux)
}

covid_sampling <- function(N1, M, U, N, u, n, eta, alpha=1.01, beta=201, alpha_q=4.5, beta_q=1.5){
  pi.samples <- q.aux.samples <- q.samples <- numeric(N1)
  pi.samples[1] <- max(max(n_s/M_s), rbeta(1, alpha, beta)); q.aux.samples[1] <- rbeta(1, alpha_q, beta_q)
  
  for(i in 1:(N1-1)){
    up <- update_pow(pi.samples[i], q.aux.samples[i], M, U, N, u, n, eta, alpha, beta, alpha_q, beta_q)
    pi.samples[i+1] <- up[1]
    q.aux.samples[i+1] <- up[2]
  }
  
  for(i in 1:N1){
    q.samples[i] <- rbeta(1, sum(n)+alpha_q, sum(pi.samples[i]*M-n)+beta_q)
  }
  list("pi" = pi.samples, "q" = q.samples)
}

# LOOCV on n
elpd_n_syn <- numeric(n.eta)
for(t in 1:n.eta){
  pred_dens <- numeric(n_data)
  for(i in 1:length(pred_dens)){
    s <- covid_sampling(N1=1000000, M_s, U_s, N_s[-i], u_s, n_s[-i], eta=eta.cand[t], 1.01, 201, 4.5, 1.5)
    post.pi <- s$pi[100*(1:10000)]
    post.q <- s$q[100*(1:10000)]
    pred_dens[i] <- mean(exp(dbinom(n_s[i], ceiling(post.pi*M_s), post.q, log = TRUE) + dbinom(N_s[i]-n_s[i], M_s-ceiling(post.pi*M_s), (N_s[i]-n_s[i])/M_s, log = TRUE)))
  }
  elpd_n_syn[t] <- mean(log(pred_dens))
}

# Fig.19(a)
plot(eta.cand, elpd_n_syn, type = "l", xlab = "eta", ylab = "elpd")

# LOOCV on u
elpd_u_syn <- numeric(n.eta)
for(t in 1:n.eta){
  pred_dens <- numeric(n_data)
  for(i in 1:length(pred_dens)){
    s <- covid_sampling(N1=1000000, M_s, U_s, N_s[-i], u_s, n_s[-i], eta=eta.cand[t], 1.01, 201, 4.5, 1.5)
    post.pi <- s$pi[100*(1:10000)]
    pred_dens[i] <- mean(dhyper(u_s[i], ceiling(post.pi*M_s), M_s-ceiling(post.pi*M_s), U_s[i]))
  }
  elpd_u_syn[t] <- mean(log(pred_dens))
}

# Fig.19(b)
plot(eta.cand, elpd_u_syn, type = "l", xlab = "eta", ylab = "elpd")