# Here is the code for MCMC and to produce figures in Appendix B.

# load data
load("example.RData")

# EDA (Fig.1)
# week number with U>0
valid_index <- (southeast_df$Nr>0)*(1:length(southeast_df$Nr)) 

plot(southeast_df$nr[valid_index]/southeast_df$Nr[valid_index], xlab = "week", ylab = "u/U", type = "l")
axis(1, at = c(0,5,10,15,20,25,30,35,40))

plot(southeast_df$nt[valid_index]/southeast_df$Nt[valid_index], xlab = "week", ylab = "n/N", type = "l")
axis(1, at = c(0,5,10,15,20,25,30,35,40))

# In section 4 and 5, we use first week data to analyze.
M <- southeast_df$M[1]
N <- southeast_df$Nt[1]
n <- southeast_df$nt[1]
U <- southeast_df$Nr[1]
u <- southeast_df$nr[1]
p.hat <- (N-n)/M

# In section 7, we use data from week 3 to week 10.
index1 <- c(5,6,9,10,11,13,14,15)
M_multi <- 9180135
N_multi <- southeast_df$Nt[index1]
n_multi <- southeast_df$nt[index1]
U_multi <- southeast_df$Nr[index1]
u_multi <- southeast_df$nr[index1]

# prior parameter
alpha = 1.01; beta = 1.01*199 # (pi prior, Eq.7)
alpha_q = 4.5; beta_q = 1.5 # (q prior, Eq.8)

# subjective prior
# log likelihood of g_eta(pi,u,n) under subjective prior (Eq.16 & Eq.29)
pi_loglik <- function(pi, M, U, N, u, n, eta, alpha, beta, alpha_q, beta_q){
  lbeta(sum(n)*eta+alpha_q, sum(ceiling(pi*M)-n)*eta+beta_q) 
  + sum(dhyper(u, ceiling(pi*M), M-ceiling(pi*M), U, log = TRUE)) 
  + eta*sum(dbinom(N-n, M-ceiling(pi*M), p.hat, log = TRUE)) 
  + eta*sum(lchoose(ceiling(pi*M),n)) + (alpha-1)*log(pi) + (beta-1)*log(1-pi)
}

# Metropolis within Gibbs sampling targeting pi and tilde q smi posterior under subjective prior (Eq.13)
update_pow <- function(pi, q.aux, M, U, N, u, n, eta, alpha, beta, alpha_q, beta_q){
  
  # update pi (MH)
  logdensity0 <- pi_loglik(pi, M, U, N, u, n, eta, alpha, beta, alpha_q, beta_q)
  pi.new <- rbeta(1, alpha, beta)
  if((pi.new >= max(n/M)) & (pi.new <= min((M-N+n)/M))){
    logdensity1 <- pi_loglik(pi.new, M, U, N, u, n, eta, alpha, beta, alpha_q, beta_q)
    r <- logdensity1 - logdensity0
    if(log(runif(1)) < r) pi <- pi.new
  }
  
  # update q
  q.aux <- rbeta(1, sum(n)*eta+alpha_q, sum(pi*M-n)*eta+beta_q)
  
  # return new samples
  c(pi, q.aux)
}

# MCMC sampling for pi, tilde q, and q (Algorithm 1) under subjective prior
covid_sampling <- function(N1, M, U, N, u, n, eta, alpha, beta, alpha_q, beta_q){
  pi.samples <- q.aux.samples <- q.samples <- numeric(N1)
  pi.samples[1] <- max(max(n/M),rbeta(1, alpha, beta)); q.aux.samples[1] <- rbeta(1, alpha_q, beta_q)
  
  for(i in 1:(N1-1)){
    pi.samples[i+1] <- update_pow(pi.samples[i], q.aux.samples[i], M, U, N, u, n, eta, alpha, beta, alpha_q, beta_q)[1]
    q.aux.samples[i+1] <- update_pow(pi.samples[i], q.aux.samples[i], M, U, N, u, n, eta, alpha, beta, alpha_q, beta_q)[2]
  }
  
  for(i in 1:N1){
    q.samples[i] <- rbeta(1, sum(n)+alpha_q, sum(pi.samples[i]*M-n)+beta_q)
  }
  list("pi" = pi.samples, "q" = q.samples)
}

# uniform prior
# log likelihood of pi posterior unnormalized under uniform prior
pi_loglik_unif <- function(pi, M, U, N, u, n, eta){
  lbeta(n*eta+1, eta*(ceiling(pi*M)-n)+1) + dhyper(u, ceiling(pi*M), M-ceiling(pi*M), U, log = TRUE) 
  + eta*dbinom(N-n, M-ceiling(pi*M), p.hat, log = TRUE) + eta*lchoose(ceiling(pi*M),n)
}

# Metropolis within Gibbs sampling targeting pi and tilde q smi posterior under uniform prior
update_pow_unif <- function(pi, q.aux, M, U, N, u, n, eta){
  
  # update pi (MH)
  logdensity0 <- pi_loglik_unif(pi, M, U, N, u, n, eta)
  pi.new <- runif(1)
  if((pi.new >= (n/M)) & (pi.new <= (M-N+n)/M)){
    logdensity1 <- pi_loglik_unif(pi.new, M, U, N, u, n, eta)
    r <- (logdensity1 - logdensity0)
    if(log(runif(1)) < r) pi <- pi.new
  }
  
  # update q
  q.aux <- rbeta(1, n*eta+1, (pi*M-n)*eta+1)
  
  # return new samples
  c(pi, q.aux)
}

# MCMC sampling for pi, tilde q, and q under uniform prior
covid_sampling_unif <- function(N1, M, U, N, u, n, eta){
  pi.samples <- q.aux.samples <- q.samples <- numeric(N1)
  pi.samples[1] <- runif(1); q.aux.samples[1] <- runif(1)
  
  for(i in 1:(N1-1)){
    pi.samples[i+1] <- update_pow_unif(pi.samples[i], q.aux.samples[i], M, U, N, u, n, eta)[1]
    q.aux.samples[i+1] <- update_pow_unif(pi.samples[i], q.aux.samples[i], M, U, N, u, n, eta)[2]
  }
  
  for(i in 1:N1){
    q.samples[i] <- rbeta(1, n+1, pi.samples[i]*M-n+1)
  }
  list("pi" = pi.samples, "q" = q.samples)
}

# Display some results

## Fig.10 (eta=1 under subjective prior)
covid_sample_1 <- covid_sampling(N1=100000, M, U, N, u, n, eta=1, alpha=1.01, beta=1.01*199, alpha_q=4.5, beta_q=1.5)
plot(covid_sample_1$pi[100*(1:1000)], type = "l", xlab = "MCMC step t=1, 2, 3,...", ylab = "MCMC state pi_t")
plot(covid_sample_1$q[100*(1:1000)], type = "l", xlab = "MCMC step t=1, 2, 3,...", ylab = "MCMC state q_t")

## Fig.11 (ACF plot)
acf(covid_sample_1$pi, lag.max = 200, main="")
acf(covid_sample_1$pi[100*(1:1000)], main="", lag.max = 200)

## Fig.12 (eta=1 under uniform prior)
covid_sample_unif_1 <- covid_sampling_unif(N1=1000000, M, U, N, u, n, eta=1)
plot(covid_sample_unif_1$pi[100*(1:10000)], type = "l", xlab = "MCMC step t=1, 2, 3,...", ylab = "MCMC state pi_t")
plot(covid_sample_unif_1$q[100*(1:10000)], type = "l", xlab = "MCMC step t=1, 2, 3,...", ylab = "MCMC state q_t")

## Fig.14 (eta=0.92 under subjective prior)
covid_sample_opt <- covid_sampling(N1=100000, M, U, N, u, n, eta=0.92, alpha=1.01, beta=1.01*199, alpha_q=4.5, beta_q=1.5)
plot(covid_sample_opt$pi[100*(1:1000)], type = "l", xlab = "MCMC step t=1, 2, 3,...", ylab = "MCMC state pi_t")
plot(covid_sample_opt$q[100*(1:1000)], type = "l", xlab = "MCMC step t=1, 2, 3,...", ylab = "MCMC state q_t")

## Fig.15 (eta=0.92 under uniform prior)
covid_sample_unif_opt <- covid_sampling_unif(N1=1000000, M, U, N, u, n, eta=0.92)
plot(covid_sample_unif_opt$pi[1000*(1:1000)], type = "l", xlab = "MCMC step t=1, 2, 3,...", ylab = "MCMC state pi_t")
plot(covid_sample_unif_opt$q[1000*(1:1000)], type = "l", xlab = "MCMC step t=1, 2, 3,...", ylab = "MCMC state q_t")

## Fig.17 (eta=1 multiple data under subjective prior)
covid_multi <- covid_sampling(N1=100000, M_multi, U_multi, N_multi, u_multi, n_multi, eta=1, alpha=1.01, beta=1.01*199, alpha_q=4.5, beta_q=1.5)
plot(covid_multi$pi[10*(2:10000)], type = "l", xlab = "MCMC step t=1, 2, 3,...", ylab = "MCMC state pi_t")
plot(covid_multi$q[10*(2:10000)], type = "l", xlab = "MCMC step t=1, 2, 3,...", ylab = "MCMC state q_t")

## Fig.18 (ACF plot)
acf(covid_multi$pi[10*(2:10000)], main="", lag.max = 100)
acf(covid_multi$q[10*(2:10000)], main="", lag.max = 100)


