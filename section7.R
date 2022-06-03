# Here is the code for section 7

# data
index1 <- c(5,6,9,10,11,13,14,15)
M_multi <- 9180135
N_multi <- southeast_df$Nt[index1]
n_multi <- southeast_df$nt[index1]
U_multi <- southeast_df$Nr[index1]
u_multi <- southeast_df$nr[index1]

n.eta <- 40
eta.cand <- seq(from=0, to=1, length.out=n.eta)

# LOOCV on n
elpd_n <- numeric(n.eta)
for(t in 1:n.eta){
  pred_dens <- numeric(8)
  for(i in 1:length(pred_dens)){
    up <- covid_sampling(N1=1000000, M_multi, U_multi, N_multi[-i], u_multi, n_multi[-i], eta.cand[t], 1.01, 201, 4.5, 1.5)
    post.pi <- up$pi[100*(1:10000)]
    post.q <- up$q[100*(1:10000)]
    pred_dens[i] <- mean(exp(dbinom(n_multi[i], ceiling(post.pi*M_multi), post.q, log = TRUE) + dbinom(N_multi[i]-n_multi[i], M_multi-ceiling(post.pi*M_multi), (N_multi[i]-n_multi[i])/M_multi, log = TRUE)))
  }
  elpd_n[t] <- mean(log(pred_dens))
}

# LOOCV on u
elpd_u <- numeric(n.eta)
for(t in 1:n.eta){
  pred_dens <- numeric(8)
  for(i in 1:length(pred_dens)){
    up <- covid_sampling(N1=1000000, M_multi, U_multi[-i], N_multi, u_multi[-i], n_multi, eta.cand[t], 1.01, 201, 4.5, 1.5)
    post.pi <- up$pi[100*(1:10000)]
    pred_dens[i] <- mean(dhyper(u_multi[i], ceiling(post.pi*M_multi), M_multi-ceiling(post.pi*M_multi), U_multi[i]))
  }
  elpd_u[t] <- mean(log(pred_dens))
}

# LOOCV on both u and n
elpd_both <- numeric(n.eta)
for(t in 1:n.eta){
  pred_dens <- numeric(8)
  for(i in 1:length(pred_dens)){
    up <- covid_sampling(N1=1000000, M_multi, U_multi[-i], N_multi[-i], u_multi[-i], n_multi[-i], eta.cand[t], 1.01, 201, 4.5, 1.5)
    post.pi <- up$pi[100*(1:10000)]
    post.q <- up$q[100*(1:10000)]
    pred_dens[i] <- mean(exp(dbinom(n_multi[i], ceiling(post.pi*M_multi), post.q, log = TRUE) + dbinom(N_multi[i]-n_multi[i], M_multi-ceiling(post.pi*M_multi), (N_multi[i]-n_multi[i])/M_multi, log = TRUE) + dhyper(u_multi[i], ceiling(post.pi*M_multi), M_multi-ceiling(post.pi*M_multi), U_multi[i], log = TRUE)))
  }
  elpd_both[t] <- mean(log(pred_dens))
}

## Fig.8(a,b,c)
plot(eta.cand, elpd_n, type="l", xlab="eta", ylab="elpd")
plot(eta.cand, elpd_u, type="l", xlab="eta", ylab="elpd")
plot(eta.cand, elpd_both, type="l", xlab="eta", ylab="elpd")

## Fig.8(d)
covid_sample_0_multi <- covid_sampling(N1=10000, M_multi, U_multi, N_multi, u_multi, n_multi, eta=0)
covid_sample_opt_multi <- covid_sampling(N1=10000, M_multi, U_multi, N_multi, u_multi, n_multi, eta=0.7)
covid_sample_1_multi <- covid_sampling(N1=10000, M_multi, U_multi, N_multi, u_multi, n_multi, eta=1)

plot(covid_sample_0_multi$pi, covid_sample_0_multi$q, pch=".", xlab="pi", ylab="q")
lines(covid_sample_opt_multi$pi, covid_sample_opt_multi$q, col="green", pch=".")
lines(covid_sample_1_multi$pi, covid_sample_1_multi$q, col="red", pch=".")
legend("topright", legend = c("eta=0 (Cut)","eta=1 (Bayes)", "eta=0.7"), col=c("black","red","green"), cex=0.6)


# prior-to-posterior probability (multiple data version)
f_m <- function(pi){
  exp(lbeta(sum(n_multi)*eta+alpha_q, sum(ceiling(pi*M_multi)-n_multi)*eta+beta_q)
      + sum(dhyper(u_multi, ceiling(pi*M_multi), M_multi-ceiling(pi*M_multi), U_multi, log = TRUE)) 
      + eta*sum(dbinom(N_multi-n_multi, M_multi-ceiling(pi*M_multi), (N_multi - n_multi)/M_multi, log = TRUE)) 
      + eta*sum(lchoose(ceiling(pi*M_multi),n_multi)) + (alpha-1)*log(pi) + (beta-1)*log(1-pi))
}

p_cut_m <- function(pi){
  exp(sum(dhyper(u_multi, ceiling(pi*M_multi), M_multi-ceiling(pi*M_multi), U_multi, log = TRUE)) + 
        dbeta(pi, alpha, beta, log=TRUE))
}

L <- max(n_multi/M_multi); K <- 1000
x <- L + delta^(-(1:K)); max_x_index <- 0
for(k in 1:K){
  if(x[k] <= min((M_multi-N_multi+n_multi)/M_multi)) max_x_index <- max_x_index+1
}
x <- x[K:(K-max_x_index+1)]; K <- max_x_index

inte_cut <- 0
for(i in 1:(K-1)){inte_cut <- inte_cut + p_cut_m(x[i])*(x[i+1]-x[i])}
p_cut1 <- numeric(K-1)
for(i in 1:(K-1)){p_cut1[i] <- p_cut_m(x[i])/inte_cut}
cdf_cut_m <- cumsum(p_cut1*(x[2:K]-x[1:(K-1)]))

f_eta_pi_n_m <- function(pi, inte, inte_cut){
  lbeta(sum(n_multi)*eta+alpha_q, sum(ceiling(pi*M_multi)-n_multi)*eta+beta_q)
  + sum(dhyper(u_multi, ceiling(pi*M_multi), M_multi-ceiling(pi*M_multi), U_multi, log = TRUE)) 
  + eta*sum(dbinom(N_multi-n_multi, M_multi-ceiling(pi*M_multi), p.hat, log = TRUE)) 
  + eta*sum(lchoose(ceiling(pi*M_multi),n_multi)) + (alpha-1)*log(pi) + (beta-1)*log(1-pi) 
  - sum(dhyper(u_multi, ceiling(pi*M_multi), M_multi-ceiling(pi*M_multi), U_multi, log = TRUE))
  - dbeta(pi, alpha, beta, log=TRUE) - log(inte) + log(inte_cut)
}

n.eta <- 101
eta.cand <- seq(from=0, to=1, length.out=n.eta)
d_nobs_eta <- numeric(n.eta)
for(i in 1:n.eta){
  eta <- eta.cand[i]
  inte <- 0; for(z in 1:(K-1)){inte <- inte + f_m(x[z])*(x[z+1]-x[z])}
  f_m1 <- numeric(K-1)
  for(k in 1:(K-1)){f_m1[k] <- f_m(x[k])/inte}
  cdf_em <- cumsum(f_m1*(x[2:K]-x[1:(K-1)]))
  pi_smi_sample <- pi_inver(N=1000) 
  d_nobs_eta[i] <- mean(f_eta_pi_n_m(pi_smi_sample, inte, inte_cut))
}

n_mcmc_m <- function(N1, pi, q){
  n_mcmc_sample <- numeric(N1)
  n_mcmc_sample[1] <- n_multi[1]
  for(i in 1:(N1-1)){
    n0 <- n_mcmc_sample[i]
    logdensity0 <- 8*dbinom(n0, ceiling(pi*M), q, log = TRUE)+sum(dbinom(N_multi-n0, M-ceiling(pi*M), (N_multi-n0)/M, log = TRUE))
    n1 <- sample(c(n0+1, n0-1),1)
    if((n1>=pi*M) | (n1 <= max(N_multi-M+pi*M))) n_mcmc_sample[i+1] <- n0
    else{
      logdensity1 <- 8*dbinom(n1, ceiling(pi*M), q, log = TRUE)+sum(dbinom(N_multi-n1, M-ceiling(pi*M), (N_multi-n1)/M, log = TRUE))
      r <- logdensity1 - logdensity0
      if(log(runif(1)) < r) n_mcmc_sample[i+1] <- n1
      else{n_mcmc_sample[i+1] <- n0}
    }
  }
  return(n_mcmc_sample)
}

pi_inver_cut_m <- function(N=10000){
  pi_sample <- numeric(N)
  for(i in 1:N){
    r <- runif(1)
    for(j in 1:(K-1)){
      if((r >= cdf_cut_m[j]) & (r < cdf_cut_m[j+1])){
        pi_sample[i] <- x[j]
      }
    }
  }
  pi_sample
}

pi_cut_sample <- pi_inver_cut_m(200) # sample pi from the cut posterior
q_sample <- rbeta(200, alpha_q, beta_q) # sample q from prior
n_sample <- matrix(NA, nrow = 200, ncol = 10)
for(i in 1:200){ # we use S = 1
  n_sample[i,] <- n_mcmc_m(N1=100000, pi_cut_sample[i], q_sample[i])[5000*(10:19)]
}

eta <- 1; d_n <- numeric(200)
for(i in 1:200){
  n_multi <- n_sample[i,3:10]
  L <- max(n_multi/M); K <- 1000; delta <- (10^8)^(1/K)
  x <- L + delta^(-(1:K)); max_x_index <- 0
  for(k in 1:K){
    if(x[k] <= min((M-N_multi+n_multi)/M)) max_x_index <- max_x_index+1
  }
  x <- x[K:(K-max_x_index+1)]; K <- max_x_index
  inte_full <- 0
  for(z in 1:(K-1)){inte_full <- inte_full + f_m(x[z])*(x[z+1]-x[z])}
  if(inte_full == 0){d_n[i] <- 0}
  else{
    for(k in 1:(K-1)){f_m1[k] <- f_m(x[k])/inte_full}
    cdf_em <- cumsum(f_m1*(x[2:K]-x[1:(K-1)]))
    pi_full_sample <- pi_inver(1000)
    inte_cut_new_n <- 0
    for(z in 1:(K-1)){inte_cut_new_n <- inte_cut_new_n + p_cut_m(x[z])*(x[z+1]-x[z])}
    d_n[i] <- mean(f_eta_pi_n_m(pi_full_sample, inte_full, inte_cut_new_n))
  }
}

d_n_valid <- d_n[d_n > 0]

pr <- numeric(n.eta)
for(i in 1:n.eta){
  pr[i] <- mean(d_n_valid > d_nobs_eta[i])
}

count <- 0
for(i in 1:n.eta){
  if(pr[i] >= 0.05) count <- count+1
}

# Fig.9
plot(eta.cand, pr, xlab = "eta", ylab = "tail probability", type = "l")
abline(v=eta.cand[16], lty=2, col="red")
abline(h=0.05, col="grey", lty=2)