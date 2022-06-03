# code: compute prior-to-posterior tail probability, and choose the optimal eta
# corresponding to section 5 in the paper

# import data
data_num <- 1 # choose first week data (this should change each time)

M <- southeast_df$M[data_num]
N <- southeast_df$Nt[data_num]
n <- southeast_df$nt[data_num]
U <- southeast_df$Nr[data_num]
u <- southeast_df$nr[data_num]
p.hat <- (N-n)/M

# MCMC to sample n 
n_mcmc <- function(N1, pi, q){
  n_mcmc_sample <- numeric(N1)
  n_mcmc_sample[1] <- n
  for(i in 1:(N1-1)){
    n0 <- n_mcmc_sample[i]
    logdensity0 <- dbinom(n0, ceiling(pi*M), q, log = TRUE)+dbinom(N-n0, M-ceiling(pi*M), (N-n0)/M, log = TRUE)
    n1 <- sample(c(n0+1, n0-1),1)
    if((n1>=pi*M) | (n1 <= N-M+pi*M)) n_mcmc_sample[i+1] <- n0
    else{
      logdensity1 <- dbinom(n1, ceiling(pi*M), q, log = TRUE)+dbinom(N-n1, M-ceiling(pi*M), (N-n1)/M, log = TRUE)
      r <- logdensity1 - logdensity0
      if(log(runif(1)) < r) n_mcmc_sample[i+1] <- n1
      else{n_mcmc_sample[i+1] <- n0}
    }
  }
  return(n_mcmc_sample)
}

# candidate eta
eta.cand <- seq(from=0, to=1, length.out=101)
n.eta <- length(eta.cand)

# unnormalized smi-posterior for pi
f <- function(pi){ 
  exp(lbeta(n*eta+alpha_q, ceiling(pi*M)*eta-n*eta+beta_q) +
 dhyper(u, ceiling(pi*M), M-ceiling(pi*M), U, log = TRUE) + 
 eta*dbinom(N-n, M-ceiling(pi*M), p.hat, log = TRUE) + eta*lchoose(ceiling(pi*M),n) 
 + (alpha-1)*log(pi) + (beta-1)*log(1-pi))
}

# set points to estimate smi-posterior for pi

L <- n/M; K <- 1000
x <- L + delta^(-(1:K)); max_x_index <- 0
for(k in 1:K){
  if(x[k] <= (M-N+n)/M) max_x_index <- max_x_index+1
}
x <- x[K:(K-max_x_index+1)]; K <- max_x_index # make sure that each x is in the support

# estimate normalizing constant for smi-posterior for pi
inte <- 0
for(i in 1:(K-1)){inte <- inte + f(x[i])*(x[i+1]-x[i])}

# unnormalized cut posterior for pi
p_cut <- function(pi){
  exp(dhyper(u, ceiling(pi*M), M-ceiling(pi*M), U, log = TRUE) + 
   dbeta(pi, alpha, beta, log=TRUE))
}

# estimate normalizing constant for cut posterior for pi
inte_cut <- 0
for(i in 1:(K-1)){inte_cut <- inte_cut + p_cut(x[i])*(x[i+1]-x[i])}

# estimate cdf of cut posterior
cdf_cut <- cumsum((p_cut(x[1:(K-1)])/inte_cut)*(x[2:K]-x[1:(K-1)]))

# inversion method to sample pi from cut posterior
pi_inver_cut <- function(N=10000){
  pi_sample <- numeric(N)
  for(i in 1:N){
    r <- runif(1)
    for(j in 1:(K-1)){
      if((r >= cdf_cut[j]) & (r < cdf_cut[j+1])){
        pi_sample[i] <- x[j]
      }
    }
  }
  pi_sample
}

# f_eta(pi;n) (Eq.26)
f_eta_pi_n <- function(pi, inte, inte_cut){
  lbeta(n*eta+alpha_q, ceiling(pi*M)*eta-n*eta+beta_q) +
    dhyper(u, ceiling(pi*M), M-ceiling(pi*M), U, log = TRUE) + 
    eta*dbinom(N-n, M-ceiling(pi*M), p.hat, log = TRUE) + eta*lchoose(ceiling(pi*M),n) 
  + (alpha-1)*log(pi) + (beta-1)*log(1-pi) - dhyper(u, ceiling(pi*M), M-ceiling(pi*M), U, log = TRUE) 
  - dbeta(pi, alpha, beta, log=TRUE) - log(inte) + log(inte_cut)
}

## Following Algorithm 3
# d(n_obs;eta) for different eta values
d_nobs_eta <- numeric(n.eta)
for(i in 1:n.eta){
  eta <- eta.cand[i]
  inte <- 0; for(z in 1:(K-1)){inte <- inte + f(x[z])*(x[z+1]-x[z])}
  cdf_em <- cumsum((f(x[1:(K-1)])/inte)*(x[2:K]-x[1:(K-1)]))
  pi_smi_sample <- pi_inver(N=1000) 
  d_nobs_eta[i] <- mean(f_eta_pi_n(pi_smi_sample, inte, inte_cut))
}

# d(n;eta=1)
pi_cut_sample <- pi_inver_cut(100) # sample pi from the cut posterior
q_sample <- rbeta(100, alpha_q, beta_q) # sample q from prior
n_sample <- numeric(100)
for(i in 1:100){ # we use S = 100
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

# epsilon_eta (Eq.22)
pr <- numeric(n.eta)
for(i in 1:n.eta){
  pr[i] <- mean(d_n > d_nobs_eta[i])
}

# choose the maximum eta with pr_eta > 0.05

## Fig.5(a) smi posterior of pi with n_obs, varied eta
eta.sub <- seq(from=0, to=1, length.out=11) # select eta=0, 0.1, 0.2, ..., 1
eta <- eta.sub[1]
h <- hist(pi_inver(), breaks=40, plot=F)
plot(h$mids, h$density, type="l", xlab = "pi", ylab = "density")
for(i in 2:11){
  eta <- eta.sub[i]
  h <- hist(pi_inver(), breaks=40, plot=F)
  lines(h$mids, h$density, type="l", col=i)
}

## Fig.5(b) smi posterior of pi with random n, eta=1
eta <- 1
pi_cut_sample <- pi_inver_cut(10) # sample pi from the cut posterior
q_sample <- rbeta(10, alpha_q, beta_q)# sample q from prior
u_sample <- numeric(10)
n_sample <- numeric(10)
for(i in 1:10){ # we use S = 1
  u_sample[i] <- rhyper(1, ceiling(pi_cut_sample[i]*M), M-ceiling(pi_cut_sample[i]*M), U)
  n_sample[i] <- n_mcmc(N1=100000, pi_cut_sample[i], q_sample[i])[100000]
}

# plot the first data
n <- n_sample[1]
u <- u_sample[1]
L <- n/M; K <- 1000
x <- L + delta^(-(1:K)); max_x_index <- 0
for(k in 1:K){
  if(x[k] <= (M-N+n)/M) max_x_index <- max_x_index+1
}
x <- x[K:(K-max_x_index+1)]; K <- max_x_index
inte_full <- 0
for(z in 1:(K-1)){inte_full <- inte_full + f(x[z])*(x[z+1]-x[z])}
cdf_em <- cumsum((f(x[1:(K-1)])/inte_full)*(x[2:K]-x[1:(K-1)]))

inte_cut <- 0
for(i in 1:(K-1)){
  inte_cut <- inte_cut + p_cut(x[i])*(x[i+1]-x[i])
}
cdf_cut <- cumsum((p_cut(x[1:(K-1)])/inte_cut)*(x[2:K]-x[1:(K-1)]))

hr <- hist(pi_inver(), breaks=10, plot=F)
hcut <- hist(pi_inver_cut(), breaks=10, plot=F)
plot(hr$mids, hr$density, type="l", xlab="pi (sampled n)", ylab="density", xlim=c(0,0.3), ylim=c(0,400))
lines(hcut$mids, hcut$density, lty=2) # the dotted line is the associated cut posterior

# plot the following data
for(i in 2:10){
  n <- n_sample[i]
  u <- u_sample[i]
  L <- n/M; K <- 1000
  x <- L + delta^(-(1:K)); max_x_index <- 0
  for(k in 1:K){
    if(x[k] <= (M-N+n)/M) max_x_index <- max_x_index+1
  }
  x <- x[K:(K-max_x_index+1)]; K <- max_x_index
  inte_full <- 0
  for(z in 1:(K-1)){inte_full <- inte_full + f(x[z])*(x[z+1]-x[z])}
  cdf_em <- cumsum((f(x[1:(K-1)])/inte_full)*(x[2:K]-x[1:(K-1)]))
  inte_cut <- 0
  for(z in 1:(K-1)){
    inte_cut <- inte_cut + p_cut(x[i])*(x[i+1]-x[i])
  }
  cdf_cut <- cumsum((p_cut(x[1:(K-1)])/inte_cut)*(x[2:K]-x[1:(K-1)]))
  hr <- hist(pi_inver(), breaks=10, plot=F)
  hcut <- hist(pi_inver_cut(), breaks=10, plot=F)
  lines(hr$mids, hr$density, type="l", col=i)
  lines(hcut$mids, hcut$density, col=i, lty=2)
}










