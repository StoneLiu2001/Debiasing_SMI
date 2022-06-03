# Here is the code for inversion method to sample pi, q from smi posterior and to produce Fig.3.
# corresponding to Appendix E in the paper

# unnormalized smi-posterior for pi (g_eta(pi,u,n), Eq.16)
f <- function(pi){ 
  exp(lbeta(n*eta+alpha_q, ceiling(pi*M)*eta-n*eta+beta_q) +
        dhyper(u, ceiling(pi*M), M-ceiling(pi*M), U, log = TRUE) + 
        eta*dbinom(N-n, M-ceiling(pi*M), p.hat, log = TRUE) + eta*lchoose(ceiling(pi*M),n) 
      + (alpha-1)*log(pi) + (beta-1)*log(1-pi))
}

# set points to estimate smi-posterior for pi
L <- n/M
K <- 1000
delta <- (10^8)^(1/K)
x <- L + delta^(-(1:K))
x <- x[K:1]

# estimate normalizing constant for smi-posterior for pi
inte <- 0
for(i in 1:(K-1)){inte <- inte + f(x[i])*(x[i+1]-x[i])}

# estimated (empirical) cdf
cdf_em <- cumsum((f(x[1:(K-1)])/inte)*(x[2:K]-x[1:(K-1)]))

# inversion method to sample pi from smi posterior
pi_inver <- function(N=10000){
  pi_sample <- numeric(N)
  for(i in 1:N){
    r <- runif(1)
    for(j in 1:(K-1)){
      if((r >= cdf_em[j]) & (r < cdf_em[j+1])){
        pi_sample[i] <- x[j]
      }
    }
  }
  pi_sample
}

# joint scatter plot of (pi,q) from smi posterior with eta=1 using both MCMC and inversion (Fig.13)
eta <- 1
inte <- 0
for(i in 1:(K-1)){inte <- inte + f(x[i])*(x[i+1]-x[i])}
cdf_em <- cumsum((f(x[1:(K-1)])/inte)*(x[2:K]-x[1:(K-1)]))

pi_inverse_sample <- pi_inver(100000)
q_sample <- numeric(100000)
for(i in 1:100000){q_sample[i] <- rbeta(1, n+alpha_q, pi_inverse_sample[i]*M-n+beta_q)}

covid_sample_1 <- covid_sampling(N1=100000, M, U, N, u, n, eta=1, alpha=1.01, beta=1.01*199, alpha_q=4.5, beta_q=1.5)
plot(covid_sample_1$pi[100*(1:1000)], covid_sample_1$q[100*(1:1000)], xlab = "pi (MCMC)", ylab = "q",pch='.', col="blue")
points(pi_inverse_sample[100*(1:1000)], q_sample[100*(1:1000)], col="red",pch='.')
legend("topright", legend = c("MCMC", "inversion"), col=c("blue", "red"), pch=1)

# Fig.16 (eta=0.92, optimal)
eta <- 0.92
inte_opt <- 0
for(i in 1:(K-1)){inte_opt <- inte_opt + f(x[i])*(x[i+1]-x[i])}
cdf_em_opt <- cumsum((f(x[1:(K-1)])/inte_opt)*(x[2:K]-x[1:(K-1)]))

pi_inverse_sample_opt <- pi_inver(100000)
q_sample_opt <- numeric(100000)
for(i in 1:100000){q_sample_opt[i] <- rbeta(1, n+alpha_q, pi_inverse_sample_opt[i]*M-n+beta_q)}

covid_sample_opt <- covid_sampling(N1=100000, M, U, N, u, n, eta=0.92, alpha=1.01, beta=1.01*199, alpha_q=4.5, beta_q=1.5)
plot(covid_sample_1$pi[100*(1:1000)], covid_sample_1$q[100*(1:1000)], xlab = "pi (MCMC)", ylab = "q",pch='.', col="blue")
points(pi_inverse_sample_opt[100*(1:1000)], q_sample_opt[100*(1:1000)], col="red",pch='.')
legend("topright", legend = c("MCMC", "inversion"), col=c("blue", "red"), pch=1)


## under uniform prior

# unnormalized smi-posterior for pi under uniform prior
f_unif <- function(pi){
  exp(lbeta(n*eta+1, ceiling(pi*M)*eta-n*eta+1) +
        dhyper(u, ceiling(pi*M), M-ceiling(pi*M), U, log = TRUE) + eta*dbinom(N-n, M-ceiling(pi*M), p.hat, log = TRUE) + eta*lchoose(ceiling(pi*M),n))
}

# integrate f_unif from 0 to 1, normalizing constant
inte_unif <- 0
for(i in 1:(K-1)){
  inte_unif <- inte_unif + f_unif(x[i])*(x[i+1]-x[i])
}

# estimated cdf under uniform prior
cdf_em_unif <- cumsum((f_unif(x[1:(K-1)])/inte_unif)*(x[2:K]-x[1:(K-1)]))

# inversion method to sample pi,q from smi-posterior under uniform prior
pi_inver_unif <- function(N=10000){
  pi_sample <- numeric(N)
  for(i in 1:N){
    r <- runif(1)
    for(j in 1:(K-1)){
      if((r >= cdf_em_unif[j]) & (r < cdf_em_unif[j+1])){
        pi_sample[i] <- x[j]
      }
    }
  }
  pi_sample
}

pi_inverse_sample_unif <- pi_inver_unif(N=100000)
q_inverse_unif <- numeric(100000)
for(i in 1:100000){q_inverse_unif[i] <- rbeta(1, n+1, pi_inverse_sample[i]*M-n+1)}

## Fig.3(a)
# MCMC + subjective prior
h <- hist(covid_sample_1$pi[100*(1:1000)], plot = F, breaks = seq(n/M, 0.005, length.out=100))

# inversion + subjective prior
h_inv <- hist(pi_inverse_sample, plot = F, breaks = seq(n/M, 0.005, length.out=100))

# prior
r1 <- rbeta(10000, 1.01, 201)
h1 <- hist(r1, plot=F, breaks = 100)

# MCMC + uniform prior
h_unif <- hist(covid_sample_unif_1$pi[500*(1:2000)], plot = F, breaks = 30)

# inversion + uniform prior
h_unif_inv <- hist(pi_inverse_sample_unif, plot = F, breaks = 30)

plot(h_inv$mids, h_inv$density, type="l", col="red", xlab = "pi", ylab = "density")
lines(h$mids, h$density, type="l")
lines(h_unif$mids, h_unif$density, col="green")
lines(seq(n/M, 0.005, length.out=50), dbeta(seq(n/M, 0.005, length.out=50), 1.01, 201), col="blue")
abline(v=n/M, lty=2,col="purple")
lines(seq(n/M, 0.005, length.out=100), f(seq(n/M, 0.005, length.out=100))/inte, col="grey")
lines(h_unif_inv$mids, h_unif_inv$density, col="pink")
legend("topright", legend = c("prior", "posterior with subjective prior (MCMC)", "posterior with subjective prior (inversion)", "posterior with uniform prior (MCMC)", "posterior with uniform prior (inversion)", "true density"), col = c("blue", "black", "red", "green", "pink", "grey"), lty=1, cex=0.6)

## Fig.3(b)
# MCMC + subjective prior
h_q <- hist(covid_sample_1$q, plot = F, breaks = 20)

# prior
r2 <- rbeta(10000, 4.5, 1.5)
h2 <- hist(r2, plot=F, breaks=20)

# inversion + subjective prior
h_q_inv <- hist(q_sample, plot = F, breaks = 20)

# MCMC + uniform prior
h_unif_q <- hist(covid_sample_unif_1$q, plot = F, breaks = 50)

# inversion + uniform prior
h_unif_q_inv <- hist(q_inverse_unif, plot = F, breaks = 80)

plot(h_unif_q$mids, h_unif_q$density, col="green", type="l", xlab = "q", ylab = "density")
lines(h_q_inv$mids, h_q_inv$density, col="red")
lines(h_q$mids, h_q$density)
lines(h2$mids, h2$density, col="blue")
lines(h_unif_q_inv$mids, h_unif_q_inv$density, col="pink")
legend("topright", legend = c("prior", "posterior with subjective prior (MCMC)", "posterior with subjective prior (inversion)", "posterior with uniform prior (MCMC)", "posterior with uniform prior (inversion)"), col = c("blue", "black", "red", "green", "pink"), lty=1, cex=0.6)



