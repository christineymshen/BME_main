slurm_id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
n_cpus <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
library(doMC)
library(rstan)
source("functions.R")

registerDoMC(cores = n_cpus)

simnum <- 3
folder <- "gglm"

# p-value to be calculated
p_label <- "ptm"

# fixed parameters
burnin <- 1000
thin <- 5; nchains <- 1
nsim <- 1000 # number of MCMC posterior samples/ posterior predictive rep datasets
nsim1 <- rep(1,nsim)

# number of test statistics/ discrepancy measures
d_labels <- c("scaled deviance", "chi-squared", "in-model score", "out-model score", "KS")
nd <- length(d_labels)

p <- 6
n <- 50

set.seed(1)
X <- matrix(rnorm(n*(p-1)), nrow=n)
X <- cbind(rep(1,n),X)
beta <- rnorm(p,0,1.5)
z <- rnorm(n)
alpha <- 2

X1 <- cbind(X,z)
beta1 <- c(beta,0.1)
m <- exp(X1 %*% beta1)

set.seed(slurm_id)
y <- rgamma(n, alpha, alpha/m)

# prior
prior_labels <- c("good", "bad mean", "bad var", "both bad")
nprior <- length(prior_labels)

# prior for mu
tau2s <- c(4, 4/p, 4, 4/p) # used 4 instead of true para 1.5^2 because will set large prior for s2
mu0s <- sqrt(tau2s)*qnorm(0.95)* c(0,1,0,1)

# IG prior para for s2
# chose 3 because want prior sample size to be roughly 6=p (with ref to Bayesian textbook)
# then chose 14 s.t. the true para value is roughly at 0.05 ptile of the distribution
a1s <- c(0,3,3,3) # 0 -> use half cauchy prior
b1s <- c(0,8,14,14) # 0 -> use half cauchy prior

# Gamma prior for alpha
a2s <- c(0,6,12,12)
b2s <- c(0,3,3,3)

# read in compilation files
# note: need to compile separately on different machine
mod <- readRDS("../../20240131/R/gglm_mod3.rds")

calptm <- function(y,X,p_m_alphas,p_m_betas){
  
  # using nsim, nsim1, nd, X from the global environment
  ptm <- numeric(nd)

  # ptm1
  p_m_m <- c(exp(X %*% p_m_betas))
  m_ratio <- y/p_m_m
  m_prop_diff <- m_ratio-1; m_log_diff <- log(m_ratio)
  
  ptm[1] <- sum(m_prop_diff - m_log_diff) * p_m_alphas
  ptm[2] <- sum(m_prop_diff^2) * p_m_alphas
  ptm[3] <- -sum(m_prop_diff - m_log_diff) + n*(log(p_m_alphas)-digamma(p_m_alphas))
  ptm[3] <- ptm[3] + p_m_alphas*sum(m_prop_diff %*% X)
  ptm[4] <- p_m_alphas*sum(m_prop_diff * z)
  cdf_ptm <- pgamma(y,p_m_alphas,p_m_alphas/p_m_m)
  ptm[5] <- ks.test(cdf_ptm, "punif")$statistic

  return(ptm)

}

runtime <- system.time({
  
  output <- matrix(0,nprior,nd,dimnames = list(prior_labels,d_labels))
  
  for (q in 1:nprior){
    
    data <- list(n=n,p=p,y=y,X=X,mu0=mu0s[q],tau2=tau2s[q],
                 a1=a1s[q],b1=b1s[q],a2=a2s[q],b2=b2s[q])
    stan_mod <- sampling(mod, data = data,
                         chains = nchains, warmup = burnin,iter = (burnin+thin*nsim/nchains),
                         cores = nchains, thin = thin, seed = slurm_id, refresh=0)

    p_alphas <- as.matrix(stan_mod, pars="alpha")
    p_m_betas <- colMeans(as.matrix(stan_mod, pars="beta"))
    p_m_alphas <- mean(p_alphas)
    p_m <- as.matrix(stan_mod, pars="m")
    p_b <- sweep(1/p_m, 1, p_alphas, "*")
    set.seed(slurm_id+1)
    yrep <- matrix(rgamma(nsim*n,p_alphas, c(p_b)),nrow=nsim)
    
    ptm_obs <- calptm(y,X,p_m_alphas,p_m_betas)
    
    ptm_rep <- foreach(i=1:nsim, .combine="rbind") %dopar% {
      data$y <- yrep[i,]
      
      stan_mod_i <- sampling(mod, data = data,
                             chains = nchains, warmup = burnin,iter = (burnin+thin*nsim/nchains),
                             cores = nchains, thin = thin, seed = slurm_id, refresh=0)
      p_m_betas_i <- colMeans(as.matrix(stan_mod_i, pars="beta"))
      p_m_alphas_i <- mean(as.matrix(stan_mod_i, pars="alpha"))
 
      ptm_rep_i <- calptm(yrep[i,],X,p_m_alphas_i,p_m_betas_i)

      return(ptm_rep_i)
    }
    
    output[q,] <- pr_ge(ptm_rep,outer(nsim1, ptm_obs))

  }
})

res <- list(runtime=runtime,output=output)

dir.create(paste0("../res/", folder))
output_name <- paste0("../res/",folder,"/sim", simnum, "_", p_label,"_",slurm_id, ".rds")
saveRDS(res, output_name)
