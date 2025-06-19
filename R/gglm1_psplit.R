slurm_id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
library(rstan)
source("functions.R")

simnum <- 1
folder <- "gglm"

# p-value to be calculated
p_label <- "psplit"

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
m <- exp(X %*% beta)

set.seed(slurm_id)
y <- rgamma(n, alpha, alpha/m)
Xobs <- X[1:(n/2),]; yobs <- y[1:(n/2)]; zobs <- z[1:(n/2)]
Xnew <- X[(n/2+1):n,]; ynew <- y[(n/2+1):n]; znew <- z[(n/2+1):n]

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

calpsplit <- function(){
  
  # using n,nsim,nd,y,yrep_spt,ynew,Xnew,Xobs,p_spt_alphas,p_spt_betas from the global environment
  # decided to keep using "_spt" in variable names for ease of debugging
  
  Dobs <- Drep <- matrix(0,nsim,nd)
  
  p_new_m <- exp(tcrossprod(p_spt_betas, Xnew))
  p_new_b <- sweep(1/p_new_m, 1, p_spt_alphas, "*")
  
  ratio_obs <- outer(nsim1,ynew)/p_new_m
  ratio_rep <- yrep_spt/p_spt_m
  log_diff_obs <- log(ratio_obs); log_diff_rep <- log(ratio_rep)
  prop_diff_obs <- ratio_obs-1; prop_diff_rep <- ratio_rep-1
  
  # deviance (without scaling by alphas since won't affect results)
  Dobs[,1] <- rowSums(prop_diff_obs - log_diff_obs)
  Drep[,1] <- rowSums(prop_diff_rep - log_diff_rep)
  
  # chi squared test stats
  Dobs[,2] <- rowSums(prop_diff_obs^2) #* p_spt_alphas
  Drep[,2] <- rowSums(prop_diff_rep^2) #* p_spt_alphas
  
  # in model score
  Dobs3_beta <- rowSums(prop_diff_obs %*% Xnew) * p_spt_alphas
  Dobs[,3] <- -Dobs[,1] + Dobs3_beta
  
  Drep3_beta <- rowSums(prop_diff_rep %*% Xobs) * p_spt_alphas
  Drep[,3] <- -Drep[,1] + Drep3_beta
  
  # out of model score
  Dobs[,4] <- rowSums(prop_diff_obs * outer(nsim1, znew)) * p_spt_alphas
  Drep[,4] <- rowSums(prop_diff_rep * outer(nsim1, zobs)) * p_spt_alphas
  
  # KS
  cdf_obs <- matrix(pgamma(rep(ynew,each=nsim),rep(p_spt_alphas,n/2),c(p_new_b)),nrow=nsim)
  cdf_rep <- matrix(pgamma(c(yrep_spt),rep(p_spt_alphas,n/2),c(p_spt_b)), nrow=nsim)
  
  Dobs[,5] <- sapply(c(1:nsim), function (i) 
    ks.test(cdf_obs[i,], "punif")$statistic )
  Drep[,5] <- sapply(c(1:nsim), function (i) 
    ks.test(cdf_rep[i,], "punif")$statistic )
  
  pr_ge(Drep,Dobs)
}

runtime <- system.time({
  
  output <- matrix(0,nprior,nd,dimnames = list(prior_labels,d_labels))
  
  for (q in 1:nprior){
    
    data_spt <- list(n=n/2,p=p,y=y[1:(n/2)],X=X[1:(n/2),],mu0=mu0s[q],tau2=tau2s[q],
                     a1=a1s[q],b1=b1s[q],a2=a2s[q],b2=b2s[q])
    stan_mod_spt <- sampling(mod, data = data_spt,
                             chains = nchains,warmup = burnin,iter = (burnin+thin*nsim/nchains),
                             cores = nchains, thin = thin, seed = slurm_id, refresh=0)

    p_spt_alphas <- as.matrix(stan_mod_spt, pars="alpha")
    p_spt_betas <- as.matrix(stan_mod_spt, pars="beta")
    p_spt_m <- as.matrix(stan_mod_spt, pars="m")
    p_spt_b <- sweep(1/p_spt_m, 1, p_spt_alphas, "*")
    
    set.seed(slurm_id+1)
    yrep_spt <- matrix(rgamma(nsim*n/2,p_spt_alphas, c(p_spt_b)),nrow=nsim)

    output[q,] <- calpsplit()

  }
})

res <- list(runtime=runtime,output=output)

dir.create(paste0("../res/", folder))
output_name <- paste0("../res/",folder,"/sim", simnum, "_", p_label,"_",slurm_id, ".rds")
saveRDS(res, output_name)
