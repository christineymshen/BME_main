slurm_id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
n_cpus <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
library(doMC)
library(rstan)
source("functions.R")

registerDoMC(cores = n_cpus)

simnum <- 1
folder <- "gglm"

# p-value to be calculated
p_label <- "pptpost"

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

calppost <- function(y,yrep,X,p_b,p_m,p_alphas){
  
  # using nsim, nsim1, nd, X from the global environment
  
  Dobs <- Drep <- matrix(0,nsim,nd)
  
  ratio_obs <- outer(nsim1,y)/p_m
  ratio_rep <- yrep/p_m
  log_diff_obs <- log(ratio_obs); log_diff_rep <- log(ratio_rep)
  prop_diff_obs <- ratio_obs-1; prop_diff_rep <- ratio_rep-1
  
  # scaled deviance
  tmp_obs <- rowSums(prop_diff_obs - log_diff_obs)
  tmp_rep <- rowSums(prop_diff_rep - log_diff_rep)
  
  Dobs[,1] <- tmp_obs * p_alphas
  Drep[,1] <- tmp_rep * p_alphas
  
  # chi squared test stats
  Dobs[,2] <- rowSums(prop_diff_obs^2) * p_alphas
  Drep[,2] <- rowSums(prop_diff_rep^2) * p_alphas
  
  # in model score
  Dobs3_beta <- rowSums(prop_diff_obs %*% X) * p_alphas
  Dobs[,3] <- -tmp_obs + n*(log(p_alphas)-digamma(p_alphas)) + Dobs3_beta
  
  Drep3_beta <- rowSums(prop_diff_rep %*% X) * p_alphas
  Drep[,3] <- -tmp_rep + n*(log(p_alphas)-digamma(p_alphas)) + Drep3_beta
  
  # out of model score
  Dobs[,4] <- rowSums(prop_diff_obs * outer(nsim1, z)) * p_alphas
  Drep[,4] <- rowSums(prop_diff_rep * outer(nsim1, z)) * p_alphas
  
  # KS
  cdf_obs <- matrix(pgamma(rep(y,each=nsim),rep(p_alphas,n),c(p_b)),nrow=nsim)
  cdf_rep <- matrix(pgamma(c(yrep),rep(p_alphas,n),c(p_b)), nrow=nsim)
  
  Dobs[,5] <- sapply(c(1:nsim), function (i) 
    ks.test(cdf_obs[i,], "punif")$statistic )
  Drep[,5] <- sapply(c(1:nsim), function (i) 
    ks.test(cdf_rep[i,], "punif")$statistic )
  
  return(pr_ge(Drep,Dobs))
  
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
    p_m <- as.matrix(stan_mod, pars="m")
    p_b <- sweep(1/p_m, 1, p_alphas, "*")
    set.seed(slurm_id+1)
    yrep <- matrix(rgamma(nsim*n,p_alphas, c(p_b)),nrow=nsim)
    
    ppost <- calppost(y,yrep,X,p_b,p_m,p_alphas)
    
    ppost_rep <- foreach(i=1:nsim, .combine="rbind") %dopar% {
      data$y <- yrep[i,]
      
      stan_mod_i <- sampling(mod, data = data,
                             chains = nchains, warmup = burnin,iter = (burnin+thin*nsim/nchains),
                             cores = nchains, thin = thin, seed = slurm_id, refresh=0)
      p_alphas_i <- as.matrix(stan_mod_i, pars="alpha")
      p_m_i <- as.matrix(stan_mod_i, pars="m")
      p_b_i <- sweep(1/p_m_i, 1, p_alphas_i, "*")
      
      set.seed(slurm_id+1)
      yrep_rep <- matrix(rgamma(nsim*n,p_alphas_i, c(p_b_i)),nrow=nsim)
      
      ppost_rep_i <- calppost(yrep[i,],yrep_rep,X,p_b_i,p_m_i,p_alphas_i)

      return(ppost_rep_i)
    }
    
    output[q,] <- pr_ge(-ppost_rep,outer(nsim1, -ppost))

  }
})

res <- list(runtime=runtime,output=output)

dir.create(paste0("../res/", folder))
output_name <- paste0("../res/",folder,"/sim", simnum, "_", p_label,"_",slurm_id, ".rds")
saveRDS(res, output_name)
