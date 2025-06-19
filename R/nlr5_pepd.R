slurm_id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
n_cpus <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
library(doMC)
source("functions.R")

registerDoMC(cores = n_cpus)

simnum <- 5
folder <- "nlr"

# p-value to be calculated
p_label <- "pepd"

dir.create(paste0(outfolder, folder))
output_name <- paste0(outfolder,folder,"/sim", simnum, "_", p_label,"_",slurm_id, ".rds")

# fixed parameters
nsim <- 1000 # number of MCMC posterior samples/ posterior predictive rep datasets
nsim1 <- rep(1,nsim)

# number of test statistics/ discrepancy measures
d_labels <- c("chi-squared","KS","sufficient")
nd <- length(d_labels)

# general parameters
n <- 50
p <- 6
sigma <- 1

# true parameters
set.seed(1)
V <- matrix(rnorm(n*p), nrow=n)
VTV <- crossprod(V)
inv_VTV <- solve(VTV)
set.seed(2)
theta <- rnorm(p,1,1.5)

# simulate observed data
set.seed(slurm_id)
# y <- rnorm(n, V %*% theta, sigma)

# power testing, lognormal
y <- c(V %*% theta) + exp(rnorm(n,0,sigma))

# prior parameter
prior_labels <- "good"
nprior <- length(prior_labels)

mus <- matrix(rep(1,p),nrow=1)
Psis <- array(0,dim=c(nprior, p,p))
Psis[1,,] <- diag(p)*4

calpepd <- function(y,pred_m,p_thetas,MLE){
  
  Dobs <- matrix(0,nsim,nd)
  
  # chi-squared test
  Dobs[,1] <- rowSums((outer(nsim1,y) - pred_m)^2)/sigma^2

  # KS test
  cdf_obs <- matrix(pnorm(rep(y,each=nsim),c(pred_m),sigma),nrow=nsim)
  Dobs[,2] <- sapply(c(1:nsim), function (i) ks.test(cdf_obs[i,], "punif")$statistic )

  Dobs[,3] <- apply(outer(nsim1,MLE)-p_thetas,1,function(x) t(x) %*% VTV %*% x / sigma^2)

  pepd <- colMeans(Dobs)
  
  return(pepd)

}

runtime <- system.time({
  
  MLE <- c(inv_VTV %*% crossprod(V,y))
  
  output <- matrix(0,nprior,nd,dimnames = list(prior_labels,d_labels))

  for (q in 1:nprior){
    
    p_thetas <- nlr_MCMC(V,y,sigma,mus[q,],Psis[q,,],nsim,seed=slurm_id)
  
    pred_m <- tcrossprod(p_thetas,V)
    
    set.seed(slurm_id+1)
    yrep <- matrix(rnorm(nsim*n, c(pred_m),sigma), nrow=nsim)
    MLE_rep <- apply(yrep,1, function(x) c(inv_VTV %*% crossprod(V,x)))
    
    pepd_obs <- calpepd(y,pred_m,p_thetas,MLE)
    
    pepd_rep <- foreach(i=1:nsim, .combine="rbind") %dopar% {
      p_thetas_i <- nlr_MCMC(V,yrep[i,],sigma,mus[q,], Psis[q,,],nsim,seed=slurm_id)
      pred_m_i <- p_thetas_i %*% t(V)
      
      pepd_rep_i <- calpepd(yrep[i,],pred_m_i,p_thetas_i,MLE_rep[,i])
      return(pepd_rep_i)
    }
    
    output[q,] <- pr_ge(pepd_rep,outer(nsim1, pepd_obs))
  }
})

res <- list(runtime=runtime,output=output)

saveRDS(res, output_name)
