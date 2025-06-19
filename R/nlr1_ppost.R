slurm_id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

source("functions.R")


simnum <- 1
folder <- "nlr"

# p-value to be calculated
p_label <- "ppost"

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
y <- rnorm(n, V %*% theta, sigma)

# prior parameter
prior_labels <- c("flat", "good", "bad", "unit info")
nprior <- length(prior_labels)

mus <- rbind(rep(0,p),
             rep(1,p),
             rep(-3,p),
             rep(0,p)) # this is placeholder for the unit information prior
Psis <- array(0,dim=c(nprior, p,p))

Psis[1,,] <- Inf
Psis[2,,] <- diag(p)*4
Psis[3,,] <- diag(p)/4


calppost <- function(y,yrep,pred_m,p_thetas,MLE,MLE_rep){
  
  Dobs <- Drep <- matrix(0,nsim,nd)
  
  # chi-squared test
  Dobs[,1] <- rowSums((outer(nsim1,y) - pred_m)^2)/sigma^2
  Drep[,1] <- rowSums((yrep - pred_m)^2)/sigma^2
  
  # KS test
  cdf_obs <- matrix(pnorm(rep(y,each=nsim),c(pred_m),sigma),nrow=nsim)
  cdf_rep <- matrix(pnorm(c(yrep),c(pred_m),sigma), nrow=nsim)
  Dobs[,2] <- sapply(c(1:nsim), function (i) ks.test(cdf_obs[i,], "punif")$statistic )
  Drep[,2] <- sapply(c(1:nsim), function (i) ks.test(cdf_rep[i,], "punif")$statistic )
  
  Dobs[,3] <- apply(outer(nsim1,MLE)-p_thetas,1,function(x) t(x) %*% VTV %*% x / sigma^2)
  Drep[,3] <- apply(t(MLE_rep)-p_thetas,1,function(x) t(x) %*% VTV %*% x / sigma^2)
  
  ppost <- pr_ge(Drep,Dobs)
  
  return(ppost)

}

runtime <- system.time({
  
  # for unit information prior
  mus[nprior,] <- inv_VTV %*% crossprod(V, y)
  Psis[nprior,,] <- n*sigma^2*inv_VTV
  MLE <- c(inv_VTV %*% crossprod(V,y))
  
  output <- matrix(0,nprior,nd,dimnames = list(prior_labels,d_labels))
  
  for (q in 1:nprior){
    
    p_thetas <- nlr_MCMC(V,y,sigma,mus[q,],Psis[q,,],nsim,seed=slurm_id)
  
    pred_m <- tcrossprod(p_thetas,V)
    
    set.seed(slurm_id+1)
    yrep <- matrix(rnorm(nsim*n, c(pred_m),sigma), nrow=nsim)
    MLE_rep <- apply(yrep,1, function(x) c(inv_VTV %*% crossprod(V,x)))
    
    output[q,] <- calppost(y,yrep,pred_m,p_thetas,MLE,MLE_rep)
  }
})

res <- list(runtime=runtime,output=output)
saveRDS(res, output_name)
