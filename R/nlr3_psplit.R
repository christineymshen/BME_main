slurm_id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
source("functions.R")

simnum <- 3
folder <- "nlr"

# p-value to be calculated
p_label <- "psplit"

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

set.seed(2)
theta <- rnorm(p,1,1.5)

# simulate observed data
set.seed(slurm_id)

# power testing
tau <- 0.5
y <- rnorm(n, V %*% theta, tau)

Vobs <- V[1:(n/2),]; yobs <- y[1:(n/2)]
Vnew <- V[(n/2+1):n,]; ynew <- y[(n/2+1):n]
VTVnew <- crossprod(Vnew); VTVobs <- crossprod(Vobs)
inv_VTVnew <- solve(VTVnew); inv_VTVobs <- solve(VTVobs)

# prior parameter
prior_labels <- "good"
nprior <- length(prior_labels)

mus <- matrix(rep(1,p),nrow=1)
Psis <- array(0,dim=c(nprior, p,p))
Psis[1,,] <- diag(p)*4

calpsplit <- function(){
  # decided to keep using "_spt" in variable names for ease of debugging
  
  Dobs <- Drep <- matrix(0,nsim,nd)
  
  # chi-squared test
  Dobs[,1] <- rowSums((outer(nsim1,ynew) - pred_new_m)^2)/sigma^2
  Drep[,1] <- rowSums((yrep_spt - pred_spt_m)^2)/sigma^2
  
  # KS test
  cdf_obs <- matrix(pnorm(rep(ynew,each=nsim),c(pred_new_m),sigma),nrow=nsim)
  cdf_rep <- matrix(pnorm(c(yrep_spt),c(pred_spt_m),sigma), nrow=nsim)
  Dobs[,2] <- sapply(c(1:nsim), function (i) ks.test(cdf_obs[i,], "punif")$statistic )
  Drep[,2] <- sapply(c(1:nsim), function (i) ks.test(cdf_rep[i,], "punif")$statistic )
  
  Dobs[,3] <- apply(outer(nsim1,MLE_new_spt)-p_spt_thetas,1,function(x) t(x) %*% VTVnew %*% x / sigma^2)
  Drep[,3] <- apply(t(MLE_rep_spt)-p_spt_thetas,1,function(x) t(x) %*% VTVobs %*% x / sigma^2)
  
  pr_ge(Drep,Dobs)
}

runtime <- system.time({

  MLE <- c(inv_VTVobs %*% crossprod(Vobs,yobs))
  MLE_new_spt <- c(inv_VTVnew %*% crossprod(Vnew,ynew))
  
  output <- matrix(0,nprior,nd,dimnames = list(prior_labels,d_labels))
  
  for (q in 1:nprior){
    
    p_spt_thetas <- nlr_MCMC(Vobs,yobs,sigma,mus[q,], Psis[q,,],nsim,seed=slurm_id)
    pred_spt_m <- tcrossprod(p_spt_thetas,Vobs)
    pred_new_m <- tcrossprod(p_spt_thetas,Vnew)
    
    set.seed(slurm_id+1)
    yrep_spt <- matrix(rnorm(nsim*n/2, c(pred_spt_m),sigma), nrow=nsim)
    MLE_rep_spt <- apply(yrep_spt,1, function(x) c(inv_VTVobs %*% crossprod(Vobs,x)))
    
    output[q,] <- calpsplit()
  }
})

res <- list(runtime=runtime,output=output)

saveRDS(res, output_name)
