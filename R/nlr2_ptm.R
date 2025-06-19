slurm_id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

source("functions.R")

simnum <- 2
folder <- "nlr"

# p-value to be calculated
p_label <- "ptm"

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
theta <- c(rnorm(p,1,1.5),3)

z <- rgamma(n,6,2)
Vz <- cbind(V,z)

# simulate observed data
set.seed(slurm_id)
# y <- rnorm(n, V %*% theta, sigma)

# power testing
y <- rnorm(n, Vz %*% theta, sigma)

# prior parameter
prior_labels <- "good"
nprior <- length(prior_labels)

mus <- matrix(rep(1,p),nrow=1)
Psis <- array(0,dim=c(nprior, p,p))
Psis[1,,] <- diag(p)*4

runtime <- system.time({

  MLE <- c(inv_VTV %*% crossprod(V,y))
  
  output <- matrix(0,nprior,nd,dimnames = list(prior_labels,d_labels))
  
  ptm_obs <- numeric(nd)
  ptm_rep <- matrix(0,nsim,nd)
  for (q in 1:nprior){
    
    p_thetas <- nlr_MCMC(V,y,sigma,mus[q,],Psis[q,,],nsim,seed=slurm_id)
    pred_m <- tcrossprod(p_thetas,V)
    
    set.seed(slurm_id+1)
    yrep <- matrix(rnorm(nsim*n, c(pred_m),sigma), nrow=nsim)
    MLE_rep <- apply(yrep,1, function(x) c(inv_VTV %*% crossprod(V,x)))
    
    if (Psis[q,1,1]==Inf){
      inv_Psi <- matrix(0,p,p)
    } else {
      inv_Psi <- solve(Psis[q,,])
    }
    
    S <- solve(VTV/sigma^2 + inv_Psi)
    m_yobs <- S %*% (crossprod(V, y)/sigma^2 + inv_Psi %*% mus[q,])
    m_yrep <- sweep(yrep %*% V/sigma^2,2,c(inv_Psi %*% mus[q,]),"+") %*% S
    
    Vm_yobs <- V %*% m_yobs
    Vm_yrep <- tcrossprod(m_yrep,V)
    
    ptm_obs[1] <- sum((y-Vm_yobs)^2)/sigma^2
    ptm_obs[2] <- ks.test(pnorm(y,Vm_yobs,sigma), "punif")$statistic
    ptm_obs[3] <- t(MLE-m_yobs) %*% VTV %*% (MLE-m_yobs) / sigma^2
    
    ptm_rep[,1] <- rowSums((yrep - Vm_yrep)^2)/sigma^2
    
    cdf_rep <- matrix(pnorm(c(yrep),c(Vm_yrep),sigma), nrow=nsim)
    ptm_rep[,2] <- sapply(c(1:nsim), function (i) ks.test(cdf_rep[i,], "punif")$statistic )
    ptm_rep[,3] <- apply(t(MLE_rep)-m_yrep,1,function(x) t(x) %*% VTV %*% x / sigma^2)

    output[q,] <- pr_ge(ptm_rep,outer(nsim1, ptm_obs))
  }
})

res <- list(runtime=runtime,output=output)

saveRDS(res, output_name)
