testing <- F
override <- T
simnum <- 3
simlabel <- "sm"
p_label <- "pepd" # p-value to be calculated

if (testing){
  # for testing
  burnin <- 200
  nsim <- 24 # number of posterior samples, also number of predictive datasets
  slurm_id <- 1; n_cpus <- 6
} else {
  # real run
  burnin <- 1000
  nsim <- 1000
  slurm_id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
  n_cpus <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))-2
}

thin <- 5
nchains <- 1 # sometimes server run has problems collecting all the results from multiple chains, therefore only using one chain
niter <- nsim*thin/nchains + burnin

library(rstan)
source("functions.R")
registerDoMC(cores = n_cpus)

output_name <- paste0(outfolder,simlabel,"/sim", simnum, "_", p_label,"_",slurm_id, ".rds")
idxlist <- readRDS(paste0("../res/",simlabel,simnum,"_idxlist.rds"))

if (!file.exists(output_name) | override){
  
# survival data true parameters
n <- 100; m<-2; p<-3
gammas <- c(1,1) # scale
alphas <- c(1/0.3,1/0.5) # shape
k <- 30; a <- 2.5; b <- 10

# to test for betas that change with t
betas0 <- matrix(c(-10,-1,-0.2,0.3,
                   -6.5,-0.4,-1.1,0.2), nrow=m, byrow=T)
betas1 <- matrix(c(0,0.016,0.024,0.028,
                   0,0.028,0.024,0.016), nrow=m, byrow=T)

### fixed parameters

# for p-values
one_nsim <- rep(1,nsim)
one_n <- rep(1,n)
zero_m <- rep(0,m)

# number of test statistics/ discrepancy measures
d_labels <- c("chi-squared","KS")
nd <- length(d_labels)

# number of priors
prior_labels <- c("good")
mus <- c(0); taus <- c(10)
nprior <- length(prior_labels)
mod <- readRDS("crsm_mgp7_mod.rds")

calpepd <- function(survdata,delta,s,p_beta,p_lambda,X){
  
  # using the following variables in global environment: n, nsim, nd, k
  n <- length(survdata$time)
  Dobs <- matrix(0,nsim,nd)

  Ls_obs <- vapply(survdata$time,findInterval,s, FUN.VALUE=integer(1))
  
  s_diff <- tail(s,k)-head(s,k)
  
  # posterior mean
  p_lambda_m <- matrix(colMeans(p_lambda),ncol=m) # k x m
  p_beta_m <- matrix(colMeans(p_beta),ncol=m) # p x m
  
  dim(p_lambda) <- c(nsim,k,m)
  p_lambda_ext <- array(0,dim=c(nsim,k+1,m))
  p_lambda_ext[,1:k,] <- p_lambda
  p_lambda_ext[,k+1,] <- p_lambda_ext[,k,]  
  p_lambda_ext <- aperm(p_lambda_ext,c(2,3,1))
  
  cumLambda_ext <- array(0,dim=c(k+1,m,nsim))
  cumLambda_ext[-1,,] <- vapply(c(1:nsim), function(i) 
    colCumsums(p_lambda[i,,]*s_diff), matrix(0,k,m))
  cumLambda_obs <- cumLambda_ext[Ls_obs,,]
  
  dim(p_beta) <- c(nsim,p,m)
  eXbeta <- exp(vapply(c(1:nsim), function(i) X %*% p_beta[i,,], matrix(0,n,m)))
  
  p_lambda_obs <- p_lambda_ext[Ls_obs,,]
  tmp <- p_lambda_obs * eXbeta
  probs_obs <- vapply(c(1:nsim), function(i) tmp[,,i]/rowSums(tmp[,,i]), matrix(0,n,m))
  
  # chi-squared test
  Dobs[,1] <- vapply(c(1:nsim), function(i) vec_get_chisqr(delta, probs_obs[,,i]), numeric(1))
  
  # KS
  Ts_obs <- survdata$time - s[Ls_obs]
  H0_obs <- cumLambda_obs + vapply(c(1:nsim), function(i) Ts_obs*p_lambda_obs[,,i],matrix(0,n,m))
  
  # the following two lines are faster than directly using apply
  tmp <- H0_obs * eXbeta
  H_obs <- vapply(c(1:nsim), function(i) rowSums(tmp[,,i]),numeric(n))
  
  Dobs[,2] <- vapply(c(1:nsim), function (i) 
    ks.test(H_obs[,i], "pexp")$statistic, numeric(1) )
  
  return(colMeans(Dobs))
  
}

output <- matrix(0,nprior,nd,dimnames = list(prior_labels,d_labels))

# fix design matrix for all simulations
set.seed(1)
x1 <- rbinom(n,1,0.5)
x2 <- rnorm(n)
x3 <- runif(n,0,2)
X <- cbind(one_n, x1,x2,x3)

data <- simdata7(alphas,gammas,X,betas0,betas1,seed=idxlist[slurm_id])
ise <- data$data$event != 0 # is event, not censored
Xe <- data$X[ise,]

max_time <- ceiling(max(data$data$time))

runtime <- system.time({
  
  for (q in 1:nprior){
    s <- seq(0,max_time,length.out=k+1)
    
    dataset <- list(n=data$n,p=data$p,k=k,m=data$m,y=data$data$time, s=s,
                    delta=data$delta, X=data$X, kappa0=0.1, a=a,b=b, mu0=mus[q],
                    tau0=taus[q])
    runres <- sampling(mod,data=dataset,chains=nchains,warmup=burnin,iter=niter,
                       cores=nchains,thin=thin,seed=slurm_id,refresh=0,
                       include=T,pars=c("beta","lambda"))
    p_beta <- as.matrix(runres,par="beta")
    p_lambda <- as.matrix(runres,par="lambda")
    
    # get replicated datasets
    data_rep <- to_sim_ppdata2(data,s,p_beta,p_lambda,seed=slurm_id)
    
    pepd_obs <- calpepd(data$data[ise,],data$delta[ise,],s,p_beta,p_lambda,Xe)
    survdata_i <- data$data[ise,]
    
    pepd_rep <- foreach(i=1:nsim, .combine="rbind") %dopar% {
      
      max_time_i <- ceiling(max(data_rep$res$time_rep[i,]))
      s_i <- seq(0,max_time_i,length.out=k+1)
      dataset$s <- s_i
      dataset$y <- data_rep$res$time_rep[i,]
      dataset$delta <- data_rep$res$delta_rep[i,,]
      
      runres <- sampling(mod,data=dataset,chains=nchains,warmup=burnin,iter=niter,
                         cores=nchains,thin=thin,seed=slurm_id,refresh=0,
                         include=T,pars=c("beta","lambda"))
      p_beta_i <- as.matrix(runres,par="beta")
      p_lambda_i <- as.matrix(runres,par="lambda")
      survdata_i$time <- data_rep$res$time_rep_e[i,]; survdata_i$event <- data_rep$res$event_rep_e[i,]
      
      pepd_i <- calpepd(survdata_i,data_rep$res$delta_rep_e[i,,],s_i,p_beta_i,p_lambda_i,Xe)
      
      return(pepd_i)
    }
    
    output[q,] <- pr_ge(pepd_rep,outer(one_nsim, pepd_obs))
    
  }

})

res <- list(runtime=runtime,output=output)

dir.create(paste0(outfolder, simlabel))
saveRDS(res, output_name)

}
