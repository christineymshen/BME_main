testing <- F
override <- T
simnum <- 2
simlabel <- "sm"
p_label <- "ptm" # p-value to be calculated

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

if (!file.exists(output_name) | override){
  
# survival data true parameters
n <- 100; m<-2; p<-3
gammas <- c(1,1) # scale
alphas <- c(1/0.3,1/0.5) # shape
# beta, m x (p+1)
betas <- matrix(c(-10,-0.6,0.4,1,-1.5,0.5,1,
                  -6.5,0.3,-0.5,-0.8,0.7,-1.5,-1), nrow=m, byrow=T)

### fixed parameters

# for p-values
one_nsim <- rep(1,nsim)
one_n <- rep(1,n)
zero_m <- rep(0,m)
k <- 30; a <- 2.5; b <- 10

# number of test statistics/ discrepancy measures
nd <- 2
d_labels <- c("chi-squared","KS")

# number of priors
nprior <- 2
prior_labels <- c("good", "bad")
mus <- c(0,10); taus <- c(10,0.1)
mod <- readRDS("crsm_mgp7_mod.rds")

calptm <- function(survdata,delta,s,p_beta,p_lambda,X){
  
  # using the following variables in global environment: n, nsim, nd, k
  n <- length(survdata$time)
  Dobs <- numeric(nd)

  Ls_obs <- vapply(survdata$time,findInterval,s, FUN.VALUE=integer(1))
  
  s_diff <- tail(s,k)-head(s,k)
  
  # posterior mean
  p_lambda_m <- matrix(colMeans(p_lambda),ncol=m) # k x m
  p_beta_m <- matrix(colMeans(p_beta),ncol=m) # p x m
  
  Ts_obs <- survdata$time - s[Ls_obs]
  
  cumLambda_m <- colCumsums(p_lambda_m*s_diff); cumLambda_m <- rbind(zero_m,cumLambda_m)
  p_lambda_m <- rbind(p_lambda_m, p_lambda_m[k,])
  p_lambda_m_obs <- p_lambda_m[Ls_obs,]

  eXbeta_m <- exp(X %*% p_beta_m)
  tmp <- p_lambda_m_obs * eXbeta_m
  probs_m <- tmp / rowSums(tmp)
  Dobs[1] <- vec_get_chisqr(delta, probs_m)

  H0_m <- cumLambda_m[Ls_obs,] + p_lambda_m_obs * Ts_obs # checked no need manually broadcast  
  H_m <- rowSums(H0_m * eXbeta_m)
  
  Dobs[2] <- ks.test(H_m, "pexp")$statistic
  
  return(Dobs)
  
}

output <- matrix(0,nprior,nd,dimnames = list(prior_labels,d_labels))

# fix design matrix for all simulations
set.seed(1)
x1 <- rbinom(n,1,0.5)
x2 <- rnorm(n)
x3 <- rnorm(n,1,1)
x4 <- x2^2
x5 <- x3^3
x6 <- x2*x3
X <- cbind(one_n, x1,x2,x3,x4,x5,x6)
X_obs <- cbind(x1,x2,x3)

data <- simdata6(alphas,gammas,X,betas,seed=slurm_id)
data$X <- X_obs; data$p <- 3; 
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
    
    ptm_obs <- calptm(data$data[ise,],data$delta[ise,],s,p_beta,p_lambda,Xe)
    
    survdata_i <- data$data[ise,]
    ptm_rep <- foreach(i=1:nsim, .combine="rbind") %dopar% {
      
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
      
      ptm_i <- calptm(survdata_i,data_rep$res$delta_rep_e[i,,],s_i,p_beta_i,p_lambda_i,Xe)
      
      return(ptm_i)
    }
    
    output[q,] <- pr_ge(ptm_rep,outer(one_nsim, ptm_obs))
    
    # # for debugging
    # survdata <- data$data[ise,]
    # delta <- data$delta[ise,]
    # dtime_rep <- data_rep$res$time_rep_e
    # delta_rep <- data_rep$res$delta_rep_e
    # X <- data$X[ise,]
  }

})

res <- list(runtime=runtime,output=output)

dir.create(paste0(outfolder, simlabel))
saveRDS(res, output_name)

}