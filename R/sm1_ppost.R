testing <- F
overwrite <- T
simnum <- 1
simlabel <- "sm"
p_label <- "ppost" # p-value to be calculated

if (testing){
  # for testing
  burnin <- 200
  nsim <- 20 # number of posterior samples, also number of predictive datasets
  slurm_id <- 1
} else {
  # real run
  burnin <- 1000
  nsim <- 1000 # number of posterior samples, also number of predictive datasets
  slurm_id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
}

thin <- 5
nchains <- 1 # sometimes server run has problems collecting all the results from multiple chains, therefore only using one chain
niter <- nsim*thin/nchains

library(cmdstanr)
source("functions.R")

output_name <- paste0(outfolder,simlabel,"/sim",simnum,"_",p_label,"_",slurm_id,".rds")

if (!file.exists(output_name)){
  
  # survival data true parameters
  n <- 100; m<-2; p<-3
  gammas <- c(1,1) # scale
  alphas <- c(1/0.3,1/0.5) # shape
  # beta, m x (p+1) if no censoring, (m+1) x (p+1) if with censoring
  betas <- matrix(c(-10,-0.6,0.4,1,
                    -6.5,0.3,-0.5,-0.8), nrow=m, byrow=T)
  
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
  mod <- cmdstan_model("crsm_mgp7.stan")
  
  calppost <- function(survdata,delta,dtime_rep,delta_rep,s,p_beta,p_lambda,X){
    
    # using the following variables in global environment: n, nsim, nd, k
    n <- length(survdata$time)
    Dobs <- Drep <- matrix(0,nsim,nd)
  
    Ls_obs <- vapply(survdata$time,findInterval,s, FUN.VALUE=integer(1))
    Ls_rep <- t(vapply(c(1:nsim), function(i) findInterval(dtime_rep[i,],s), integer(n)))
    
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
    cumLambda_rep <- vapply(c(1:nsim), function(i) cumLambda_ext[Ls_rep[i,],,i], matrix(0,n,m))
    
    dim(p_beta) <- c(nsim,p,m)
    eXbeta <- exp(vapply(c(1:nsim), function(i) X %*% p_beta[i,,], matrix(0,n,m)))
    
    p_lambda_obs <- p_lambda_ext[Ls_obs,,]
    tmp <- p_lambda_obs * eXbeta
    probs_obs <- vapply(c(1:nsim), function(i) tmp[,,i]/rowSums(tmp[,,i]), matrix(0,n,m))
    
    p_lambda_rep <- vapply(c(1:nsim), function(i) p_lambda_ext[Ls_rep[i,],,i], matrix(0,n,m))
    tmp <- p_lambda_rep * eXbeta
    probs_rep <- vapply(c(1:nsim), function(i) tmp[,,i]/rowSums(tmp[,,i]), matrix(0,n,m))
    
    # chi-squared test
    Dobs[,1] <- vapply(c(1:nsim), function(i) vec_get_chisqr(delta, probs_obs[,,i]), numeric(1))
    Drep[,1] <- vapply(c(1:nsim), function(i) vec_get_chisqr(delta_rep[i,,], probs_rep[,,i]), numeric(1))
    
    # KS
    Ts_obs <- survdata$time - s[Ls_obs]
    H0_obs <- cumLambda_obs + vapply(c(1:nsim), function(i) Ts_obs*p_lambda_obs[,,i],matrix(0,n,m))
    
    # the following two lines are faster than directly using apply
    tmp <- H0_obs * eXbeta
    H_obs <- vapply(c(1:nsim), function(i) rowSums(tmp[,,i]),numeric(n))
  
    Ts_rep <- dtime_rep - array(s[Ls_rep],dim=dim(Ls_rep))
    H0_rep <- cumLambda_rep + vapply(c(1:nsim), function(i) Ts_rep[i,]*p_lambda_rep[,,i],matrix(0,n,m))
    tmp <- H0_rep * eXbeta
    H_rep <- vapply(c(1:nsim), function(i) rowSums(tmp[,,i]),numeric(n))
    
    Dobs[,2] <- vapply(c(1:nsim), function (i) 
      ks.test(H_obs[,i], "pexp")$statistic, numeric(1) )
    Drep[,2] <- vapply(c(1:nsim), function (i) 
      ks.test(H_rep[,i], "pexp")$statistic, numeric(1) )
    
    return(pr_ge(Drep,Dobs))
    
  }
  
  output <- matrix(0,nprior,nd,dimnames = list(prior_labels,d_labels))
  
  # fix design matrix for all simulations
  set.seed(1)
  x1 <- rbinom(n,1,0.5)
  x2 <- rnorm(n)
  x3 <- rnorm(n,1,1)
  
  X <- cbind(one_n, x1,x2,x3)
  
  data <- simdata6(alphas,gammas,X,betas,seed=slurm_id)
  ise <- data$data$event != 0 # is event, not censored
  Xe <- data$X[ise,]
  
  max_time <- ceiling(max(data$data$time))
  
  runtime <- system.time({
    
    for (q in 1:nprior){
      s <- seq(0,max_time,length.out=k+1)
      
      dataset <- list(n=data$n,p=data$p,k=k,m=data$m,y=data$data$time, s=s,
                      delta=data$delta, X=data$X, kappa0=0.1, a=a,b=b, mu0=mus[q],
                      tau0=taus[q])
      cmdrun <- mod$sample(data=dataset,chains=nchains,parallel_chains=nchains,
                           iter_warmup=burnin,iter_sampling=niter,thin=thin,init=0.5,
                           seed=slurm_id,refresh=0,show_messages=F,show_exceptions=F)
      p_beta <- matrix(as_draws(cmdrun,variable="beta"),nrow=nsim)
      p_lambda <- matrix(as_draws(cmdrun,variable="lambda"),nrow=nsim) # faster than apply(,identity)
      
      # get replicated datasets
      data_rep <- to_sim_ppdata2(data,s,p_beta,p_lambda,seed=slurm_id)
      
      output[q,] <- calppost(data$data[ise,],data$delta[ise,],data_rep$res$time_rep_e,
                             data_rep$res$delta_rep_e,s,p_beta,p_lambda,Xe)
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