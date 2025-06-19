testing <- F
overwrite <- T
simnum <- 1
simlabel <- "sm"
p_label <- "psplit" # p-value to be calculated

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

if (!file.exists(output_name) | overwrite){

  # survival data true parameters
  n <- 100; m<-2; p<-3
  gammas <- c(1,1) # scale
  alphas <- c(1/0.3,1/0.5) # shape
  # beta, m x (p+1) if no censoring, (m+1) x (p+1) if with censoring
  betas <- matrix(c(-10,-0.6,0.4,1,
                    -6.5,0.3,-0.5,-0.8), nrow=m, byrow=T)
  
  ### fixed parameters
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
  # k=30 was tested to be better, smoother curve compared to k=10
  # kappa1=10 was tested to be a pretty good choice, therefore the parameters for IG
  mus <- c(0,10); taus <- c(10,0.1)
  mod <- cmdstan_model("crsm_mgp7.stan")
  
  calpsplit <- function(){
    
    Dobs <- Drep <- matrix(0,nsim,nd)
    # by our setup, number of events in the hold-out set will always be either same as the events in the training set, or smaller by 1
    # for consistency, we will use this smaller number of events for p-value calculation
    nho <- sum(ise_ho)
    
    # drop one sample if needed so that training and hold-out set have the same number of event
    dtime_rep_spt <- data_rep_spt$res$time_rep_e[,1:nho] # nsim x nho
    X_tr_e <- data_spt_tr$X[ise_tr,][1:nho,] # nho x p
    delta_rep_spt <- data_rep_spt$res$delta_rep_e[,1:nho,] # nsim x nho x m
    
    Ls_ho_e <- vapply(time_ho_e,findInterval,s, FUN.VALUE=integer(1)) # "events" of the "hold-out" set, hence ho_e, length = nho
    Ls_rep_spt <- t(vapply(c(1:nsim), function(i) findInterval(dtime_rep_spt[i,],s), integer(nho))) # nsim x nho
    
    s_diff <- tail(s,k)-head(s,k)
    
    dim(p_lambda_spt) <- c(nsim,k,m)
    p_lambda_spt_ext <- array(0,dim=c(nsim,k+1,m))
    p_lambda_spt_ext[,1:k,] <- p_lambda_spt
    p_lambda_spt_ext[,k+1,] <- p_lambda_spt_ext[,k,]  
    p_lambda_spt_ext <- aperm(p_lambda_spt_ext,c(2,3,1))
    cumLambda_spt_ext <- array(0,dim=c(k+1,m,nsim))
    cumLambda_spt_ext[-1,,] <- vapply(c(1:nsim), function(i) colCumsums(p_lambda_spt[i,,]*s_diff), matrix(0,k,m))
    cumLambda_spt_rep <- vapply(c(1:nsim), function(i) cumLambda_spt_ext[Ls_rep_spt[i,],,i], matrix(0,nho,m))
    
    dim(p_beta_spt) <- c(nsim,p,m)
    
    eXbeta_ho_e <- exp(vapply(c(1:nsim), function(i) X_ho_e %*% p_beta_spt[i,,], matrix(0,nho,m))) # nho x m x nsim
    eXbeta_tr_e <- exp(vapply(c(1:nsim), function(i) X_tr_e %*% p_beta_spt[i,,], matrix(0,nho,m))) # nho x m x nsim
    
    p_lambda_ho_e <- p_lambda_spt_ext[Ls_ho_e,,]
    tmp <- p_lambda_ho_e * eXbeta_ho_e
    probs_ho_e <- vapply(c(1:nsim), function(i) tmp[,,i]/rowSums(tmp[,,i]), matrix(0,nho,m)) # nho x m x nsim
    
    p_lambda_rep_spt <- vapply(c(1:nsim), function(i) p_lambda_spt_ext[Ls_rep_spt[i,],,i], matrix(0,nho,m)) # nho x m x nsim
    tmp <- p_lambda_rep_spt * eXbeta_tr_e
    probs_rep_spt <- vapply(c(1:nsim), function(i) tmp[,,i]/rowSums(tmp[,,i]), matrix(0,nho,m)) # nho x m x nsim
    
    # chi-squared test
    Dobs[,1] <- vapply(c(1:nsim), function(i) vec_get_chisqr(delta_ho_e, probs_ho_e[,,i]), numeric(1))
    Drep[,1] <- vapply(c(1:nsim), function(i) vec_get_chisqr(delta_rep_spt[i,,], probs_rep_spt[,,i]), numeric(1))
    
    # KS
    Ts_ho_e <- time_ho_e - s[Ls_ho_e]
    H0_ho_e <- cumLambda_spt_ext[Ls_ho_e,,] + vapply(c(1:nsim), function(i) Ts_ho_e*p_lambda_ho_e[,,i],matrix(0,nho,m))
    tmp <- H0_ho_e * eXbeta_ho_e
    H_ho_e <- vapply(c(1:nsim), function(i) rowSums(tmp[,,i]),numeric(nho)) # nho x nsim
    
    Ts_rep_spt <- dtime_rep_spt - array(s[Ls_rep_spt],dim=dim(Ls_rep_spt))
    H0_rep_spt <- cumLambda_spt_rep + vapply(c(1:nsim), function(i) Ts_rep_spt[i,]*p_lambda_rep_spt[,,i],matrix(0,nho,m))
    tmp <- H0_rep_spt * eXbeta_tr_e
    H_rep_spt <- vapply(c(1:nsim), function(i) rowSums(tmp[,,i]),numeric(nho)) # nho x nsim
    
    Dobs[,2] <- vapply(c(1:nsim), function (i) 
      ks.test(H_ho_e[,i], "pexp")$statistic, numeric(1) )
    Drep[,2] <- vapply(c(1:nsim), function (i) 
      ks.test(H_rep_spt[,i], "pexp")$statistic, numeric(1) )
    
    colMeans(Drep >= Dobs)
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
  
  max_time <- ceiling(max(data$data$time))
  
  # for psplit
  # training vs hold-out sets
  
  # try to evenly split based on number of events
  # if total number of events is odd, have 1 more event in the training set, later drop it for comparability
  idx_tr_e <- base::sample(c(1:n)[ise],ceiling(sum(ise)/2))
  idx_tr_c <- base::sample(c(1:n)[!ise],n/2-ceiling(sum(ise)/2))
  idx_tr <- rep(F,n)
  idx_tr[c(idx_tr_e,idx_tr_c)] <- T
  
  # training data, with censoring
  # created data_spt_tr for later calling simdata6 function
  data_spt_tr <- data
  data_spt_tr$data <- data_spt_tr$data[idx_tr,]
  data_spt_tr$X <- data_spt_tr$X[idx_tr,]
  data_spt_tr$delta <- data_spt_tr$delta[idx_tr,]
  data_spt_tr$n <- n/2
  
  # hold out, we only need the events for chi-squared test
  ise_ho <- data$data$event[!idx_tr] != 0
  ise_tr <- data$data$event[idx_tr] != 0
  time_ho_e <- data$data$time[!idx_tr][ise_ho]
  X_ho_e <- data$X[!idx_tr,][ise_ho,]
  delta_ho_e <- data$delta[!idx_tr,][ise_ho,]
  
  runtime <- system.time({
    
    for (q in 1:nprior){

      s <- seq(0,max_time,length.out=k+1)
      
      dataset_spt <- list(n=data_spt_tr$n,p=data_spt_tr$p,k=k,m=data_spt_tr$m,y=data_spt_tr$data$time, s=s,
                          delta=data_spt_tr$delta, X=data_spt_tr$X, kappa0=0.1, a=a,b=b, mu0=mus[q],
                          tau0=taus[q])
      cmdrun_spt <- mod$sample(data=dataset_spt,chains=nchains,parallel_chains = nchains,
                               iter_warmup=burnin,iter_sampling=niter,
                               thin=thin,seed=slurm_id, refresh=0, show_messages = F, show_exceptions = F)
      
      p_beta_spt <- matrix(as_draws(cmdrun_spt,variable="beta"),nrow=nsim)
      p_lambda_spt <- matrix(as_draws(cmdrun_spt,variable="lambda"),nrow=nsim) # this is faster than apply(,identity)
      
      # get replicated datasets
      data_rep_spt <- to_sim_ppdata2(data_spt_tr,s,p_beta_spt,p_lambda_spt,seed=slurm_id)
      
      output[q,] <- calpsplit()
  
    }
  
  })
  
  res <- list(runtime=runtime,output=output)
  
  dir.create(paste0(outfolder, simlabel))
  saveRDS(res, output_name)

}