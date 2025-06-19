library(doMC)
library(Matrix) # create block diagonal matrix
library(tidyverse)
library(abind)
library(wesanderson)
library(bde)
library(amen)
library(matrixStats)

outfolder <- "/work/ys288/BME/"

### Survival Model ====

# this function works for n x m matrices
vec_get_chisqr <- function(delta,prob,tol=1e-8){
  
  # delta: n x m
  # prob: n x m
  # tol: threshold to decide whether the probability is 0
  n <- dim(delta)[1]
  idx <- prob>tol
  
  # removes everything that's too small, after that, remove the first element
  # in this way, if a prob vector is almost like (0,1,0), it won't be counted at all
  idx <- t(vapply(c(1:n), function(i) {
    idx[i,][idx[i,]][1] <- F
    return(idx[i,])}, logical(m) ))
  
  sum_p <- rowSums(prob*idx)
  r_p <- 1- sum_p # remainder
  dfs <- rowSums(idx) # degree of freedoms?
  
  vec_pi <- t(prob)[t(idx)] # ensure row major order
  vec_Vipi <- rep(1+sum_p/r_p, times=dfs)
  vec_ei <- t(delta)[t(idx)] # ensure row major order
  
  tmp <- (1/prob + 1/r_p)*delta
  tmp[!idx] <- 0
  eiViei <- sum(tmp)
  eiVipi <- sum(vec_ei * vec_Vipi)
  piVipi <- sum(vec_pi * vec_Vipi)
  
  (eiViei - 2*eiVipi + piVipi)/n/(m-1)
  
}

simdata6 <- function(alphas,gammas,X,betas,seed=1,pc=0.2){
  
  # change from simdata3:
  # I realized that if I were to condition on the censoring events, I better make sure I have enough other events in the data set for analysis
  # in v5, I don't have a good control on the proportion of censored data. For n=100, some of the simulated dataset have more than 60 censored data points. Then later when I simulate posterior predictive datasets, some of them essentially have too few events for one of the risk types that it's even not enough to estimate the regression coefficients
  # therefore in this version, I'm going to change the way I simulate data by separately simulating censored data
  
  # added this new variable pc, proportion of censored data points, default at 0.2
  
  # number of risk types
  m <- length(alphas)
  n <- dim(X)[1]; p <- dim(X)[2]
  nc <- floor(pc*n)
  
  set.seed(seed)
  Us <- matrix(runif(n*m),n,m)
  tmp <- sweep(-log(Us)/exp(tcrossprod(X,betas)),2,gammas,"/")
  T_latent <- sweep(tmp,2,1/alphas,"^")
  
  # get censored data
  idx_c <- sample(1:n, nc)
  
  ctime_upper <- rowMins(T_latent[idx_c,])
  ctime <- runif(nc,0,ctime_upper)
  
  time <- rowMins(T_latent)
  time[idx_c] <- ctime
  event <- apply(T_latent,1,which.min)
  event[idx_c] <- 0
  data <- data.frame(cbind(time,event))
  
  # remove intercepts
  for (i in 1:p){
    if (all(X[,i]==rep(1,n))){
      X <- X[,-i]
      p <- p-1
      betas_intercept <- betas[,i]
      betas <- betas[,-i]
      break
    }
  }

  data <- cbind(data,X)
  delta <- matrix(0,n,m)
  idx <- cbind(c(1:n), event)[event!=0,]
  delta[idx] <- 1
  
  list(n=n,p=p,m=m,X=X,data=data,delta=delta,
       betas_intercept=betas_intercept,
       betas=betas,alphas=alphas,gammas=gammas)
}

# adapted from simdata6 above to allow beta's to change linearly w.r.t. time
simdata7 <- function(alphas,gammas,X,betas0,betas1,seed=1,pc=0.2){

  # number of risk types
  m <- length(alphas)
  n <- dim(X)[1]; p <- dim(X)[2]
  nc <- floor(pc*n)
  
  set.seed(seed)
  Us <- matrix(runif(n*m),n,m)
  
  # solve equation ax+blog(x)=c
  cs <- log(-log(Us)) - tcrossprod(X,betas0) - outer(rep(1,n),log(gammas))
  as <- tcrossprod(X,betas1)
  
  f <- function(x,as,bs,cs) {
    as*x + bs*log(x) - cs
  }
  
  T_latent <- matrix(0,n,m)
  for (i in 1:n){
    for (j in 1:m){
      T_latent[i,j] <- uniroot(f,interval=c(1E-10,100),
                     as=as[i,j], bs=alphas[j], cs=cs[i,j])$root
    }
  }

  # get censored data
  idx_c <- sample(1:n, nc)
  
  ctime_upper <- rowMins(T_latent[idx_c,])
  ctime <- runif(nc,0,ctime_upper)
  
  time <- rowMins(T_latent)
  time[idx_c] <- ctime
  event <- apply(T_latent,1,which.min)
  event[idx_c] <- 0
  data <- data.frame(cbind(time,event))
  
  # remove intercepts
  for (i in 1:p){
    if (all(X[,i]==rep(1,n))){
      X <- X[,-i]
      p <- p-1
      break
    }
  }
  
  data <- cbind(data,X)
  delta <- matrix(0,n,m)
  idx <- cbind(c(1:n), event)[event!=0,]
  delta[idx] <- 1
  
  list(n=n,p=p,m=m,X=X,data=data,delta=delta)
}

# I realized that sometimes extreme data can be simulated. E.g., a dataset with very few events 
# of one risk type, even fewer than the number of coefficients to be estimated
# hence creating issues in the model fitting
# this is a wrapper function to keep simulating data until this tragedy doesn't happen

# v2 of this function also return the num_fail variable
to_sim_ppdata2 <- function(data,s,p_beta,p_lambda,seed=1){
  
  # count returns the number of iterations needed to get the replicated dataset
  # num_fail returns the number of replicated dataset that fails the check when first simulated
  
  res <- sim_ppdata(data,s,p_beta,p_lambda,seed)
  
  count_e <- vapply(c(1:nsim), function(i) colSums(res$delta_rep_e[i,,]), numeric(2))
  idx <- which(colSums(count_e <= data$p+1) >0)
  num_fail <- length(idx)
  
  count <- 0
  while(length(idx)!=0 & count<=1000){
    count <- count +1
    res2 <- sim_ppdata(data,s,p_beta[idx,,drop=F], p_lambda[idx,,drop=F], seed=seed+count)
    count_e <- vapply(c(1:length(idx)), function(i) colSums(res2$delta_rep_e[i,,]), numeric(2))
    idx2 <- which(colSums(count_e <= data$p+1) >0)
    
    if (length(idx2)==0){
      res$time_rep[idx,] <- res2$time_rep
      res$event_rep[idx,] <- res2$event_rep
      res$delta_rep[idx,,] <- res2$delta_rep
      res$time_rep_e[idx,] <- res2$time_rep_e
      res$event_rep_e[idx,] <- res2$event_rep_e
      res$delta_rep_e[idx,,] <- res2$delta_rep_e
    } else {
      res$time_rep[idx[-idx2],] <- res2$time_rep[-idx2,]
      res$event_rep[idx[-idx2],] <- res2$event_rep[-idx2,]
      res$delta_rep[idx[-idx2],,] <- res2$delta_rep[-idx2,,]
      res$time_rep_e[idx[-idx2],] <- res2$time_rep_e[-idx2,]
      res$event_rep_e[idx[-idx2],] <- res2$event_rep_e[-idx2,]
      res$delta_rep_e[idx[-idx2],,] <- res2$delta_rep_e[-idx2,,]
    }
    idx <- idx[idx2]
  }
 
  return(list(res=res,num_fail=num_fail,count=count))
}

to_sim_ppdata <- function(data,s,p_beta,p_lambda,seed=1){
  
  res <- sim_ppdata(data,s,p_beta,p_lambda,seed)
  
  count_e <- vapply(c(1:nsim), function(i) colSums(res$delta_rep_e[i,,]), numeric(2))
  idx <- which(colSums(count_e <= data$p+1) >0)
  
  count <- 0
  while(length(idx)!=0 & count<=1000){
    count <- count +1
    res2 <- sim_ppdata(data,s,p_beta[idx,,drop=F], p_lambda[idx,,drop=F], seed=seed+count)
    count_e <- vapply(c(1:length(idx)), function(i) colSums(res2$delta_rep_e[i,,]), numeric(2))
    idx2 <- which(colSums(count_e <= data$p+1) >0)
    
    if (length(idx2)==0){
      res$time_rep[idx,] <- res2$time_rep
      res$event_rep[idx,] <- res2$event_rep
      res$delta_rep[idx,,] <- res2$delta_rep
      res$time_rep_e[idx,] <- res2$time_rep_e
      res$event_rep_e[idx,] <- res2$event_rep_e
      res$delta_rep_e[idx,,] <- res2$delta_rep_e
    } else {
      res$time_rep[idx[-idx2],] <- res2$time_rep[-idx2,]
      res$event_rep[idx[-idx2],] <- res2$event_rep[-idx2,]
      res$delta_rep[idx[-idx2],,] <- res2$delta_rep[-idx2,,]
      res$time_rep_e[idx[-idx2],] <- res2$time_rep_e[-idx2,]
      res$event_rep_e[idx[-idx2],] <- res2$event_rep_e[-idx2,]
      res$delta_rep_e[idx[-idx2],,] <- res2$delta_rep_e[-idx2,,]
    }
    idx <- idx[idx2]
  }
  
  return(res)
}

# simulate posterior predictive datasets
# condition on censoring, i.e., censoring events are held unchanged for all posterior predictive datasets
# using latent failure time method
# changed name to sim_ppdata, pp for posterior predictive
sim_ppdata <- function(data,s,p_beta,p_lambda,seed=1){
  
  # p_beta: nsim x ? x ? 
  
  set.seed(seed)
  nsim <- dim(p_beta)[1]
  k <- length(s)-1; m <- data$m; p <- data$p; n <- data$n
  s_diff <- tail(s,k)-head(s,k)
  
  # only get information for those that need to be simulated
  ise <- data$data$event != 0 # is event, not censored
  ne <- sum(ise) # number of events, not censored
  Xe <- data$X[ise,] # design matrix for the events
  
  Us <- array(runif(nsim*ne*m), dim=c(ne,m,nsim))
  
  dim(p_beta) <- c(nsim,p,m)
  dim(p_lambda) <- c(nsim,k,m)
  p_lambda_ext <- array(0,dim=c(nsim,k+1,m))
  p_lambda_ext[,1:k,] <- p_lambda
  p_lambda_ext[,k+1,] <- p_lambda_ext[,k,]  
  
  cumLambda <- array(0,dim=c(k+1,m,nsim))
  cumLambda[-1,,] <- vapply(c(1:nsim), function(i) colCumsums(p_lambda[i,,]*s_diff), 
                            matrix(0,k,m))
  
  # ne x m x nsim
  eXbeta <- exp(vapply(c(1:nsim), function(i) Xe %*% p_beta[i,,], matrix(0,ne,m)))
  H_rep <- -log(Us)/eXbeta
  
  x1 <- rep(c(1:nsim),each=m)
  x2 <- rep(c(1:m),nsim)
  x3 <- cbind(x1,x2)
  
  x1 <- rep(1:(nsim*m),each=ne)
  x4 <- x3[x1,]
  
  L <- vapply(c(1:(nsim*m)), function(i) 
    findInterval(H_rep[,x3[i,2],x3[i,1]], cumLambda[,x3[i,2],x3[i,1]]), integer(ne))
  
  x1 <- cbind(c(L),x4[,2],x4[,1])
  x2 <- cbind(x4[,1],c(L),x4[,2])
  
  cumlambda_l <- array(cumLambda[x1], dim=c(ne,m,nsim))
  lambda_l <- array(p_lambda_ext[x2], dim=c(ne,m,nsim))
  T_latent <- (H_rep - cumlambda_l)/lambda_l + array(s[L], dim=c(ne,m,nsim))
  
  # # 20241222 checking
  # res <- array(0,dim=c(ne,m,nsim))
  # for (i in 1:nsim){
  #   for (j in 1:m){
  #     for (l in 1:ne){
  #       tmp <- findInterval(H_rep[l,j,i],cumLambda[,j,i])
  #       res[l,j,i] <- (H_rep[l,j,i]-cumLambda[tmp,j,i])/p_lambda_ext[i,tmp,j] + s[tmp]
  #     }
  #   }
  # }
  # all.equal(T_latent,res)
  
  # "e" for events
  time_rep_e <- vapply(c(1:nsim), function(i) rowMins(T_latent[,,i]), numeric(ne)) 
  event_rep_e <- vapply(c(1:nsim), function(i) max.col(-T_latent[,,i], ties="first"), 
                        numeric(ne))
  
  delta_rep_e <- array(0,dim=c(nsim,ne,m))
  x1 <- cbind(rep(c(1:nsim),each=ne),rep(c(1:ne),nsim),c(event_rep_e))
  delta_rep_e[x1] <- 1
  
  if (ne==n){
    time_rep <- time_rep_e
    event_rep <- event_rep_e
    delta_rep <- delta_rep_e
  } else {
    time_rep <- event_rep <- matrix(0,n,nsim)
    delta_rep <- array(0,dim=c(nsim,n,m))
    time_rep[!ise,] <- data$data$time[!ise]
    time_rep[ise,] <- time_rep_e
    event_rep[ise,] <- event_rep_e
    delta_rep[,ise,] <- delta_rep_e
  }
  
  list(time_rep=t(time_rep), event_rep=t(event_rep), delta_rep=delta_rep,
       time_rep_e=t(time_rep_e), event_rep_e=t(event_rep_e), delta_rep_e=delta_rep_e)
}


### MCMC ====

# normal linear model
nlr_MCMC <- function(V,y,sigma,mu,Psi,nsim,seed=1){
  set.seed(seed)
  
  if (Psi[1,1]==Inf){
    S <- solve(crossprod(V))*sigma^2
    m <- S %*% (crossprod(V, y)/sigma^2)
  } else {
    S <- solve(crossprod(V)/sigma^2+solve(Psi))
    m <- S %*% (crossprod(V, y)/sigma^2+ solve(Psi) %*% mu)
  }
  
  rmvnorm(nsim, mu=m, Sigma=S)
  
}

# normal mean model
normal_MCMC <- function(y,sigma,mu0,s0,nsim,seed=1){
  set.seed(seed)
  n <- length(y)
  
  v2 <- 1/(n/sigma^2+1/s0^2)
  m <- v2*(sum(y)/sigma^2+mu0/s0^2)
  rnorm(nsim, m, sqrt(v2))
}

# social mobility data log linear model
mob_MCMC <- function(y,X,a0s,b0s,burn=2000,nscan=10000,odens=10,seed=1){
  
  # n: number of observations per country
  # q: number of coefficients
  
  # y: n x 1, counts
  # X: n x q, design matrix
  
  # a0s, b0s: q-dim vector, priors
  
  set.seed(seed)
  n <- length(y); q <- dim(X)[2]
  
  stopifnot(length(a0s)==q)
  stopifnot(length(b0s)==q)
  
  # store results
  n_op <- floor(nscan / odens)
  LAMBDA <- matrix(0, nrow=n_op, ncol=q)
  j <- 0
  
  # initialize
  lambda <- c(mean(y),rep(1,q-1))
  
  a1s <- colSums(X * outer(y, rep(1,q))) + a0s
  idx <- X>0; ic <- colSums(idx)
  
  for (i in 1:(nscan+burn)){
    
    for (k in 1:q){
      
      # tested this isn't slow
      tmp <- t(replicate(ic[k],lambda[-k]))^X[idx[,k],-k]
      b1k <- sum(apply(tmp,1,prod)) + b0s[k]
      
      lambda[k] <- rgamma(1, a1s[k], b1k) 
    }
    
    # store results  
    if ((i > burn) && (i %% odens == 0)){
      j <- j+1
      LAMBDA[j,] <- lambda
    }
  }
  
  LAMBDA
  
}

### plot p ====
plotp <- function(res, idx=NULL){
  
  pval <- names(res); np <- length(pval)
  titles <- colnames(res[[1]])
  
  if (is.null(idx)) idx<-c(1:dim(res[[1]])[2])
  
  par(mfrow=c(np,length(idx)),mar=c(0,2.5,1,0), mgp=c(1.75,0.75,0), oma=c(2,2,1,1))
  cex_main <- 0.9
  
  for (i in seq_along(pval)){
    for (j in seq_along(idx)){
      
      if (i==1){
        
        plot(res[[i]][,idx[j]], ylim=c(0,1), ylab=ifelse(j==1, pval[i],""), pch=20, cex=0.5,
             font.main = 1, cex.main=cex_main, main=titles[j], xlab="", xaxt="n", yaxt="n")
        
        
      } else {
        plot(res[[i]][,idx[j]], ylim=c(0,1), ylab=ifelse(j==1, pval[i],""), pch=20, cex=0.5,
             font.main = 1, cex.main=cex_main, main="", xlab="", xaxt="n", yaxt="n")
      }
      if (j==1) axis(side = 2, at = seq(0,1,0.2), labels = T)
      #if (i==length(pval)) axis(side = 1, at = seq(0,1000,200), labels = T)
      
    }
  }
  
}

# in this version, we plot all discrepancy measures in the same plot
# 20240806, we fix the color mapping for the p-values
# 20240806, we change the legends to be using latex
plotp_density2 <- function(res, spec=NULL){
  
  # things can be specified in spec:
  # op_density: 1 for usual density; 2 (default) for bde; 3 for the density codes online
  # - d_idx: numeric vector, discrepancies to be plotted
  # - fs: font size, default at 0.8
  # - prior_idx: numeric vector, priors to be plotted
  # - p_idx: numeric vector, p-values to be plotted
  # - prior_label: character vector, labels for the priors, the order of this should follow the original order in the data rather than prior_idx
  # - d_label: character vector, labels for the discrepancy measures
  # - ymax: to control the plot ymax
  # - na.rm: if T, remove all na's. Default to be F so that I can get warnings for NAs
  
  # these labels are in the order I always want to be presenting (based on the structural comparison figure in the draft)
  p_labels_all <- c("ppost","pdpost","psplit","pdsplit","ptm","pepd","pptpost","pdptpost")
  p_labels_all_latex <- c(expression(p[post]),
                          expression(p[dpost]),
                          expression(p[split]),
                          expression(p[dsplit]),
                          expression(p[tm]),
                          expression(p[epd]),
                          expression(p[ptpost]),
                          expression(p[dptpost]))
  colors_all <- c(1,1,4,4,3,2,6,6)
  color_mapping <- data.frame(p_label=p_labels_all,color=colors_all)
  
  # note cannot change to use ifelse when the TRUE/FALSE values are not same length
  if (is.null(spec$d_idx)){
    d_idx <- c(1:length(res))
  } else {
    d_idx <- spec$d_idx
  }
  
  if (is.null(spec$fs)){
    fs <- 0.8
  } else {
    fs <- spec$fs
  }
  
  if (is.null(spec$prior_idx)){
    prior_idx <- c(1:dim(res[[1]][[1]])[2])
  } else {
    prior_idx <- spec$prior_idx
  }
  
  if (is.null(spec$p_idx)){
    p_idx <- c(1:length(res[[1]]))
  } else {
    p_idx <- spec$p_idx
  }
  
  if (is.null(spec$prior_label)){
    prior_label <- colnames(res[[1]][[1]])
  } else {
    prior_label <- spec$prior_label
  }
  
  if (is.null(spec$p_label)){
    p_label <- names(res[[1]])
  } else {
    p_label <- spec$p_label
  }
  
  if (is.null(spec$d_label)){
    d_label <- names(res)
  } else {
    d_label <- spec$d_label
  }
  
  if (is.null(spec$op_density)) {
    op_density <- 2
  } else {
    op_density <- spec$op_density
  }
  
  if (!is.null(spec$ymax)) ymax <- spec$ymax

  na.rm <- ifelse(is.null(spec$na.rm), FALSE, spec$na.rm)

  nprior <- length(prior_idx); np <- length(p_idx); nd<-length(d_idx)
  
  # map colors -- this is for plotting
  #mycols <- brewer.pal(np, "Spectral")
  mycols_idx <- data.frame(p_label=p_label) %>%
    left_join(color_mapping, join_by(p_label)) %>%
    select(color) %>% pull()
  mycols <- wes_palette("Zissou1", 6, type = "continuous")[mycols_idx]
  
  # map latex legend labels -- this is for legends
  p_label_idx <- which(p_labels_all %in% p_label)
  p_labels_latex <- p_labels_all_latex[p_label_idx]
  
  p_label_ordered <- p_labels_all[p_labels_all %in% p_label]
  mycols_idx_ordered <- data.frame(p_label=p_label_ordered) %>%
    left_join(color_mapping, join_by(p_label)) %>%
    select(color) %>% pull()
  mycols_ordered <- wes_palette("Zissou1", 6, type = "continuous")[mycols_idx_ordered]
  
  # plotting
  par(mfrow=c(nd,nprior),mar=c(0,0,1,1), mgp=c(1.75,0.75,0), oma=c(4,4,1,1))
  for (i in 1:nd){
    dm <- d_idx[i]
    for (j in 1:nprior){
      prior <- prior_idx[j]
      if (is.null(spec$ymax)){
        tmp <- sapply(res[[dm]], function(x) x[,prior])
        ymax <- round(getrange(tmp,op_density,na.rm),2)
      }
      
      for (k in 1:np){
        res_i <- res[[i]][[k]][,prior]
        if (na.rm) res_i <- res_i[!is.na(res_i)]
        
        if (op_density==1){
          est_density <- density(res_i)
        } else if (op_density==2) {
          # tested different options of mu under the boundary kernel
          # the default (mu=1) seems to be the most smooth
          est_density <- bde(res_i, estimator="boundarykernel", options = list(mu=1)) 
        } else if (op_density==3){
          est_density <- bd_density(res_i)
        }

        if (k==1){
          plot(est_density, xlim=c(0,1), type="l",xlab="",lwd=2, main="",
               col=mycols[1],ylim=c(0,ymax),ylab="",
               xaxt="n",yaxt="n",cex.axis=fs, cex.lab=fs)
        } else {
          lines(est_density, col=mycols[k],lwd=2)
        }
      }

      #if (i==1 & j==nprior)
      #  legend(legend_loc, legend = p_labels_latex, y.intersp=0.75,
      #         col = mycols_ordered, lwd=rep(2,np), bty="n", cex = fs)
      if (i==1) title(main=list(prior_label[prior],cex=fs*1.1),line=0.25)
      if (j==1) mtext(side=2, text=d_label[i], line=1.75, cex=fs*0.8)
      if (i==nd) axis(side = 1, at = seq(0,1,0.2), labels = T, cex.axis=fs)
      if (j==1) axis(side=2, labels=T, cex.axis=fs)
    }
  }
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend('bottom', legend=p_labels_latex, col=mycols_ordered, lwd=rep(2,np), 
         xpd = TRUE, horiz = TRUE, cex = fs*1.1, seg.len=2, bty = 'n')
  
}

plotp_hist2 <- function(res, spec=NULL){
  
  nd <- length(res); np <- length(res[[1]]); nprior <- dim(res[[1]][[1]])[2]; nsim <- dim(res[[1]][[1]])[1]
  
  p_label <- names(res[[1]])
  p_df <- rep(p_label, each=nsim)
  d_label <- names(res)
  d_df <- rep(d_label, each=nsim*np)
  
  df <- data.frame(do.call("rbind",res[[1]]))
  for (d in 2:nd){
    df <- rbind(df, data.frame(do.call("rbind",res[[d]])))
  }
  colnames(df) <- colnames(res[[1]][[1]])
  df$p <- rep(p_df,nd)
  df$d <- d_df
  
  df <- pivot_longer(df, cols=!c("p","d"), names_to = "prior", values_to = "value")
  
  df %>% ggplot()+
    geom_histogram(aes(x=value,y=..density..,fill=p),alpha = 0.5, position = "identity") + 
    facet_grid(rows=vars(d),cols=vars(prior),scales = "free_y") +
    scale_fill_manual(values=mycols) +
    theme_bw()
}

plotp_density <- function(res, spec=NULL){
  
  if (is.null(spec$fs)){
    fs <- 0.8
  } else {
    fs <- spec$fs
  }
  
  if (is.null(spec$prior_idx)){
    prior_idx <- c(1:dim(res[[1]])[2])
  } else {
    prior_idx <- spec$prior_idx
  }
      
  if (is.null(spec$p_idx)){
    p_idx <- c(1:length(res))
  } else {
    p_idx <- spec$p_idx
  }
  
  if (is.null(spec$prior_label)){
    prior_label <- colnames(res[[1]])
  } else {
    prior_label <- spec$prior_label
  }
  
  if (is.null(spec$p_label)){
    p_label <- names(res)
  } else {
    p_label <- spec$p_label
  }
  
  nprior <- length(prior_idx); np <- length(p_idx)
  
  #mycols <- brewer.pal(np, "Spectral")
  mycols <- wes_palette("Zissou1", np, type = "continuous")

  par(mfrow=c(1,nprior),mar=c(0,3,2,0), mgp=c(1.75,0.75,0), oma=c(2,2,1,1))
  for (i in seq_along(prior_idx)){
    prior <- prior_idx[i]
    tmp <- sapply(res, function(x) x[,prior])
    ymax <- round(getrange(tmp),2)
    
    if (i>1){
      plot(density(res[[1]][,prior]), xlim=c(0,1), main=prior_label[prior], 
           xlab="",col=mycols[1],lwd=2, ylim=c(0,ymax),ylab="", 
           cex.main=fs, cex.axis=fs)
    } else{
      plot(density(res[[1]][,prior]), xlim=c(0,1), main=prior_label[prior], 
          xlab="",col=mycols[1],lwd=2, ylim=c(0,ymax), 
          cex.main=fs, cex.axis=fs)
    }
    for (j in 2:np){
      lines(density(res[[j]][,prior]), col=mycols[j],lwd=2)
    }
    if (i==1)
      legend("topleft", legend = p_label,
             col = mycols, lwd=rep(2,np), bty="n", cex = fs)
  }
  
}

combine <- function(labels,...){
  res <- list(...)
  nr <- length(res)
  nd <- length(res[[1]]); np <- length(res[[1]][[1]])
  
  for (i in 1:nd){
    for (j in 1:np){
      for (k in 2:nr){
        res[[1]][[i]][[j]] <- cbind(res[[1]][[i]][[j]], res[[k]][[i]][[j]])
      }
      colnames(res[[1]][[i]][[j]]) <- labels
    }
  }
  res[[1]]
}

bd_density <- function(x){
  
  # Construct histograms
  his = hist(x, breaks = 60, plot = F)
  X = his$counts
  u = his$mids
  
  # Prepare basis (I-mat) and penalty (1st difference)
  B = diag(length(X))
  D1 = diff(B, diff = 1)
  lambda = 1e6 # fixed but can be selected (e.g. AIC)
  P = lambda * t(D1) %*% D1
  
  # Smooth
  tol = 1e-8
  eta = log(X + 1)
  for (it in 1:20){
    mu = exp(eta)
    z = X - mu + mu * eta
    a = solve(t(B) %*% (c(mu) * B) + P, t(B) %*% z)
    etnew = B %*% a
    de = max(abs(etnew - eta))
    # cat('Crit', it, de, '\n')
    if(de < tol) break
    eta = etnew
  }
  
  smooth.spline(u, his$density)
  
}

getrange <- function(res, op_density=1, na.rm=F){
  p <- dim(res)[2]
  upper <- numeric(p)
  
  for (i in 1:dim(res)[2]){
    
    if (na.rm){
      res_i <- res[,i][!is.na(res[,i])]
    } else {
      res_i <- res[,i]
    }
    
    if (op_density==1){
      upper[i] <- max(density(res_i)$y)
    } else if (op_density==2) {
      upper[i] <- max(bde(res_i, estimator="boundarykernel")@densityCache )
    } else if (op_density==3) {
      upper[i] <- max(bd_density(res_i[[1]][,prior])$y)
    }
  }
  
  # if the largest is too large compared to the others, return the second largest
  if (max(max(upper) / upper) >3 && max(upper) > 5){
    max(upper[-which.max(upper)])
  } else {
    max(upper)
  }
}


### aggregation functions ====

to_agg <- function(label,folder){
  
  p_vals <- c("ppost","pptpost","pepd","ptm","psplit")
  
  for (i in seq_along(p_vals)){
    agg(label,folder,p_vals[i])
  }
  
}

# this version is used for runs in the 20240804 folder, with run time
agg <- function(runnum,folder,p_label){
  
  fdir <- paste0(outfolder,folder,"/")
  
  # Files to aggregate 
  f <- list.files(fdir)
  f <- f[stringr::str_detect(f, paste0(runnum,"_",p_label))]
  
  res <- readRDS(paste0(fdir,f[1]))
  res <- res$output
  
  nd <- dim(res)[2]
  labels <- dimnames(res)
  discrepancy <- labels[[2]]
  priors <- labels[[1]]
  
  # Read into memory
  output <- vector(mode="list", length=length(nd))
  runtime <- vector(mode="list",length=length(f))
  for (d in 1:nd){
    
    out <- NULL
    for (i in seq_along(f)){
      
      idx <- as.integer(str_extract(f[i],"(?<=_)\\d+(?=\\.rds)"))

      res <- readRDS(paste0(fdir,f[i]))
      runtime[[idx]] <- res$runtime
      out <- rbind(out, (res$output[,d]))

    }
    
    colnames(out) <- priors
    output[[d]] <- out
  }
  names(output) <- discrepancy
  res <- list(output=output,runtime=runtime)
  saveRDS(res, paste0("../res/",folder, "_",runnum,"_",p_label,".rds"))

}

wanted <- function(runnum,folder,p_label,nrun=1000){
  
  fdir <- paste0(outfolder,folder,"/")
  check <- logical(nrun)
  
  f <- list.files(fdir)
  f <- f[stringr::str_detect(f, paste0(runnum,"_",p_label,"_"))]
  for (i in seq_along(f)){
    idx <- as.integer(str_extract(f[i],"(?<=_)\\d+(?=\\.rds)"))
    check[idx] <- T
  }
  list(total=sum(check==F), list=paste(which(check==F), collapse = ","))
  
}

check_log <- function(simname,simnum,nrun){
  check <- logical(nrun)
  for (i in 1:nrun){
    log_i <- read.csv(file=paste0("../out/",simname,simnum,"_",i))
    check[i] <- any(sapply(tail(log_i,5), str_detect, "elapsed"))
  }
  paste(which(check==F), collapse = ",")
}

# Misc ----

# I use this function to calculate the MC results of the p-values because for some special 
# runs where all the p-values are the same, if I directly check number of rep values >= obs value
# it might not give 1, but instead gives some seemingly meaningful numbers due to rounding errors

pr_ge <- function(x,y,tol=1.5e-8){
  
  # probability of greater or equal to 0 (ge from the latex greater than or equal to)
  # x: matrix, m x n
  # y: matrix, m x n
  
  diff <- x-y
  diff[abs(diff)<tol] <- 0
  colMeans(diff >= 0)

}

# function to combine p value results into previous format
combine_ps <- function(labels_out,...,res=NULL){
  
  if (is.null(res)) res <- list(...)
  
  np <- length(res)
  nd <- length(res[[1]]); nprior <- dim(res[[1]][[1]])[2]
  
  d_labels <- names(res[[1]])
  
  res_out <- vector(mode="list",length=nd)
  for (i in 1:nd){
    for (j in 1:np){
      res_out[[i]][[j]] <- res[[j]][[i]]
    }
    names(res_out[[i]]) <- labels_out
  }
  names(res_out) <- d_labels
  return(res_out)
}

# function to combine runtime results
combine_rts <- function(labels_out,...,res=NULL){
  
  if (is.null(res)) res <- list(...)
  np <- length(res)
  
  res_out <- vector(mode="list",length=np)
  for (i in 1:np){
    res_out[[i]] <- map_dbl(res[[i]], `[`, "elapsed")
  }
  
  df <- do.call(cbind,res_out)
  colnames(df) <- labels_out
  
  return(df)
  
}


# using two sets of labels, because of the thing ppost vs pdpost
read_nlr_res <- function(runname,num,labels_in,labels_out){
  
  np <- length(labels_in)
  res <- vector(mode="list",length=np)
  
  for (i in 1:np){
    res[[i]] <- readRDS(paste0("res/", runname, "_",num,"_",labels_in[[i]],".rds"))
  }
  
  res_output <- sapply(res,`[`,"output")
  res_runtime <- sapply(res,`[`,"runtime")
  
  output <- combine_ps(labels_out, res=c(res_output))
  runtime <- combine_rts(labels_out, res=c(res_runtime))
  
  list(output=output,runtime=runtime)  
}

