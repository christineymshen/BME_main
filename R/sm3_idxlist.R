source("functions.R")

# survival data true parameters
n <- 100; m<-2; p<-3
gammas <- c(1,1) # scale
alphas <- c(1/0.3,1/0.5) # shape

# to test for betas that change with t
betas0 <- matrix(c(-10,-1,-0.2,0.3,
                   -6.5,-0.4,-1.1,0.2), nrow=m, byrow=T)
betas1 <- matrix(c(0,0.016,0.024,0.028,
                   0,0.028,0.024,0.016), nrow=m, byrow=T)

### fixed parameters

# for p-values
one_n <- rep(1,n)

# fix design matrix for all simulations
set.seed(1)
x1 <- rbinom(n,1,0.5)
x2 <- rnorm(n)
x3 <- runif(n,0,2)
X <- cbind(one_n, x1,x2,x3)

data <- vector("list",length=1000)
idxlist <- numeric(1000)
idx <- 0
for (slurm_id in 1:2000){
  tmp <- tryCatch({simdata7(alphas,gammas,X,betas0,betas1,seed=slurm_id)},
                   error=function(cond) "error")
  if (sum(tmp=="error")!=1){
    idx <- idx + 1
    data[[idx]] <- tmp
    idxlist[idx] <- slurm_id
    
    if (idx==1000) break
  }
}
saveRDS(idxlist,"res/sm3_idxlist.rds")
