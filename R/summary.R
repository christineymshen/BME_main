
source("functions.R")

labels_in <- c("ppost","pptpost","pepd","ptm","psplit")
labels_out <- c("pdpost","pdptpost","pepd","ptm","pdsplit")

# NLR


res1_ppost <- readRDS("res/nlr_1_ppost.rds")
res1_pptpost <- readRDS("res/nlr_1_pptpost.rds")
res1_pepd <- readRDS("res/nlr_1_pepd.rds")
res1_ptm <- readRDS("res/nlr_1_ptm.rds")
res1_psplit <- readRDS("res/nlr_1_psplit.rds")

res2 <- read_nlr_res("nlr",2,labels_in,labels_out)
res3 <- read_nlr_res("nlr",3,labels_in,labels_out)
res4 <- read_nlr_res("nlr",4,labels_in,labels_out)
res5 <- read_nlr_res("nlr",5,labels_in,labels_out)

# level
output_level <- combine_ps(labels_out, res1_ppost$output, res1_pptpost$output, res1_pepd$output, res1_ptm$output, res1_psplit$output)

spec <- vector(mode="list")
spec$ymax <- 5
spec$prior_idx <- c(1,4,2,3)
spec$fs <- 1

pdf(file = paste0("../draft/fig/nlr_level.pdf"), width=8, height=5)
plotp_density2(output_level, spec=spec)
dev.off()

# power
res_power <- combine(labels=c("overdispersion", "underdispersion","lognormal","missing covariates"), 
                     res4$output,res3$output,res5$output,res2$output)

spec <- vector(mode="list")
spec$ymax <- 5
spec$fs <- 1

plotp_density2(res_power, spec=spec)

# GGLM

# level
res_1 <- read_nlr_res("gglm",1,labels_in,labels_out)
res_2 <- read_nlr_res("gglm",2,labels_in,labels_out)
# power
res_3 <- read_nlr_res("gglm",3,labels_in,labels_out)
res_4 <- read_nlr_res("gglm",4,labels_in,labels_out)
res_5 <- read_nlr_res("gglm",5,labels_in,labels_out)

# level
spec <- vector(mode="list")
spec$ymax <- 6
spec$fs <- 1
spec$prior_idx <- c(1,3,2)

pdf(file = paste0("../draft/fig/gglm_fixed_level.pdf"), width=6, height=7)
plotp_density2(res_1$output, spec=spec)
dev.off()

# level joint model
spec <- vector(mode="list")
spec$ymax <- 6
spec$fs <- 1
spec$prior_idx <- c(1,3,2)

pdf(file = paste0("../draft/fig/gglm_joint_level.pdf"), width=6, height=7)
plotp_density2(res_2$output, spec=spec)
dev.off()

# power

res_3_1 <- lapply(res_3$output, function(x) lapply(x, function(y) y[,1]))
res_4_1 <- lapply(res_4$output, function(x) lapply(x, function(y) y[,1]))
res_5_1 <- lapply(res_5$output, function(x) lapply(x, function(y) y[,1]))

res_345 <- combine(labels=c("missing covariates","wrong link","lognormal"), res_3_1,res_4_1,res_5_1)

spec <- vector(mode="list")
spec$ymax <- 6
spec$fs <- 1

pdf(file = paste0("../draft/fig/gglm_fixed_power.pdf"), width=6, height=7)
plotp_density2(res_345, spec=spec)
dev.off()

## Survival Model


res_1 <- read_nlr_res("sm",1,labels_in,labels_out)
res_2 <- read_nlr_res("sm",2,labels_in,labels_out)
res_3 <- read_nlr_res("sm",3,labels_in,labels_out)

res_2_1 <- lapply(res_2$output, function(x) lapply(x, function(y) y[,1]))
res_3_1 <- lapply(res_3$output, function(x) lapply(x, function(y) y[,1]))

res_23 <- combine(labels=c("missing covariates","time-varying coefficients"), res_2_1,res_3_1)


spec <- vector(mode="list")
spec$ymax <- 6
spec$fs <- 1

pdf(file = paste0("../draft/fig/sm_level.pdf"), width=5, height=4)
plotp_density2(res_1$output, spec=spec)
dev.off()

pdf(file = paste0("../draft/fig/sm_power.pdf"), width=5, height=4)
plotp_density2(res_23, spec=spec)
dev.off()