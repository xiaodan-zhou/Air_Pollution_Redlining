###########################################################################
######################### Generate Data ###################################
###########################################################################
### we use geometry of redlining data
### to reduce runtime, we take a subset
###########################################################################
library(tidyr)
library(dplyr)
library(ggplot2)

###########################################################################
mimic_redlining_u1 <- function(out.folder, id, input.path, outcome, 
                               ratio.u = Inf, ratio.uy = .6, 
                               n.min = 500, seed = 128) {
  source("./src/data_gen.R")
  source("./src/spline.R")
  
  set.seed(seed)
  
  load(input.path) # include model and dt
  
  list2env(lapply(colMeans(model$out[model$out$post==1,]), as.vector), .GlobalEnv)
  argg <- as.list(colMeans(model$out[model$out$post==1,]))
  SPL.copy <- dt$SPL
  
  ## select so that sample size is n
  sel_SPL <- c()
  total_rows <- 0
  while(total_rows < n.min) {
    new_spl <- sample(setdiff(SPL.copy, sel_SPL), 1)
    sel_SPL <- c(sel_SPL, new_spl)
    total_rows <- sum(SPL.copy %in% sel_SPL)
  } 
  ## sort!!! 
  sel_SPL <- sort(sel_SPL)
  
  idx <- dt$SPL %in% sel_SPL
  SPL <- dt$SPL[idx]
  n <- sum(idx)
  tobit.pct <- mean(dt$V == min(dt$V)) 
  n.city <- length(sel_SPL)
  
  splineBU.temp <- readRDS(paste0("data/nhgis0008_csv/only 1940/RedliningSplineMat_R", ratio.u*100, ".rds"))
  splineBZ.temp <- readRDS(paste0("data/nhgis0008_csv/only 1940/RedliningSplineMat_R", ratio.uy*100, ".rds"))
  splineBU <- splineBU.temp[sel_SPL]; rm(splineBU.temp)
  splineBZ <- splineBZ.temp[sel_SPL]; rm(splineBZ.temp)
  
  ncol.BU <- unlist(lapply(splineBU, function(x) if(is.matrix(x)) ncol(x) else 1))
  ncol.BZ <- unlist(lapply(splineBZ, function(x) if(is.matrix(x)) ncol(x) else 1))
  
  sigma.u <- sqrt(sigma2.u)
  sigma.z <- sqrt(sig2.Z)
  
  ZZ <- rnorm(sum(ncol.BZ), 0, sd=sigma.z)  
  BZ <- compute.spline(splineBZ, ZZ)
  
  UU <- list(); BU <- c()
  for (i in 1:n.city) {
    UU[[i]] <- matrix(rnorm(ncol.BU[i], 0, sd=sigma.u), ncol=1)
    BU <- rbind(BU, splineBU[[i]] %*% UU[[i]]) 
  }
  
  W1 <- alpha.w1 + alpha.w1u*BU + rnorm(n, 0, sd=sqrt(sigma2.w1))
  W2 <- alpha.w2 + alpha.w2u*BU + rnorm(n, 0, sd=sqrt(sigma2.w2))
  V <- alpha.v + alpha.vu*BU + rnorm(n, 0, sd=sqrt(sigma2.v))
  tobit.cut <- quantile(V, tobit.pct)
  V[V < tobit.cut] <- tobit.cut
  
  prob.AB <- exp(eta.1 + eta.2*BU + BZ)/(1+exp(eta.1 + eta.2*BU + BZ))
  
  X <- rep(NA, n)
  for(i in 1:n) X[i] <- which.max(rmultinom(1, size=1, prob=c(prob.AB[i], 1-prob.AB[i])))
  
  Y <- theta.1*(X==1) + theta.2*(X==2) + alpha.yu*BU + omega*BZ + rnorm(n, 0, sd=sqrt(sigma2.y))
  
  dt_new <- data.frame(Y=Y, X=X, W1=W1, W2=W2, V=V, BU=BU,BZ=BZ)
  dt_new$SPL <- as.numeric(factor(SPL, levels = unique(SPL)))
  # plot(SPL, dt_new$SPL) 
  
  idxSPL <- sel_SPL ## sorted. 
  
  dt <- dt_new 
  
  dt <- dt %>% group_by(SPL) %>% mutate(Y = Y - mean(Y, na.rm=T)) %>% data.frame(.) 
  dt$W1 <- dt$W1/sd(dt$W1,na.rm=T)
  dt$W2 <- dt$W2/sd(dt$W2,na.rm=T)
  dt$V <- dt$V/sd(dt$V,na.rm=T)
  
  if (!is.na(out.folder)) {
    outpath <- paste0(out.folder, "/data", id, "_", outcome, ".RData")
    save(dt, idxSPL, argg, file = outpath)
    return()
  }
  
  return(list(dt=dt, idxSPL=idxSPL, argg=argg))
}


############### Set-up ################
out.folder <- "output/simured1_randUZ_data"
dir.create(out.folder, recursive = TRUE)
set.seed(12)

n.min <- 500 
nrep <- 500
seeds <- runif(nrep, min=0, max=nrep*100)

############### PM25 ################
input.path <- "output/redlining_mcmc_latent_u1w3/pm25_1u_3w_sInf_s60.RData"
outcome <- "pm25"

for (i in 1:nrep) {
  mimic_redlining_u1(out.folder, i, input.path, outcome, n.min = n.min, seed = seeds[i])
}

############### PM25 ################
input.path <- "output/redlining_mcmc_latent_u1w3/no2_1u_3w_sInf_s60.RData"
outcome <- "no2"

for (i in 1:nrep) {
  mimic_redlining_u1(out.folder, i, input.path, outcome, n.min = n.min, seed = seeds[i])
}
