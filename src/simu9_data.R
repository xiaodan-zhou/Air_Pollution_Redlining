###########################################################################
######################### Generate Data ###################################
###########################################################################
library(tidyr)
library(dplyr)
library(ggplot2)

simu_data <- function(outpath, n = 490, seed = 128, 
                      ratio.u = Inf, ratio.uy = .6, return.sp = F) {
  
  source("./src/data_gen.R")
  source("./src/spline.R")
  
  ## store parameters
  argg <- as.list(environment())
  
  L1 <- L2 <- 20 
  
  sigma.u <- 1; sigma.z <- 1
  eta.1 <- -1; eta.2 <- -2 
  sigma2.w1 <- 1; alpha.w1u <- 2
  sigma2.w2 <- 1; alpha.w2u <- -.5
  sigma2.v <- 1; alpha.vu <- -1
  sigma2.y <- .25; alpha.yu <- -1 
  theta.1 <- 0; theta.2 <- 0.2
  omega <- -1
  tobit.pct <- .05
  
  #### simu9 noisy proxy originally this is 3 or 10 for all 
  sigma2.w1 <- .5
  sigma2.w2 <- .5
  sigma2.v <- .5 
  #### 
  
  ## use 49 regions in each city
  n.city <- as.integer(n / 49)
  n.region <- ceiling(n / n.city)
  SPL <- rep(1:n.city, each=n.region)[1:n] ## spatial level for cities
  
  geometry <- get.geometry(SPL)
  
  splineBU <- get.spline(SPL, geometry, L1=L1, L2=L2, ratio=ratio.u, npin=100)
  ncol.BU <- unlist(lapply(splineBU, function(x) if(is.matrix(x)) ncol(x) else 1))
  
  splineBZ <- get.spline(SPL, geometry, L1=L1, L2=L2, ratio=ratio.uy, npin=100) 
  ncol.BZ <- unlist(lapply(splineBZ, function(x) if(is.matrix(x)) ncol(x) else 1))
  
  ## seed ## 
  set.seed(seed)
  
  ################## -> Uy ##################
  ZZ <- rnorm(sum(ncol.BZ), 0, sd=sigma.z) # Uy has same scale as U
  BZ <- compute.spline(splineBZ, ZZ)
  
  UU <- list(); BU <- c()
  for (i in 1:n.city) {
    UU[[i]] <- matrix(rnorm(ncol.BU[i], 0, sd=sigma.u), ncol=1)
    BU <- rbind(BU, splineBU[[i]] %*% UU[[i]]) 
  }
  
  ## seed ## 
  
  argg$ZZ <- ZZ
  argg$UU <- UU
  
  W1 <- alpha.w1u*BU + rnorm(n, 0, sd=sqrt(sigma2.w1))
  W1 <- W1 - mean(W1, na.rm=T)
  W1 <- W1/sd(W1, na.rm=T) 
  
  W2 <- alpha.w2u*BU + rnorm(n, 0, sd=sqrt(sigma2.w2))
  W2 <- W2 - mean(W2, na.rm=T)
  W2 <- W2/sd(W2, na.rm=T)
  
  V <- alpha.vu*BU + rnorm(n, 0, sd=sqrt(sigma2.v))
  tobit.cut <- quantile(V, tobit.pct); V[V < tobit.cut] <- tobit.cut
  V <- V/sd(V,na.rm=T)
  
  prob.AB <- exp(eta.1 + eta.2*BU + BZ)/(1+exp(eta.1 + eta.2*BU + BZ))
  
  X <- rep(NA, n)
  for(i in 1:n) X[i] <- which.max(rmultinom(1, size=1, prob=c(prob.AB[i], 1-prob.AB[i])))
  
  Y <- theta.1*(X==1) + theta.2*(X==2) + alpha.yu*BU + omega*BZ + rnorm(n, 0, sd=sqrt(sigma2.y))
  
  dt <- data.frame(Y=Y, X=X, W1=W1, W2=W2, V=V, BU=BU, BZ=BZ, SPL=SPL)
  dt$geometry <- geometry
  
  dt <- dt %>% group_by(SPL) %>% mutate(Y = Y - mean(Y, na.rm=T)) %>% data.frame(.) 
  
  if(return.sp) {
    save(dt, argg, splineBU, splineBZ, file = outpath) 
  } else {
    save(dt, argg, file = outpath)
  }
  
  return()
  
}


############### Generate dataset ################
out.folder <- "output/simu9_randUZ_data"
dir.create(out.folder, recursive = TRUE)
set.seed(19)

n <- 490 
nrep <- 500
seeds <- runif(nrep, min=0, max=nrep*100)

for (i in 1:nrep) {
  outpath <- paste0(out.folder, "/data", i, ".Rdata")
  simu_data(outpath, n = n, seed = seeds[i], 
            ratio.u = Inf, ratio.uy = .4, 
            return.sp = F)
}

############### Generate data with different number of splines with the same grid ################
i <- 1
ratio.uy <- 0 # not matter

for (ratio.u in c(seq(0, 1, 0.1), Inf)) {
  outpath <- file.path(out.folder, paste0("data", i, "_s", ratio.u*100, ".Rdata"))
  simu_data(outpath, n = 490, seed = 1, 
            ratio.u = ratio.u, ratio.uy = ratio.uy,
            return.sp = T)
}
