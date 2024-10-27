library(MASS)
library(LaplacesDemon)
library(zoo)
library(mvtnorm)
library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(furrr)

source("./src/helper.R")
source("./src/spline.R")

gibbs_sample <- function(Y, X, W1, W2, V, SPL, 
                         splineBZ, BZZ=NA, 
                         iter=2000, nburn=1000, nthin=1, 
                         test.idx=NA,
                         knowZ=F) {
  n <- length(Y)

  n.X <- length(unique(X))
  n.theta <- n.X
  
  n.w <- 3
  n.u <- 1
  ########################## set up prior ########################## 
  a <- b <- 0.01
  
  ########################## parse data ##########################
  thetaMat <- matrix(data=0, nrow=n, ncol=n.theta)
  thetaMat[cbind(1:n, X)] <- 1
  
  ## no need to center W1, W2, V though  
  W1 <- W1 - mean(W1, na.rm=T)
  W2 <- W2 - mean(W2, na.rm=T)
  W1[is.na(W1)] <- mean(W1, na.rm=T)
  W2[is.na(W2)] <- mean(W2, na.rm=T)
  WWW <- cbind(W1, W2, V)
  
  ########################## set up outputs ##########################
  theta <- matrix(NA, nrow=iter, ncol=n.theta)
  theta[1, ] <- rep(mean(Y, na.rm=T), n.theta) 
  
  sigma2.y <- matrix(NA, nrow=iter, ncol=1)
  sigma2.y[1] <- sd(Y, na.rm=T)^2
  
  tau2.y <- matrix(NA, nrow=iter, ncol=1)
  tau2.y[1] <- 10^4
  
  alpha.yw <- matrix(NA, nrow=iter, ncol=n.w)
  alpha.yw[1, ] <- 0
  
  ## spatial latent 
  omega <- rep(NA, iter)
  omega[1] <- -1 # estimated from results, do not change. 
  
  sig2.Z <- matrix(NA, nrow=iter, ncol=1)
  sig2.Z[1] <- .1
  
  tau2.omega <- matrix(NA, nrow=iter, ncol=1)
  tau2.omega[1] <- 10
  
  ########################## set up spatial ########################## 
  v.city <- table(SPL)
  n.city <- length(unique(SPL))
  
  ncol.By <- unlist(lapply(splineBZ, function(x) if(is.matrix(x)) ncol(x) else 1))
  ZZ <- matrix(NA, nrow=iter, ncol=sum(ncol.By))
  ZZ[1, ] <- runif(sum(ncol.By))*100 # do not change 
  
  ########################## set up initial (okay to use data information) ########################## 
  
  acc.Z <- rep(0, n.city)
  att.Z <- rep(0, n.city)
  MH.Z <- rep(0.1, n.city)
  
  keep.MH.Z <- matrix(NA, nrow=iter, ncol=n.city)
  accept.Z <- matrix(NA, nrow=iter, ncol=n.city)
  accept.Z[1, ] <- TRUE
  
  acc.omega <- 0
  att.omega <- 0
  MH.omega <- 0.1
  
  keep.MH.omega <- matrix(NA, nrow=iter, ncol=1)
  accept.omega <- matrix(NA, nrow=iter, ncol=1)
  accept.omega[1] <- TRUE
  
  ########################## initial ########################## 
  if (!knowZ) { 
    ## init Z
    ## in a way s.t. fits X ~ spatial(Z) well  
    ## this is bad when #spline >= #region 
    idx0 <- 1; idy0 <- 1
    for (kk in 1:n.city) {
      idx1 <- sum(ncol.By[1:kk])
      idy1 <- sum(v.city[1:kk])
      if (idy1 == idy0) { idx0 <- idx1 + 1; idy0 <- idy1 + 1; next }
      try(
        ZZ[1, idx0:idx1] <- lm(X[idy0:idy1] ~ splineBZ[[kk]])$coefficients[-1]
      )
      
      idx0 <- idx1 + 1; idy0 <- idy1 + 1 
    }
    ZZ[1, is.na(ZZ[1, ])] <- 0 
    BZZ <- compute.spline(splineBZ, ZZ[1, ])
  }
  
  ZZ.mean <- ZZ.M2 <- BZ.mean <- BZ.M2 <- 0
  
  ## mask test data by TRUE
  mask <- rep(FALSE, n)
  if(length(test.idx)>1) {mask[test.idx] <- TRUE}
  idy0 <- 1; mask.city <- list()
  for (kk in 1:n.city) {
    idy1 <- sum(v.city[1:kk])
    mask.city[[kk]] <- mask[idy0:idy1]
    idy0 <- idy1 + 1 }
  
  logpdf1 <- matrix(0, nrow=sum(!mask), ncol=5)
  logpdf2 <- matrix(0, nrow=sum(!mask), ncol=5)
  
  log.score <- matrix(0, nrow=n, ncol=5) # The masked part is true log-score of test data 
  
  ########################################### Loop of Sampling ############################################
  for (i in 2:iter) {
    if(!knowZ) BZZ <- compute.spline(splineBZ, ZZ[i-1, ])

    curll <- matrix(0, nrow=n, ncol=5) # four spatial process proxy1/proxy2/trt/outcome 
    
    if(i%%1000==0) {print(paste(i, timestamp()))}
    
    ##################################################  sampling for X, U -> Y ############################################
    ################ update non-spline part in standard regression ################
    YPart <- Y - omega[i-1]*c(BZZ)
    XMat <- cbind(thetaMat, WWW)
    
    txx <- t(XMat[!mask,])%*%XMat[!mask,]
    txy <- t(XMat[!mask,])%*%YPart[!mask]
    p1 <- n.theta + n.w
    
    VVV  <- solve(txx/sigma2.y[i-1] + diag(p1)/tau2.y[i-1])
    MMM  <- txy/sigma2.y[i-1] + 0/tau2.y[i-1]
    beta1 <- as.vector(VVV%*%MMM+t(chol(VVV))%*%rnorm(p1))
    
    theta[i, ] <- beta1[1:(n.theta)]
    
    alpha.yw[i, ] <- beta1[(n.theta+1):(n.theta+n.w)]
  
    YPart <- Y - XMat%*%beta1 # effect of treatment, latent 
    
    ################ update spline part in standard regression ################
    omega.cand <- omega[i-1] + rnorm(1, 0, MH.omega)
    
    TEMP.prev <- log.dnorm(YPart - omega[i-1]*BZZ, 0, sigma2.y[i-1])
    TEMP.cand <- log.dnorm(YPart - omega.cand*BZZ, 0, sigma2.y[i-1])
    
    lp.prev <- sum(TEMP.prev[!mask]) - 0.5*omega[i-1]^2/tau2.omega[i-1]
    lp.cand <- sum(TEMP.cand[!mask]) - 0.5*omega.cand^2/tau2.omega[i-1]
    
    accept <- (min(1, exp(lp.cand-lp.prev)) > runif(1))
    
    if(i<100) {
      omega[i] <- omega[i-1]
    } else {
      if(accept) {
        omega[i] <- omega.cand 
      } else {
        omega[i] <- omega[i-1]}
    }
    
    accept.omega[i] <- accept 
    keep.MH.omega[i] <- MH.omega
    
    temp <- updateMH(accept, acc.omega, att.omega, MH.omega, i, nburn)
    acc.omega <- temp[[1]]; att.omega <- temp[[2]]; MH.omega <- temp[[3]]
    
    tau2.omega[i] <- 1/rgamma(1,1/2+a,omega[i]^2/2+b)
    
    resid <- YPart-omega[i]*BZZ
    SSS <- sum(resid^2)
    sigma2.y[i] <- 1/rgamma(1, n/2 + a, SSS/2+b)
    tau2.y[i] <- 1/rgamma(1,p1/2+a,sum(beta1^2)/2+b) 
    
    if(accept) {
      curll[, 4] <- TEMP.cand
    } else {
      curll[, 4] <- TEMP.prev }

    ##################################################  sampling for Z -> X ##################################################  
    out <- matrix(NA, nrow=length(BZZ), ncol=n.X)
    l1 <- BZZ
    l1[l1 < -10] = -10; l1[l1 > 10] = 10
    out[, 1] <- exp(l1) / (exp(l1) + 1); out[, 2] <- (1 - out[, 1])
    curll[, 3] <- out[cbind(1:n, X)]
    
    ################################################## sampling Z ##################################################
    if (!knowZ) {
      idx0 <- 1
      idy0 <- 1
      for(j in 1:n.city) {
        Bk <- splineBZ[[j]]
        
        idx1 <- idx0 + ncol(Bk) - 1 # index of splines 
        idy1 <- sum(v.city[1:j]) # index of samples 
        
        xk <- X[idy0:idy1]
        yk <- Y[idy0:idy1]
        thetamatk <- thetaMat[idy0:idy1, ]
        www <- WWW[idy0:idy1, ] 
        q <- nrow(Bk)
        
        ZZ.curr <- ZZ[i-1, idx0:idx1]
        ZZ.cand <- oneCol(ZZ.curr + rnorm(ncol(Bk), 0, MH.Z[j]))
        
        ####### lp.prev, lp.cand, treatment part and outcome part ####### 
        bz      <- Bk %*% ZZ.curr
        bz.cand <- Bk %*% ZZ.cand
        
        out <- matrix(NA, nrow=length(bz), ncol=n.X)
        l1 <- bz
        l1[l1 < -10] = -10; l1[l1 > 10] = 10
        out[, 1] <- exp(l1) / (exp(l1) + 1); out[, 2] <- (1 - out[, 1])
        temp.prev <- out[cbind(1:q, xk)]
        
        out <- matrix(NA, nrow=length(bz.cand), ncol=n.X)
        l1 <- bz.cand
        l1[l1 < -10] = -10; l1[l1 > 10] = 10
        out[, 1] <- exp(l1) / (exp(l1) + 1); out[, 2] <- (1 - out[, 1])
        temp.cand <- out[cbind(1:q, xk)]
        
        yTEMP <- thetamatk %*% oneCol(theta[i, ]) + www %*% oneCol(alpha.yw[i, ])
        
        temp.prev <- temp.prev + log.dnorm(yk, yTEMP + omega[i]*c(bz     ), sigma2.y[i])
        temp.cand <- temp.cand + log.dnorm(yk, yTEMP + omega[i]*c(bz.cand), sigma2.y[i])
        
        lp.prev <- sum(temp.prev[!mask.city[[j]]]) - 0.5*sum(ZZ.curr^2)/sig2.Z[i-1]
        lp.cand <- sum(temp.cand[!mask.city[[j]]]) - 0.5*sum(ZZ.cand^2)/sig2.Z[i-1]
        
        ############################
        accept <- (min(1, exp(lp.cand-lp.prev)) > runif(1))
        
        if(accept) {
          ZZ[i, idx0:idx1] <- ZZ.cand 
        } else {
          ZZ[i, idx0:idx1] <- ZZ.curr }
        
        accept.Z[i, j] <- accept 
        keep.MH.Z[i, j] <- MH.Z[j] 
        
        temp <- updateMH(accept, acc.Z[j], att.Z[j], MH.Z[j], i, nburn)
        acc.Z[j] <- temp[[1]]; att.Z[j] <- temp[[2]]; MH.Z[j] <- temp[[3]]
        
        idx0 <- idx1 + 1
        idy0 <- idy1 + 1
      } # end for(j in 1:n.city)
      
      sig2.Z[i] <- 1/rgamma(1, sum(ncol.By)/2+a, sum(ZZ[i, ]^2)/2+b)
      if(i<100) sig2.Z[i] <- sig2.Z[i-1] # fix it for 100 times for stability 
    }
    
    ######################################## post-burn ######################################################
    if(i>nburn){
      
      ZZ.mean <- ZZ.mean + ZZ[i,]/(iter-nburn)
      ZZ.M2   <- ZZ.M2 + ZZ[i,]^2/(iter-nburn)
      
      BZ.mean <- BZ.mean + compute.spline(splineBZ, ZZ[i,])/(iter-nburn)
      BZ.M2   <- BZ.M2 + compute.spline(splineBZ, ZZ[i,])^2/(iter-nburn)
      
      curll[, 5] <- rowSums(curll[, 1:4])
      
      logpdf1 <- logpdf1 + (curll[!mask,]  )/(iter-nburn) # MEAN of curll after burn-in, n-by-1
      logpdf2 <- logpdf2 + (curll[!mask,]^2)/(iter-nburn) # VAR  of curll after burn-in, n-by-1
      
      log.score <- log.score + curll/(iter-nburn)
    }
  } 
  
  # WAIC computations
  # the last one is for the whole mode 
  mn_logpdf  <- logpdf1
  var_logpdf <- logpdf2 - logpdf1^2
  pW         <- colSums(var_logpdf)
  WAIC       <- list(WAIC=-2*colSums(mn_logpdf)+2*pW,pW=pW)
  
  keep <- seq(1, iter, by=nthin)
  
  out <- list(iteration=keep, 
              post=1*(keep>nburn),
              alpha.yw=alpha.yw[keep, ], 
              sigma2.y=sigma2.y[keep], 
              theta=theta[keep, 1:n.theta])
  
  # out[['ZZ']] <- ZZ[keep, ] # to save storage
  out[['sig2.Z']] <- sig2.Z[keep]
  out[['omega']] <- omega[keep]

  ZZ.var <- ZZ.M2 - ZZ.mean^2
  BZ.var <- BZ.M2 - BZ.mean^2
  
  return(list(out=as.data.frame(out), 
              WAIC=WAIC, log.score=log.score, 
              mask=mask, mask.city=mask.city, 
              mn_logpdf=mn_logpdf, var_logpdf=var_logpdf, 
              ZZ.mean=ZZ.mean, ZZ.var=ZZ.var, 
              BZ.mean=BZ.mean, BZ.var=BZ.var)) 
}
