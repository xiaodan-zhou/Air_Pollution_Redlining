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
                         splineBU, splineBZ, n.u, 
                         BUU=NA, BZZ=NA, 
                         iter=2000, nburn=1000, nthin=1, 
                         no.proxy=F, # no.proxy override n.u, splineBU etc
                         test.idx=NA, # for CV
                         knowU=F, knowZ=F # for testing
) {
  ########################## note ################################# 
  # when splineBU == I, non-spatial model 
  # when splineBU == 0, latent not exist 
  # splineBU %*% UU -> BUU
  # splineBZ %*% ZZ -> BZZ
  # to force no latent, set coefficients to zero
  #####
  
  n <- length(Y)
  n.theta <- n.X <- length(unique(X))
  if(no.proxy) n.u <- 1 
  
  ## TODO force n.u == 1 here? yes. 
  
  a <- b <- 0.01
  
  ########################## parse data ##########################
  thetaMat <- matrix(data=0, nrow=n, ncol=n.theta)
  thetaMat[cbind(1:n, X)] <- 1
  
  ## fill a few missing by mean values 
  W1 <- W1 - mean(W1, na.rm=T)
  W2 <- W2 - mean(W2, na.rm=T)
  W1[is.na(W1)] <- mean(W1, na.rm=T)
  W2[is.na(W2)] <- mean(W2, na.rm=T)
  
  ########################## set up outputs ##########################
  ### Latent ### 
  sigma2.u <- matrix(NA, nrow=iter, ncol=n.u)
  sigma2.u[1, ] <- .5 # matter! 
  
  ### Proxy Normal ###
  alpha.ww1 <- matrix(NA, nrow=iter, ncol=1+n.u)
  sigma2.w1 <- matrix(NA, nrow=iter, ncol=1)
  tau2.w1 <- matrix(NA, nrow=iter, ncol=1)
  
  alpha.ww2 <- matrix(NA, nrow=iter, ncol=1+n.u)
  sigma2.w2 <- matrix(NA, nrow=iter, ncol=1)
  tau2.w2 <- matrix(NA, nrow=iter, ncol=1)
  
  alpha.ww1[1, ] <- 0
  alpha.ww1[1, 2] <- 1 # don't change this, reparam 
  sigma2.w1[1] <- sd(W1, na.rm=T)^2
  tau2.w1[1] <- 10^4
  
  alpha.ww2[1, ] <- 0
  sigma2.w2[1] <- sd(W2, na.rm=T)^2
  tau2.w2[1] <- 10^4
  
  ### Proxy TOBIT ###
  alpha.vv <- matrix(NA, nrow=iter, ncol=1+n.u)
  sigma2.v <- matrix(NA, nrow=iter, ncol=1)
  tau2.v <- matrix(NA, nrow=iter, ncol=1)

  alpha.vv[1, ] <- c(0, rep(1, n.u))
  sigma2.v[1] <- sd(V, na.rm=T)^2
  tau2.v[1] <- 10^4
  
  ### Outcome ###
  theta <- matrix(NA, nrow=iter, ncol=n.theta)
  theta[1, ] <- rep(mean(Y, na.rm=T), n.theta) 
  
  alpha.yu <- matrix(NA, nrow=iter, ncol=n.u)
  alpha.yu[1, ] <- 0
  
  sigma2.y <- matrix(NA, nrow=iter, ncol=1)
  sigma2.y[1] <- sd(Y, na.rm=T)^2
  
  tau2.y <- matrix(NA, nrow=iter, ncol=1)
  tau2.y[1] <- 10^4
  
  ## spatial latent 
  omega <- rep(NA, iter)
  omega[1] <- -1 # sign estimated from results, do not change, otherwise stable. 
  
  sig2.Z <- matrix(NA, nrow=iter, ncol=1)
  sig2.Z[1] <- .1
  
  tau2.omega <- matrix(NA, nrow=iter, ncol=1)
  tau2.omega[1] <- 10
  
  ### Treatment ###
  eta <- matrix(NA, nrow=iter, ncol=(1+n.u)*(n.X-1))
  sig2.eta <- matrix(NA, nrow=iter, ncol=1)
  
  eta[1, ] <- rep(0, (n.u+1)*(n.X-1)) # do not change this 
  sig2.eta[1] <- .1
  
  ########################## set up spatial ########################## 
  spline.u <- (is.list(splineBU)) # always TRUE 
  spline.z <- (is.list(splineBZ)) # always TRUE 
  ## if no latent, up a all-zero splineBU and aplineBZ 
  
  v.city <- table(SPL)
  n.city <- length(unique(SPL))

  ncol.B <- unlist(lapply(splineBU, function(x) if(is.matrix(x)) ncol(x) else 1)) # the number of columns in each spline basis
  UU.full <- matrix(rnorm(sum(ncol.B)*n.u, 0, sd=1), ncol=n.u)
  UU <- as.list(split.data.frame(UU.full, rep(1:n.city, ncol.B))) 
  
  ncol.By <- unlist(lapply(splineBZ, function(x) if(is.matrix(x)) ncol(x) else 1))
  ZZ <- matrix(NA, nrow=iter, ncol=sum(ncol.By))
  ZZ[1, ] <- runif(sum(ncol.By))*100 # do not change 
  
  ########################## set up MH ########################## 
  ## initial is not important for those with conjugate prior, proxy and outcome model 
  ## but important for the assignment model.
  
  ### Latent ### 
  acc.u <- matrix(0, nrow=n.u, ncol=n.city)
  att.u <- matrix(0, nrow=n.u, ncol=n.city)
  MH.u <- matrix(0.1, nrow=n.u, ncol=n.city)
  
  acc.sig2u <- rep(0, n.u)
  att.sig2u <- rep(0, n.u)
  MH.sig2u <- rep(.1, n.u)
  
  ### Proxy ### 
  acc.v <- rep(0, 1+n.u)
  att.v <- rep(0, 1+n.u)
  MH.v <- rep(.1, 1+n.u)
  
  keep.MH.v <- matrix(NA, nrow=iter, ncol=1+n.u)
  accept.v <- matrix(NA, nrow=iter, ncol=1+n.u)
  accept.v[1, ] <- TRUE
  
  ### TOBIT ###
  thr.v <- min(V, na.rm=T) # threshold where data get truncated 
  idx.miss <- which(is.na(V))
  idx.obs <- which(!is.na(V)) # idx.obs idx.o
  idx.above <- which((!is.na(V)) & (V > thr.v))
  idx.below <- which((!is.na(V)) & (V <= thr.v))
  ## check 
  # sum(idx.miss, idx.above, idx.below) - sum(1:n) == 0
  
  ### Treatment ### 
  acc.eta <- rep(0, (1+n.u)*(n.X-1))
  att.eta <- rep(0, (1+n.u)*(n.X-1))
  MH.eta <- rep(0.1, (1+n.u)*(n.X-1))
  
  keep.MH.eta <- matrix(NA, nrow=iter, ncol=(1+n.u)*(n.X-1))
  accept.eta <- matrix(NA, nrow=iter, ncol=(1+n.u)*(n.X-1))
  accept.eta[1, ] <- TRUE 
  
  ### Outcome ### 
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
  if (!knowU) { 
    if (spline.u) {
      ## init UU
      ## in a way s.t. fits W1 ~ U1 well
      idy0 <- 1
      for (kk in 1:n.city) {
        idy1 <- sum(v.city[1:kk])
        temp <- lm(W1[idy0:idy1] ~ splineBU[[kk]])$coefficients
        UU[[kk]][,1] <- temp[-1]
        idy0 <- idy1 + 1 }
      
      ### ignore this 
      if (n.u==2) {
        ## init UU 
        ## in a way s.t. fits V ~ U2 well
        idy0 <- 1
        for (kk in 1:n.city) {
          idy1 <- sum(v.city[1:kk])
          temp <- lm(V[idy0:idy1] ~ splineBU[[kk]])$coefficients
          UU[[kk]][,2] <- temp[-1]
          idy0 <- idy1 + 1 }
      }
    }

    for (kk in 1:n.city) {
      for (jj in 1:n.u) {
        temp <- is.na(UU[[kk]][,jj])
        UU[[kk]][temp,jj] <- rnorm(n=sum(temp), 0, 1) }}
    
    UU.full <- do.call(rbind, UU)
    
    BUU <- matrix(data=NA, nrow=n, ncol=n.u)
      for (jj in 1:n.u) {
        BUU[,jj] <- compute.spline(splineBU, UU.full[,jj]) 
      }
    
    if(no.proxy) BUU[] <- 0 
  }
  if (!knowZ) { 
    ## init ZZ
    ## in a way s.t. fits X ~ spatial(U + ZZ) well  
    ## this is bad when #spline >= #region 
    idx0 <- 1; idy0 <- 1
    for (kk in 1:n.city) {
      idx1 <- sum(ncol.By[1:kk])
      idy1 <- sum(v.city[1:kk])
      if (idy1 == idy0) { idx0 <- idx1 + 1; idy0 <- idy1 + 1; next }

      try(
        if(spline.u) {
          temp <- glm(X[idy0:idy1] ~ splineBU[[kk]], family="poisson")$residuals
          ZZ[1, idx0:idx1] <- lm(temp ~ splineBZ[[kk]])$coefficients[-1]
        } else {
          ZZ[1, idx0:idx1] <- lm(X[idy0:idy1] ~ splineBZ[[kk]])$coefficients[-1]
        }
      )
      
      idx0 <- idx1 + 1; idy0 <- idy1 + 1 
    }
    ZZ[1, is.na(ZZ[1, ])] <- 0 
    BZZ <- compute.spline(splineBZ, ZZ[1, ])
  
  } ## end knowU
  
  ## init eta 
  ## such that treatment model fits well
  temp <- cbind(BUU, BZZ)
  eta[1, ] <- - glm(X ~ temp, family="poisson")$coefficients[-(ncol(temp)+1)]  
  
  if(no.proxy) eta[1, ] <- 0
  
  UU.mean <- UU.M2 <- BU.mean <- BU.M2 <- 0
  ZZ.mean <- ZZ.M2 <- BZ.mean <- BZ.M2 <- 0
  
  ## Cross Validation: mask test data by TRUE
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
    
    curll <- matrix(0, nrow=n, ncol=5) # four spatial process proxy/proxy/trt/outcome 
    
    if(i%%1000==0) {print(paste(i, timestamp()))
      # print(ggplot() +
      #         geom_histogram(aes(BZZ[X==1]),col="red",alpha=.5) +
      #         geom_histogram(aes(BZZ[X==2]),col="blue",alpha=.5))
      # plot(theta[,2]-theta[,1])
    }
    
    ########################## sampling UU ##########################
    if (!knowU) {
      ## U loop 
      for (j in 1:n.u) {
        ## get all U.cand 
        UU.cand  <- UU
        BUU.cand <- BUU
        
        idy0 <- 1
        for(kk in 1:n.city) {
          idy1 <- sum(v.city[1:kk]) # index of samples 
      
          UU.cand[[kk]][,j] <- UU[[kk]][,j] + rnorm(ncol.B[kk], 0, MH.u[j, kk]) 
          BUU.cand[idy0:idy1,j] <- compute.spline(splineBU[kk], UU.cand[[kk]][,j])
    
          idy0 <- idy1 + 1
        }
        
        ## P(Y|U) 
        yTEMP <- thetaMat %*% oneCol(theta[i-1, ]) + omega[i-1]*BZZ
        lp.prev.y <- log.dnorm(Y, yTEMP + BUU      %*% oneCol(alpha.yu[i-1, ]), sigma2.y[i-1])
        lp.cand.y <- log.dnorm(Y, yTEMP + BUU.cand %*% oneCol(alpha.yu[i-1, ]), sigma2.y[i-1])
        
        ## P(X|U) 
        lp.prev.x <- eta2lp(eta[i-1, ], BUU,      n.X, n.u, extra=BZZ)[cbind(1:n, X)]
        lp.cand.x <- eta2lp(eta[i-1, ], BUU.cand, n.X, n.u, extra=BZZ)[cbind(1:n, X)]
        
        ## P(W|U)
        lp.prev.w1 <- log.dnorm(W1, alpha.ww1[i-1, 1] + BUU      %*% oneCol(alpha.ww1[i-1, -1]), sigma2.w1[i-1]) 
        lp.cand.w1 <- log.dnorm(W1, alpha.ww1[i-1, 1] + BUU.cand %*% oneCol(alpha.ww1[i-1, -1]), sigma2.w1[i-1])    
        lp.prev.w2 <- log.dnorm(W2, alpha.ww2[i-1, 1] + BUU      %*% oneCol(alpha.ww2[i-1, -1]), sigma2.w2[i-1]) 
        lp.cand.w2 <- log.dnorm(W2, alpha.ww2[i-1, 1] + BUU.cand %*% oneCol(alpha.ww2[i-1, -1]), sigma2.w2[i-1])    
        
        ## P(V|U) 
        lp.prev.v <- tobit.loglik(BUU,      V, thr.v, alpha.vv[i-1, ], sigma2.v[i-1], n, idx.obs, idx.above, idx.below)
        lp.cand.v <- tobit.loglik(BUU.cand, V, thr.v, alpha.vv[i-1, ], sigma2.v[i-1], n, idx.obs, idx.above, idx.below)
        
        ## sampling U, we don't record curll nor log.score
        ## avoid using test data to determine the model 
        ## simply mask test data with zero 
        lp.prev.y[mask] <- lp.cand.y[mask] <- 0
        lp.prev.x[mask] <- lp.cand.x[mask] <- 0
        lp.prev.w1[mask] <- lp.cand.w1[mask] <- 0
        lp.prev.w2[mask] <- lp.cand.w2[mask] <- 0
        lp.prev.v[mask] <- lp.cand.v[mask] <- 0
        
        idy0 <- 1
        for(kk in 1:n.city) {
          idy1 <- sum(v.city[1:kk]) # index of samples 
          
          ## P(U|prior)
          lp.prev <- sum(lp.prev.y[idy0:idy1] + lp.prev.x[idy0:idy1] + lp.prev.w1[idy0:idy1] + lp.prev.w2[idy0:idy1] + lp.prev.v[idy0:idy1]) - 0.5*sum(UU[[kk]][,j]     ^2)/sigma2.u[i-1, j]
          lp.cand <- sum(lp.cand.y[idy0:idy1] + lp.cand.x[idy0:idy1] + lp.cand.w1[idy0:idy1] + lp.cand.w2[idy0:idy1] + lp.cand.v[idy0:idy1]) - 0.5*sum(UU.cand[[kk]][,j]^2)/sigma2.u[i-1, j]
          
          ######### accept for i, j, kk 
          accept <- (min(1, exp(lp.cand-lp.prev)) > runif(1))
          if(accept) { UU[[kk]][,j] <- UU.cand[[kk]][,j] }
          
          temp <- updateMH(accept, acc.u[j, kk], att.u[j, kk], MH.u[j, kk], i, nburn)
          acc.u[j, kk] <- temp[[1]]; att.u[j, kk] <- temp[[2]]; MH.u[j, kk] <- temp[[3]]
          
          idy0 <- idy1 + 1
        } # end for(kk in 1:n.city) 
        
        UU.full <- do.call(rbind, UU)
        if(spline.u) {
          BUU[,j] <- compute.spline(splineBU, UU.full[,j])
        } else {
          BUU[,j] <- UU.full[,j]
        }
        sigma2.u[i, j] <- 1/rgamma(1,sum(ncol.B)/2+a,sum(UU.full[,j]^2)/2+b)
        
      } # end for(j in 1:n.u)
    } # end if (knowU)
    ## end U loop 
    
    if(no.proxy) BUU[] <- 0 # key 
    
    #################################################### sampling for U -> W normal ####################################################
    ## todo need to add if(no.proxy) alpha.ww[] <- alpha.ww.cand[] <- 0 
    if(n.u == 1) {
      ############ W1 ############
      ### intercept 
      WTEMP <- W1 - BUU %*% oneCol(alpha.ww1[i-1,-1])
      
      VVV   <- n/sigma2.w1[i-1] + 1/tau2.w1[i-1]  
      MMM   <- sum(WTEMP)/sigma2.w1[i-1] + 0/tau2.w1[i-1]
      alpha.ww1[i, 1] <- rnorm(1,MMM/VVV,1/sqrt(VVV))
      
      alpha.ww1[i, 2] <- 1 ## reparam
      
      ### variance
      resid <- W1 - alpha.ww1[i,1] - BUU %*% oneCol(alpha.ww1[i,-1])
      SSE <- sum(resid^2)
      sigma2.w1[i] <- 1/rgamma(1, n/2+a, SSE/2+b)
      tau2.w1[i] <- 1/rgamma(1, 1/2+a, sum(alpha.ww1[i, 1]^2)/2+b)
      
      curll[, 1] <- curll[, 1] + log.dnorm(resid,0,sigma2.w1[i])
      
      ############ W2 ############
      XTEMP <- cbind(1, BUU)
      WTEMP <- W2

      txx <- t(XTEMP)%*%XTEMP
      txy <- t(XTEMP)%*%WTEMP 
      p <- ncol(XTEMP)
      
      VVV  <- solve(txx/sigma2.w2[i-1] + diag(p)/tau2.w2[i-1])
      MMM  <- txy/sigma2.w2[i-1] + 0/tau2.w2[i-1] 
      beta <- as.vector(VVV%*%MMM+t(chol(VVV))%*%rnorm(p))

      alpha.ww2[i, ] <- beta
      
      resid <- WTEMP - XTEMP%*%beta
      SSE <- sum(resid^2)
      sigma2.w2[i] <- 1/rgamma(1, n/2+a, SSE/2+b)  
      tau2.w2[i] <- 1/rgamma(1, p/2+a, sum(alpha.ww2[i, ]^2)/2+b)
      
      curll[, 1] <- curll[, 1] + log.dnorm(resid,0,sigma2.w2[i]) 
    
      ## ignore this 
    } else if (n.u == 2) {
      alpha.ww1[i, 2] <- 1 ## reparam
      
      XTEMP <- cbind(1, BUU[,-1])
      WTEMP <- W1 - BUU[,1]*alpha.ww1[i,2]
      
      txx <- t(XTEMP)%*%XTEMP
      txy <- t(XTEMP)%*%WTEMP 
      p <- ncol(XTEMP)
      
      VVV  <- solve(txx/sigma2.w1[i-1] + diag(p)/tau2.w1[i-1])
      MMM  <- txy/sigma2.w1[i-1] + 0/tau2.w1[i-1] 
      beta <- as.vector(VVV%*%MMM+t(chol(VVV))%*%rnorm(p))
      
      alpha.ww1[i,-2] <- beta
      
      resid <- WTEMP - XTEMP%*%beta
      SSE <- sum(resid^2)
      sigma2.w1[i] <- 1/rgamma(1, n/2+a, SSE/2+b)  
      tau2.w1[i] <- 1/rgamma(1, p/2+a, sum(beta^2)/2+b)
      
      curll[, 1] <- curll[, 1] + log.dnorm(resid,0,sigma2.w1[i])
      
      ############
      XTEMP <- cbind(1, BUU)
      WTEMP <- W2
      
      txx <- t(XTEMP)%*%XTEMP
      txy <- t(XTEMP)%*%WTEMP 
      p <- ncol(XTEMP)
      
      VVV  <- solve(txx/sigma2.w2[i-1] + diag(p)/tau2.w2[i-1])
      MMM  <- txy/sigma2.w2[i-1] + 0/tau2.w2[i-1] 
      beta <- as.vector(VVV%*%MMM+t(chol(VVV))%*%rnorm(p))
      
      alpha.ww2[i, ] <- beta
      
      resid <- WTEMP - XTEMP%*%beta
      SSE <- sum(resid^2)
      sigma2.w2[i] <- 1/rgamma(1, n/2+a, SSE/2+b)  
      tau2.w2[i] <- 1/rgamma(1, p/2+a, sum(beta^2)/2+b)
      
      curll[, 1] <- curll[, 1] + log.dnorm(resid,0,sigma2.w2[i]) 
    }
    
    #################################################### sampling for U -> V TOBIT ####################################################
    
    ################ update beta part ################ 
    ## use att.v[1] acc.v[1] MH.v[1] accept.v[1] keep.MH.v[1]
    alpha.vv.cand <- alpha.vv[i-1, ] + rnorm(length(alpha.vv[i-1, ]), 0, MH.v[1])
    
    ## ignore this 
    if(n.u == 2) {alpha.vv.cand[3] <- 1 } ## reparam
    
    TEMP.prev <- tobit.loglik(BUU, V, thr.v, alpha.vv[i-1, ], sigma2.v[i-1], n, idx.obs, idx.above, idx.below)
    TEMP.cand <- tobit.loglik(BUU, V, thr.v, alpha.vv.cand,   sigma2.v[i-1], n, idx.obs, idx.above, idx.below)
    
    lp.prev <- sum(TEMP.prev[!mask]) - 0.5*sum(alpha.vv[i-1, ]^2)/sigma2.v[i-1]  
    lp.cand <- sum(TEMP.cand[!mask]) - 0.5*sum(alpha.vv.cand  ^2)/sigma2.v[i-1] 
    
    accept <- (min(1, exp(lp.cand-lp.prev)) > runif(1))
    accept.v[i, 1] <- accept
    
    keep.MH.v[i, 1] <- MH.v[1]
    
    temp <- updateMH(accept, acc.v[1], att.v[1], MH.v[1], i, nburn)
    acc.v[1] <- temp[[1]]; att.v[1] <- temp[[2]];  MH.v[1] <- temp[[3]]
    
    if(accept) {
      alpha.vv[i, ] <- alpha.vv.cand
    } else {
      alpha.vv[i, ] <- alpha.vv[i-1, ]}
    
    ################ update variance part ################   
    ## use att.v[2] acc.v[2] MH.v accept.v[2] keep.MH.v[2]
    
    sigma2.v.cand <- rlnorm(1, meanlog = log(sigma2.v[i-1]), sdlog = MH.v[2])
    
    TEMP.prev <- tobit.loglik(BUU, V, thr.v, alpha.vv[i, ], sigma2.v[i-1], n, idx.obs, idx.above, idx.below)
    TEMP.cand <- tobit.loglik(BUU, V, thr.v, alpha.vv[i, ], sigma2.v.cand, n, idx.obs, idx.above, idx.below)
    
    lp.prev <- sum(TEMP.prev[!mask]) + dgamma(1/sigma2.v[i-1], a, b, log=T)
    lp.cand <- sum(TEMP.cand[!mask]) + dgamma(1/sigma2.v.cand, a, b, log=T)
    
    accept <- (min(1, exp(lp.cand-lp.prev)) > runif(1))
    accept.v[i, 2] <- accept
    
    keep.MH.v[i, 2] <- MH.v[2]
    
    temp <- updateMH(accept, acc.v[2], att.v[2], MH.v[2], i, nburn)
    acc.v[2] <- temp[[1]]; att.v[2] <- temp[[2]];  MH.v[2] <- temp[[3]]
    
    if(accept) {
      sigma2.v[i] <- sigma2.v.cand
    } else {
      sigma2.v[i] <- sigma2.v[i-1]}
    
    if(accept) {
      curll[, 2] <- TEMP.cand
    } else {
      curll[, 2] <- TEMP.prev }
    
    
    ##################################################  sampling for X, U -> Y ############################################
    if(spline.z) {
      ################ update non-spline part in standard regression ###############
      YPart <- Y - omega[i-1]*c(BZZ)
      XMat <- cbind(thetaMat, BUU)
      
      txx <- t(XMat[!mask,])%*%XMat[!mask,]
      txy <- t(XMat[!mask,])%*%YPart[!mask]
      p1 <- n.theta + n.u
      
      VVV  <- solve(txx/sigma2.y[i-1] + diag(p1)/tau2.y[i-1])
      MMM  <- txy/sigma2.y[i-1] + 0/tau2.y[i-1]
      beta1 <- as.vector(VVV%*%MMM+t(chol(VVV))%*%rnorm(p1))
      
      theta[i, ] <- beta1[1:(n.theta)]
      
      alpha.yu[i, ] <- beta1[(n.theta+1):(n.theta+n.u)]
      if(no.proxy) alpha.yu[i, ] <- 0 # just in case
      
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
      tau2.y[i] <- 1/rgamma(1,p1/2+a,sum(beta1^2)/2+b) ## maybe b1, beta1 should but cut when no.proxy... 
      
      if(accept) {
        curll[, 4] <- TEMP.cand
      } else {
        curll[, 4] <- TEMP.prev }
      
      ## ignore this 
    } else {
      ################ update non-spline part in standard regression ################
      XMat <- cbind(thetaMat, BUU)
      
      txx <- t(XMat)%*%XMat
      txy <- t(XMat)%*%Y
      p <- n.theta + n.u
      
      VVV  <- solve(txx/sigma2.y[i-1] + diag(p)/tau2.y[i-1])
      MMM  <- txy/sigma2.y[i-1] + 0/tau2.y[i-1]
      beta <- as.vector(VVV%*%MMM+t(chol(VVV))%*%rnorm(p))
      
      theta[i, ] <- beta[1:(n.theta)]
      
      ## just in case 
      alpha.yu[i, ] <- beta1[(n.theta+1):(n.theta+n.u)]
      if(no.proxy) alpha.yu[i, ] <- 0

      ################ update variance ################
      Yfit <- XMat%*%beta
      resid <- Y- Yfit
      SSS <- sum(resid^2)
      sigma2.y[i] <- 1/rgamma(1, n/2+a, SSS/2+b)
      tau2.y[i] <- 1/rgamma(1,p/2+a,sum(beta^2)/2+b)
      
      # curll[, 4] <- log.dnorm(resid,0,sigma2.y[i]) ## maybe b1, beta1 should but cut when no.proxy... 
    }  
    
    ##################################################  sampling for U -> X ##################################################  
    eta[i, ] <- eta[i-1, ] # eta init as zero and BUU are all zero if(no.proxy)
    
    ################ update non-spline part ################
    TEMP <- eta2lp(eta[i-1, ], BUU, n.X, n.u, extra=BZZ)[cbind(1:n, X)]
    
    lp.prev <- sum(TEMP[!mask]) - 0.5*sum(eta[i-1, ]^2)/sig2.eta[i-1]
    
    for(j in 1:((1+n.u)*(n.X-1))){
      eta.cand    <- eta[i, ]
      eta.cand[j] <- eta[i, j] + rnorm(1, 0, MH.eta[j])
      
      TEMP.cand <- eta2lp(eta.cand, BUU, n.X, n.u, extra=BZZ)[cbind(1:n, X)]
      lp.cand <- sum(TEMP.cand[!mask]) - 0.5*sum(eta.cand^2)/sig2.eta[i-1] 
      
      accept <- (min(1, exp(lp.cand-lp.prev)) > runif(1))
      
      if(accept) {
        eta[i, j] <- eta.cand[j]
        lp.prev <- lp.cand 
        TEMP <- TEMP.cand }
      
      accept.eta[i, j] <- accept
      keep.MH.eta[i, j] <- MH.eta[j]
      
      temp <- updateMH(accept, acc.eta[j], att.eta[j], MH.eta[j], i, nburn)
      acc.eta[j] <- temp[[1]]; att.eta[j] <- temp[[2]]; MH.eta[j] <- temp[[3]]
    }
    
    if(no.proxy) eta[i, 2:(1+n.u)] <- 0 ## just in case
    
    ################ update variance ################
    sig2.eta[i] <- 1/rgamma(1,length(eta[i, ])/2+a,sum(eta[i, ]^2)/2+b)
    if(i<100) sig2.eta[i] <- sig2.eta[i-1] # fix it for 100 times for stability 
    
    curll[, 3] <- eta2lp(eta[i, ], BUU, n.X, n.u, extra=BZZ)[cbind(1:n, X)] 
    
    ################################################## sampling ZZ ##################################################
    ## current model does not allow spline.z=F 
    if (!knowZ) {
      if (spline.z) {
        idx0 <- 1
        idy0 <- 1
        for(j in 1:n.city) {
          Bk <- splineBZ[[j]]
          
          idx1 <- idx0 + ncol(Bk) - 1 # index of splines 
          idy1 <- sum(v.city[1:j]) # index of samples 
          
          buu <- BUU[idy0:idy1, ] 
          xk <- X[idy0:idy1]
          yk <- Y[idy0:idy1]
          thetamatk <- thetaMat[idy0:idy1, ]
          
          q <- nrow(Bk)
          
          ZZ.curr <- ZZ[i-1, idx0:idx1]
          ZZ.cand <- oneCol(ZZ.curr + rnorm(ncol(Bk), 0, MH.Z[j]))
          
          ####### lp.prev, lp.cand, treatment part and outcome part ####### 
          bz      <- Bk %*% ZZ.curr
          bz.cand <- Bk %*% ZZ.cand
          
          temp.prev <- eta2lp(eta[i, ], buu, n.X, n.u, extra=bz     )[cbind(1:q, xk)]
          temp.cand <- eta2lp(eta[i, ], buu, n.X, n.u, extra=bz.cand)[cbind(1:q, xk)]
          
          yTEMP <- thetamatk %*% oneCol(theta[i, ]) + buu %*% oneCol(alpha.yu[i, ])
          
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
        
      } # end if (spline.u)
    }
    
    ######################################## post-burn ######################################################
    if(i>nburn){
      
      UU.mean <- UU.mean + UU.full/(iter-nburn)
      UU.M2   <- UU.M2 + UU.full^2/(iter-nburn)
      
      BU.mean <- BU.mean + BUU/(iter-nburn)
      BU.M2   <- BU.M2 + BUU^2/(iter-nburn)
      
      ZZ.mean <- ZZ.mean + ZZ[i,]/(iter-nburn)
      ZZ.M2   <- ZZ.M2 + ZZ[i,]^2/(iter-nburn)
      
      BZ.mean <- BZ.mean + compute.spline(splineBZ, ZZ[i,])/(iter-nburn)
      BZ.M2   <- BZ.M2 + compute.spline(splineBZ, ZZ[i,])^2/(iter-nburn)
      
      curll[, 5] <- rowSums(curll[, 1:4])
      if(no.proxy) curll[, 5] <- rowSums(curll[, 3:4])
      
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
              alpha.w1=alpha.ww1[keep, 1], 
              alpha.w1u=alpha.ww1[keep, -1], 
              sigma2.w1=sigma2.w1[keep], 
              alpha.w2=alpha.ww2[keep, 1], 
              alpha.w2u=alpha.ww2[keep, -1], 
              sigma2.w2=sigma2.w2[keep], 
              eta=eta[keep, ], 
              sig2.eta=sig2.eta[keep], 
              alpha.yu=alpha.yu[keep, ], 
              sigma2.y=sigma2.y[keep], 
              theta=theta[keep, 1:n.theta],
              sigma2.u=sigma2.u[keep, ],
              alpha.v = alpha.vv[keep, 1],
              alpha.vu=alpha.vv[keep, -1],
              sigma2.v=sigma2.v[keep])
  
  if(spline.u) {
    # out[['ZZ']] <- ZZ[keep, ] # to save storage
    out[['sig2.Z']] <- sig2.Z[keep]
    out[['omega']] <- omega[keep]
  }
  
  UU.var <- UU.M2 - UU.mean^2
  BU.var <- BU.M2 - BU.mean^2
  ZZ.var <- ZZ.M2 - ZZ.mean^2
  BZ.var <- BZ.M2 - BZ.mean^2
  
  return(list(out=as.data.frame(out), 
              WAIC=WAIC, log.score=log.score, 
              mask=mask, mask.city=mask.city, 
              mn_logpdf=mn_logpdf, var_logpdf=var_logpdf, 
              UU.mean=UU.mean, UU.var=UU.var, 
              BU.mean=BU.mean, BU.var=BU.var, 
              ZZ.mean=ZZ.mean, ZZ.var=ZZ.var, 
              BZ.mean=BZ.mean, BZ.var=BZ.var)) 
}
