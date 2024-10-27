library(stringr)
library(ggplot2)
library(sf)

LOG2PI <- log(2*pi)

oneCol <- function(x) {return(matrix(x, ncol=1))}
oneRow <- function(x) {return(matrix(x, nrow=1))}

log.dnorm <- function(x, mu, sd2) { -.5*(LOG2PI+(x-mu)^2/sd2+log(sd2))}

etaVtoM <- function(et, n.X) {matrix(et, nrow=n.X-1, byrow=T)}

eta2lp <- function(et, u, n.X, n.u, extra=NULL) {
  ## convert eta to the probability vector for multinomial distribution given U
  ## input et is the six eta parameters
  ## input U can be scalar or vector, n-by-?
  ## output is of n-by-?
  
  et <- etaVtoM(et, n.X)
  out <- matrix(NA, nrow=length(u), ncol=n.X)
  
  if(is.null(extra)) {
    extra <- 0 
  }
  
  if(n.u==1) {
    l1 <- et[1, 1] + et[1, 2] * u + extra
  } else if(n.u == 2) {
    l1 <- et[1, 1] + u %*% oneCol(et[1, -1]) + extra
  }
  l1[l1 < -10] = -10; l1[l1 > 10] = 10
  
  out[, 1] <- exp(l1) / (exp(l1) + 1)
  out[, 2] <- (1 - out[, 1])
  
  ## more efficient way is to put log() outside this function
  ## however, it turns out R could pass 9e-14 out and becomes 0
  ## with log() here we avoid this problem
  return(log(out))
}

updateMH <- function(accept, acc, att, MH, i, nburn) {
  acc <- acc + accept
  att <- att + 1 
  if(i < nburn) {
    if(att>50){
      if(acc/att < 0.3) {
        MH <- MH * 0.8 }
      if(acc/att > 0.5) {
        MH <- MH * 1.2 }
      acc <- att <- 0   }  }
  return(list(acc, att, MH))
}

tobit.loglik <- function(latentU, V, Vthreshold, beta, err.var, n, 
                         idx.obs, idx.above, idx.below) {
  ## return a vector of n-by-1 
  
  beta0 <- beta[1]
  beta1 <- oneCol(beta[-1])
  
  log.prob <- rep(0, n)
  Vestimate <- rep(NA, n)
  
  if(length(beta) == 2) {
    Vestimate[idx.obs] <- beta0 + latentU[idx.obs] * beta1
  } else {
    Vestimate[idx.obs] <- beta0 + latentU[idx.obs, ] %*% beta1 }
  
  log.prob[idx.above] <- log.prob[idx.above] +
    log.dnorm((V[idx.above] - Vestimate[idx.above]) / sqrt(err.var), 0, err.var) + 
    -log(sqrt(err.var))
  
  log.prob[idx.below] <- log.prob[idx.below] + 
    log(1 - pnorm( (Vestimate[idx.below] - Vthreshold) / sqrt(err.var), 
                   mean=0, sd=sqrt(err.var), log=F) )
  
  log.prob[log.prob == -Inf] <- -20 # log(.00000001)
  log.prob
}


## logit function 
## implemented for a few scenarios only 
eta2pi <- function(et, u, n.X, n.u, extra=0) {
  if(length(et) != (n.u+1)*(n.X-1)) 
    return("error in et")
  
  ## u for one observation each time 
  et <- matrix(et, nrow=n.X-1, byrow=T) #fixed 
  out <- matrix(NA, nrow=1, ncol=n.X)

  if(is.null(extra)) extra <- 0
  
  if(n.X == 4) {
    if(n.u==1) {
      out[, 1] <- 1 / (1 + exp(  -(et[1, 1] + et[1, 2] * u + extra)))
      ratio2 <- 1 / (1 + exp(  -(et[2, 1] + et[2, 2] * u + extra)))
      ratio3 <- 1 / (1 + exp(  -(et[3, 1] + et[3, 2] * u + extra)))
    } else if(n.u == 2) {
      out[, 1] <- 1 / (1 + exp(  -(et[1, 1] + oneRow(et[1, -1]) %*% oneCol(u) + extra)))
      ratio2 <- 1 / (1 + exp(  -(et[2, 1] + oneRow(et[2, -1]) %*% oneCol(u) + extra)))
      ratio3 <- 1 / (1 + exp(  -(et[3, 1] + oneRow(et[3, -1]) %*% oneCol(u) + extra)))
    }
    out[, 2] <- (1 - out[, 1]) * ratio2
    out[, 3] <- (1 - out[, 1] - out[, 2]) * ratio3
    out[, 4] <- (1 - out[, 1] - out[, 2] - out[, 3])
    
  } else if (n.X == 2) {
    if(n.u==1) {
      out[, 1] <- 1 / (1 + exp(  -(et[1, 1] + et[1, 2] * u + extra)))
    } else if(n.u == 2) {
      out[, 1] <- 1 / (1 + exp(  -(et[1, 1] + oneRow(et[1, -1]) %*% oneCol(u) + extra)))
    }
    out[, 2] <- (1 - out[, 1])  
  }

  return(out) 
}
 
# Function to create a square polygon
make_square <- function(x, y, dx, dy) {
  matrix(c(x, y, 
           x + dx, y, 
           x + dx, y + dy, 
           x, y + dy, 
           x, y), 
         ncol=2, byrow=TRUE)
}


get.geometry <- function(SPL) {
  ## set up up-right corner of each city 
  nc <- ceiling(sqrt(length(unique(SPL))))
  s1 <- seq(from=-120, to=-100, length.out = nc)
  s2 <- seq(from=25, to=40, length.out = nc)
  s <- as.matrix(expand.grid(s1,s2))
  s <- s[1:length(unique(SPL)),]
  
  dx <- 0.05; dy <- 0.05
  
  geometry <- st_sfc(lapply(integer(0), function(x) { NULL })) 
  for (i in 1:nrow(s)) {
    ng <- ceiling(sqrt(sum(SPL == unique(SPL)[i])))
    grid_list <- lapply(seq(from=s[i, 1], by=dx, length.out=ng), function(x) {
      lapply(seq(from=s[i, 2], by=dy, length.out=ng), function(y) {
        st_polygon(list(make_square(x, y, dx, dy)))
      })
    })
    
    grid_sf <- st_sfc(do.call(c, grid_list))
    grid_sf <- grid_sf[1:sum(SPL == unique(SPL)[i])]
    geometry <- c(st_geometry(geometry), st_geometry(grid_sf))
  }
  # plot(geometry, col = rgb(0, 0, 1, alpha=0.5), border = 'darkblue')
  return(geometry)
}

output.create <- function(id, id2="") {
  ### create a folder and a sets of output files 
  output.folder <- file.path("output", id) 
  dir.create(output.folder, recursive = TRUE)

  output.jags <- file.path(output.folder, paste0("jags", id2, ".txt"))
  output.post <- file.path(output.folder, paste0("post", id2, ".txt"))
  output.mcmc <- file.path(output.folder, paste0("mcmc", id2, ".pdf"))
  output.boxplot <- file.path(output.folder, paste0("boxplot", id2, ".pdf"))
  output.mixing <- file.path(output.folder, paste0("mixing", id2, ".pdf"))
  output.autocor <- file.path(output.folder, paste0("autocorrelation", id2, ".pdf"))
  output.rda <- file.path(output.folder, paste0("data", id2, ".RData"))
  output.rda.2 <- file.path(output.folder, paste0("data.2", id2, ".RData"))
  return(list(output.folder = output.folder, 
              output.jags = output.jags, 
              output.post = output.post, 
              output.mcmc = output.mcmc, 
              output.boxplot = output.boxplot,
              output.mixing = output.mixing, 
              output.autocor = output.autocor, 
              output.rda = output.rda, 
              output.rda.2=output.rda.2))
}


