###########################################################################
############################# SETUP in HPC ################################
###########################################################################
args = commandArgs(TRUE)

source("./src/helper.R")
source("./src/bayes_mcmc_latent_u1w3.R")
source("./src/spline.R")
library(tidyr)
library(caret)

cases <- rbind(
  as.matrix(expand.grid(c(0,2,4,6,8)/10, # ratio.uy
                        c(Inf), # ratio.u 
                        c(1), # n.u  
                        c(1), # nothing 
                        c(FALSE, TRUE), # no.proxy 
                        c(1:100) # data 
  )))

## if no.proxy=TRUE, model without latent,  
## this ignore value of n.u, ratio.u

iter=10*10^4; nburn=5*10^4; nthin=10

## Test
# iter=500; nburn=300; nthin=10
# args = 1

i <- as.integer(args[1])
ratio.uy <- cases[i, 1]
ratio.u <- cases[i, 2]
n.u <- cases[i, 3]
# if(cases[i, 4] == 1) {outcome="pm25"} else {outcome="no2"} 
no.proxy = cases[i, 5]
i.data <- cases[i, 6]

ratio.u <- Inf
############################# load redlining data ############################# 
input.folder <- "output/simu9_randUZ_data"
dir.create("output/simu9_randUZ_mcmc_noproxy", recursive = TRUE)
dir.create("output/simu9_randUZ_mcmc_latent_u1w3", recursive = TRUE)

## load data: dt idxSPL, argg 
load(paste0(input.folder, "/data", i.data, ".Rdata"))

## load splines, need to take subset 
envtemp <- new.env()
load(file.path(input.folder, paste0("data1_s", ratio.u*100, ".Rdata")), envir = envtemp)
splineBU <- envtemp$splineBU

envtemp <- new.env()
load(file.path(input.folder, paste0("data1_s", ratio.uy*100, ".Rdata")), envir = envtemp)
splineBZ <- envtemp$splineBU

########################### execute ####################################
list2env(lapply(as.list(dt[, c("Y", "W1", "W2", "V", "X", "SPL")]), as.vector), .GlobalEnv)

if (no.proxy) {
  output.folder <- "output/simu9_randUZ_mcmc_noproxy"
  output.file <- file.path(output.folder, paste0("noproxy_s", ratio.u*100, "_s", ratio.uy*100, "_data", i.data, ".Rdata"))
} else {
  output.folder <- "output/simu9_randUZ_mcmc_latent_u1w3"
  output.file <- file.path(output.folder, paste0("1u_3w_s", ratio.u*100, "_s", ratio.uy*100, "_data", i.data, ".Rdata"))
}

timestamp()
## Test
# knowU=F; censor.w=F; two.proxy=T; tobit.v=T; test.idx=NA

model <- gibbs_sample(Y=Y, X=X, W1=W1, W2=W2, V=V, SPL=SPL, 
                      splineBU=splineBU, splineBZ=splineBZ, n.u=n.u, 
                      iter=iter, nburn=nburn, nthin=nthin,
                      no.proxy=no.proxy) 

## save memory 
outputs <- list(WAIC=model$WAIC,
                estimate=(model$out$theta.2-model$out$theta.1)[model$out$post == 1])

if(i.data==1) {
  save(dt, model, outputs, file= output.file)
} else {
  save(dt, outputs, file= output.file) 
}





