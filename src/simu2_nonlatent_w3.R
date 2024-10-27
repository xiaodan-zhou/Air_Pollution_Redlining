###########################################################################
############################# SETUP in HPC ################################
###########################################################################
input.folder <- "output/simu2_randUZ_data"
output.folder <- "output/simu2_randUZ_mcmc_nonlatent_w3"
dir.create("output/simu2_randUZ_mcmc_nonlatent_w3", recursive = TRUE)

# < 10min per 

args = commandArgs(TRUE)

source("./src/helper.R")
source("./src/bayes_mcmc_nonlatent_3w.R")
source("./src/spline.R")
library(tidyr)
library(caret)

cases <- rbind(as.matrix(expand.grid(c(0,2,4,6,8)/10, # ratio.uy
                                     c(1), # nothing 
                                     c(1:100) # data 
))) 

iter=10*10^4; nburn=5*10^4; nthin=10

## Test
# iter=500; nburn=300; nthin=10
# args = 1

i <- as.integer(args[1])
ratio.uy <- cases[i, 1]
# if(cases[i, 2] == 1) {outcome="pm25"} else {outcome="no2"} 
i.data <- cases[i, 3]

ratio.u <- Inf 
############################# load redlining data ############################# 
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
output.file <- file.path(output.folder, paste0("3w_s", ratio.uy*100, "_data", i.data, ".Rdata"))

timestamp()
## Test
# knowU=F; censor.w=F; two.proxy=T; tobit.v=T; test.idx=NA

model <- gibbs_sample(Y=Y, X=X, W1=W1, W2=W2, V=V, SPL=SPL,
                      splineBZ=splineBZ, 
                      iter=iter, nburn=nburn, nthin=nthin) 
## save memory 
outputs <- list(WAIC=model$WAIC,
                estimate=(model$out$theta.2-model$out$theta.1)[model$out$post == 1])

if(i.data==1) {
  save(dt, model, outputs, file= output.file)
} else {
  save(dt, outputs, file= output.file) 
}


