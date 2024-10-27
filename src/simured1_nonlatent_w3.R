###########################################################################
############################# SETUP in HPC ################################
###########################################################################
# < 10min per 

args = commandArgs(TRUE)

source("./src/data_gen.R")
source("./src/bayes_mcmc_nonlatent_3w.R")
source("./src/spline.R")
library(tidyr)
library(caret)

cases <- rbind(as.matrix(expand.grid(c(0,2,4,6,8)/10, # ratio.uy
                                     c(1,2), # pm25 or no2 
                                     c(101:500) # data 
))) 

iter=10*10^4; nburn=5*10^4; nthin=10

## Test
# iter=500; nburn=300; nthin=10
# args = 1

args <- as.integer(args[1])
for (i in c(args, args+1000, args+2000, args+3000)) {
  ratio.uy <- cases[i, 1]
  if(cases[i, 2] == 1) {outcome="pm25"} else {outcome="no2"} 
  i.data <- cases[i, 3]
  
  ############################# load redlining data ############################# 
  input.folder <- "output/simured1_randUZ_data"
  output.folder <- "output/simured1_randUZ_mcmc_nonlatent_w3"
  dir.create("output/simured1_randUZ_mcmc_nonlatent_w3", recursive = TRUE)
  
  ## load data: dt idxSPL, argg 
  load(paste0(input.folder, "/data", i.data, "_", outcome, ".RData"))
  
  ## load splines, need to take subset 
  splineBZ.temp <- readRDS(paste0("data/nhgis0008_csv/only 1940/RedliningSplineMat_R", ratio.uy*100, ".rds"))
  splineBZ <- splineBZ.temp[idxSPL]
  rm(splineBZ.temp)
  
  ########################### execute ####################################
  list2env(lapply(as.list(dt[, c("Y", "W1", "W2", "V", "X", "SPL")]), as.vector), .GlobalEnv)
  output.file <- file.path(output.folder, paste0(outcome, "_3w_s", ratio.uy*100, "_data", i.data, ".RData"))
  
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
}

