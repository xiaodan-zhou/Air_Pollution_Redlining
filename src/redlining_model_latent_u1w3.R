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
  as.matrix(expand.grid(0:9/10, # ratio.uy
                        c(Inf), # ratio.u 
                        c(1), # n.u  
                        c(1,2), # pm25 or no2 
                        c(TRUE) # c(FALSE, TRUE) # no.proxy 
                        )))

## if no.proxy=TRUE, model without latent,  
## this ignore value of n.u, ratio.u

iter=30*10^4; nburn=15*10^4; nthin=10; nc=1 ## < 1min per 1000 iter, < 5h total 

## Test 
# iter=30; nburn=15; nthin=1; # args=2; outcome="pm25"
# iter=2000; nburn=1000; nthin=1; nc=1
# iter=3*10^4; nburn=1.5*10^4; nthin=10; nc=1

i <- as.integer(args[1])
ratio.uy <- cases[i, 1]
ratio.u <- cases[i, 2]
n.u <- cases[i, 3]
if(cases[i, 4] == 1) {outcome="pm25"} else {outcome="no2"} 
no.proxy = cases[i, 5]


############################# load redlining data ############################# 
splineBU <- readRDS(paste0("data/cleaned/RedliningSplineMat_R", ratio.u*100, ".rds"))
splineBZ <- readRDS(paste0("data/cleaned/RedliningSplineMat_R", ratio.uy*100, ".rds"))

dt <- readRDS("data/cleaned/final_census1940.rds")

dt <- dt %>% group_by(SPL) %>% mutate(
  wt.mean.rent.center = wt.mean.rent - mean(wt.mean.rent, na.rm=T), 
  pm25 = pm25 - mean(pm25, na.rm=T),
  no2 = no2 - mean(no2, na.rm=T)
  ) %>% data.frame(.) 

dt$W1 <- dt$bc.rate.employed
dt$W2 <- dt$wt.mean.rent.center
dt$V <- dt$pct.black.pop40.rank

########################## proxy, treatment, outcome ############################# 
dt$W1 <- dt$W1/sd(dt$W1,na.rm=T)
dt$W2 <- dt$W2/sd(dt$W2,na.rm=T)
dt$V <- dt$V/sd(dt$V,na.rm=T)

dt$Y <- dt[[outcome]]

dt$X <- dt$holc
dt$X[dt$X <= 2] = 1; dt$X[dt$X >= 3] = 2 ## binary treatment 

########################### execute ####################################
list2env(lapply(as.list(dt[, c("Y", "W1", "W2", "V", "X", "SPL")]), as.vector), .GlobalEnv)
dir.create("output/redlining_mcmc_noproxy", recursive = TRUE)
dir.create("output/redlining_mcmc_latent_u1w3", recursive = TRUE)

if (no.proxy) {
  output.folder <- "output/redlining_mcmc_noproxy"
  output.file <- file.path(output.folder, paste0(outcome, "_noproxy_s", ratio.u*100, "_s", ratio.uy*100, ".Rdata"))
} else {
  output.folder <- "output/redlining_mcmc_latent_u1w3"
  output.file <- file.path(output.folder, paste0(outcome, "_1u_3w_s", ratio.u*100, "_s", ratio.uy*100, ".Rdata"))
}

timestamp()
## Test
# knowU=F; censor.w=F; two.proxy=T; tobit.v=T; test.idx=NA

model <- gibbs_sample(Y=Y, X=X, W1=W1, W2=W2, V=V, SPL=SPL, 
                      splineBU=splineBU, splineBZ=splineBZ, n.u=n.u, 
                      iter=iter, nburn=nburn, nthin=nthin,
                      no.proxy=no.proxy) 

save(dt, model, file= output.file)

# plot(model$out$theta.2-model$out$theta.1)

















################################################################
################################################################
################################################################
library(ggplot2)
library(viridis)
library(tidyr)
library(dplyr)

input <- output <- "output/redlining_mcmc_latent_u1w3" 

u <- 1
for (outcome in c("pm25", "no2")) {
  tb1 <- c()
  for (ratio.u in c(Inf,.1,.15,.2,.25,.3,.35,.4,.45,.5,.6,.7,.8,.9)) {
    for (ratio.uy in c(0,.1,.15,.2,.25,.3,.35,.4,.45,.5,.6,.7,.8,.9)) {
      file_path <- file.path(input, paste0(outcome, "_1u_3w_s", ratio.u*100, "_s", ratio.uy*100, ".Rdata"))
      info <- file.info(file_path)
      
      if (is.na(info$size)) { warning(paste("File not exist:", file_path))
      } else if (info$size < .1) { warning(paste("File empty:", file_path))
      } else { continue_read <- TRUE
      
      tryCatch(load(file_path), error = function(e) { continue_read <<- FALSE})
      if(continue_read) {
        
        values <- (model$out$theta.2-model$out$theta.1)[model$out$post==1]
        q025 <- quantile(values, .025)
        q975 <- quantile(values, .975)
        
        tb1 <- rbind(tb1, c(u, ratio.u*100, ratio.uy*100, mean(values), sd(values),
                            q025, q975, model$WAIC$WAIC[4], sum(model$WAIC$WAIC[1:4])))
        rm(model)
      }
      }}}
  
  tb1 <- data.frame(tb1)
  colnames(tb1) = c("u", "ratio.u", "ratio.uy", "estimate", "sd", 
                    "low", "high", "waic", "waic.full")
  print(head(tb1))
  saveRDS(tb1, file = paste0(output, "/", outcome, "_u1w3.rds"))
}

################################################################
input <- output <- "output/redlining_mcmc_latent_u1w3" 
u <- 1
for (outcome in c("pm25", "no2")) {
  tb1 <- c()
  for (ratio.u in c(Inf,.1,.15,.2,.25,.3,.35,.4,.45,.5,.6,.7,.8,.9)) {
    for (ratio.uy in c(0,.1,.15,.2,.25,.3,.35,.4,.45,.5,.6,.7,.8,.9)) {
      file_path <- file.path(input, paste0(outcome, "_1u_3w_s", ratio.u*100, "_s", ratio.uy*100, ".Rdata"))
      info <- file.info(file_path)
      
      if (is.na(info$size)) { warning(paste("File not exist:", file_path))
      } else if (info$size < .1) { warning(paste("File empty:", file_path))
      } else { continue_read <- TRUE
      
      tryCatch(load(file_path), error = function(e) { continue_read <<- FALSE})
      if(continue_read) {
        tb1 <- rbind(tb1, cbind(u, ratio.u*100, ratio.uy*100, dt$neighborho, 
                                model$BU.mean, model$BZ.mean))
        rm(model) }}}}
  
  tb1 <- data.frame(tb1)
  colnames(tb1) = c("u", "ratio.u", "ratio.uy", "location", "BU", "BZ")
  print(head(tb1))
  saveRDS(tb1, file = paste0(output, "/", outcome, "_BZ.rds"))
}

################################################################
input <- output <- "output/redlining_mcmc_latent_u1w3" 
u <- 1
for (outcome in c("pm25", "no2")) {
  tb1 <- c()
  for (ratio.u in c(Inf,.1,.15,.2,.25,.3,.35,.4,.45,.5,.6,.7,.8,.9)) {
    for (ratio.uy in c(0,.1,.15,.2,.25,.3,.35,.4,.45,.5,.6,.7,.8,.9)) {
      file_path <- file.path(input, paste0(outcome, "_1u_3w_s", ratio.u*100, "_s", ratio.uy*100, ".Rdata"))
      info <- file.info(file_path)
      
      if (is.na(info$size)) { warning(paste("File not exist:", file_path))
      } else if (info$size < .1) { warning(paste("File empty:", file_path))
      } else { continue_read <- TRUE
      
      tryCatch(load(file_path), error = function(e) { continue_read <<- FALSE})
      if(continue_read) {
        tb1 <- rbind(tb1, c(u, ratio.u*100, ratio.uy*100, colMeans(model$out[model$out$post==1,])))
        # rm(model) 
      }}}}
  
  tb1 <- data.frame(tb1)
  colnames(tb1) = c("u", "ratio.u", "ratio.uy", names(model$out))
  print(head(tb1))
  saveRDS(tb1, file = paste0(output, "/", outcome, "_params.rds"))
}

################################################################
################################################################
################################################################
input <- output <- "output/redlining_mcmc_noproxy"

u <- 1
for (outcome in c("pm25", "no2")) {
  tb1 <- c()
  for (ratio.u in c(Inf,.1,.15,.2,.25,.3,.35,.4,.45,.5,.6,.7,.8,.9)) {
    for (ratio.uy in c(0,.1,.15,.2,.25,.3,.35,.4,.45,.5,.6,.7,.8,.9)) {
      file_path <- file.path(input, paste0(outcome, "_noproxy_s", ratio.u*100, "_s", ratio.uy*100, ".Rdata"))
      info <- file.info(file_path)
      
      if (is.na(info$size)) { warning(paste("File not exist:", file_path))
      } else if (info$size < .1) { warning(paste("File empty:", file_path))
      } else { continue_read <- TRUE
      
      tryCatch(load(file_path), error = function(e) { continue_read <<- FALSE})
      if(continue_read) {
        
        values <- (model$out$theta.2-model$out$theta.1)[model$out$post==1]
        q025 <- quantile(values, .025)
        q975 <- quantile(values, .975)
        
        tb1 <- rbind(tb1, c(u, ratio.u*100, ratio.uy*100, mean(values), sd(values),
                            q025, q975, model$WAIC$WAIC[4], sum(model$WAIC$WAIC[1:4])))
        rm(model)
      }
      }}}
  
  tb1 <- data.frame(tb1)
  colnames(tb1) = c("u", "ratio.u", "ratio.uy", "estimate", "sd", 
                    "low", "high", "waic", "waic.full")
  print(head(tb1))
  saveRDS(tb1, file = paste0(output, "/", outcome, "_noproxy.rds"))
}

################################################################
input <- output <- "output/redlining_mcmc_noproxy"

u <- 1
for (outcome in c("pm25", "no2")) {
  tb1 <- c()
  for (ratio.u in c(Inf)) { 
    for (ratio.uy in c(0,.1,.15,.2,.25,.3,.35,.4,.45,.5,.6,.7,.8,.9)) {
      file_path <- file.path(input, paste0(outcome, "_noproxy_s", ratio.u*100, "_s", ratio.uy*100, ".Rdata"))
      info <- file.info(file_path)
      
      if (is.na(info$size)) { warning(paste("File not exist:", file_path))
      } else if (info$size < .1) { warning(paste("File empty:", file_path))
      } else { continue_read <- TRUE
      
      tryCatch(load(file_path), error = function(e) { continue_read <<- FALSE})
      if(continue_read) {
        model$BU.mean <- NA 
        tb1 <- rbind(tb1, cbind(u, ratio.u*100, ratio.uy*100, dt$neighborho, model$BU.mean, model$BZ.mean))
        rm(model) }}}}
  
  tb1 <- data.frame(tb1)
  colnames(tb1) = c("u", "ratio.u", "ratio.uy", "location", "BU", "BZ")
  print(head(tb1))
  saveRDS(tb1, file = paste0(output, "/", outcome, "_BZ.rds"))
}

################################################################
input <- output <- "output/redlining_mcmc_noproxy"

u <- 1
for (outcome in c("pm25", "no2")) {
  tb1 <- c()
  for (ratio.u in c(Inf,.1,.15,.2,.25,.3,.35,.4,.45,.5,.6,.7,.8,.9)) {
    for (ratio.uy in c(0,.1,.15,.2,.25,.3,.35,.4,.45,.5,.6,.7,.8,.9)) {
      file_path <- file.path(input, paste0(outcome, "_noproxy_s", ratio.u*100, "_s", ratio.uy*100, ".Rdata"))
      info <- file.info(file_path)
      
      if (is.na(info$size)) { warning(paste("File not exist:", file_path))
      } else if (info$size < .1) { warning(paste("File empty:", file_path))
      } else { continue_read <- TRUE
      
      tryCatch(load(file_path), error = function(e) { continue_read <<- FALSE})
      if(continue_read) {
        tb1 <- rbind(tb1, c(u, ratio.u*100, ratio.uy*100, colMeans(model$out[model$out$post==1,])))
        # rm(model) 
      }}}}
  
  tb1 <- data.frame(tb1)
  colnames(tb1) = c("u", "ratio.u", "ratio.uy", names(model$out))
  print(head(tb1))
  saveRDS(tb1, file = paste0(output, "/", outcome, "_params.rds"))
}
