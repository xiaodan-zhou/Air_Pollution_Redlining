###########################################################################
############################# SETUP in HPC ################################
###########################################################################
args = commandArgs(TRUE)

source("./src/helper.R")
source("./src/bayes_mcmc_latent_u1w3_random_effect.R")
source("./src/spline.R")
library(tidyr)
library(caret)

cases <- rbind(
  as.matrix(expand.grid(5:7/10, # ratio.uy
                        c(Inf), # ratio.u 
                        c(1), # n.u  
                        c(1,2), # pm25 or no2 
                        c(FALSE) # c(FALSE, TRUE) # no.proxy 
  )))

## if no.proxy=TRUE, model without latent,  
## this ignore value of n.u, ratio.u

iter=30*10^4; nburn=15*10^4; nthin=10; nc=1 ## < 1min per 1000 iter, < 5h total 

## Test 
# iter=1000; nburn=600; nthin=1; 
# args=6; outcome="no2"
# iter=3000; nburn=2000; nthin=1 
# iter=3*10^4; nburn=1.5*10^4; nthin=10; nc=1 # 10min

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

####################### load initial values ###############################
init_env <- new.env()
load(paste0("output/redlining_mcmc_latent_u1w3/", outcome, "_1u_3w_sInf_s", ratio.uy*100, ".RData"), 
     envir = init_env)
model <- init_env$model
init <- c(model, colMeans(model$out[model$out$post==1,]))
init <- init[!names(init) %in% c("out")]
rm(model, init_env)
###########################################################################

########################### execute ####################################
list2env(lapply(as.list(dt[, c("Y", "W1", "W2", "V", "X", "SPL")]), as.vector), .GlobalEnv)
dir.create("output/redlining_mcmc_latent_u1w3_random_effect", recursive = TRUE)

output.folder <- "output/redlining_mcmc_latent_u1w3_random_effect"
output.file <- file.path(output.folder, paste0(outcome, "_1u_3w_s", ratio.u*100, "_s", ratio.uy*100, ".Rdata"))

timestamp()
## Test
# knowU=F; censor.w=F; two.proxy=T; tobit.v=T; test.idx=NA

model <- gibbs_sample(Y=Y, X=X, W1=W1, W2=W2, V=V, SPL=SPL, 
                      splineBU=splineBU, splineBZ=splineBZ, n.u=n.u, 
                      iter=iter, nburn=nburn, nthin=nthin,
                      no.proxy=no.proxy, init=init)

save(dt, model, file= output.file)

# model <- copy
# plot(model$out$mu.theta-model$out$theta0)
# plot(model$out$theta0)
# plot(model$out$mu.theta)
# mean((model$out$mu.theta-model$out$theta0)[model$out$post==1]) # 0.211 
# 
# model$out[,paste0("theta.", 1:3)] %>% 
#   pivot_longer(cols = everything(), names_to = "theta", values_to = "value") %>% 
#   mutate(xaxis=rep(1:1000, 3)) %>%
#   as.data.frame() %>% 
#   ggplot(aes(x = xaxis, y = value, col=theta)) + 
#   geom_smooth() + theme(legend.position="none")
# 
# plot(model$out$sig2.theta)
# hist(model$effect)
# sum(model$effect<0)
# hist(model$effect.higher)
# sum(model$effect.higher>.95)
# hist(model$effect.var)
# colMeans(model$out[model$out$post==1,])
# 
# copy <- model 
# 
# load("/Users/mac/Documents/GitHub/Redlining_Air_Pollution/output/redlining_mcmc_latent_u1w3/no2_1u_3w_sInf_s50.RData")
# colMeans(model$out[model$out$post==1,])
# plot(model$BU.mean, copy$BU.mean)
# plot(model$BZ.mean, copy$BZ.mean)
# plot(model$BU.var, copy$BU.var) 
# plot(model$BZ.var, copy$BZ.var)

################################################################
################################################################
################################################################
library(ggplot2)
library(viridis)
library(tidyr)
library(dplyr)

input <- output <- "output/redlining_mcmc_latent_u1w3_random_effect" 

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
input <- output <- "output/redlining_mcmc_latent_u1w3_random_effect" 
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
input <- output <- "output/redlining_mcmc_latent_u1w3_random_effect" 
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
