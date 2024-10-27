################################################################
################################################################
################################################################
library(ggplot2)
library(viridis)
library(tidyr)
library(dplyr)

input <- output <- "output/simu9_latent_u1w3" 
ratio.u <- Inf 
u <- 1
for (outcome in c("Y")) {
  tb1 <- c()
  for (i.data in 1:100) {
    for (ratio.uy in c(0:9)/10) {
      file_path <- file.path(input, paste0("1u_3w_s", ratio.u*100, "_s", ratio.uy*100, "_data", i.data, ".Rdata"))
      info <- file.info(file_path)
      
      if (is.na(info$size)) { warning(paste("File not exist:", file_path))
      } else if (info$size < .1) { warning(paste("File empty:", file_path))
      } else { continue_read <- TRUE
      
      tryCatch(load(file_path), error = function(e) { continue_read <<- FALSE})
      if(continue_read) {
        
        values <- outputs$estimate
        q025 <- quantile(values, .025)
        q975 <- quantile(values, .975)
        
        tb1 <- rbind(tb1, c(u, ratio.u*100, ratio.uy*100, mean(values), sd(values),
                            q025, q975, outputs$WAIC$WAIC[4], sum(outputs$WAIC$WAIC[1:4]), i.data))
        rm(outputs)
      }
      }}}
  
  tb1 <- data.frame(tb1)
  colnames(tb1) = c("u", "ratio.u", "ratio.uy", "estimate", "sd", 
                    "low", "high", "waic", "waic.full", "data")
  print(head(tb1))
  saveRDS(tb1, file = paste0(output, "/u1w3.rds"))
}



################################################################
################################################################
################################################################
input <- output <- "output/simu9_nonlatent_w3"
u <- 0; ratio.u <- 0

for (outcome in c("Y")) {
  tb1 <- c()
  for (i.data in 1:100) {
    for (ratio.uy in c(0:9)/10) {
      file_path <- file.path(input, paste0("3w_s", ratio.uy*100, "_data", i.data, ".Rdata"))
      info <- file.info(file_path)
      
      if (is.na(info$size)) { warning(paste("File not exist:", file_path))
      } else if (info$size < .1) { warning(paste("File empty:", file_path))
      } else { continue_read <- TRUE
      
      tryCatch(load(file_path), error = function(e) { continue_read <<- FALSE})
      if(continue_read) {
        
        values <- outputs$estimate
        q025 <- quantile(values, .025)
        q975 <- quantile(values, .975)
        
        tb1 <- rbind(tb1, c(u, ratio.u*100, ratio.uy*100, mean(values), sd(values),
                            q025, q975, outputs$WAIC$WAIC[4], sum(outputs$WAIC$WAIC[1:4]), i.data))
        rm(outputs)
      }}}
  }
  
  tb1 <- data.frame(tb1)
  colnames(tb1) = c("u", "ratio.u", "ratio.uy", "estimate", "sd",
                    "low", "high", "waic", "waic.full", "data")
  print(head(tb1))
  saveRDS(tb1, file = paste0(output, "/w3.rds"))
}

################################################################
################################################################
################################################################
input <- output <- "output/simu9_no_adjustment"
u <- 0; ratio.u <- 0

for (outcome in c("Y")) {
  tb1 <- c()
  for (i.data in 1:100) {
    for (ratio.uy in c(0:9)/10) {
      file_path <- file.path(input, paste0("noproxy_sInf_s", ratio.uy*100, "_data", i.data, ".Rdata"))
      info <- file.info(file_path)
      
      if (is.na(info$size)) { warning(paste("File not exist:", file_path))
      } else if (info$size < .1) { warning(paste("File empty:", file_path))
      } else { continue_read <- TRUE
      
      tryCatch(load(file_path), error = function(e) { continue_read <<- FALSE})
      if(continue_read) {
        
        values <- outputs$estimate
        q025 <- quantile(values, .025)
        q975 <- quantile(values, .975)
        
        tb1 <- rbind(tb1, c(u, ratio.u*100, ratio.uy*100, mean(values), sd(values),
                            q025, q975, outputs$WAIC$WAIC[4], sum(outputs$WAIC$WAIC[1:4]), i.data))
        rm(outputs)
      }}}
  }
  
  tb1 <- data.frame(tb1)
  colnames(tb1) = c("u", "ratio.u", "ratio.uy", "estimate", "sd",
                    "low", "high", "waic", "waic.full", "data")
  print(head(tb1))
  saveRDS(tb1, file = paste0(output, "/noproxy.rds"))
}


################################################################
###################### estimate ################################
################################################################
## get the true effect
# load("output/simu9_data/data1.RData") forget to save 
true_effect <- .2 # as.numeric(argg["theta.2"]-argg["theta.1"])

t2 <- readRDS("output/simu9_nonlatent_w3/w3.rds")
t2$method <- "Outcome Adjustment"
t2$withZ <- "Yes"; t2$withZ[t2$ratio.uy == 0] <- "No"

t4 <- readRDS("output/simu9_latent_u1w3/u1w3.rds")
t4 <- t4[t4$ratio.u==Inf,]
t4$method <- "Latent model"
t4$withZ <- "Yes"; t4$withZ[t4$ratio.uy == 0] <- "No"

t1 <- readRDS("output/simu9_no_adjustment/noproxy.rds")
t1$method <- "No Adjustment"
t1$withZ <- "Yes"; t1$withZ[t1$ratio.uy == 0] <- "No"

tb2 <- rbind(t2,t4,t1)

tb2$true_effect <- true_effect
tb2$cover <- NA 
tb2$cover <- (tb2$low<=tb2$true_effect)&(tb2$high>=tb2$true_effect)

saveRDS(tb2, file = "output/simu9_latent_u1w3/combined.rds")


