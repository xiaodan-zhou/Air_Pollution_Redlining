################################################################
################################################################
################################################################
library(ggplot2)
library(viridis)
library(tidyr)
library(dplyr)

input <- output <- "output/simured1_latent_u1w3" 
ratio.u <- Inf 
u <- 1
for (outcome in c("pm25", "no2")) {
  tb1 <- c()
  for (i.data in 1:500) {
    for (ratio.uy in c(0:9)/10) {
      file_path <- file.path(input, paste0(outcome, "_1u_3w_s", ratio.u*100, "_s", ratio.uy*100, "_data", i.data, ".Rdata"))
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
  saveRDS(tb1, file = paste0(output, "/", outcome, "_u1w3.rds"))
}



################################################################
################################################################
################################################################
input <- output <- "output/simured1_nonlatent_w3"
u <- 0; ratio.u <- 0

for (outcome in c("pm25", "no2")) {
  tb1 <- c()
  for (i.data in 1:500) {
    for (ratio.uy in c(0:9)/10) {
      file_path <- file.path(input, paste0(outcome, "_3w_s", ratio.uy*100, "_data", i.data, ".Rdata"))
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
  saveRDS(tb1, file = paste0(output, "/", outcome, "_w3.rds"))
}

################################################################
################################################################
################################################################
input <- output <- "output/simured1_no_adjustment"
u <- 0; ratio.u <- 0

for (outcome in c("pm25", "no2")) {
  tb1 <- c()
  for (i.data in 1:500) {
    for (ratio.uy in c(0:9)/10) {
      file_path <- file.path(input, paste0(outcome, "_noproxy_sInf_s", ratio.uy*100, "_data", i.data, ".Rdata"))
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
  saveRDS(tb1, file = paste0(output, "/", outcome, "_noproxy.rds"))
}


################################################################
###################### estimate ################################
################################################################
## get the true effect
load("output/simured1_randUZ_data/data1_no2.RData")
no2_effect <- as.numeric(argg$theta.2-argg$theta.1)

load("output/simured1_randUZ_data/data1_pm25.RData")
pm25_effect <- as.numeric(argg$theta.2-argg$theta.1)
# argg["theta.2"]-argg["theta.1"])

t2 <- readRDS("output/simured1_nonlatent_w3/no2_w3.rds")
t2$outcome <- "no2"
t2$method <- "Outcome Adjustment"
t2$withZ <- "Yes"; t2$withZ[t2$ratio.uy == 0] <- "No"

t3 <- readRDS("output/simured1_nonlatent_w3/pm25_w3.rds")
t3$outcome <- "pm25"
t3$method <- "Outcome Adjustment"
t3$withZ <- "Yes"; t3$withZ[t3$ratio.uy == 0] <- "No"

t4 <- readRDS("output/simured1_latent_u1w3/no2_u1w3.rds")
t4 <- t4[t4$ratio.u==Inf,]
t4$outcome <- "no2"
t4$method <- "Latent model"
t4$withZ <- "Yes"; t4$withZ[t4$ratio.uy == 0] <- "No"

t5 <- readRDS("output/simured1_latent_u1w3/pm25_u1w3.rds")
t5 <- t5[t5$ratio.u==Inf,]
t5$outcome <- "pm25"
t5$method <- "Latent model"
t5$withZ <- "Yes"; t5$withZ[t5$ratio.uy == 0] <- "No"

t1 <- readRDS("output/simured1_noproxy/no2_noproxy.rds")
t1$outcome <- "no2"
t1$method <- "No Adjustment"
t1$withZ <- "Yes"; t1$withZ[t1$ratio.uy == 0] <- "No"

t6 <- readRDS("output/simured1_noproxy/pm25_noproxy.rds")
t6$outcome <- "pm25"
t6$method <- "No Adjustment"
t6$withZ <- "Yes"; t6$withZ[t6$ratio.uy == 0] <- "No"

tb2 <- rbind(t1,t2,t3,t4,t5,t6)

tb2$true_effect <- NA
tb2$true_effect[tb2$outcome=="pm25"] <- pm25_effect
tb2$true_effect[tb2$outcome=="no2"] <- no2_effect
tb2$cover <- NA 
tb2$cover <- (tb2$low<=tb2$true_effect)&(tb2$high>=tb2$true_effect)


#####################
tb_selected <- tb2 %>%
  group_by(outcome, data) %>%
  filter(waic == min(waic)) %>%
  ungroup() %>% 
  group_by(outcome) %>%
  summarise(true=mean(true_effect), 
            estimate=mean(estimate), 
            bias=mean(estimate-true_effect), 
            abs.bias=mean(abs(estimate-true_effect)), 
            mse = mean((estimate-true_effect)^2), 
            sd=mean(sd), 
            cover=mean(cover)*100, 
            waic=mean(waic),
            n=n(),
            u=u[1], 
            ratio.u=ratio.u[1],
            withZ=withZ[1]) %>% 
  pivot_longer(cols = c("bias", "mse", "waic", "cover"), names_to = "metric", values_to = "value_selected") %>% 
  select(outcome, metric, value_selected) 


## compute the coverage probability by method
## tb2$cp <- ...
## group by method, compute cp etc, point: cancel datacolumn so that we can make plot??? 
## do some careful check
## ignore some of them, make no sense to compute average 
tb3 <- tb2 %>% 
  group_by(method, outcome, ratio.uy) %>% 
  summarise(true=mean(true_effect), 
            estimate=mean(estimate), 
            bias=mean(estimate-true_effect), 
            abs.bias=mean(abs(estimate-true_effect)), 
            mse = mean((estimate-true_effect)^2), 
            sd=mean(sd), 
            cover=mean(cover)*100, 
            waic=mean(waic),
            n=n(),
            u=u[1], 
            ratio.u=ratio.u[1],
            withZ=withZ[1]) %>% 
  pivot_longer(cols = c("bias", "mse", "waic", "cover"), names_to = "metric", values_to = "value") %>% 
  mutate(metric = factor(metric, levels = c("bias", "mse", "cover", "waic")),
         xaxis = as.numeric(ratio.uy),
         yintp = case_when(
           metric == "bias" ~ 0,
           metric == "mse" ~ 0,
           metric == "cover" ~ 95,
           metric == "waic" ~ NA
         )) %>% 
  merge(tb_selected, by = c("outcome", "metric")) %>% 
  mutate(metric_label = factor(metric, levels = c("bias", "mse", "waic", "cover"), 
                               labels = c("Bias", "MSE", "WAIC", "Coverage Probability")))

tb3$xaxis <- as.numeric(tb3$ratio.uy)
# tb3$xaxis[tb3$method=="No Adjustment"] <- tb3$xaxis[tb3$method=="No Adjustment"] - 1
# tb3$xaxis[tb3$method=="Outcome Adjustment"] <- tb3$xaxis[tb3$method=="Outcome Adjustment"] + 1

tb3$outcome <- factor(tb3$outcome, levels = c("pm25", "no2"))


tb4 <- tb2 %>% 
  mutate(method = factor(method, levels = c("Latent model", "Outcome Adjustment", "No Adjustment"))) %>% 
  group_by(method, outcome, data) %>% 
  filter(waic == min(waic)) %>%
  ungroup() %>% 
  group_by(outcome, method) %>% 
  summarise(bias=mean(estimate-true_effect), 
            mse = mean((estimate-true_effect)^2), 
            cover=mean(cover)*100, 
            waic=mean(waic),
            ratio.uy=round(mean(ratio.uy),0)
  )
library(xtable)
print(xtable(tb4), include.rownames=FALSE)

            
## if we directly use true ratio.uy
tb4 <- tb2 %>% 
  mutate(method = factor(method, levels = c("Latent model", "Outcome Adjustment", "No Adjustment"))) %>% 
  group_by(method, outcome, data) %>% 
  filter(ratio.uy == 60) %>% # or 60 
  ungroup() %>% 
  group_by(outcome, method) %>% 
  summarise(bias=mean(estimate-true_effect), 
            mse = mean((estimate-true_effect)^2), 
            cover=mean(cover)*100, 
            waic=mean(waic),
            ratio.uy=round(mean(ratio.uy),0)
  )
tb4


## get S.D. 
tb5 <- tb2 %>% 
  mutate(method = factor(method, levels = c("Latent model", "Outcome Adjustment", "No Adjustment"))) %>% 
  group_by(method, outcome, data) %>% 
  filter(waic == min(waic)) %>%
  ungroup() %>% 
  group_by(outcome, method) %>% 
  summarise(bias=sd(estimate-true_effect), 
            mse = sd((estimate-true_effect)^2)
  )
print(xtable(tb5, digits = 3), include.rownames=FALSE)
