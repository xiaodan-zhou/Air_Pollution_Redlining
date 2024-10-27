#####################################################################################################
# This script generate simulation result table. 
#####################################################################################################

library(tidyverse)
library(xtable)

tb2 <- rbind(
  cbind(readRDS(file = "output/simu1_latent_u1w3/combined.rds"), simu = '(1)'),
  cbind(readRDS(file = "output/simu9_latent_u1w3/combined.rds"), simu = '(2)'),
  cbind(readRDS(file = "output/simu0_latent_u1w3/combined.rds"), simu = '(3)'),
  cbind(readRDS(file = "output/simu2_latent_u1w3/combined.rds"), simu = '(4)'),
  cbind(readRDS(file = "output/simu7_latent_u1w3/combined.rds"), simu = '(5)'),
  cbind(readRDS(file = "output/simu8_latent_u1w3/combined.rds"), simu = '(6)'),
  cbind(readRDS(file = "output/simu3_latent_u1w3/combined.rds"), simu = '(7)')
) 


tb_selected <- tb2 %>%
  filter(simu=="(1)") %>%
  group_by(simu, method, data) %>%
  filter(waic == min(waic)) %>%
  ungroup() %>% 
  group_by(simu) %>%
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
  select(simu, metric, value_selected) 

tb3 <- tb2 %>%
  filter(method %in% c("Latent model", "No Adjustment", "Outcome Adjustment")) %>% 
  group_by(simu, method, ratio.uy) %>% 
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
  merge(tb_selected, by = c("simu", "metric")) %>% 
  mutate(metric_label = factor(metric, levels = c("bias", "mse", "waic", "cover"), 
                               labels = c("Bias", "MSE", "WAIC", "Coverage Probability")))
  

tb4 <- tb2 %>% 
  mutate(method = factor(method, levels = c("Latent model", "Outcome Adjustment", "No Adjustment"))) %>% 
  group_by(simu, method, data) %>%  
  filter(waic == min(waic)) %>%
  ungroup() %>% 
  group_by(simu, method) %>%
  summarise(bias=mean(estimate-true_effect), 
            mse = mean((estimate-true_effect)^2), 
            cover=mean(cover)*100, 
            waic=mean(waic),
            ratio.uy=round(mean(ratio.uy),0),
            )

print(xtable(tb4, digits = 3), include.rownames=FALSE)
print(xtable(tb4, digits = 2), include.rownames=FALSE)

## get S.D. 
tb5 <- tb2 %>% 
  mutate(method = factor(method, levels = c("Latent model", "Outcome Adjustment", "No Adjustment"))) %>% 
  group_by(simu, method, data) %>%  
  filter(waic == min(waic)) %>%
  ungroup() %>% 
  group_by(simu, method) %>%
  summarise(bias=sd(estimate-true_effect), 
            mse = sd((estimate-true_effect)^2)
  )
print(xtable(tb5, digits = 3), include.rownames=FALSE)


## if we directly use true ratio.uy
tb4 <- tb2 %>% 
  mutate(method = factor(method, levels = c("Latent model", "Outcome Adjustment", "No Adjustment"))) %>% 
  group_by(simu, method, data) %>%  
  filter(ratio.uy == 60) %>% # or 60 
  ungroup() %>% 
  group_by(simu, method) %>%
  summarise(bias=mean(estimate-true_effect), 
            mse = mean((estimate-true_effect)^2), 
            cover=mean(cover)*100, 
            waic=mean(waic),
            ratio.uy=round(mean(ratio.uy),0),
  )
tb4
