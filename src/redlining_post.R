library(ggplot2)
library(viridis)
library(tidyr)
library(dplyr)
library(xtable)

################################################################
###################### estimate ################################
################################################################
t2 <- readRDS("output/redlining_mcmc_nonlatent_w3/no2_w3.rds")
t2$outcome <- "no2"
t2$method <- "Outcome Reg with Proxy"
t2$withZ <- "Yes"; t2$withZ[t2$ratio.uy == 0] <- "No"

t3 <- readRDS("output/redlining_mcmc_nonlatent_w3/pm25_w3.rds")
t3$outcome <- "pm25"
t3$method <- "Outcome Reg with Proxy"
t3$withZ <- "Yes"; t3$withZ[t3$ratio.uy == 0] <- "No"

t4 <- readRDS("output/redlining_mcmc_latent_u1w3/no2_u1w3.rds")
t4 <- t4[t4$ratio.u==Inf,]
t4$outcome <- "no2"
t4$method <- "Latent Adjustment"
t4$withZ <- "Yes"; t4$withZ[t4$ratio.uy == 0] <- "No"

t5 <- readRDS("output/redlining_mcmc_latent_u1w3/pm25_u1w3.rds")
t5 <- t5[t5$ratio.u==Inf,]
t5$outcome <- "pm25"
t5$method <- "Latent Adjustment"
t5$withZ <- "Yes"; t5$withZ[t5$ratio.uy == 0] <- "No"

t1 <- readRDS("output/redlining_mcmc_noproxy/no2_noproxy.rds")
t1$outcome <- "no2"
t1$method <- "No Adjustment"
t1$withZ <- "Yes"; t1$withZ[t1$ratio.uy == 0] <- "No"

t6 <- readRDS("output/redlining_mcmc_noproxy/pm25_noproxy.rds")
t6$outcome <- "pm25"
t6$method <- "No Adjustment"
t6$withZ <- "Yes"; t6$withZ[t6$ratio.uy == 0] <- "No"

tb2 <- rbind(t1,t2,t3,t4,t5,t6)
tb2$xaxis <- as.numeric(tb2$ratio.uy)
tb2$xaxis[tb2$method=="No Adjustment"] <- tb2$xaxis[tb2$method=="No Adjustment"] - 1
tb2$xaxis[tb2$method=="Outcome Reg with Proxy"] <- tb2$xaxis[tb2$method=="Outcome Reg with Proxy"] + 1


### main results, constant effect estimate
output <- "output/redlining_plots"

for (outc in c('pm25', 'no2')) {
  pdf(paste0(output, "/", outc, "_estimate_nolegend.pdf"), width = 3, height = 3)
  p1 <- ggplot(tb2[tb2$outcome==outc,]) +
    geom_point(aes(x = xaxis, y = estimate, color = as.factor(method)), size = 2) +  
    geom_hline(yintercept = 0, color="black") + 
    geom_errorbar(aes(x = xaxis, y = estimate, ymin = low, ymax = high, color = as.factor(method)), width = 2) +
    theme_bw() + 
    scale_color_discrete(name = "Method") + 
    theme(legend.position = "none")
  
  if (outc == 'pm25') {
    p1 <- p1 + labs(title = bquote(PM[2.5]), x = "Spline Ratio", y = "Effect Estimate")
  } else {
    p1 <- p1 + labs(title = bquote(NO[2]), x = "Spline Ratio", y = "Effect Estimate")
  }   
  
  print(p1)
  dev.off()
}

 
for (outc in c('pm25', 'no2')) {
  pdf(paste0(output, "/", outc, "_estimate_legend.pdf"), width = 5, height = 3)
  p1 <- ggplot(tb2[tb2$outcome==outc,]) +
    geom_point(aes(x = xaxis, y = estimate, color = as.factor(method)), size = 2) +  
    geom_hline(yintercept = 0, color="black") +  
    geom_errorbar(aes(x = xaxis, y = estimate, ymin = low, ymax = high, color = as.factor(method)), width = 2) +
    theme_bw() + scale_color_discrete(name = "Method")  
  
  if (outc == 'pm25') {
    p1 <- p1 + labs(title = bquote(PM[2.5]), x = "Spline Ratio", y = "")
  } else {
    p1 <- p1 + labs(title = bquote(NO[2]), x = "Spline Ratio", y = "")
  } 
  print(p1)
  dev.off()
}


for (outc in c('pm25', 'no2')) {
  
  pdf(paste0(output, "/", outc, "_waic_nolengend.pdf"), width = 3, height = 3)
  p1 <- ggplot(tb2[tb2$outcome==outc,]) +
    geom_point(aes(x = xaxis, y = waic, color = as.factor(method))) +  
    geom_line(aes(x = xaxis, y = waic, color = as.factor(method))) +  
    theme_bw() + 
    scale_color_discrete(name = "Method") + 
    theme(legend.position = "none") 
  
  if (outc == 'pm25') {
    p1 <- p1 + labs(title = bquote(PM[2.5]), x = "Spline Ratio", y = "WAIC")
  } else {
    p1 <- p1 + labs(title = bquote(NO[2]), x = "Spline Ratio", y = "WAIC")
  } 
  
  print(p1)
  dev.off()

}

for (outc in c('pm25', 'no2')) {

  pdf(paste0(output, "/", outc, "_waic_legend.pdf"), width = 5, height = 3)
  p1 <- ggplot(tb2[tb2$outcome=="no2",]) +
    geom_point(aes(x = xaxis, y = waic, color = as.factor(method))) +  
    geom_line(aes(x = xaxis, y = waic, color = as.factor(method))) +  
    theme_bw() + 
    scale_color_discrete(name = "Method")  
  
  if (outc == 'pm25') {
    p1 <- p1 + labs(title = bquote(PM[2.5]), x = "Spline Ratio", y = "")
  } else {
    p1 <- p1 + labs(title = bquote(NO[2]), x = "Spline Ratio", y = "")
  } 
  print(p1)
  dev.off()
}


#### print more stats for paper writing 
print(xtable(tb2[(tb2$ratio.uy==60)&(tb2$ratio.u==Inf)&(tb2$method=="Latent Adjustment"), c("estimate", "low", "high")]))





################################################################
###################### overlap #################################
################################################################
t4 <- readRDS("output/redlining_mcmc_latent_u1w3/no2_BZ.rds")
t4 <- t4[t4$ratio.u==Inf,]
t4$outcome <- "no2"
t4$method <- "Latent U model, with Z(s)"
t4$method[t4$ratio.uy == 0] <- "Latent U model, without Z(s)"

t5 <- readRDS("output/redlining_mcmc_latent_u1w3/pm25_BZ.rds")
t5 <- t5[t5$ratio.u==Inf,]
t5$outcome <- "pm25"
t5$method <- "Latent U model, with Z(s)"
t5$method[t5$ratio.uy == 0] <- "Latent U model, without Z(s)"

t0 <- rbind(t4,t5)
t0$xaxis <- as.numeric(t0$ratio.uy)
t0$xaxis[t0$u==0] <- t0$xaxis[t0$u==0] - 1
t0$xaxis[t0$method=="no W anywhere"] <- t0$xaxis[t0$method=="no W anywhere"] + 1

############ original data ############
dt <- readRDS("data/nhgis0008_csv/only 1940/final_census1940.rds")
dt <- dt %>% group_by(SPL) %>% mutate(
  wt.mean.rent.center = wt.mean.rent - mean(wt.mean.rent, na.rm=T), 
  pm25 = pm25 - mean(pm25, na.rm=T),
  no2 = no2 - mean(no2, na.rm=T)
) %>% data.frame(.) 
dt$X <- dt$holc
dt$X[dt$X <= 2] = 1; dt$X[dt$X >= 3] = 2 

t0$W1 <- t0$bc.rate.employed <- dt$bc.rate.employed
t0$W2 <- t0$wt.mean.rent.center <- dt$wt.mean.rent.center
t0$V <- t0$pct.black.pop40.rank <- dt$pct.black.pop40.rank
t0$pm25 <- dt$pm25
t0$no2 <- dt$no2
t0$X <- dt$X

summary_BU <- t0  %>% 
  filter(!is.na(BU)) %>% 
  group_by(outcome, method, ratio.uy) %>% summarise(
  n = n(),
  mean = mean(BU, na.rm=T),
  sd = sd(BU, na.rm=T)
)

### we flip U in paper writing so that the higher the worse
t0$BU <- -t0$BU

output <- "output/redlining_info"

####################################
for (outc in c("pm25", "no2")) {
  pdf(file = paste0(output, "/", outc, "_overlap_of_latent_example_nolegend.pdf"), width = 2.5, height = 2)
  p1 <- t0 %>% 
    filter(outcome==outc, ratio.uy == 60, !is.na(BU)) %>% 
    ggplot() + geom_density(aes(BU, col=as.factor(X))) +
    theme(strip.text = element_text(size = 4)) +  
    scale_color_manual(
      values=c("#00B4F0", "#F8766D"), name = "", 
      labels = c("Non-redlined", "Redlined")) + 
    theme_bw() + labs(x="Latent U", y="Density") + theme(legend.position = "none") 
  
  if (outc == 'pm25') {
    p1 <- p1 + ggtitle(bquote(PM[2.5]))
  } else {
    p1 <- p1 + ggtitle(bquote(NO[2])) 
  } 
  print(p1)
  
  dev.off()
  
  pdf(file = paste0(output, "/", outc, "_overlap_of_latent_example_legend.pdf"), width = 4, height = 2)
  p1 <- t0 %>% 
    filter(outcome==outc, ratio.uy == 60, !is.na(BU)) %>% 
    ggplot() + geom_density(aes(BU, col=as.factor(X))) +
    theme(strip.text = element_text(size = 4)) +  
    scale_color_manual(
      values=c("#00B4F0", "#F8766D"), name = "", 
      labels = c("Non-redlined", "Redlined")) + 
    theme_bw() + labs(x="Latent U", y="Density") + theme(legend.title=element_blank()) 
  if (outc == 'pm25') {
    p1 <- p1 + ggtitle(bquote(PM[2.5]))
  } else {
    p1 <- p1 + ggtitle(bquote(NO[2])) 
  } 
  print(p1)
  dev.off()
}


# cor(t4$BU[t4$ratio.uy == 60], t5$BU[t5$ratio.uy==60])

for (y in c("pm25", "no2")) {
  pdf(file = paste0("output/redlining_info/", y, "_overlap_of_latent.pdf"), width = 16, height = 4)
  print(
    t0 %>% 
      filter(outcome==y, !is.na(BU)) %>% 
      ggplot() + geom_density(aes(BU, col=as.factor(X))) +
      theme(strip.text = element_text(size = 4)) + 
      theme_bw() + labs(title="", x="Latent U") + theme(legend.title=element_blank()) + 
      facet_wrap(. ~ ratio.uy, scales = "free", ncol = length(unique(t0$ratio.uy)))
  )
  dev.off()
}

for (y in c("pm25", "no2")) {
  pdf(file = paste0("output/redlining_info/", y, "_latent_and_proxy.pdf"), width = 16, height = 4)
  print(
    t0 %>% 
      filter(outcome==y, !is.na(BU)) %>% 
      ggplot() + geom_point(aes(BU, bc.rate.employed, col=as.factor(X))) +
      theme(strip.text = element_text(size = 4)) + 
      theme_bw() + labs(title="", x="Latent U", y="bc.rate.employed") + theme(legend.title=element_blank()) + 
      facet_wrap(. ~ ratio.uy, scales = "free", ncol = length(unique(t0$ratio.uy)))
  )
  print(
    t0 %>% 
      filter(outcome==y, !is.na(BU)) %>% 
      ggplot() + geom_point(aes(BU, wt.mean.rent.center, col=as.factor(X))) +
      theme(strip.text = element_text(size = 4)) + 
      theme_bw() + labs(title="", x="Latent U", y="wt.mean.rent.center") + theme(legend.title=element_blank()) + 
      facet_wrap(. ~ ratio.uy, scales = "free", ncol = length(unique(t0$ratio.uy)))
  )
  print(
    t0 %>% 
      filter(outcome==y, !is.na(BU)) %>% 
      ggplot() + geom_point(aes(BU, pct.black.pop40.rank, col=as.factor(X))) +
      theme(strip.text = element_text(size = 4)) + 
      theme_bw() + labs(title="", x="Latent U", y="pct.black.pop40.rank") + theme(legend.title=element_blank()) + 
      facet_wrap(. ~ ratio.uy, scales = "free", ncol = length(unique(t0$ratio.uy)))
  )
  print(
    t0 %>% 
      filter(outcome==y, !is.na(BU)) %>% 
      ggplot() + geom_point(aes(BU, BZ, col=as.factor(X))) +
      theme(strip.text = element_text(size = 4)) + 
      theme_bw() + labs(title="", x="Latent U", y="Spatial adjustment Z(s)") + theme(legend.title=element_blank()) + 
      facet_wrap(. ~ ratio.uy, scales = "free", ncol = length(unique(t0$ratio.uy)))
  )
  dev.off()
}


for (y in c("pm25", "no2")) {
  pdf(file = paste0("output/redlining_info/", y, "_spatial_latent.pdf"), width = 16, height = 4)
  print(
    t0 %>% 
      filter(outcome==y, !is.na(BZ)) %>% 
      ggplot() + geom_point(aes(BZ, bc.rate.employed, col=as.factor(X))) +
      theme(strip.text = element_text(size = 4)) + 
      theme_bw() + labs(title="", x="Spatial adjustment Z(s)", y="bc.rate.employed") + theme(legend.title=element_blank()) + 
      facet_wrap(. ~ ratio.uy, scales = "free", ncol = length(unique(t0$ratio.uy)))
  )
  print(
    t0 %>% 
      filter(outcome==y, !is.na(BZ)) %>% 
      ggplot() + geom_point(aes(BZ, wt.mean.rent.center, col=as.factor(X))) +
      theme(strip.text = element_text(size = 4)) + 
      theme_bw() + labs(title="", x="Spatial adjustment Z(s)", y="wt.mean.rent.center") + theme(legend.title=element_blank()) + 
      facet_wrap(. ~ ratio.uy, scales = "free", ncol = length(unique(t0$ratio.uy)))
  )
  print(
    t0 %>% 
      filter(outcome==y, !is.na(BZ)) %>% 
      ggplot() + geom_point(aes(BZ, pct.black.pop40.rank, col=as.factor(X))) +
      theme(strip.text = element_text(size = 4)) + 
      theme_bw() + labs(title="", x="Spatial adjustment Z(s)", y="pct.black.pop40.rank") + theme(legend.title=element_blank()) + 
      facet_wrap(. ~ ratio.uy, scales = "free", ncol = length(unique(t0$ratio.uy)))
  )
  dev.off()
}


################################################################
###################### Propensity Score ########################
################################################################
pm25_params <- readRDS("output/redlining_mcmc_latent_u1w3/pm25_params.rds")
pm25_params <- rbind(pm25_params, readRDS("output/redlining_sensitivity1_latent_u1w3/pm25_params.rds"))
pm25_params <- rbind(pm25_params, readRDS("output/redlining_sensitivity2_latent_u1w3/pm25_params.rds"))
pm25_params$outcome <- "pm25"
pm25_params$treatment <- c("AB/CD", "ABC/D", "A/BCD")

no2_params <- readRDS("output/redlining_mcmc_latent_u1w3/no2_params.rds")
no2_params <- rbind(no2_params, readRDS("output/redlining_sensitivity1_latent_u1w3/no2_params.rds"))
no2_params <- rbind(no2_params, readRDS("output/redlining_sensitivity2_latent_u1w3/no2_params.rds"))
no2_params$outcome <- "no2"
no2_params$treatment <- c("AB/CD", "ABC/D", "A/BCD")

params <- rbind(pm25_params, no2_params)
params <- params[(params$ratio.u == Inf)&(params$ratio.uy == 70), ]

qt <- c(0, .025, .25, .5, .75, .975, 1)
ps <- c()
for (i in 1:6) {
  
  load(paste0("output/redlining_mcmc_latent_u1w3/", params$outcome[i], "_1u_3w_sInf_s", params$ratio.uy[i], ".RData"))
  list2env(lapply(as.list(colMeans(model$out[model$out$post==1,])), as.vector), .GlobalEnv)

  U <- quantile(model$BU.mean, qt) # seq(-2, 3, length.out=11)
  
  l1 <- params$eta.1[i] + params$eta.2[i] * U
  prob.AB <- exp(l1) / (exp(l1) + 1)
  prob.CD <- 1 - prob.AB
  ps <- rbind(ps, unname(prob.CD))
}

ps <- cbind(params[, c("outcome", "treatment")], as.data.frame(ps))
names(ps)[3:ncol(ps)] <- paste0("Q", qt*100)

library(xtable)
print(xtable(ps), include.rownames = F)

################################################################
###################### Interpretation ##########################
################################################################
### example case is enough 
tb2[tb2$outcome=="pm25", c("ratio.u", "ratio.uy", "alpha.w1u", "alpha.w2u", "alpha.vu", "alpha.yu", "eta.2", "omega")]
tb2[tb2$outcome=="no2", c("ratio.u", "ratio.uy", "alpha.w1u", "alpha.w2u", "alpha.vu", "alpha.yu", "eta.2", "omega")]

# load("output/redlining_mcmc_latent_u1w3/pm25_1u_3w_sInf_s70.RData")
# load("output/redlining_mcmc_latent_u1w3/no2_1u_3w_sInf_s70.RData")
# list2env(lapply(as.list(colMeans(model$out[model$out$post==1,])), as.vector), .GlobalEnv)
# apply(model$out[model$out$post==1,], FUN=sd, MARGIN = 2)


###################################################################
########################## Sample Maps ############################
###################################################################
library(sf)
library(scales)
library(viridisLite)
library(patchwork)
library(RColorBrewer)
library(dplyr)
library(ggplot2)

load("data/nhgis0008_csv/only 1940/final_census1940_sp.Rdata") # dt 
geometry <- dt$geometry

dt$X <- dt$holc
dt$X[dt$X <= 2] = 1; dt$X[dt$X >= 3] = 2 ## binary treatment 

n <- nrow(dt)

for (id in c(52, 3, 35)) { # unique(dt$SPL)

  idx <- dt$SPL == id

  if(sum(idx)>30) {
    location <- dt$location[dt$SPL == id][1]
    city <- dt$city[dt$SPL == id][1]
    state <- dt$state[dt$SPL == id][1]

    # temp <- c(dt$Y[idx], dt$y.spline[idx], dt$y.theta[idx], dt$y.u[idx])
    # ylim <- c(-max(abs(min(temp)), max(temp)), max(abs(min(temp)), max(temp)))
    # col <- scale_fill_gradientn(colors = c(viridis(n = 5), rev(viridis(n = 5, option = "inferno"))), limit=ylim)

    p1 <- ggplot(data=dt$geometry[idx]) + ggtitle("Redlined Regions") +
      geom_sf(aes(fill=as.factor(dt$X[idx])), color="black") +
      theme_void() + # theme(legend.position = "None") + 
      scale_fill_manual(
        values=c("#00B4F0", "#F8766D"), name = "", 
        labels = c("Non-redlined", "Redlined"), 
        guide = guide_legend(keyheight = unit(.4, "cm"), keywidth = unit(0.5, "cm"))) +  
      theme(legend.key.size = unit(0.2, "cm"),  
        legend.text = element_text(size = 6))  


    p2 <- ggplot(data=dt$geometry[idx]) + ggtitle(bquote(PM[2.5] ~ "("~mu*"g/m"^3~")")) +             
      geom_sf(aes(fill=dt$pm25[idx]), color="black") +
      theme_void() + 
      scale_fill_gradientn(colors = rev(viridis(n = 10, option = "inferno")), 
                           guide = guide_colorbar(barwidth = 0.3, barheight = unit(4, "cm"))) + 
      theme(legend.key.size = unit(0.5, "cm"),  
            legend.text = element_text(size = 6), 
            legend.title = element_blank()) 

    p3 <- ggplot(data=dt$geometry[idx]) + ggtitle(bquote(NO[2] ~ "(ppb)")) +
      geom_sf(aes(fill=dt$no2[idx]), color="black") +
      theme_void() + 
      scale_fill_gradientn(colors = rev(viridis(n = 10, option = "inferno")), 
                           guide = guide_colorbar(barwidth = 0.3, barheight = unit(4, "cm"))) + 
      theme(legend.key.size = unit(0.5, "cm"),  
            legend.text = element_text(size = 6),  
            legend.title = element_blank()) 
  
    # TODO 
    # p3 <- ggplot(data=dt$geometry[idx]) + ggtitle(bquote(NO[2])) +
    #   geom_sf(aes(fill=dt$wt.mean.rent[idx]), color="black") +
    #   theme_void() + theme(legend.position = "None") +
    #   scale_fill_gradientn(colors = c(viridis(n = 5), rev(viridis(n = 5, option = "inferno"))))
    
    
    layout <- (p1|p2|p3) + plot_layout(widths = c(1, 1, 1)) # Adjust the widths

    outpath <- paste0("./output/redlining_info/maps_SPL", id, "_", city, "_", state, ".pdf")
    ggsave(outpath, plot = layout, width = 6, height = 2.5)
  }
}


for (id in c(52, 3, 35)) { # unique(dt$SPL)
  
  idx <- dt$SPL == id
  
  if(sum(idx)>30) {
    location <- dt$location[dt$SPL == id][1]
    city <- dt$city[dt$SPL == id][1]
    state <- dt$state[dt$SPL == id][1]
    
    p1 <- ggplot(data=dt$geometry[idx]) + ggtitle("Unemployment (%)") +
      geom_sf(aes(fill=100-dt$rate.employed[idx]*100), color="black") +
      theme_void() + 
      scale_fill_gradientn(colors = rev(viridis(n = 10, option = "inferno")), 
                           guide = guide_colorbar(barwidth = 0.3, barheight = unit(4, "cm"))) + 
      theme(legend.key.size = unit(0.5, "cm"),   # Adjust this value to make the legend smaller
            legend.text = element_text(size = 6), # Adjust the text size of the legend
            legend.title = element_blank()) 
    
    p2 <- ggplot(data=dt$geometry[idx]) + ggtitle('Mean House Rent ($)') +
      geom_sf(aes(fill=dt$wt.mean.rent[idx]), color="black") +
      theme_void() +  
      scale_fill_gradientn(colors = rev(viridis(n = 10, option = "inferno")), 
                           guide = guide_colorbar(barwidth = 0.3, barheight = unit(4, "cm"))) + 
      theme(legend.key.size = unit(0.5, "cm"),  
            legend.text = element_text(size = 6), 
            legend.title = element_blank()) 
    
    p3 <- ggplot(data=dt$geometry[idx]) + ggtitle('Black Population (%)') +
      geom_sf(aes(fill=dt$pct.black.pop40[idx]*100), color="black") +
      theme_void() + 
      scale_fill_gradientn(colors = rev(viridis(n = 10, option = "inferno")), 
                           guide = guide_colorbar(barwidth = 0.3, barheight = unit(4, "cm"))) + 
      theme(legend.key.size = unit(0.5, "cm"),   
            legend.text = element_text(size = 6), 
            legend.title = element_blank()) 
    
    layout <- (p1|p2|p3) + plot_layout(widths = c(1, 1, 1)) # Adjust the widths
    
    outpath <- paste0("./output/redlining_info/maps_SPL", id, "_", city, "_", state, "_proxy.pdf")
    ggsave(outpath, plot = layout, width = 6, height = 2.5)
  }
}


########################################
ratio.u <- Inf
for (outcome in c("pm25", "no2")) {
  for (ratio.uy in c(0, .1, .7)) {
    
    ## load spatial data !!! be careful using these values in dt for visualization since they are not transformed 
    load("data/nhgis0008_csv/only 1940/final_census1940_sp.Rdata") # dt 
    geometry <- dt$geometry; rm(dt) 
    
    ## load parameter set 
    if (ratio.uy == 0) {
      load(paste0("output/redlining_mcmc_noproxy/", outcome, "_noproxy_s", ratio.u*100, "_s", ratio.uy*100, ".RData")) # !!! use dt here       
    } else {
      load(paste0("output/redlining_mcmc_latent_u1w3/", outcome, "_1u_3w_s", ratio.u*100, "_s", ratio.uy*100, ".RData")) # !!! use dt here 
    }

    dt$geometry <- geometry
    
    n <- nrow(dt)
    
    list2env(as.list(colMeans(model$out[model$out$post==1, ])), envir = .GlobalEnv)
    
    thetaMat <- matrix(data=0, nrow=n, ncol=2)
    thetaMat[cbind(1:n, dt$X)] <- 1
    
    dt$y.spline <- omega*model$BZ.mean
    dt$y.theta <- thetaMat %*% t(t(c(theta.1, theta.2))) ## make sure we get the right residual 
    dt$y.u <- alpha.yu*model$BU.mean
    dt$y.fit <- dt$y.spline + dt$y.theta + dt$y.u
    dt$resid <- dt$Y - dt$y.fit
    
    dt$y.theta <- (theta.2-theta.1)*(dt$X==2) ## recalculate for visualization 
    
    U <- model$BU.mean
    
    for (id in c(52, 3, 35)) { # unique(dt$SPL) 52 is Philadelphia
      
      idx <- dt$SPL == id
      
      if(sum(idx)>30) {
        location <- dt$location[dt$SPL == id][1]
        city <- dt$city[dt$SPL == id][1]
        state <- dt$state[dt$SPL == id][1]
        
        ss1 <- data.frame(dt$y.spline, dt$y.u, dt$y.theta, dt$resid, dt$X)[idx, ]
        colnames(ss1) <- c("spline", "latent", "effect", "resid", "A")
        ss1$id <- as.factor(1:nrow(ss1))
        
        ss1_long <- tidyr::gather(ss1, key = "Source", value = "Value", -id, -A)
        ordered_cities <- ss1_long %>% dplyr::filter(Source == "spline") %>% dplyr::arrange(A, Value) %>% dplyr::pull(id)
        ss1_long$id <- factor(ss1_long$id, levels = ordered_cities)
        ss1_long$Source <- factor(ss1_long$Source, levels = c("latent", "resid", "effect", "spline"))
        
        if (outcome=="pm25") {
          p1 <- ggplot(ss1_long, aes(x = id, y = Value, fill = Source)) +
            geom_bar(stat = "identity", position = "stack") +
            labs(title = bquote(PM[2.5] ~ "(" ~ r ~ "=" ~ .(ratio.uy) ~ ")"))
            # labs(title = bquote(PM[2.5] ~" Decomposition, "~ .(location)))
            # ~ " (" ~ mu * "g/m"^3 ~ "), "
        } else {
          p1 <- ggplot(ss1_long, aes(x = id, y = Value, fill = Source)) +
            geom_bar(stat = "identity", position = "stack") +
            labs(title = bquote(NO[2] ~ "(" ~ r ~ "=" ~ .(ratio.uy) ~ ")"))
        } 
        
        p1.nolegend <- p1 +  theme_minimal() +
          scale_fill_manual(
            values=c("#00B4F0", "#B79F00", "#F8766D", "#7CAE00"), name = "",
            labels = c("Latent U", "Residual", "Redlining Effect", "Spatial Adjustment")) +
          geom_vline(aes(xintercept = sum(ss1$A==1)+.5), col="red") +
          theme(text = element_text(size = 12),
                legend.position = "none",
                axis.ticks = element_blank(),
                axis.text.x = element_blank(),
                strip.text = element_text(face="bold", size=12),
                strip.background = element_blank()) + ylab("") + xlab("")

        outpath <- paste0("./output/redlining_info/map_", outcome, "_SPL", id, "_", city, "_", state, "_s", ratio.uy*100, "_nolegend.pdf")
        pdf(outpath, width = 3.8, height = 4)
        print(p1.nolegend)
        dev.off()
        
        p1.legend <- p1 + theme_minimal() + 
          scale_fill_manual(
            values=c("#00B4F0", "#B79F00", "#F8766D", "#7CAE00"), name = "", 
            labels = c("Latent U", "Residual", "Redlining Effect", "Spatial Adjustment")) + 
          geom_vline(aes(xintercept = sum(ss1$A==1)+.5), col="red") + 
          theme(text = element_text(size = 12),
                legend.title = element_blank(),
                axis.ticks = element_blank(),
                axis.text.x = element_blank(),
                strip.text = element_text(face="bold", size=12),
                strip.background = element_blank()) + ylab("") + xlab("")
  
        outpath <- paste0("./output/redlining_info/map_", outcome, "_SPL", id, "_", city, "_", state, "_s", ratio.uy*100, "_legend.pdf")
        pdf(outpath, width = 6, height = 4)
        print(p1.legend)
        dev.off()
      }
    }
  }
}




###################################################################
########################## Map Results ############################
###################################################################
ratio.u <- Inf

for (outcome in c("pm25", "no2")) {
  for (ratio.uy in c(.1, .3, .5)) {
    
    ## load spatial data !!! be careful using these values in dt for visualization since they are not transformed 
    load("data/nhgis0008_csv/only 1940/final_census1940_sp.Rdata") # dt 
    geometry <- dt$geometry; rm(dt) 
    
    ## load parameter set 
    load(paste0("output/redlining_mcmc_latent_u1w3/", outcome, "_1u_3w_s", ratio.u*100, "_s", ratio.uy*100, ".RData")) # !!! use dt here 
    dt$geometry <- geometry
    
    n <- nrow(dt)
    
    list2env(as.list(colMeans(model$out[model$out$post==1, ])), envir = .GlobalEnv)
    
    thetaMat <- matrix(data=0, nrow=n, ncol=2)
    thetaMat[cbind(1:n, dt$X)] <- 1
    
    dt$y.spline <- omega*model$BZ.mean
    dt$y.theta <- thetaMat %*% t(t(c(theta.1, theta.2))) ## make sure we get the right residual 
    dt$y.u <- alpha.yu*model$BU.mean
    dt$y.fit <- dt$y.spline + dt$y.theta + dt$y.u
    dt$resid <- dt$Y - dt$y.fit
    
    dt$y.theta <- (theta.2-theta.1)*(dt$X==2) ## recalculate for visualization 
    
    U <- model$BU.mean
    
    for (id in unique(dt$SPL)) {
      
      idx <- dt$SPL == id
      
      if(sum(idx)>30) {
        location <- dt$location[dt$SPL == id][1]
        city <- dt$city[dt$SPL == id][1]
        state <- dt$state[dt$SPL == id][1]
        
        ss1 <- data.frame(dt$y.spline, dt$y.u, dt$y.theta, dt$resid, dt$X)[idx, ]
        colnames(ss1) <- c("spline", "latent", "effect", "resid", "A")
        ss1$id <- as.factor(1:nrow(ss1))
        
        ss1_long <- tidyr::gather(ss1, key = "Source", value = "Value", -id, -A)
        ordered_cities <- ss1_long %>% dplyr::filter(Source == "spline") %>% dplyr::arrange(A, Value) %>% dplyr::pull(id)
        ss1_long$id <- factor(ss1_long$id, levels = ordered_cities)
        ss1_long$Source <- factor(ss1_long$Source, levels = c("latent", "resid", "effect", "spline"))
        
        p1 <- ggplot(ss1_long, aes(x = id, y = Value, fill = Source)) +
          geom_bar(stat = "identity", position = "stack") +
          labs(title = bquote(PM[2.5] ~" Decomposition, "~ .(location))) +  # ~ " (" ~ mu * "g/m"^3 ~ "), "
          theme_minimal() + 
          scale_fill_manual(
            values=c("#00B4F0", "#B79F00", "#F8766D", "#7CAE00"), name = "", 
            labels = c("Latent U", "Residual", "Redlining Effect", "Spatial Adjustment")) + 
          geom_vline(aes(xintercept = sum(ss1$A==1)+.5), col="red") +
          theme(text = element_text(size = 12),
                legend.title = element_blank(),
                # panel.grid.major = element_blank(),
                # legend.position = "none",
                axis.ticks = element_blank(),
                axis.text.x = element_blank(),
                strip.text = element_text(face="bold", size=12),
                strip.background = element_blank()) + ylab("") + xlab("") 
        
        # geom_segment(aes(x = sum(ss1$A==1)+.5, y = max(ss1_long$Value), xend = sum(ss1$A==1)+12, yend = max(ss1_long$Value)),
        #              arrow = arrow(type = "closed", length = unit(0.04, "inches")), color = "gray") + 
        # annotate("text", x = sum(ss1$A==1)+12, y = max(ss1_long$Value), label = "Redlined", hjust = 0, vjust = 0, col="gray")
        
        temp <- c(dt$Y[idx], dt$y.spline[idx], dt$y.theta[idx], dt$y.u[idx])
        ylim <- c(-max(abs(min(temp)), max(temp)), max(abs(min(temp)), max(temp)))
        col <- scale_fill_gradientn(colors = c(viridis(n = 5), rev(viridis(n = 5, option = "inferno"))), limit=ylim)
        
        # library(viridis); viridis(9)
        p2 <- ggplot(data=dt$geometry[idx]) + ggtitle(bquote(PM[2.5] * " Explained by Spatial Adjustment")) + 
          geom_sf(aes(fill=dt$y.spline[idx]), color="black") + 
          theme_void() + col + theme(legend.position = "None")
        
        p3 <- ggplot(data=dt$geometry[idx]) + ggtitle(bquote(PM[2.5])) + 
          geom_sf(aes(fill=dt$Y[idx]), color="black") + 
          theme_void() + col + theme(legend.position = "None") # TODO 
        
        # c("#fb6964", "#6ea524", "#00b8bc", "#c26efa")
        p4 <- ggplot(data=dt$geometry[idx]) + ggtitle("Redlined Regions") + 
          geom_sf(aes(fill=as.factor(dt$X[idx])), color="black") + 
          theme_void() + theme(legend.position = "None")
        
        ### 
        col2 <- scale_fill_gradientn(colors = c(viridis(n = 5), rev(viridis(n = 5, option = "inferno"))))
        
        p5 <- ggplot(data=dt$geometry[idx]) + ggtitle(bquote(PM[2.5] * " Explained by Spatial Adjustment")) + 
          geom_sf(aes(fill=U[idx,1]), color="black") + 
          theme_void() + col2 + theme(legend.position = "None")
        
        layout <- (p4|p3|p1) +  # p2|
          plot_layout(widths = c(1, 1, 1)) # Adjust the widths
        
        outpath <- paste0("/output/redlining_info/redlining_mcmc_latent_u1w3/", outcome, "_s", ratio.u*100, "_s", ratio.uy*100, 
                          "_SPL", id, "_", city, "_", state, ".pdf")
        ggsave(outpath, plot = layout, width = 10, height = 4)
        
        # layout <- (p4+p3+p2+p1+p5) # Adjust the widths
        # outpath <- paste0("./output/redlining_info/temp/region", id, "_s", ratio*100, ".pdf")
        # ggsave(outpath, plot = layout, width = 10, height = 4) 
      }
    }
  }
}


################################################################
######################### trace plot ###########################
################################################################
for (outcome in c("pm25", "no2")) {
  pdf(paste0("output/redlining_plots/", outcome, "_trace.pdf"), width = 6, height = 4)
  for (ratio.uy in c(10,40,70)) {
    load(paste0("output/redlining_mcmc_latent_u1w3/", outcome, "_1u_3w_sInf_s", ratio.uy, ".RData"))
    trace_value <- (model$out$theta.2 - model$out$theta.1)[model$out$post==1]
    plot(trace_value, xlab="", ylab="", type="l")
  }
  dev.off()
}


####################################################################################
################################# trace plot #######################################
####################################################################################
# load("output/redlining_mcmc_latent_u1w3/pm25_1u_3w_sInf_s90.RData")
# out <- model$out[model$out$post==1,]
# n <- nrow(out)
# params <- names(out)[2:ncol(out)]
# out <- out %>% tidyr::gather(key = "parameter", value, -iteration)
# 
# plot.list <- list()
# for (i in 1:length(params)) {
#   plot.list[[i]] <- out %>% 
#     dplyr::filter(parameter == params[i]) %>% 
#     ggplot(.) + 
#     geom_line(aes(1:n, value), size=.3) +
#     theme_classic() + 
#     ggtitle(params[i]) + 
#     theme(panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank(),
#           axis.title.x = element_blank(),  
#           axis.title.y = element_blank(),
#           legend.position = "none") 
# }
# 
# outpath <- "output/redlining_info/test_redlining_mix.pdf"
# pdf(outpath, width = 6, height = .5 * length(params))
# do.call('grid.arrange',c(plot.list, ncol = 2))
# dev.off()

pm25_params <- readRDS("output/redlining_mcmc_latent_u1w3/pm25_params.rds")
pm25_params <- rbind(pm25_params, readRDS("output/redlining_sensitivity1_latent_u1w3/pm25_params.rds"))
pm25_params <- rbind(pm25_params, readRDS("output/redlining_sensitivity2_latent_u1w3/pm25_params.rds"))
pm25_params$outcome <- "pm25"
pm25_params$treatment <- c("AB/CD", "ABC/D", "A/BCD")

no2_params <- readRDS("output/redlining_mcmc_latent_u1w3/no2_params.rds")
no2_params <- rbind(no2_params, readRDS("output/redlining_sensitivity1_latent_u1w3/no2_params.rds"))
no2_params <- rbind(no2_params, readRDS("output/redlining_sensitivity2_latent_u1w3/no2_params.rds"))
no2_params$outcome <- "no2"
no2_params$treatment <- c("AB/CD", "ABC/D", "A/BCD")

params <- rbind(pm25_params, no2_params)
params <- params[(params$ratio.u == Inf)&(params$ratio.uy == 70), ]

xtable(params[,c("outcome", "treatment", "alpha.w1u", "alpha.w2u", "alpha.vu", "eta.2", "alpha.yu")], include.rownames = F)



################################################################
######################### diagnostic ###########################
################################################################
load("data/nhgis0008_csv/only 1940/final_census1940_sp.Rdata") # dt
geometry <- dt$geometry

## load one data
load("output/redlining_mcmc_latent_u1w3/pm25_1u_3w_sInf_s70.RData")
dt$geometry <- geometry

params <- readRDS("output/redlining_mcmc_latent_u1w3/pm25_params.rds")
list2env(lapply(as.list(params[(params$ratio.u == Inf)&(params$ratio.uy == 70), ]), as.vector), .GlobalEnv)

# forget to center
dt$W1 <- dt$W1 - mean(dt$W1, na.rm=T)
dt$W2[is.na(dt$W2)] <- mean(dt$W2, na.rm=T)
mean(dt$W1)
sd(dt$W1)
mean(dt$W2)
sd(dt$W2)

## get fitted values of proxy, treatment, outcome?
BU <- model$BU.mean
BZ <- model$BZ.mean

dt$W1.fit <- alpha.w1 + alpha.w1u * BU; dt$W1.resid <- dt$W1 - dt$W1.fit
dt$W2.fit <- alpha.w2 + alpha.w2u * BU; dt$W2.resid <- dt$W2 - dt$W2.fit

## which one?
dt$V.fit <- alpha.v + alpha.vu * BU; dt$V.resid <- dt$V - dt$V.fit
dt$V.fit <- pmax(alpha.v + alpha.vu * BU, min(dt$V)); dt$V.resid <- dt$V - dt$V.fit

dt$BU <- BU
dt$BZ <- BZ

dt$Y.fit <- theta.1*(dt$X==1) + theta.2*(dt$X==2) + alpha.yu * BU + omega * BZ
dt$X.probAB <- exp(eta.1 + eta.2 * BU + BZ) / (1 + exp(eta.1 + eta.2 * BU + BZ))
dt$X.probCD <- 1 - dt$X.probAB

plot(dt$W1, dt$W1.fit)
plot(dt$W2, dt$W2.fit)
plot(dt$V, dt$V.fit)

plot(dt$W1.resid, dt$W2.resid)
plot(dt$W1.resid, dt$V.resid)
plot(dt$V.resid, dt$W2.resid)

## If we only use two proxy
## we hope they are correlated only through U
##
cor.test(dt$W1.resid, dt$W2.resid) # -0.6540171
cor.test(dt$W1.resid, dt$V.resid) # -0.1396984
cor.test(dt$V.resid, dt$W2.resid) # 0.2700661 # when using these two proxy, quite biased

cor.test(dt$W1, dt$W2) # 0.747508
cor.test(dt$W1, dt$V) # -0.2161492
cor.test(dt$V, dt$W2) # -0.04474902 # <-- this too week, may find a week U

eta.2 # 2.02

plot(dt$W1, dt$W1.resid)
plot(dt$W2, dt$W2.resid)
plot(dt$V, dt$V.resid)

idx <- dt$SPL == 62 # 52

idx <- dt$location == "Seattle WA"

ggplot(data=dt$geometry[idx]) + ggtitle(bquote(PM[2.5])) +
  geom_sf(aes(fill=dt$W1[idx]), color="black") +
  theme_void() + theme(legend.position = "None") +
  scale_fill_gradientn(colors = c(viridis(n = 5), rev(viridis(n = 5, option = "inferno"))))

ggplot(data=dt$geometry[idx]) + ggtitle(bquote(PM[2.5])) +
  geom_sf(aes(fill=dt$W1.resid[idx]-mean(dt$W1.resid[idx])), color="black") +
  theme_void() + theme(legend.position = "None") +
  scale_fill_gradientn(colors = c(viridis(n = 5), rev(viridis(n = 5, option = "inferno"))))

ggplot(data=dt$geometry[idx]) + ggtitle(bquote(PM[2.5])) +
  geom_sf(aes(fill=dt$W2.resid[idx]-mean(dt$W2.resid[idx])), color="black") +
  theme_void() + theme(legend.position = "None") +
  scale_fill_gradientn(colors = c(viridis(n = 5), rev(viridis(n = 5, option = "inferno"))))

ggplot(data=dt$geometry[idx]) + ggtitle(bquote(PM[2.5])) +
  geom_sf(aes(fill=dt$V.resid[idx]), color="black") +
  theme_void() + # theme(legend.position = "None") +
  scale_fill_gradientn(colors = c(viridis(n = 5), rev(viridis(n = 5, option = "inferno"))))

ggplot(data=dt$geometry[idx]) + ggtitle(bquote(PM[2.5])) +
  geom_sf(aes(fill=dt$BZ[idx]), color="black") +
  theme_void() + # theme(legend.position = "None") +
  scale_fill_gradientn(colors = c(viridis(n = 5), rev(viridis(n = 5, option = "inferno"))))

ggplot(data=dt$geometry[idx]) + ggtitle(bquote(PM[2.5])) +
  geom_sf(aes(fill=dt$BU[idx]), color="black") +
  theme_void() + # theme(legend.position = "None") +
  scale_fill_gradientn(colors = c(viridis(n = 5), rev(viridis(n = 5, option = "inferno"))))


