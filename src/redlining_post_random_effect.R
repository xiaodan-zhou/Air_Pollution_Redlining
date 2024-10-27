library(ggplot2)
library(dplyr)

boxxy <- function(x) {
  r <- quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  return(r)
}

outside <- function(x) {
  out.x <- subset(x, x < hpd(x)[1] | x > hpd(x)[2])
  out.x[seq(1, length(out.x), by = 15)] # thinning to length 1000
}

hpd <- function(x, alpha = 0.05){ 
  # Highest Posterior Density intervals
  n <- length(x)
  m <- round(n * alpha)
  x <- sort(x)
  y <- x[(n - m + 1):n] - x[1:m]
  z <- min(y)
  k <- which(y == z)[1]
  c(x[k], x[n - m + k])
}

### redlining data 
dt <- readRDS("data/nhgis0008_csv/cleaned/final_census1940.rds")
cities <- unique(dt$city.unique)

n.u <- 1 
ratio.u <- Inf 
stt <- c()

trtname <- list("redlining_mcmc_latent_u1w3_random_effect" = "AB.CD", 
                "redlining_sensitivity1_latent_u1w3_random_effect" = "ABC.D", 
                "redlining_sensitivity2_latent_u1w3_random_effect" = "A.BCD")

for (subfolder in c("redlining_mcmc_latent_u1w3_random_effect", 
               "redlining_sensitivity1_latent_u1w3_random_effect", 
               "redlining_sensitivity2_latent_u1w3_random_effect")
     ) {
  for (outcome in c("pm25", "no2")) {
    for (ratio.uy in c(.6)) {
      file_path <- paste0("output/", subfolder, "/", outcome, "_1u_3w_sInf_s", ratio.uy*100, ".Rdata")
      
      info <- file.info(file_path)
      
      if (is.na(info$size)) { warning(paste0("File does not exist:", file_path)) 
      } else if (info$size < .1) { print(paste0("The file is empty", file_path)) 
      } else { skip_to_next <- TRUE 
      tryCatch(load(file_path), error = function(e) { skip_to_next <<- FALSE})
      if(skip_to_next) {
        outputs <- model
        
        pooled.effects <- data.frame(
          Pooled = model$out$mu.theta[model$out$post==1] - model$out$theta0[model$out$post==1]) 
        
        city.effects <- data.frame(
          model$out[model$out$post==1, paste0("theta.", 1:69)] - model$out$theta0[model$out$post==1])
        
        names(city.effects) <- unique(dt$location)
        
        all.effects <- cbind(city.effects, pooled.effects)
        all.effects <- pivot_longer(all.effects, cols = everything(), names_to = "location", values_to = "value")
        
        x.effect <- data.frame(location=c(names(city.effects), "Pooled"), 
                               x=c(unname(colMeans(city.effects)), mean(pooled.effects$Pooled)))
        
        all.effects <- merge(all.effects, x.effect, by="location")
        
        all.effects$location <- reorder(all.effects$location, -all.effects$x)
        
        label_colors <- ifelse(levels(reorder(all.effects$location, -all.effects$x)) == "Pooled", "red", "black")
        
        pdf(paste0("output/redlining_plots/", outcome, "_", trtname[subfolder], "_s", ratio.uy*100, "_random_effects.pdf"), 
            width = 6, height = 10)
        
        print(
          ggplot() + 
            stat_summary(fun.data = boxxy, geom = "boxplot", 
                         data = mutate(all.effects, location = reorder(all.effects$location, -x), 
                                       box_col = location=="Pooled"),
                         mapping = aes(x = location, y = value, fill = box_col), lwd = 0.4, width = 0.8) + 
            stat_summary(fun = outside, geom = "point",
                         data = mutate(all.effects, location = reorder(all.effects$location, -x), 
                                       box_col = location=="Pooled"),
                         mapping = aes(x = location, y = value, fill = box_col), size = 0.6) +
            geom_hline(yintercept = 0, color = "blue") + 
            theme_bw() + 
            coord_flip() + 
            theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                  axis.title.y=element_blank(),
                  axis.text.y = element_text(color = label_colors),
                  legend.position = "none") + 
            labs(title="", x="", y="") + 
            scale_fill_manual(values = c("#00BA38", "red")) +
            ggtitle("")
        )
        
        dev.off()
        
      }}}}}




for (subfolder in c("redlining_mcmc_latent_u1w3_random_effect", 
                    "redlining_sensitivity1_latent_u1w3_random_effect", 
                    "redlining_sensitivity2_latent_u1w3_random_effect")
) {
  for (outcome in c("pm25", "no2")) {
    for (ratio.uy in c(.5,.6,.7)) {
      file_path <- paste0("output/", subfolder, "/", outcome, "_1u_3w_sInf_s", ratio.uy*100, ".Rdata")
      
      info <- file.info(file_path)
      
      if (is.na(info$size)) { warning(paste0("File does not exist:", file_path)) 
      } else if (info$size < .1) { print(paste0("The file is empty", file_path)) 
      } else { skip_to_next <- TRUE 
      tryCatch(load(file_path), error = function(e) { skip_to_next <<- FALSE})
      if(skip_to_next) {
        outputs <- model
        
        ## pooled effect
        pooled.effects <- model$out$mu.theta[model$out$post==1] - model$out$theta0[model$out$post==1]
        pooled <- mean(pooled.effects)
        pooled.high <- quantile(pooled.effects, .975)
        pooled.low <- quantile(pooled.effects, .025)
        pooled.q25 <- quantile(pooled.effects, .25)
        pooled.q75 <- quantile(pooled.effects, .75)
        
        ## city effect
        city.effects <- model$out[model$out$post==1, paste0("theta.", 1:69)] - model$out$theta0[model$out$post==1]
        city.effect <- unname(colMeans(city.effects))
        city.high <- unname(unlist(lapply(city.effects, quantile, probs = 0.975)))
        city.low <- unname(unlist(lapply(city.effects, quantile, probs = 0.025)))
        city.q25 <- unname(unlist(lapply(city.effects, quantile, probs = 0.25)))
        city.q75 <- unname(unlist(lapply(city.effects, quantile, probs = 0.75)))
        
        # matrix(c(1,2,3,4,5,6), nrow=2) - c(1,-1)
        
        ## city effect and pooled effect
        tb <- data.frame(high = c(city.high, pooled.high), 
                         low  = c(city.low, pooled.low),
                         q25 = c(city.q25, pooled.q25),
                         q75 = c(city.q75, pooled.q75),
                         estimate = c(city.effect, pooled),
                         city = c(cities, ""), 
                         pooled= c(rep(F, length(cities)), T)
        )
        
        tb$treatment=unname(trtname[subfolder])[[1]]
        tb$outcome=outcome
        tb$ratio.uy=ratio.uy*100
        tb$n.high=sum(model$effect.higher >= .95)
        tb$n.low=sum((1-model$effect.higher) >= .95)
        tb$waic=outputs$WAIC$WAIC[4]
        
        stt <- rbind(stt, tb)
        
      }}}}}

stt <- data.frame(stt)
saveRDS(stt,file="output/redlining_mcmc_latent_u1w3_random_effect/combined.rds")

###################################################################################################
for (outc in c("pm25", "no2")) {
  for (trt in c("AB.CD", "ABC.D", "A.BCD")) {
    
    stt <- readRDS("output/redlining_mcmc_latent_u1w3_random_effect/combined.rds")
    stt$city[stt$city == ""] <- "Population Mean"
    nx <- length(cities)+1
    rz <- 60
    
    ### select the case that decide the order of cities
    temp <- stt %>% filter(treatment==trt, outcome==outc, ratio.uy==rz)
    temp <- temp[order(temp$estimate),]
    city_order <- unique(temp$city)

    stt$city <- factor(stt$city, levels = city_order)
    stt <- stt %>% arrange(city)

    sum(is.na(stt$city))
    
    ## color texts     
    x_colors <- ifelse(unique(stt$city) == "Population Mean", "red", "black")
    names(x_colors) <- unique(stt$city)

    pdf(paste0("output/redlining_plots/", outc, "_", trt, "_s", rz, "_random_effects_v2.pdf"),
        width = 5, height = 10)
    
    p1 <- stt %>%
      filter(treatment==trt, outcome==outc, ratio.uy==rz) %>%
      ggplot() + geom_point(aes(x = city, y = estimate, col = x_colors)) +
      geom_hline(yintercept = 0, linetype = 2) +
      geom_errorbar(aes(x = city, ymin = low, ymax = high), width = 0.2, col = x_colors) +
      theme_bw() + 
      scale_color_identity() + # use x_colors as color name directly 
      coord_flip() + 
      theme(legend.position = "None",
            axis.text.y = element_text(colour = x_colors)) # after coord_flip, use y instead of x 
      
    
    if (outc == 'pm25') {
      p1 <- p1 + labs(title = bquote(PM[2.5]), y = "Estimate and 95% Credible Interval", x = "")
    } else {
      p1 <- p1 + labs(title = bquote(NO[2]), y = "Estimate and 95% Credible Interval", x = "")
    }   
    
    print(p1) 
    
    dev.off()
  }
}



####################################################################################################
library(sf)
library(ggplot2)
library(rnaturalearth) # 'naturalearth' package to get the state boundaries
library(rnaturalearthdata)
library(dplyr)

states <- ne_states(country = "united states of america", returnclass = "sf")
states <- states[!(states$name %in% c("Alaska", "Hawaii")),]

rz <- 60

for (outc in c("pm25", "no2")) {
  for (trt in c("AB.CD", "ABC.D", "A.BCD")) {
    load("data/nhgis0008_csv/cleaned/final_census1940.RData")
    df <- st_centroid(dt$geometry[!duplicated(dt$location)])
    df <- as.data.frame(df)
    df$city <- dt$city.unique[!duplicated(dt$city.unique)]
    
    stt <- readRDS("output/redlining_mcmc_latent_u1w3_random_effect/combined.rds")
    stt <- stt %>% dplyr::filter(treatment==trt, outcome==outc, ratio.uy==rz)
    stt <- stt[stt$city != "", ]
    df <- merge(df, stt, by="city"); rm(stt)
    
    df$size <- (abs(df$estimate)/max(abs(df$estimate)))^2*10
    # df$alpha <- (!(df$q25 < 0 & df$q75 < 0))/4+.7
    df$alpha <- (!(df$low < 0 & df$high < 0))/4+.75
    
    if(length(unique(df$alpha))==1) df$alpha <- 1
    pdf(paste0("output/redlining_plots/map_", outc, "_", trt, "_s", rz, "_random_effects.pdf"), width = 5, height = 3)
    
    p1 <- ggplot() +
      geom_sf(data = states, fill = NA, color = "gray") + 
      geom_sf(data = df[df$estimate<=0,], 
              aes(geometry = geometry, size = size), color = "blue", shape = 1) + 
      geom_sf(data = df[df$estimate>0,], 
              aes(geometry = geometry, size = size), color = "red", shape = 1) + 
      theme_minimal() + 
      theme(axis.line = element_blank(),
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.text = element_blank(), 
            axis.title = element_blank(), 
            axis.ticks = element_blank(), 
            legend.position = "None"
      ) 
    if (outc == 'pm25') {
      p1 <- p1 + labs(title = bquote(PM[2.5]), x = "", y = "")
    } else {
      p1 <- p1 + labs(title = bquote(NO[2]), x = "", y = "")
    }
    
    print(p1)
    
    dev.off()

    
  }
}

##### some data for paper writing 
trt <- "ABC.D"
stt <- readRDS("output/redlining_mcmc_latent_u1w3_random_effect/combined.rds")
stt$city[stt$city == ""] <- "Pooled"

stt[(stt$ratio.uy==60)&(stt$treatment==trt)&(stt$city=="Pooled"), c("estimate", "low", "high", "outcome")]

stt[(stt$ratio.uy==60)&(stt$treatment==trt)&(stt$city!="Pooled"),] %>% 
  group_by(outcome) %>% 
  summarise(sig.pos = sum(low>0), sig.neg = sum(high<0))

stt[(stt$ratio.uy==60)&(stt$treatment==trt)&(stt$city!="Pooled"),] %>% 
  mutate(sig.pos = low>0, sig.neg = high<0) %>% 
  filter((sig.pos==1)|(sig.neg==1)) %>% 
  select(city, outcome, estimate, sig.pos, sig.neg)










########################################
ratio.u <- Inf
for (outcome in c("pm25", "no2")) {
  for (ratio.uy in c(0.6)) {
    
    ## load spatial data !!! be careful using these values in dt for visualization since they are not transformed 
    load("data/nhgis0008_csv/cleaned/final_census1940_sp.Rdata") # dt 
    geometry <- dt$geometry; rm(dt) 
    
    ## load parameter set 
    if (ratio.uy == 0) {
      load(paste0("output/redlining_mcmc_noproxy/", outcome, "_noproxy_s", ratio.u*100, "_s", ratio.uy*100, ".RData")) # !!! use dt here       
    } else {
      load(paste0("output/redlining_mcmc_latent_u1w3_random_effect/", outcome, "_1u_3w_s", ratio.u*100, "_s", ratio.uy*100, ".RData")) # !!! use dt here 
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
            labels = c("Latent Process U", "Residual", "Redlining Effect", "Spatial Latent Z")) +
          geom_vline(aes(xintercept = sum(ss1$A==1)+.5), col="red") +
          theme(text = element_text(size = 12),
                legend.position = "none",
                axis.ticks = element_blank(),
                axis.text.x = element_blank(),
                strip.text = element_text(face="bold", size=12),
                strip.background = element_blank()) + ylab("") + xlab("")
        
        outpath <- paste0("./output/redlining_info/map_", outcome, "_SPL", id, "_", city, "_", state, "_s", ratio.uy*100, "_nolegend_random_effect.pdf")
        pdf(outpath, width = 3.8, height = 4)
        print(p1.nolegend)
        dev.off()
        
        p1.legend <- p1 + theme_minimal() + 
          scale_fill_manual(
            values=c("#00B4F0", "#B79F00", "#F8766D", "#7CAE00"), name = "", 
            labels = c("Latent Process U", "Residual", "Redlining Effect", "Spatial Latent Z")) +
          geom_vline(aes(xintercept = sum(ss1$A==1)+.5), col="red") + 
          theme(text = element_text(size = 12),
                legend.title = element_blank(),
                axis.ticks = element_blank(),
                axis.text.x = element_blank(),
                strip.text = element_text(face="bold", size=12),
                strip.background = element_blank()) + ylab("") + xlab("")
        
        outpath <- paste0("./output/redlining_info/map_", outcome, "_SPL", id, "_", city, "_", state, "_s", ratio.uy*100, "_legend_random_effect.pdf")
        pdf(outpath, width = 6, height = 4)
        print(p1.legend)
        dev.off()
      }
    }
  }
}


