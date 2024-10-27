library(maps)
library(splines)
library(ggplot2)
library(tidyr)
library(sf)

is_integer_zero <- function(x) is.integer(x) && length(x) == 0

replace_empty_with_NA <- function(x) {
  if (length(x) == 0) return(NA)
  return(x)
}

replace_len2_with_first <- function(x) {
  if (length(x) == 2) return(x[1])
  return(x)
}

get.spline <- function(SPL, geometry, L1=NA, L2=NA, ratio=.25, npin=100) {
  # 2D b-spline basis
  ################## test ################## 
  # load("/Users/mac/Documents/GitHub/Redlining_Air_Pollution/output/simulation_v3/simu.data1.Rdata")
  # geometry <- dt$geometry
  # SPL <- dt$SPL
  # L1=5; L2=5; npin=100
  ################## test ################## 

  uq.city <- unique(SPL)
  
  splineMat <- list()
  
  if(ratio==Inf) {
    for(i in 1:length(uq.city)){
      dsub <- geometry[SPL == uq.city[i]]
      splineMat[[i]] <- diag(rep(1, length(dsub)))
    }
    return(splineMat)   
  }
  
  if(ratio==0) {
    for(i in 1:length(uq.city)){
      dsub <- geometry[SPL == uq.city[i]]
      splineMat[[i]] <- diag(rep(1, length(dsub)))*0
    }
    return(splineMat)   
  }

  ## "although coordinates are longitude/latitude, st_intersects assumes that they are planar"
  sf_use_s2(FALSE) 
  
  custL <- is.na(L1) | is.na(L2)
  
  ## if all city has same geometry
  i <- 1 # for(i in 1:length(uq.city)){
  dsub <- geometry[SPL == uq.city[i]]
  
  if(custL) {
    L1  <- floor(sqrt(length(dsub)))
    L2 <- ceiling(length(dsub) / L1)
  }

  ## Define a grid of points, create point geometry
  s1 <- seq(st_bbox(dsub)[1], st_bbox(dsub)[3], length.out = npin)
  s2 <- seq(st_bbox(dsub)[2], st_bbox(dsub)[4], length.out = npin) 
  s <- expand.grid(s1,s2)
  rm(s1, s2)
  
  ## Figure out which region each grid point belongs to
  s.sp <- st_as_sf(s, coords = c("Var1", "Var2"), crs = st_crs(dt)$epsg) # <---------------------------
  s2poly <- apply(st_intersects(s.sp, dsub, sparse = FALSE), 1, which)

  ## some points belong to multiple/no regions
  s2poly <- lapply(s2poly, FUN=function(x) ifelse(identical(x, integer(0)), NA, x))

  # lengths <- unlist(lapply(s2poly, length))
  # summary(lengths)
  # which(lengths==2)
  # which(lengths==0)
  # s2poly[[1101]]
  # s2poly[[6795]]
  # summary(unlist(s2poly))
  # length(unlist(s2poly)) == npin*npin 
  
  ## visualization
  # plot(s.sp,col="red",cex=.1)
  # plot(dsub$geometry, col = rgb(0, 0, 1, alpha=0.5), border = 'darkblue', add=TRUE)
  # 
  ## check whether some regions has no pin points 
  # print(paste(i, sum(!(dsub$id %in% s2poly)))) ## should be zero all the time 
  # 
  ## check coverage areas 
  # (sum(st_area(dsub)) / st_area(st_as_sfc(st_bbox(dsub))))
  # length(unlist(apply(st_intersects(s.sp, dsub$geometry, sparse = FALSE), 1, which))) / 100 / 100
  
  # Create a 2D spline basis (outer product of univariate bases) at the points
  B1 <- bs(s[,1], df=L1)
  B2 <- bs(s[,2], df=L2)
  Bs  <- NULL
  for(l1 in 1:L1){for(l2 in 1:L2){
    Bs <- cbind(Bs, B1[,l1]*B2[,l2])
  }}
  rm(B1, B2)

  ## check whether there are full-zero columns 
  # sum(colSums(Bs) == 0)
  
  ## visualization 
  # plot(s,cex=2*Bs[,64]+.01) # problem, many splines are wasted
  # plot(dsub$geometry, col = rgb(0, 0, 1, alpha=0.5), border = 'darkblue', add=TRUE)
  
  # Aggregate Bs by polygons
  B <- matrix(NA, length(dsub), ncol(Bs))
  for(k in 1:length(dsub)){
    these <- which(s2poly == k)
    if(length(these)==1){
      B[k, ] <- Bs[these, ] }
    if(length(these)>1){
      B[k, ] <- colMeans(Bs[these, ]) } 
  }

  ################### visualize ################### 
  # S <- matrix(0, length(dsub), 2)
  # for(k in 1:length(dsub)){
  #   these <- which(s2poly == k)
  #   S[k, ] <- colMeans(s[these, ])
  # }
  # pdf("./output/redlining_spline1.pdf", width = 6, height = 6) 
  # for (j in 1:ncol(B)) {
  #   plot(s,cex=2*Bs[,j]+.01, col="black",
  #        xlim=c(st_bbox(dsub)[1], st_bbox(dsub)[3]), ylim=c(st_bbox(dsub)[2], st_bbox(dsub)[4]))
  #   plot(dsub, col = rgb(0, 0, 1, alpha=0.2), border = 'darkblue', add=TRUE)
  #   plot(S,cex=20*B[,j]+.01, col="red",
  #        xlim=c(st_bbox(dsub)[1], st_bbox(dsub)[3]), ylim=c(st_bbox(dsub)[2], st_bbox(dsub)[4]))
  #   plot(dsub, col = rgb(0, 0, 1, alpha=0.2), border = 'darkblue', add=TRUE)
  # }
  # dev.off()
  ################### visualize ################### 

  ### TODO this factor is set for redlining data 
  n.spline <- max(floor(length(dsub) * ratio), 1)
  select <- order(colSums(B), decreasing=T)[1:n.spline]
  Bs <- Bs[, select]
  B <- B[, select]
  print(paste("there are ", ncol(B), " splines for ", length(dsub), " regions ||."))
  
  ## for areal data, this could be tricky and have no splines left
  select <- which(colSums(B) > 10^(-8))
  Bs <- Bs[, select]
  B <- B[, select]
  print(paste("there are ", ncol(B), " splines for ", length(dsub), " regions //."))

  ### TODO this factor is set for WAIC test 
  # n.spline <- L1*L2
  # select <- order(colSums(B), decreasing=T)[1:min(n.spline, ncol(B))]
  # Bs <- Bs[, select]
  # B <- B[, select]
  # print(paste("there are ", ncol(B), " splines for ", length(dsub), " regions ||."))
  
  ################### visualize ################### 
  # pdf("./output/redlining_spline2.pdf", width = 6, height = 6)
  # for (j in 1:ncol(B)) {
  #   plot(S,cex=20*B[,j]+.01, col="red",
  #        xlim=c(st_bbox(dsub)[1], st_bbox(dsub)[3]), ylim=c(st_bbox(dsub)[2], st_bbox(dsub)[4]))
  #   plot(dsub, col = rgb(0, 0, 1, alpha=0.2), border = 'darkblue', add=TRUE)
  # }
  # dev.off()
  ################### visualize ################### 
  
  ## scale in a way that in each spline, the maximum value is one 
  B <- sweep(B, 2, apply(B, 2, FUN=max), FUN="/")

  # splineMat <- c(splineMat, list(B))
  # } # end for(i in 1:length(uq.city))
  
  for(i in 1:length(uq.city)) {
    splineMat <- c(splineMat, list(B))
  }
  return(splineMat)
}


compute.spline <- function(m, b) {
  ## m is a list of matrix by city 
  ## b has length equal to the total number of columns in m
  result <- list()
  
  v.spline <- unlist(lapply(m, function(x) if(is.matrix(x)) ncol(x) else 1)) 
  
  for(i in 1:length(m)) {
    idx0 <- ifelse(i == 1, 1, sum(v.spline[1:(i-1)])+1)
    idx1 <- sum(v.spline[1:i])
    idx1 - idx0 + 1
    
    mat <- m[[i]]
    
    res <- mat %*% b[idx0:idx1]
    result[[i]] <- res
  }

  result <- matrix(unlist(result), byrow=F, ncol=1)
  return(result)
}

# m <- as.list(split.data.frame(matrix(runif(100000000), ncol=10), rep(1:100, 10)))
# v.spline <- unlist(lapply(m, function(x) if(is.matrix(x)) ncol(x) else 1))
# b <- runif(ncol(m[[1]])*length(m))
# b2 <- split(b, rep(1:100, v.spline)) # careful about this rep()
# 
# t0 <- Sys.time()
# aa <- compute.spline(m, b)
# print(Sys.time()-t0)
# 
# t0 <- Sys.time()
# bb <- unlist(mapply(FUN = function(x, y) x %*% y, m, b2, SIMPLIFY = T))
# print(Sys.time()-t0)
# cor(aa,bb)

redlining.spline <- function(ratio = .1) {
  
  load("data/nhgis0008_csv/only 1940/final_census1940_sp.Rdata")
  
  L1 <- L2 <- NA; npin <- 100
  
  uq.city <- unique(dt$SPL)
  
  splineMat <- list()
  
  if(ratio==Inf) {
    for(i in 1:length(uq.city)){
      dsub <- dt$geometry[dt$SPL == uq.city[i]]
      splineMat[[i]] <- diag(rep(1, length(dsub)))
    }
    saveRDS(splineMat, file=paste0("data/nhgis0008_csv/only 1940/RedliningSplineMat_RInf.rds"))
    return()   
  }
  
  if(ratio==0) {
    for(i in 1:length(uq.city)){
      dsub <- dt$geometry[dt$SPL == uq.city[i]]
      splineMat[[i]] <- diag(rep(1, length(dsub)))*0
    }
    saveRDS(splineMat, file=paste0("data/nhgis0008_csv/only 1940/RedliningSplineMat_R0.rds"))
    return()   
  }
  
  ## "although coordinates are longitude/latitude, st_intersects assumes that they are planar"
  sf_use_s2(FALSE) 
  
  custL <- is.na(L1) | is.na(L2)
  
  # temp <- dt %>% group_by(SPL) %>% summarise(n())
  # temp$`n()`

  for(i in 1:length(uq.city)){
    dsub <- dt[dt$SPL == uq.city[i], ] 
    
    if(custL) {
      L1  <- floor(sqrt(nrow(dsub)))
      L2 <- ceiling(nrow(dsub) / L1)
    }
    L1 <- max(L1, 4)
    L2 <- max(L2, 4)
    
    ## Define a grid of points, create point geometry
    s1 <- seq(st_bbox(dsub)[1], st_bbox(dsub)[3], length.out = npin)
    s2 <- seq(st_bbox(dsub)[2], st_bbox(dsub)[4], length.out = npin) 
    s <- expand.grid(s1,s2)
    rm(s1, s2)
    
    ## Figure out which region each grid point belongs to
    s.sp <- st_as_sf(s, coords = c("Var1", "Var2"), crs = st_crs(dt)$epsg)
    temp <- apply(st_intersects(s.sp, dsub, sparse = FALSE), 1, which)
    s2poly <- unlist(lapply(temp, replace_empty_with_NA))      # <-------------------- TODO change this????? 
    rm(temp)
    
    ## TODO tricky that some point belongs in both, need to improve this
    if(! (length(s2poly) == (npin*npin))) {
      print(paste0("s2poly truncated...", dsub$location[1]))
      s2poly <- s2poly[1:(npin*npin)]
    }

    # Create a 2D spline basis (outer product of univariate bases) at the points
    
    # if(L1 < 4) {
    #   temp <- (s[,1] - min(s[,1])) / (max(s[,1]) - min(s[,1]))
    #   if(L1 == 1) {
    #     B1 <- matrix(temp, ncol=1)
    #   } else if (L1 == 2) {
    #     B1 <- cbind(temp, temp^2)
    #   } else if (L1 == 3) {
    #     B1 <- cbind(temp, temp^2, temp^3)
    #   }
    # } else {
    #   B1 <- bs(s[,1], df=L1)
    # }
    # 
    # if(L2 < 4) {
    #   temp <- (s[,2] - min(s[,2])) / (max(s[,2]) - min(s[,2]))
    #   if(L2 == 1) {
    #     B2 <- matrix(temp, ncol=1)
    #   } else if (L2 == 2) {
    #     B2 <- cbind(temp, temp^2)
    #   } else if (L2 == 3) {
    #     B2 <- cbind(temp, temp^2, temp^3)
    #   }
    # } else {
    #   B2 <- bs(s[,2], df=L2)
    # }
    
    B1 <- bs(s[,1], df=L1)
    B2 <- bs(s[,2], df=L2)
    
    Bs  <- NULL
    for(l1 in 1:L1){for(l2 in 1:L2){
      Bs <- cbind(Bs, B1[,l1]*B2[,l2])
    }}
    rm(B1, B2)
    
    # Aggregate Bs by polygons
    B <- matrix(NA, nrow(dsub), ncol(Bs))
    for(k in 1:nrow(dsub)){
      these <- which(s2poly == k)

      # k=9
      if(is_integer_zero(these)) {
        B[k, ] <- 0 # no points fall into this region at all... 
      } else if(length(these)==1){
        B[k, ] <- Bs[these, ] 
      } else if(!is.matrix(Bs)) {
        B[k, ] <- Bs[these, ]
      } else {
        B[k, ] <- colMeans(Bs[these, ])
      } 
    }
    
    S <- matrix(0, nrow(dsub), 2)
    for(k in 1:nrow(dsub)){
      these <- which(s2poly == k)
      S[k, ] <- colMeans(s[these, ])
    }
    
    ################### visualize ################### 
    # pdf(paste0("data/nhgis0008_csv/only 1940/", dsub$location[1], "_temp.pdf"), width = 4, height = 4)
    # for (j in 1:ncol(B)) {
    #   plot(s,cex=2*Bs[,j]+.01, col="black",
    #        xlim=c(st_bbox(dsub)[1], st_bbox(dsub)[3]), ylim=c(st_bbox(dsub)[2], st_bbox(dsub)[4]))
    #   plot(dsub$geometry, col = rgb(0, 0, 1, alpha=0.2), border = 'darkblue', add=TRUE)
    #   plot(dsub$geometry, col = scales::col_bin("Reds", domain = B[,j])(B[,j]), main = "")
    # }
    # dev.off()
    ################### visualize ################### 

    n.spline <- max(floor(nrow(dsub) * ratio), 1)
    if(is.matrix(B)) {
      select <- order(colSums(B), decreasing=T)[1:n.spline]
      Bs <- Bs[, select]
      B <- B[, select]
      # print(paste0("there are ", ncol(B), " splines for ", nrow(dsub), " region ", dsub$location[1], i))
    } else {
      B <- matrix(B)
      # print(paste0("there are 0 or 1 splines for ", nrow(dsub), " region ", dsub$location[1], i))
    }

    ## for areal data, this could be tricky and have no splines left
    if(is.matrix(B)) {
      select <- which(colSums(B) > 10^(-8))
      Bs <- Bs[, select]
      B <- B[, select]
      # print(paste0(ncol(B), " splines for ", nrow(dsub), " region ", dsub$location[1], i))
      
      ################### visualize ################### 
      pdf(paste0("data/nhgis0008_csv/only 1940/", dsub$location[1], "_", ratio, ".pdf"), width = 6, height = 6)
      for (j in 1:ncol(B)) {
        plot(dsub$geometry, col = scales::col_bin("Reds", domain = B[,j])(B[,j]), main = "")
      }
      dev.off()
      ################### visualize ################### 
      
    } else {
      B <- matrix(B)
      # print(paste0("0 or 1 splines for ", nrow(dsub), " region ", dsub$location[1], i))

      ################### visualize ################### 
      pdf(paste0("data/nhgis0008_csv/only 1940/", dsub$location[1], "_", ratio, ".pdf"), width = 6, height = 6)
      plot(dsub$geometry, col = scales::col_bin("Reds", domain = B)(B), main = "")
      dev.off()
      ################### visualize ################### 
      
    }
    
    ## scale in a way that in each spline, the maximum value is one 
    if(is.matrix(B)) {
      B <- sweep(B, 2, apply(B, 2, FUN=max), FUN="/")
    } else {
      B <- matrix(B / max(B))
    }
    
    # print(kappa(B)) # condition number higher than 1000 is bad 
    
    ## B 1 and 0 fails 
    splineMat <- c(splineMat, list(B))
  }
  
  v.spline <- unlist(lapply(splineMat, function(x) if(is.matrix(x)) ncol(x) else 1))
  n.spline <- sum(v.spline)
  
  print(paste0("total spline count is ", n.spline))

  saveRDS(splineMat, file=paste0("data/nhgis0008_csv/only 1940/RedliningSplineMat_R", ratio*100, ".rds"))
  saveRDS(splineMat, file=paste0("data/nhgis0008_csv/only 1940/RedliningSplineMat", n.spline, ".rds"))
  return()
}

# redlining.spline(Inf)
# redlining.spline(.10)
# redlining.spline(.20)
# redlining.spline(.30)
# redlining.spline(.40)
# redlining.spline(.50)
# redlining.spline(.60)
# redlining.spline(.70)
# redlining.spline(.80)
# redlining.spline(.90)
# redlining.spline(1)
# redlining.spline(1.1)
# redlining.spline(1.2)
# redlining.spline(1.3)
# redlining.spline(1.4)
# redlining.spline(1.5)
# redlining.spline(1.6)
# redlining.spline(1.7)
# redlining.spline(1.8)
# redlining.spline(1.9)
# redlining.spline(2)

