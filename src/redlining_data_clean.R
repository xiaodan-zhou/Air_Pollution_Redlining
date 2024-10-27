##################################################################################
## redlnining in 1937 
## map view and download https://dsl.richmond.edu/panorama/redlining/ma
##################################################################################
## census tract 1940 mapped into 2010 census
## map view https://umn.maps.arcgis.com/apps/webappviewer/index.html?id=ef554b42643141829f9e8c4b8001f93a 
## donwload https://data2.nhgis.org/main
  ## GEOGRAPHIC LEVELS: TRACT
  ## YERS: 1940
  ## BASIS: 2008 TIGER/Line +
##################################################################################
## air pollution in 2010 mapped into 2010 census using area as weight
##################################################################################

library(readr)
library(dplyr)
library(ggplot2)
library(sf)
library(spdep)

##################################################################################
######################## clean up census 1940 #################################### 
##################################################################################

####################### get the census map in 1940 #######################
tract40 <- read_sf("data/nhgis0008_shape/nhgis0008_shapefile_tl2000_us_tract_1940/US_tract_1940.shp")

## some area has multipolygon in this data
tract40 <- tract40 %>% 
  st_transform(crs=st_crs("NAD83")) %>% 
  st_make_valid() %>% 
  mutate(tract_total_area=as.numeric(st_area(.))) %>% # before st_cast()
  filter(!st_is_empty(geometry))
dim(tract40) # 7563

tract40 <- tract40 %>% 
  st_cast("POLYGON") %>% 
  mutate(tract_partial_area=as.numeric(st_area(.))) # after st_cast()
dim(tract40) # 7921

## how many has multipolygon? looks like not much 
# summary(tract40$tract_partial_area / tract40$tract_total_area)

## get the census1940 table 
df <- read.csv("data/nhgis0008_csv/nhgis0008_ds76_1940_tract.csv")
df <- df %>% rename("total.pop40"="BUB001", "black.pop40"="BVG001", 
                    "male.employed"="BUD001", "female.employed"="BUD002", 
                    "male.labor"="BUW001", "male.not.labor"="BUW002", 
                    "female.labor"="BUW003", "female.not.labor"="BUW004", 
                    "median.rent"="BVJ001")  

# dim(df[df$COUNTY =="Richmond City",])

## jot the census1940 table with the census1940 map, one-for-one
tract40 <- inner_join(x=tract40, y=df, by="GISJOIN")
dim(tract40) # 7708

## *************** split the population in those multipolygons using area weight *************** 
pop.col <- c("black.pop40", "total.pop40", "male.employed", "female.employed", 
             "male.labor", "male.not.labor", "female.labor", "female.not.labor")
tract40 <- tract40 %>% mutate(across(all_of(pop.col), ~ .x *tract_partial_area/tract_total_area))
## *** tract_partial_area is the true area of each polygon in this data *** 

## check 
# t <- data.frame(tract40) %>% group_by(GISJOIN) %>% summarise(sum(tract_partial_area)/tract_total_area)
# summary(t[,2])

save(tract40, file="data/cleaned/tract40.Rdata")

####################### stats of census 1940 #######################
sum(tract40$black.pop40) # 3337510 the total of black population covered in census 
sum(tract40$total.pop40) # 37541977 the total of population covered in census 
sum(tract40$black.pop40)/sum(tract40$total.pop40)*100 # 8.8% black 
sum(tract40$tract_total_area) # 45156395948 (11)
####################### stats of census 1940 #######################

## check 
# a <- tract40 %>% filter(tract40$GISJOIN == "G55007900019") 
# a %>% ggplot() + geom_sf(col="red")
# round(a$tract_partial_area,2)


##################################################################################
############# join census 1940 into redlining maps ############################### 
##################################################################################

####################### get the holc map #######################
holc30 <- read_sf("/Users/mac/Documents/GitHub/Redlining-and-Air-Pollution-SDOH/mapping_ineq_shapefile/fullshpfile/shapefile/holc_ad_data.shp")
holc30 <- holc30 %>% 
  st_transform(crs=st_crs("NAD83")) %>% 
  st_make_valid() %>% 
  mutate(holc_total_area=as.numeric(st_area(.))) %>%   # before st_cast() 
  filter(!st_is_empty(geometry))
dim(holc30) # 8875

holc30 <- holc30 %>% st_cast("POLYGON")
dim(holc30) # 8875
## this data has no multipolygon, great

# holc30 %>% filter(name == "Norwood, Druid Hill, Walnut Hill")

####################### stats of holc 1930 #######################
sum(holc30$holc_total_area) # 12457087773 (11)
####################### stats of holc 1930 #######################

holc.intesect <- st_intersection(holc30, tract40)
#############################################
# save(holc.intesect, file="data/holc.intesect.Rdata")
# holc.intesect, not directly useful!!! 
# holc.intesect, not directly useful!!! 
# holc.intesect, not directly useful!!! 
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
rm(list = ls())
load("data/holc.intesect.Rdata") 
#############################################
rent.cat <- c(paste0("BVI00", 1:9), paste0("BVI0", 10:13)) 

sum(holc.intesect$black.pop40) # 6483550 there is duplicate, need to recalculate fields 

## bad, some redlinied areas are rarely census-ed in 1940, may need to exclude this 
ggplot(holc.intesect[holc.intesect$city == "Richmond",]) + geom_sf()

holc.intesect.area <- holc.intesect %>%
  dplyr::select(state, city, holc_id, holc_grade, neighborho, 
                 holc_total_area, tract_total_area, tract_partial_area, 
                 total.pop40, black.pop40, male.employed, female.employed, 
                 male.labor, male.not.labor, female.labor, female.not.labor, median.rent, 
                 geometry, rent.cat) %>%
  mutate(int_area = st_area(.) %>% as.numeric())
## *** int_area is the true area of each polygon in this data ***

## check 
# summary(holc.intesect.area$int_area/holc.intesect.area$holc_total_area)
# summary(holc.intesect.area$int_area/holc.intesect.area$tract_partial_area)
# t <- data.frame(holc.intesect.area) %>% group_by(neighborho) %>% summarise(sum(int_area)/holc_total_area[1])
# t[t[,2] > 1.01,]
# holc.intesect.area[holc.intesect.area$neighborho == 8678,]
# holc30[holc30$neighborho == 8678,]

holc.intesect.area <- holc.intesect.area %>% filter(neighborho != 8678)
dim(holc.intesect.area)
sum(holc.intesect.area$total.pop40)
sum(holc.intesect.area$black.pop40) # 6468381 

pop.col <- c("total.pop40", "black.pop40", "male.employed", "female.employed", 
             "male.labor", "male.not.labor", "female.labor", "female.not.labor",
             rent.cat)

holc.intesect.area <- holc.intesect.area %>% 
  mutate(pct_census_holc_area = int_area/holc_total_area) %>%                 # this is the split of holc region 
  mutate(across(all_of(pop.col), ~ .x * int_area/tract_partial_area))  # this is the split of tract region
## *** int_area is the true area of each polygon in this data ***
## *** pct_census_holc_area show this polygon occupies what percent in that holc *** 

########### summary of intersection of hocl and census 1940 ##############
sum(holc.intesect.area$total.pop40) # 25742598
sum(holc.intesect.area$black.pop40) # 2122270 < 3337646
########### summary of intersection of hocl and census 1940 ##############

dt <- holc.intesect.area %>%
  tidyr::drop_na() %>%
  group_by(neighborho) %>% 
  summarize(n=n(), state=state[1], city=city[1], holc_id=holc_id[1], holc_grade=holc_grade[1],
            geometry = st_union(geometry), 
            pct_census_holc_area = sum(pct_census_holc_area),
            # 
            across(all_of(pop.col), ~ sum(.x)), # group and sum 
            # 
            median.rent=weighted.mean(median.rent, int_area), 
            # 
            pct.black40 = sum(black.pop40)/sum(total.pop40)*100,
            pct.male.employed = sum(male.employed)/sum(total.pop40)*100,
            pct.female.employed = sum(female.employed)/sum(total.pop40)*100,
            pct.male.labor = sum(male.labor)/sum(total.pop40)*100,
            pct.male.not.labor = sum(male.not.labor)/sum(total.pop40)*100,
            pct.female.labor = sum(female.labor)/sum(total.pop40)*100,
            pct.female.not.labor = sum(female.not.labor)/sum(total.pop40)*100)
            
## check 
# summary(dt)
# hist(dt$n)
# hist(dt$sd.median.rent)
# hist(dt$median.rent)
# hist(holc.intesect.area$median.rent)
# hist(dt$pct.black40)
# holc.intesect.area$median.rent[holc.intesect.area$neighborho==holc.intesect.area$neighborho[1]]
# ggplot(dt) + geom_histogram(aes(pct.black40, fill=holc_grade),alpha=.5)
# ggplot(dt) + geom_histogram(aes(median.rent, fill=holc_grade),alpha=.5)
# ggplot(holc.intesect.area) + geom_histogram(aes(median.rent, fill=holc_grade),alpha=.5)

ggplot(dt[dt$city=="Boston",]) + ggtitle("") +
  geom_sf() + theme_void() +
  theme(legend.position = "None")

dt$holc <- dt$holc_grade
dt$holc[dt$holc=="A"] <- 1
dt$holc[dt$holc=="B"] <- 2
dt$holc[dt$holc=="C"] <- 3
dt$holc[dt$holc=="D"] <- 4
dt$holc <- as.numeric(as.character(dt$holc)) 
dt <- dt[!is.na(dt$holc),]

dt$city.unique <- paste(dt$state, dt$city)

dt[duplicated(dt[c("state","city","neighborho")]),]

#############################################
# write.csv(dt, file="data/HOLC_census40.csv", row.names=FALSE)
#############################################

class(dt)
save(dt, file="data/cleaned/census40_sp.RData")

dt[, !(names(dt) %in% "geometry")] %>% 
  as.data.frame() %>% 
  saveRDS(file="data/cleaned/census40.rds")
#############################################
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
rm(list = ls())

library(readr)
library(dplyr)
library(ggplot2)
library(sf)
library(spdep)
##################################################################################
########## join air pollution 2010 into [census40_sp.RData] ######################
##################################################################################
## to joint 
load("data/cleaned/census40_sp.RData") # dt

## load air pollution data 
load("/Users/mac/Documents/GitHub/Redlining-and-Air-Pollution-SDOH/data/merged_data.Rdata") # dat
# sum(duplicated(dat %>% select(-grade))) # disaster... 
# sum(duplicated(dat$geometry)) # disaster... 
dat <- dat[, names(dat) != "grade"]
dat <- dat[!duplicated(dat),]
sum(duplicated(dat))

dat <- sf::st_cast(dat, "POLYGON") # none

names(dat)[4] <- c("city2010")

dat$state2010 <- unlist(lapply(dat$city2010, FUN=function(x) substr(x, nchar(x) - 1, nchar(x))))
dat$city2010 <- unlist(lapply(dat$city, FUN=function(x) paste(strsplit(x, ",")[[1]][1:(length(strsplit(x, ",")[[1]])-1)], collapse=",")))

names(dat)

# use census map to reduce the air pollution map size
tract40 <- read_sf("data/nhgis0008_shape/nhgis0008_shapefile_tl2000_us_tract_1940/US_tract_1940.shp")
tract40$dig5 <- paste0(substr(tract40$NHGISST, 1, 2), substr(tract40$NHGISCTY, 1, 3))
sel.dig5 <- unique(tract40$dig5); rm(tract40)
# 
dat$dig5 <- substr(dat$GEOID, 1, 5)
dat <- dat %>% filter(dat$dig5 %in% sel.dig5)
names(dat)

######### st_intersection(holc30, dat) ######### 
ss <- NULL
for (city in unique(dt$city)) {
  for (state in unique(dt$state[dt$city == city])) {
    matches <- dat[dat$state2010 == state & dat$city2010 == city, ]
    ss <- rbind(ss, c(city, state, sum(dt$city==city), dim(matches)[1]))
    rm(matches)
  }
}
ss <- as.data.frame(ss)
ss$V4 <- as.numeric(as.character(ss$V4))
ss$V3 <- as.numeric(as.character(ss$V3))
sum(ss$V3*ss$V4) # 52700008

# plot(ss$V3, ss$V4)
# summary(ss$V4)
# a <- dat[!(dat$city2010 %in% unique(dt$city)), ]
# unique(a$city2010)

# ctt = "San Francisco"
# stt = "CA"
# which(unique(holc30$city)=="Indianapolis")

progress <- 0
for (ctt in unique(dt$city)) {
  for (stt in unique(dt$state[dt$city == ctt])) {
    
    match1 <- dt %>% filter(state == stt, city == ctt) # population spread by this 
    match2 <- dat %>% filter(state2010 == stt, city2010 == ctt)
    
    progress <- progress + dim(match1)[1]*dim(match2)[1]/sum(ss$V3*ss$V4)*100
    print(paste(ctt, stt, dim(match1)[1], dim(match2)[1], "progress", progress, "%..."))
    print(timestamp())
    
    intersection <- sf::st_intersection(match1, match2)
    save(intersection, file=paste0("data/cleaned/holc_pm_", ctt, "_", stt, ".Rdata"))
    rm(match1, match2, intersection)
  }
}
## this dt is not directly useful! we need to spread the population into the intersection!! 
## this dt is not directly useful! we need to spread the population into the intersection!! 
## this dt is not directly useful! we need to spread the population into the intersection!! 
#
#
#
#
#
#
#
#
#
#
#
#
#
######################################################################### 
rm(list = ls())

library(readr)
library(dplyr)
library(ggplot2)
library(sf)
library(spdep)

load("data/cleaned/census40_sp.RData") # dt
dt <- st_transform(dt, crs=st_crs("NAD83"))
t1 <- dt[,c("neighborho", "pct_census_holc_area")] %>% as.data.frame() %>% dplyr::select(-geometry)

holc30 <- read_sf("/Users/mac/Documents/GitHub/Redlining-and-Air-Pollution-SDOH/mapping_ineq_shapefile/fullshpfile/shapefile/holc_ad_data.shp")
t2 <- holc30 %>% st_transform(crs=st_crs("NAD83")) %>% st_make_valid() %>% mutate(holc_area=as.numeric(st_area(.))) %>% 
  dplyr::select(-geometry)  %>% as.data.frame() %>% dplyr::select(neighborho, holc_area)  %>% as.data.frame()
rm(holc30)

pollution <- NULL 
for (city in unique(dt$city)) {
  for (state in unique(dt$state[dt$city == city])) {
    print(paste(city, state))
    
    load(paste0("data/cleaned/holc_pm_", city, "_", state, ".Rdata")) # intersection 

    temp <- intersection %>% 
      left_join(t1, by="neighborho") %>% 
      left_join(t2, by="neighborho") %>% 
      st_transform(crs=st_crs("NAD83")) %>%
      st_make_valid() %>%
      filter(st_geometry_type(geometry) != "GEOMETRYCOLLECTION") %>%
      mutate(pm_area=as.numeric(st_area(.))) %>% # before st_cast()
      mutate(pct_pm_holc_area = pm_area/holc_area) %>% 
      filter(pm_area > 0) %>%
      filter(!st_is_empty(geometry)) %>% 
      group_by(neighborho) %>%
      summarise(
        city=city[1], state=state[1],
        pm25=weighted.mean(pm25, pm_area),
        no2=weighted.mean(no2, pm_area), 
        pm_area=sum(pm_area),
        pct_pm_holc_area=sum(pct_pm_holc_area),
        geometry=st_union(geometry)
        )
    
    pollution <- rbind(pollution, temp)
  }
}

pollution <- pollution %>% as.data.frame() # %>% select(-geometry) %>% as.data.frame()

########################################################################
save(pollution, file="data/cleaned/pollution.RData")
########################################################################

load("data/cleaned/census40_sp.RData") # dt
dt <- st_transform(dt, crs=st_crs("NAD83"))
dt <- left_join(dt, pollution %>% dplyr::select(neighborho, pm25, no2, pm_area, pct_pm_holc_area) %>% as.data.frame(), by="neighborho")
dt <- dt %>% mutate(census_area=as.numeric(st_area(geometry)))

plot(dt$census_area, dt$pm_area)

########################################################################
save(dt, file="data/cleaned/accurate_census1940_sp.Rdata")
########################################################################
#
#
#
#
#
#
#
#
#
#
#
##################################################################################
############# join census 1940 and air pollution 2010 ############################
##################################################################################
rm(list = ls())

library(readr)
library(dplyr)
library(ggplot2)
library(sf)
library(spdep)

## get back the original shape 
## get the holc map 
holc30 <- read_sf("/Users/mac/Documents/GitHub/Redlining-and-Air-Pollution-SDOH/mapping_ineq_shapefile/fullshpfile/shapefile/holc_ad_data.shp")
holc30 <- st_transform(holc30, crs=st_crs("NAD83"))
holc30 <- holc30 %>% st_make_valid() %>% st_cast("POLYGON")

load("data/cleaned/accurate_census1940_sp.Rdata") # dt 

## this is pretty bad... 
# ggplot(dt[dt$city=="Richmond",]) + geom_sf()
# ggplot(holc30[holc30$city=="Richmond",]) + geom_sf()
# 
# ggplot(dt[dt$city=="Los Angeles",]) + geom_sf()
# ggplot(holc30[holc30$city=="Los Angeles",]) + geom_sf()
# 
# holc.city <- unique(holc30$city)
# holc.city[!(holc.city %in% unique(dt$city))]

##############################
dt <- dt %>% as.data.frame() %>% dplyr::select(-geometry) %>% as.data.frame()

dt <- dt %>% arrange(city, neighborho)
holc30 <- holc30 %>% arrange(city, neighborho)

dt <- inner_join(holc30[,"neighborho"], dt, by="neighborho")

class(dt)

sum(dt$black.pop40) # 2062599
sum(dt$total.pop40) # 25535678
sum(dt$black.pop40)/sum(dt$total.pop40)*100 # 8.1%

## check data meaning 
# ints <- dt
# summary(ints$male.employed/ints$male.labor) # <= 1 
# summary(ints$male.employed/(ints$male.labor+ints$male.not.labor)) # <= 1 
# summary((ints$male.employed+ints$female.employed)/ints$total.pop40) # <= 1 
# summary((ints$male.labor+ints$male.not.labor+ints$female.labor+ints$female.not.labor)/ints$total.pop40) # <= 1

dt <- dt %>% 
  mutate(
    pct.black.pop40=black.pop40/total.pop40, 
    rate.male.employed=male.employed/male.labor, 
    rate.female.employed=female.employed/female.labor, 
    rate.employed=(male.employed+female.employed)/(male.labor+female.labor),
    rate.male.labor=male.labor/(male.labor+male.not.labor), 
    rate.female.labor=female.labor/(female.labor+female.not.labor)) 

dt <- dt %>% mutate(holc_area=as.numeric(st_area(geometry)))

sum(duplicated(dt$neighborho))

# a <- dt$census_area/dt$holc_area
# b <- dt$pm_area/dt$holc_area 
# summary(a)
# summary(b)
# hist(b,breaks=100)
# dt[is.na(b), ]
# dt[b>1, ]
# dt[b<.90, ]

dt$location <- paste(dt$city, dt$state)
dt$SPL <- as.numeric(factor(dt$location, levels = unique(dt$location)))
dt <- dt[order(dt$SPL), ]

##### calculate the median rent using the categorical columns 
rent.col <- c(paste0("BVI00", 1:9), paste0("BVI0", 10:13)) 
rent.val <- c(2.5, 5.5, 8, 12, 17, 22, 27, 34.5, 44.5, 54.5, 67, 87, 107)

dt <- dt %>% rowwise() %>% mutate(median.recal = median(rep(rent.val, times = c_across(all_of(rent.col))))) %>% ungroup()

short <- dt[,rent.col]
short <- short %>% as.data.frame() %>% dplyr::select(-geometry) %>% as.data.frame()
# dt$wt.mean.rent <- rowSums(t(t(dt[,rent.col]) * rent.val)) / rowSums(dt[,rent.col])
dt$wt.mean.rent <- rowSums(t(t(short) * rent.val)) / rowSums(short)

hist(dt$median.rent,breaks=100)
hist(dt$wt.mean.rent,breaks=100)
hist(dt$median.recal,breaks=100)
plot(dt$median.rent, dt$wt.mean.rent)
plot(dt$median.rent, dt$median.recal)

cor.test(dt$median.rent, dt$wt.mean.rent,na.rm=T)
cor.test(dt$median.rent, dt$median.recal,na.rm=T)
cor.test(dt$median.recal, dt$wt.mean.rent,na.rm=T)

# hist(rowSums(dt[,rent.col]),breaks=100)
# plot(rowSums(dt[,rent.col]), dt$total.pop40)

dt <- dt[!is.na(dt$pm25), ]
##################################################################################
save(dt, file="data/cleaned/final_census1940_sp.Rdata")
##################################################################################

library(RNOmni)
dt$pct.black.pop40.rank <- RNOmni::RankNorm(dt$pct.black.pop40) ## latent normal assumption for TOBIT model

library(MASS)
transformed_data <- boxcox(dt$rate.employed ~ 1, lambda = seq(0, 10, by = 0.01), plotit=F)
best_lambda <- transformed_data$x[which.max(transformed_data$y)]; print(best_lambda)
dt$bc.rate.employed <- ((dt$rate.employed^best_lambda - 1) / best_lambda) # box-cox transformed

dt %>% 
  as.data.frame() %>% 
  dplyr::select(-geometry) %>% 
  as.data.frame() %>% 
  saveRDS(file="data/cleaned/final_census1940.rds")
##################################################################################

sum(dt$black.pop40) # 2062599
sum(dt$total.pop40) # 25535678
sum(dt$black.pop40)/sum(dt$total.pop40)*100 # 8.1%

length(unique(dt$city.unique)) # 70
unique(dt$city.unique)
table(dt$SPL)




# others: 
# These boundaries are derived from the U.S. Census Bureau's 2008 TIGER/Line files, or they have been 
# conflated to align with 2008 TIGER/Line features as consistently as possible.
# The original NHGIS boundaries were derived primarily from the 2000 TIGER/Line files. 
# For 1990 and 2000 areas, NHGIS modified the 2000 TIGER/Line definitions by erasing coastal water areas. Because the 2000 
# TIGER/Line files contain no identifiers for census areas from 1980 and earlier, NHGIS researchers obtained 
# boundary definitions for those years by consulting other sources, including 1992 TIGER/Line data for 1980 
# census tracts; maps from printed census reports for 1910-1980 census tracts and other small areas; and the 
# Map Guide to the U.S. Federal Censuses, 1790-1920, by William Thorndale and William Dollarhide (Genealogical
# Publishing Co., Baltimore, MD, 1987), for counties and states back to 1790. Where the historical boundaries 
# follow 2000 TIGER/Line features, the original NHGIS boundary files re-use those TIGER/Line features. 
# Elsewhere, NHGIS researchers digitized new boundaries.
# Because the Census Bureau made major accuracy improvements to TIGER/Line features between the 2000 and 
# 2008 TIGER/Line releases, the original NHGIS shapefiles, being based on 2000 TIGER/Line features, are not 
# comparable with newer TIGER/Line data. We therefore generated new 2008-based boundary files for counties 
# and census tracts by systematically conflating the 2000-based NHGIS boundaries to fit with 2008 TIGER/Line 
# features.
# The Census Bureau subsequently made additional improvements to TIGER/Line features, so 
# the 2008-based files are not consistently comparable with 2010 and later TIGER/Line files. 
# In general, most 2008-based boundaries 
# align better than 2000-based boundaries with more recent TIGER/Line files, but the 2008-based boundaries 
# also include occasional gross inaccuracies. 
# For users who have no need to compare historical boundaries with boundaries from 2010 or later, we recommend 
# using the original 2000-based NHGIS boundary files. 
# For users who do wish to compare or overlay historical 
# boundaries with boundaries from 2010 or later, we recommend downloading and examining both the 2000- and 
# 2008-based versions of historical boundaries in order to determine which is more suitable for your study 
# area and analysis.
##################################################################################
