library(sf)
library(st)
library(tidyr)
library(dplyr)
library(xtable)

state_names <- c("Alabama", "Alaska", "Arizona", "Arkansas", "California", "Colorado", "Connecticut", 
                 "Delaware", "Florida", "Georgia", "Hawaii", "Idaho", "Illinois", "Indiana", 
                 "Iowa", "Kansas", "Kentucky", "Louisiana", "Maine", "Maryland", "Massachusetts", 
                 "Michigan", "Minnesota", "Mississippi", "Missouri", "Montana", "Nebraska", "Nevada", 
                 "New Hampshire", "New Jersey", "New Mexico", "New York", "North Carolina", 
                 "North Dakota", "Ohio", "Oklahoma", "Oregon", "Pennsylvania", "Rhode Island", 
                 "South Carolina", "South Dakota", "Tennessee", "Texas", "Utah", "Vermont", 
                 "Virginia", "Washington", "West Virginia", "Wisconsin", "Wyoming", "District Of Columbia")

state_abbr <- c("AL", "AK", "AZ", "AR", "CA", "CO", "CT", 
                "DE", "FL", "GA", "HI", "ID", "IL", "IN", 
                "IA", "KS", "KY", "LA", "ME", "MD", "MA", 
                "MI", "MN", "MS", "MO", "MT", "NE", "NV", 
                "NH", "NJ", "NM", "NY", "NC", 
                "ND", "OH", "OK", "OR", "PA", "RI", 
                "SC", "SD", "TN", "TX", "UT", "VT", 
                "VA", "WA", "WV", "WI", "WY", "DC")

us_states <- data.frame(state_full = state_names, state_abbr = state_abbr)


load("/Users/mac/Documents/GitHub/Redlining-and-Air-Pollution-SDOH/data/merged_data.Rdata") # dat, pm
dat$LON <- as.numeric(dat$LON); dat$LAT <- as.numeric(dat$LAT)

load("data/cleaned/tract40.Rdata") # tract40, census
holc30 <- read_sf("/Users/mac/Documents/GitHub/Redlining-and-Air-Pollution-SDOH/mapping_ineq_shapefile/fullshpfile/shapefile/holc_ad_data.shp") # redlining 
holc30$city.unique <- paste(holc30$state, holc30$city)

load("data/cleaned/final_census1940.RData") # 4089 sp
# dt <- readRDS("data/cleaned/final_census1940.rds") # our final data non-sp 4079 
dt <- st_as_sf(dt)


final.city <- unique(dt$city.unique)
holc30.city <- unique(holc30$city.unique)

nrow(holc30[holc30$city.unique %in% final.city,]) # 5264 in those 70 cities, finally there are 4089

### the relining data 
length(unique(holc30$city.unique)) # 202
nrow(holc30) # 8878
length(unique(holc30$state)) # 38

### the census 1940 data # not sure how many cities 
length(unique(tract40$STATE)) # 29

## redlining state that are in census state 
holc30$state_abbr <- holc30$state
holc30 <- merge(holc30, us_states, by.x="state_abbr", by.y="state_abbr", all.x=T)

unique(holc30$state_full)
unique(tract40$STATE)

sum(holc30$state_full %in% unique(tract40$STATE)) # 8219 of redlining are in census-covered state
sum(tract40$total.pop40) # 37541977
sum(tract40$total.pop40[tract40$STATE %in% unique(holc30$state_full)]) # 36878886 are in redlining-covered state 

## the final data 
length(unique(dt$city.unique)) # 70 
nrow(dt) # 4089

sum(holc30$city.unique %in% unique(dt$city.unique)) # 5264 areas area in (partially census covered) 70 cities

# ## get centrid point of holc30 polygon 
# holc30 <- holc30 %>% st_transform(crs=st_crs("NAD83")) %>% st_make_valid()
# tract40 <- tract40 %>% st_transform(crs=st_crs("NAD83")) %>% st_make_valid()
# 
# ## check if each of holc30$centroid full in any of tract40$geometry
# temp30 <- holc30[!duplicated(holc30$city.unique),]
# temp30$centroid <- st_centroid(temp30$geometry)
# dim(temp30) # 202 
# 
# temp30$in_tract40 <- NA
# for (i in 1:nrow(temp30)) {
#   print(i)
#   temp30$in_tract40[i] <- sum(st_within(temp30$centroid[i], tract40$geometry, sparse = FALSE)) > 0
# } 
# sum(temp30$in_tract40)

##################################################################################
## some holc areas are not fully observed in census tract
hist(dt$pct_census_holc_area, breaks=100)
mean(dt$pct_census_holc_area>=.98) # 80.3% of HOLC has less than 98% areas covered 
mean(dt$pct_census_holc_area>=.9) # 84.5% of HOLC has less than 90% areas covered 
mean(dt$pct_census_holc_area>=.8) # 86.2% of HOLC has less than 80% areas covered 
## some holc areas are not fully observed in air pollution ????????????????????????????????
hist(dt$pct_pm_holc_area, breaks=100)
mean(dt$pct_pm_holc_area>=.98, na.rm=T) # 23.1% has at least 98% area monitored by air pollution
mean(dt$pct_pm_holc_area>=.9, na.rm=T) # 64.4% has at least 90% area monitored by air pollution
mean(dt$pct_pm_holc_area>=.8, na.rm=T) # 82.2% has at least 90% area monitored by air pollution


##################################################################################
## check one city 
a <- dt[dt$city == "Richmond",]
a <- dt[dt$city == "Los Angeles",]
ggplot(a$geometry) + geom_sf() + theme_bw()

city <- "Richmond"; st=51; county=760
city <- "Los Angeles"; st=06; county=037  
edge <- .1

cord <- st_coordinates(st_centroid(dt$geometry[dt$city==city])[1])

# holc
holc30[holc30$city==city,] %>% ggplot() + geom_sf() + theme_bw()
# census 
tract40 %>% filter(STATEA == st, COUNTYA == county) %>% ggplot() + geom_sf() + theme_bw()
# pm
dat %>% filter(LON > cord[1]-edge & LON < cord[1]+edge & LAT > cord[2]-edge & LAT < cord[2]+edge) %>% 
  ggplot() + geom_sf() + theme_bw()
## final 
dt[dt$city==city,] %>% ggplot() + geom_sf() + theme_bw()

plot(log(dt$total.pop40), log(dt$black.pop40))

##################################################################################
# the new data set has 300 cities 
# a <- read_sf("/Users/mac/Downloads/mappinginequality.gpkg")
# length(unique(a$city))


##################################################################################
dt2 <- dt %>% as.data.frame() %>% select(-geometry) 
dt2$holc[dt2$holc<=2]=1; dt2$holc[dt2$holc>=3]=2


dt2 %>% group_by(holc) %>% 
  summarise(n=n(), 
            sd.pct.black40=sd(pct.black40, na.rm=T),
            sd.median.rent=sd(median.rent, na.rm=T),
            sd.pm25=sd(pm25, na.rm=T),
            
            miss.pct.black40=sum(is.na(pct.black40)), 
            miss.median.rent=sum(is.na(median.rent)),
            miss.pm25=sum(is.na(pm25)), 
            
            pct.black40=mean(pct.black40, na.rm=T),
            median.rent=mean(median.rent, na.rm=T),
            pm25=mean(pm25, na.rm=T),

  )



library(tidyr)
library(tidyverse)
################################################################
############### Redlining Data Summary Table ###################
################################################################
dt <- readRDS("data/cleaned/final_census1940.rds")

# rate.employed
# wt.mean.rent
# pct.black.pop40

dt$X <- dt$holc
dt$X[dt$X <= 2] = 1; dt$X[dt$X >= 3] = 2 ## binary treatment 
dt$X <- as.factor(dt$X)
levels(dt$X) <- c("AB", "CD")

tb1 <- dt %>% 
  # group_by(holc_grade) %>%
  group_by(X) %>%
  summarise(n=n(), 
            mean(pm25), sd(pm25),
            mean(no2), sd(no2),
            mean(rate.employed*100, na.rm=T), sd(rate.employed*100, na.rm=T),
            mean(wt.mean.rent, na.rm=T), sd(wt.mean.rent, na.rm=T),
            mean(pct.black.pop40*100, na.rm=T), sd(pct.black.pop40*100, na.rm=T))

dim(tb1)

dt %>% 
  # group_by(holc_grade) %>%
  group_by(X) %>%
  summarise(sum(is.na(rate.employed)), # 0 
            sum(is.na(wt.mean.rent)), # 4
            sum(is.na(pct.black.pop40))) # 0

total <- dt %>% 
  summarise(n=n(), 
            mean(pm25), sd(pm25),
            mean(no2), sd(no2),
            mean(rate.employed*100, na.rm=T), sd(rate.employed*100, na.rm=T),
            mean(wt.mean.rent, na.rm=T), sd(wt.mean.rent, na.rm=T),
            mean(pct.black.pop40*100, na.rm=T), sd(pct.black.pop40*100, na.rm=T))

temp1 <- rbind(tb1[, 2:ncol(tb1)], total)
# temp2 <- matrix(c(tb1$holc_grade, "Total"), ncol=1)
temp2 <- matrix(c(tb1$X, "Total"), ncol=1)
tb1 <- cbind(temp2, temp1)

new_names <- tb1[,1]
tb1 <- as.data.frame(as.matrix(t(tb1[,2:12])))

names(tb1) <- new_names
rownames(tb1) <- c(# "Grade", 
  "Count", 
  "mean(PM2.5)", "sd(PM2.5)",
  "mean(NO2)", "sd(NO2)",
  "Rate of Employed (%)", "sd(Rate of Employed)", 
  "Weighted Mean of House Rent ($)", "sd(Weighted Mean of House Rent)",
  "Percentage of Black Population (%)", "sd(Percentage of Black Population)")
write.table(tb1, file="output/redlining_info/redlining_data_summary.txt", row.names=F, quote=F, sep="\t")

library(xtable)
print(xtable(tb1, type = "latex", digits=2), include.rownames = T)


##################################################################################
load("data/cleaned/final_census1940.RData")
dt_centroids <- st_centroid(dt$geometry[!duplicated(dt$location)])

library(sf)
library(ggplot2)
library(rnaturalearth) # 'naturalearth' package to get the state boundaries
library(rnaturalearthdata)

states <- ne_states(country = "united states of america", returnclass = "sf")
states <- states[!(states$name %in% c("Alaska", "Hawaii")),]

pdf("output/redlining_info/map.pdf", width = 6, height = 4)
ggplot() +
  geom_sf(data = states, fill = NA, color = "gray") + 
  geom_sf(data = dt_centroids, aes(geometry = geometry), color = "red", shape = 2, size = 2) + 
  theme_minimal() + 
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
  ) + 
  labs(title = "", x = "", y = "")
dev.off()

pdf("output/redlining_info/map_in_black.pdf", width = 6, height = 4)
ggplot() +
  geom_sf(data = states, fill = NA, color = "gray") + 
  geom_sf(data = dt_centroids, aes(geometry = geometry), color="gray", fill = "black", shape = 21, size = 2, ) + 
  theme_minimal() + 
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
  ) + 
  labs(title = "", x = "", y = "")
dev.off()


##################################################################################
dt <- readRDS("data/cleaned/final_census1940.rds")

raw_diff <- c()
## AB/CD
dt$X <- dt$holc
dt$X[dt$X <= 2] = 1; dt$X[dt$X >= 3] = 2 
temp <- dt %>% group_by(X) %>% summarise(n=n(), pm25=mean(pm25, na.rm=T), no2=mean(no2, na.rm=T)) 
raw_diff <- rbind(raw_diff, c(diff(temp$pm25),diff(temp$no2)))

dt$X <- dt$holc
dt$X[dt$X <= 3] = 1; dt$X[dt$X >= 4] = 2
temp <- dt %>% group_by(X) %>% summarise(n=n(), pm25=mean(pm25, na.rm=T), no2=mean(no2, na.rm=T)) 
raw_diff <- rbind(raw_diff, c(diff(temp$pm25),diff(temp$no2)))

dt$X <- dt$holc
dt$X[dt$X <= 1] = 1; dt$X[dt$X >= 2] = 2
temp <- dt %>% group_by(X) %>% summarise(n=n(), pm25=mean(pm25, na.rm=T), no2=mean(no2, na.rm=T)) 
raw_diff <- rbind(raw_diff, c(diff(temp$pm25),diff(temp$no2)))

print(xtable(data.frame(raw_diff)), include.rownames = T)

# [1,] 0.2592244 2.478640
# [2,] 0.2286647 1.702625
# [3,] 0.3387443 3.123163
