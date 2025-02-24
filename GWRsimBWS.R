# SIMULATION

# packages
library(dplyr)
library(maps)
library(ggthemes)
library(viridis)
library(ggplot2)
library(stringr)
library(zipcodeR)
library(flexsurv)
library(survival)
library(performance)
library(fitdistrplus)
library(numDeriv)
library(parallel)
library(doParallel)
library(foreach)

cl = makeCluster(detectCores()-8)
cl = makeCluster(10)
registerDoParallel(cl)

# set your working directory for the data files
setwd("C:\\Users\\mzgra\\OneDrive\\Dissertation Research\\Mike Grabel\\Data and Code")

load('gwrdata.Rdata') # load CF GWR data 

# assuming maps package works --------------------------
states = map_data('state') 

a = ggplot(data = states, mapping = aes(x=long, y = lat, group = group)) + coord_fixed(1.3)+
  geom_polygon(color = 'black', fill = NA)

# get state centroids
uscentroid = states %>% group_by(region) %>% summarise(centroid_long = mean(range(long)), centroid_lat = mean(range(lat)))
# adjust certain centroids
uscentroid[uscentroid$region == 'idaho',3] <- 44
uscentroid[uscentroid$region == 'louisiana',2] <- -92 
uscentroid[uscentroid$region == 'florida',2] <- -81.5
uscentroid[uscentroid$region == 'maryland',2] <- -76.8
uscentroid[uscentroid$region == 'new jersey',2] <- -74.4
uscentroid[uscentroid$region == 'massachusetts',3]<-42.4
uscentroid[uscentroid$region == 'michigan',2] <- -85
uscentroid <- uscentroid[-8,] # remove DC

# plot them
a+geom_point(data = uscentroid, aes(centroid_long, centroid_lat), color = 'red', inherit.aes = FALSE, alpha = 1, size = 1)

# ------------------------------------

load("uscentroid.Rdata") # load data if cannot do above process

# ------------------------------------
# for manual testing purposes
i = 33
centroid = uscentroid[i, 2:3] # use i = 33 for Ohio
target = centroid
# ------------------------------------

# function for computing weight matrices
computeweights <- function(data, target, bandwidth){
  distances <- sqrt((data$lng - target[[1]])^2 + (data$lat - target[[2]])^2)
  weights <- case_when(distances < bandwidth ~ (1-(distances/bandwidth)^2)^2, .default = 0) 
  return(weights)
}

# ------------------------------------
# load sampled zipcode data from real data
load("uszcta.Rdata") # sampled 40K zipcodes from the real data to get the spatial dim and added the lat and lng

# visualize location data
a + geom_point(data = uszcta, aes(lng, lat), color = 'red', inherit.aes = FALSE, alpha = 0.5, size = 1)

testsim = gendata(10000)
simd2 = testsim %>% left_join(uszcta, by = 'id') # join sampled zipcodes to the sim data

# Bandwidth Selection -------------- Will take multiple days! ----------------------- 

ComputeBrierBandwidth <- function(data, target, j){ # function that works for simulated data in computing bandwidth
  Bs = c()
  if(j == 24 | j == 29){bandwidth = seq(15, 28, by = 1)}
  else{  bandwidth = seq(4, 28, by = 1)} # range of bandwidths considered
  Bs <- foreach(i = 1:length(bandwidth), .combine = rbind, .packages = c("rstpm2", 'survival', 'dplyr')) %dopar% {
    computeweights <- function(data, target, bandwidth){
      distances <- sqrt((data$lng - target[[1]])^2 + (data$lat - target[[2]])^2)
      weights <- case_when(distances < bandwidth ~ (1-(distances/bandwidth)^2)^2, .default = 0) 
      return(weights)
    }
    W <- computeweights(data, target, bandwidth[i]) # compute weights
    ind = which(W != 0) # find nonzero values
    wdata = data[ind,]
    W = W[ind]
    wdata$W = W
    # run model
    s = stpm2(Surv(encounterage, event) ~ FEV1ps + FEV1Slp + Sex + hof508 + hef508 + SESlow + PA + MRSA + isOnEnzymes + smoke +  
                ozone_conc_n + pm_conc_n + pct_green_n + dep_index_n, data = wdata, df = 3, weights = W)# tvc = list(FEV1ps=2, FEV1Slp=2, SESlow=2, PA=2, MRSA=2, CFRD=2, isOnEnzymes=2, smoke=2), df = 3, weights = W)
    
    # get predictions
    preds = predict(s, type = 'surv')
    # compute Brier Score
    Bs[i] <- mean((preds - wdata$event)^2)
    
    mean((preds - wdata$event)^2)
    
  }
  plot(bandwidth, Bs) 
  idx = which(Bs == min(Bs)) # find which index has the min score
  h = bandwidth[idx] # return optimal bandwidth
  return(h)
}

h = c() # bandwidth parameters
for (i in 1:dim(uscentroid)[1]){ # get bandwidth paramters for simulated data based on states
  centroid = uscentroid[i, 2:3] # select centroid and iterate
  print(i)
  h[i] <- ComputeBrierBandwidth(data, centroid, i) # compute optimal bandwidth for each centroid
}

stopCluster(cl)
save(h, file = 'bandwidthssim2525.Rdata')
