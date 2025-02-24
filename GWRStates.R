#GWR States

# packages
library(dplyr)
library(maps)
library(ggthemes)
library(viridis)
library(ggplot2)
library(stringr)
library(zipcodeR)
library(survival)
library(rstpm2)
library(splines)
library(parallel)
library(doParallel)
library(foreach)

cl <- makeCluster(detectCores() - 1)  # Use most available cores
cl <- makeCluster(10)
registerDoParallel(cl)

setwd("C:\\Users\\")

states = map_data('state') # get US states
# plot
a = ggplot(data = states, mapping = aes(x=long, y = lat, group = group)) + coord_fixed(1.3)+
  geom_polygon(color = 'black', fill = NA)

zp = zip_code_db # may not work with new versions of R zipcodeR package

load('zipcodedb.Rdata') # load the Rdata file if zipcode db doesnt work 

usdata = survcf[survcf$Zipcode %in% zp$zipcode,]

ll = geocode_zip(usdata$Zipcode) # this function is also from zipcodeR and may not work

colnames(ohiocfrdata)[8] <- "zipcode"
colnames(usdata)[8] <- "zipcode"
ghcodata = left_join(usdata, ll, by = 'zipcode')

data <- ghcodata[,-c(2,3,4,5,6,7,9:41, 44, 45:48, 51)]
data = na.omit(data)

save(data, file = "gwrdata.Rdata") # save CF GWR data

load('gwrdata.Rdata') # load CF GWR data if zipcodeR does not work

# extract centroids for each state
uscentroid = states %>% group_by(region) %>% summarise(centroid_long = mean(range(long)), centroid_lat = mean(range(lat)))
# manually adjust centroids for certain states
uscentroid[uscentroid$region == 'idaho',3] <- 44
uscentroid[uscentroid$region == 'louisiana',2] <- -92 
uscentroid[uscentroid$region == 'florida',2] <- -81.5
uscentroid[uscentroid$region == 'maryland',2] <- -76.8
uscentroid[uscentroid$region == 'new jersey',2] <- -74.4
uscentroid[uscentroid$region == 'massachusetts',3]<-42.4
uscentroid[uscentroid$region == 'michigan',2] <- -85
uscentroid <- uscentroid[-8,] # remove DC
# plot centroids
a+geom_point(data = uscentroid, aes(centroid_long, centroid_lat), color = 'red', inherit.aes = FALSE, alpha = 1, size = 1)

save(uscentroid, file = "uscentroid.Rdata") # save centroid data

load("uscentroid.Rdata") # load centroid data


computeweights <- function(data, target, bandwidth){ # function for computing weights based on bi-square kernel
  distances <- sqrt((data$lng - target[[1]])^2 + (data$lat - target[[2]])^2)
  weights <- case_when(distances < bandwidth ~ (1-(distances/bandwidth)^2)^2, .default = 0) 
  return(weights)
}

load('bandwidthssim2525.Rdata') # load optimal h values from sim data brier score

# --------------- will take 12+ hours to run --------------------------------
# HRs = matrix(ncol = 16, nrow = 48) # space for the HR of each area model
# LCIs = matrix(ncol = 16, nrow = 48)
# UCIs = matrix(ncol = 16, nrow = 48)

res <- matrix(NA, nrow = 48, ncol = 56)
colnames(res) <- c("FEV1psHR", "FEV1SlpHR", "SexHR", "hof508HR", "hef508HR", "SESlowHR", "PAHR", "MRSAHR", "isOnEnzymesHR", "smokeHR",
                   "ozone_conc_nHR", "pm_conc_nHR", "pct_green_nHR", "dep_index_nHR", "FEV1psLCI", "FEV1SlpLCI", "SexLCI", "hof508LCI",
                   "hef508LCI", "SESlowLCI", "PALCI", "MRSALCI", "isOnEnzymesLCI", "smokeLCI", "ozone_conc_nLCI", "pm_conc_nLCI", "pct_green_nLCI", "dep_index_nLCI",
                   "FEV1psUCI", "FEV1SlpUCI", "SexUCI", "hof508UCI", "hef508UCI", "SESlowUCI", "PAUCI", "MRSAUCI", "isOnEnzymesUCI", "smokeUCI", "ozone_conc_nUCI", "pm_conc_nUCI",
                   "pct_green_nUCI", "dep_index_nUCI", "FEV1ps", "FEV1Slp", "Sex", "hof508", "hef508", "SESlow", "PA", "MRSA", "isOnEnzymes", "smoke", "ozone_conc_n", "pm_conc_n",
                   "pct_green_n", "dep_index_n")
res = foreach (i = 1:dim(uscentroid)[1], .combine = rbind, .packages = c("survival", "rstpm2", 'dplyr')) %dopar% {
  centroid = uscentroid[i, 2:3] # get target centroid for model
  print(i) # print iteration
  
  W <- computeweights(data, centroid, h[i]) # compute weights
  
  ind = which(W!=0) # find nonzero weights
  W = W[ind]
  wdata = data[ind,] # and corresponding data with nonzero weight
  wdata$W = W # add weights to dataset for stpm2 to work

  s = stpm2(Surv(encounterage, event) ~ FEV1ps + FEV1Slp + Sex + hof508 + hef508 + SESlow + PA + MRSA  + isOnEnzymes + smoke +  
              ozone_conc_n + pm_conc_n + pct_green_n + dep_index_n, data = wdata, df = 3, weights = W) #tvc = list(FEV1ps=2, FEV1Slp=2, SESlow=2, PA=2, MRSA=2, CFRD=2, isOnEnzymes=2, smoke=2), df = 3, weights = W)
  c = summary(s)
  ps = c@coef[,4][2:15]
  result = eform(s)[2:15,] # get HRs and CIs
  
  res[i,] = c(t(result[,1]), t(result[,2]), t(result[,3]), ps)
  
}

stopCluster(cl)

HRs <- res[,1:14] # HRs
LCIs <- res[,15:28] # LCIs
UCIs <- res[,29:42] # UCIs
pss <- res[,43:56]

colnames(HRs) <- c("FEV1psHR", "FEV1SlpHR", "SexHR", "hof508HR", "hef508HR", "SESlowHR", "PAHR", "MRSAHR", "isOnEnzymesHR", "smokeHR",
                   "ozone_conc_nHR", "pm_conc_nHR", "pct_green_nHR", "dep_index_nHR")
colnames(LCIs) <- c("FEV1psLCI", "FEV1SlpLCI", "SexLCI", "hof508LCI",
                    "hef508LCI", "SESlowLCI", "PALCI", "MRSALCI", "isOnEnzymesLCI", "smokeLCI", "ozone_conc_nLCI", "pm_conc_nLCI", "pct_green_nLCI", "dep_index_nLCI")
colnames(UCIs) <- c("FEV1psUCI", "FEV1SlpUCI", "SexUCI", "hof508UCI", "hef508UCI", "SESlowUCI", "PAUCI", "MRSAUCI", "isOnEnzymesUCI", "smokeUCI", "ozone_conc_nUCI", "pm_conc_nUCI",
                     "pct_green_nUCI", "dep_index_nUCI")
#colnames(pss) <- c("FEV1ps", "FEV1Slp", "Sex", "hof508", "hef508", "SESlow", "PA", "MRSA", "isOnEnzymes", "smoke", "ozone_conc_n", "pm_conc_n",
#                    "pct_green_n", "dep_index_n")

save(HRs, file = 'HRs2.Rdata')
save(LCIs, file = 'LCIs2.Rdata')
save(UCIs, file = 'UCIs2.Rdata')
save(pss, file = 'ps2.Rdata')

# ------------------------------------------------------------------------------

fd = cbind(HRs,uscentroid) # combine the HR with the centroid data

fd2 = states %>% left_join(fd, by = 'region') # left join to the us states data

save(fd2, file = 'finaldataHR5.Rdata') # save data

load('finaldataHR4.Rdata')
load('HRs.Rdata')
load('LCIs.Rdata')
load('UCIs.Rdata')

library(writexl) #create xl file
write_xlsx(as.data.frame(HRs), 'HRstates4.xlsx')
write_xlsx(as.data.frame(uscentroid$region), 'states.xlsx')
write_xlsx(as.data.frame(LCIs), 'LCIstates2.xlsx')
write_xlsx(as.data.frame(UCIs), 'UCIstates2.xlsx')
write_xlsx(as.data.frame(pss), 'pvalues2.xlsx')


median = apply(as.matrix(HRs), 2, median)

# CI found by quantile of all HRs

quantile(HRs[,14], 0.025)
quantile(HRs[,14], 0.975)

as.data.frame(pss)

sig_per = colMeans(pss < 0.05)*100


ggplot() +
  geom_polygon(data = fd2, aes(x = long, y = lat, group = group, fill = dep_index_nHR), color = 'black') + coord_fixed(1)+
  geom_text(data = fd2, aes(x = centroid_long, y = centroid_lat, label = round(dep_index_nHR,2)), color = 'black', size = 3, check_overlap = TRUE) +
  scale_fill_gradient(name = 'HR', low = 'lightblue', high = 'red') + 
  labs(title = "Hazard Ratio for Deprivation Index") +
  theme_minimal()

ggplot() +
  geom_polygon(data = fd2, aes(x = long, y = lat, group = group, fill = pct_green_nHR), color = 'black') + coord_fixed(1)+
  geom_text(data = fd2, aes(x = centroid_long, y = centroid_lat, label = round(pct_green_nHR,2)), color = 'black', size = 3, check_overlap = TRUE) +
  scale_fill_gradient(name = 'HR', low = 'lightblue', high = 'red') + 
  labs(title = "Hazard Ratio for Percent Greenspace") +
  theme_minimal()

ggplot() +
  geom_polygon(data = fd2, aes(x = long, y = lat, group = group, fill = pm_conc_nHR), color = 'black') + coord_fixed(1)+
  geom_text(data = fd2, aes(x = centroid_long, y = centroid_lat, label = round(pm_conc_nHR,2)), color = 'black', size = 3, check_overlap = TRUE) +
  scale_fill_gradient(name = 'HR', low = 'lightblue', high = 'red') + 
  labs(title = "Hazard Ratio for PM Concentration") +
  theme_minimal()
ggplot() +
  geom_polygon(data = fd2, aes(x = long, y = lat, group = group, fill = ozone_conc_nHR), color = 'black') + coord_fixed(1)+
  geom_text(data = fd2, aes(x = centroid_long, y = centroid_lat, label = round(ozone_conc_nHR,2)), color = 'black', size = 3, check_overlap = TRUE) +
  scale_fill_gradient(name = 'HR', low = 'lightblue', high = 'red') + 
  labs(title = "Hazard Ratio for Ozone Concentration") +
  theme_minimal()

# 
# ggplot() +
#   geom_polygon(data = fd2, aes(x = long, y = lat, group = group, fill = median_income_n), color = 'black') + coord_fixed(1)+
#   geom_text(data = fd2, aes(x = centroid_long, y = centroid_lat, label = round(median_income_n,2)), color = 'black', size = 3, check_overlap = TRUE) +
#   scale_fill_gradient(name = 'HR', low = 'lightblue', high = 'red') + 
#   labs(title = "Hazard Ratio for Median Income") +
#   theme_minimal()


ggplot() +
  geom_polygon(data = fd2, aes(x = long, y = lat, group = group, fill = FEV1psHR), color = 'black') + coord_fixed(1)+
  geom_text(data = fd2, aes(x = centroid_long, y = centroid_lat, label = round(FEV1psHR,2)), color = 'black', size = 3, check_overlap = TRUE) +
  scale_fill_gradient(name = 'HR', low = 'lightblue', high = 'red') + 
  labs(title = "Hazard Ratio for FEV1ps") +
  theme_minimal()

ggplot() +
  geom_polygon(data = fd2, aes(x = long, y = lat, group = group, fill = FEV1SlpHR), color = 'black') + coord_fixed(1)+
  geom_text(data = fd2, aes(x = centroid_long, y = centroid_lat, label = round(FEV1SlpHR,2)), color = 'black', size = 3, check_overlap = TRUE) +
  scale_fill_gradient(name = 'HR', low = 'lightblue', high = 'red') + 
  labs(title = "Hazard Ratio for FEV1Slp") +
  theme_minimal()

ggplot() +
  geom_polygon(data = fd2, aes(x = long, y = lat, group = group, fill = SexHR), color = 'black') + coord_fixed(1)+
  geom_text(data = fd2, aes(x = centroid_long, y = centroid_lat, label = round(SexHR,2)), color = 'black', size = 3, check_overlap = TRUE) +
  scale_fill_gradient(name = 'HR', low = 'lightblue', high = 'red') + 
  labs(title = "Hazard Ratio for Sex") +
  theme_minimal()

ggplot() +
  geom_polygon(data = fd2, aes(x = long, y = lat, group = group, fill = hof508HR), color = 'black') + coord_fixed(1)+
  geom_text(data = fd2, aes(x = centroid_long, y = centroid_lat, label = round(hof508HR,2)), color = 'black', size = 3, check_overlap = TRUE) +
  scale_fill_gradient(name = 'HR', low = 'lightblue', high = 'red') + 
  labs(title = "Hazard Ratio for Homozygous F508") +
  theme_minimal()

ggplot() +
  geom_polygon(data = fd2, aes(x = long, y = lat, group = group, fill = hef508HR), color = 'black') + coord_fixed(1)+
  geom_text(data = fd2, aes(x = centroid_long, y = centroid_lat, label = round(hef508HR,2)), color = 'black', size = 3, check_overlap = TRUE) +
  scale_fill_gradient(name = 'HR', low = 'lightblue', high = 'red') + 
  labs(title = "Hazard Ratio for Hetrozygous F508") +
  theme_minimal()

ggplot() +
  geom_polygon(data = fd2, aes(x = long, y = lat, group = group, fill = SESlowHR), color = 'black') + coord_fixed(1)+
  geom_text(data = fd2, aes(x = centroid_long, y = centroid_lat, label = round(SESlowHR,2)), color = 'black', size = 3, check_overlap = TRUE) +
  scale_fill_gradient(name = 'HR', low = 'lightblue', high = 'red') + 
  labs(title = "Hazard Ratio for Low Socioeconomic Status") +
  theme_minimal()

ggplot() +
  geom_polygon(data = fd2, aes(x = long, y = lat, group = group, fill = PAHR), color = 'black') + coord_fixed(1)+
  geom_text(data = fd2, aes(x = centroid_long, y = centroid_lat, label = round(PAHR,2)), color = 'black', size = 3, check_overlap = TRUE) +
  scale_fill_gradient(name = 'HR', low = 'lightblue', high = 'red') + 
  labs(title = "Hazard Ratio for PA") +
  theme_minimal()

ggplot() +
  geom_polygon(data = fd2, aes(x = long, y = lat, group = group, fill = MRSAHR), color = 'black') + coord_fixed(1)+
  geom_text(data = fd2, aes(x = centroid_long, y = centroid_lat, label = round(MRSAHR,2)), color = 'black', size = 3, check_overlap = TRUE) +
  scale_fill_gradient(name = 'HR', low = 'lightblue', high = 'red') + 
  labs(title = "Hazard Ratio for MRSA") +
  theme_minimal()

# ggplot() +
#   geom_polygon(data = fd2, aes(x = long, y = lat, group = group, fill = CFRD), color = 'black') + coord_fixed(1)+
#   geom_text(data = fd2, aes(x = centroid_long, y = centroid_lat, label = round(CFRD,2)), color = 'black', size = 3, check_overlap = TRUE) +
#   scale_fill_gradient(name = 'HR', low = 'lightblue', high = 'red') + 
#   labs(title = "Hazard Ratio for CFRD") +
#   theme_minimal()

ggplot() +
  geom_polygon(data = fd2, aes(x = long, y = lat, group = group, fill = isOnEnzymesHR), color = 'black') + coord_fixed(1)+
  geom_text(data = fd2, aes(x = centroid_long, y = centroid_lat, label = round(isOnEnzymesHR,2)), color = 'black', size = 3, check_overlap = TRUE) +
  scale_fill_gradient(name = 'HR', low = 'lightblue', high = 'red') + 
  labs(title = "Hazard Ratio for On Enzymes") +
  theme_minimal()

ggplot() +
  geom_polygon(data = fd2, aes(x = long, y = lat, group = group, fill = smokeHR), color = 'black') + coord_fixed(1)+
  geom_text(data = fd2, aes(x = centroid_long, y = centroid_lat, label = round(smokeHR,2)), color = 'black', size = 3, check_overlap = TRUE) +
  scale_fill_gradient(name = 'HR', low = 'lightblue', high = 'red') + 
  labs(title = "Hazard Ratio for Smoke Exposure") +
  theme_minimal()

