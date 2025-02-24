library(rstpm2)
library(simsurv)
library(survival)
library(dplyr)
library(ggplot2)
library(MASS)
library(sp)
library(gstat)
library(splines)
library(pracma)
library(statmod)
library(parallel)
library(doParallel)
library(foreach)

# set up parallel process. If the number of cores you have is > 20 use 10.
# if there isnt enough memory to use with the multiple cores it will crash the computer.
cl <- makeCluster(detectCores() - 1)
cl = makeCluster(10)
registerDoParallel(cl)

set.seed(111)

# data generation --------------------------------------

gendata <- function(){
  
  lo = 50 # length out of grid # 50 = 2500 | 35 = 1225 | 25 = 625
  nlo = lo*lo # size of data number of grid points
  
  x <- seq(0, 10, length.out = lo)  # X coordinates
  y <- seq(0, 10, length.out = lo)  # Y coordinates
  grid <- expand.grid(x = x, y = y)  # Create a grid of all combinations of x and y
  coordinates(grid) <- ~x + y        # Convert to a spatial object
  proj4string(grid) <- CRS("+proj=longlat +datum=WGS84")  # Set coordinate reference system
  variogram_model <- vgm(psill = 1, model = "Exp", range = 300, nugget = 0.0000001)
  simulated_data <- gstat(formula = dummy ~ 1, locations = grid, dummy = TRUE, beta = 0, model = variogram_model, nmax = 20)
  simulated_values <- predict(simulated_data, newdata = grid, nsim = 1)
  grid$sim = simulated_values$sim1
  
  gsimd = as.data.frame(grid)
  
  n = nlo
  
  X = gsimd$sim
  
  time_max = 5
  beta = log(3)
  alpha = log(2)
  tde = 3
  censoring_rate = 0.45 
  lambda1 = 0.1
  lambda2 = 0.1
  gamma1 = 3
  gamma2 = 1.6
  p = 0.8
  
  X2 = rbinom(n, 1,0.5)
  times = numeric(n)
  event = numeric(n)
  event_times = numeric(n)
  censoring_times = numeric(n)
  
  for(i in 1:n){
    beta1 <- function(t){
      alpha + tde*t
    }
    
    haz <- function(t){
      baseline = (lambda1*gamma1*t^(gamma1-1)*p*exp(-lambda1*t^gamma1) + lambda2*gamma2*t^(gamma2-1)*(1-p)*exp(-lambda2*t^gamma2))/(p*exp(-lambda1*t^gamma1)+(1-p)*exp(-lambda2*t^gamma2))
      ht = baseline * exp(beta*X[i] + beta1(t)*X2[i])
    }

    cumhaz <- function(t){
      integrate(haz, lower = 0, upper = t)$value
    }
    
    u = runif(1, min = 1e-6, max = 0.9999)
    
    target = function(t){
      cumhaz(t) + log(u)
    }
  
    lower = target(1e-6)
    upper = target(time_max)
    
    if (lower * upper < 0){  
    sol = brentDekker(target, a = 1e-6, b = time_max)
    times[i] = sol$root    
    } else {
      times[i] = max(time_max, -lower / (upper - lower) * time_max)
    }
  }
  
  # adjust censoring rate 
  upper_bound = quantile(times, 1 - censoring_rate)
  for (i in 1:500){
    censoring_times = runif(n, 0, upper_bound)
    event = ifelse(times < censoring_times, 1, 0)
    
    actual = 1 - mean(event)
    
    if (abs(actual - censoring_rate) < 0.01){
      break
    }
    upper_bound = upper_bound*(actual/censoring_rate)
  }
  
  data = data.frame(survtime = times, event = event, geo = X, lat = gsimd$y, lng = gsimd$x, X2 = X2)
  return(data)
}

true = log(3) # true parameter
true2 = log(2)
centroid = c(5,5) # centroid

# plot for spatial data if you run process inside gendata

ggplot(gsimd, aes(x = x, y = y, fill = sim)) +
  geom_tile() +
  scale_fill_viridis_c() +
  labs(title = "Simulated Spatial Data", x = "Longitude", y = "Latitude")

# Weight function ----------------------------------

computeweights <- function(data, target, bandwidth){
  distances <- sqrt((data$lng - target[[1]])^2 + (data$lat - target[[2]])^2)
  weights <- case_when(distances < bandwidth ~ (1-(distances/bandwidth)^2)^2, .default = 0) 
  return(weights)
}


# AIC test for degrees of freedom for RP model (switch out i in to the tvc arguement after the df arguement)

data = gendata()
W = computeweights(data, centroid, h)
data$W = W
As = c()
for(i in 1:15){
  rp = stpm2(Surv(survtime, event) ~ geo+X2, tvc = list(geo = 5), data = data, df = i, weights = W)
  As[i] = AIC(rp)
}
idx = which(As == min(As)) # find which index has the min score
idx

# Simulation ---------------------------------------

h = 10 # bandwidth parameter
df = 15 # degrees of freedom for RP model
tdf = 5 # degrees of freedom for time dependent covariate
nsim = 1000 # number of simulations

res <- matrix(NA, nrow = nsim, ncol = 16)

# do parallel
res = foreach(i = 1:nsim, .combine = rbind, .packages = c("sp", "survival", "rstpm2", 'simsurv', 'statmod', 'splines', 'pracma', 'gstat', 'MASS', 'dplyr')) %dopar%{ # simulate 1000 times
  print(i)
  

  data = gendata()

  W <- computeweights(data, centroid, h) # compute weights
  
  data$W = W
  
  # GWR RP
  rp = stpm2(Surv(survtime, event) ~ geo+ X2, tvc = list(X2 = tdf), data = data, df = df, weights = W)

  # GWR Cox
  c = coxph(Surv(survtime, event) ~ geo + X2, data = data, weights = W)

  # Global RP
  grp = stpm2(Surv(survtime, event) ~ geo + X2, tvc = list(X2 = tdf), data = data, df = df)

  # Global Cox
  gc = coxph(Surv(survtime, event) ~ geo+X2, data =  data)

  
  res[i,] = c(coef(rp)[[2]], sqrt(vcov(rp)[2,2]), coef(c)[[1]], sqrt(vcov(c)[1,1]), coef(grp)[[2]], sqrt(vcov(grp)[2,2]), coef(gc)[[1]], sqrt(vcov(gc)[1,1]), coef(rp)[[3]], sqrt(vcov(rp)[3,3]), coef(c)[[2]], sqrt(vcov(c)[2,2]), coef(grp)[[3]], sqrt(vcov(grp)[3,3]), coef(gc)[[2]], sqrt(vcov(gc)[2,2]))
}

MAB = colMeans(abs(res[,c(1,3,5,7)] - true))
MMSE = colMeans((res[,c(1,3,5,7)] - true)^2)
MCP = colMeans(abs(res[,c(1,3,5,7)] - true) <= 1.96 * res[,c(2,4,6,8)])

MAB2 = colMeans(abs(res[,c(9,11,13,15)] - true2))
MMSE2 = colMeans((res[,c(9,11,13,15)] - true2)^2)
MCP2 = colMeans(abs(res[,c(9,11,13,15)] - true2) <= 1.96 * res[,c(10,12,14,16)])

results = data.frame(MAB = MAB, MMSE = MMSE, MCP = MCP, MAB2 = MAB2, MMSE2 = MMSE2, MCP2 = MCP2)
results

stopCluster(cl)
