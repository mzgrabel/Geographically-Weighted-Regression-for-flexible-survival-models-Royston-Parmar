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
# cl = makeCluster(10)
registerDoParallel(cl)

set.seed(111)

beta1_t <- function(t) {
  0.7 + 0.5 * log(t)
}
gendata = function(){
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
  
  beta1_t <- function(t) {
    0.7 + 0.5 * log(t)
  }
  
  X2 = data.frame(rbinom(n, 1,0.33))
  colnames(X2) <- 'X2'
  cov = data.frame(X = X, X2 = X2)
  # data = simsurv(dist = 'weibull', lambdas = 0.1, gammas = 1.3, x = cov, 
  #                betas = c(X = log(3), X2 = log(2)), tde = c(X2 = 1), tdefunction = beta1_t, maxt = 5)
  
  # data = simsurv(loghazard = function(t, x, betas) -18 + 7.3*t-11.5*t^0.5*log(t) + 9.5*t^0.5 + x$X*betas['X']+ x$X2*beta1_t(t), x = cov, 
                 # betas = c(X = log(3)), maxt = 5)
  data = simsurv(loghazard = function(t, x, betas) 0.0130 * exp(0.2043 * t) + x$X*betas['X']+ x$X2*beta1_t(t), x = cov, 
                 betas = c(X = log(3)), maxt = 5)  
  dat = cbind(data, X, X2, lat = gsimd$y, lng = gsimd$x)
  return(dat)
}

computeweights <- function(data, target, bandwidth){
  distances <- sqrt((data$lng - target[[1]])^2 + (data$lat - target[[2]])^2)
  weights <- case_when(distances < bandwidth ~ (1-(distances/bandwidth)^2)^2, .default = 0) 
  return(weights)
}

centroid = c(5,5)

h = 10 # bandwidth parameter
df = 2 # degrees of freedom for RP model
tdf = 2 # degrees of freedom for time dependent covariate
nsim = 1000 # number of simulations

# do parallel
res = foreach(i = 1:nsim, .combine = rbind, .packages = c("sp", "survival", "rstpm2", 'simsurv', 'statmod', 'splines', 'pracma', 'gstat', 'MASS', 'dplyr')) %dopar%{ 

  data = gendata()
  
  W <- computeweights(data, centroid, h)
  data$W = W
  
  # GWR-RP
  rp = stpm2(Surv(eventtime, status) ~ X + X2, tvc = list(X2 = tdf), data = data, df = df, weights = W)
  
  beta = coef(rp)[2]
  se = sqrt(rp@vcov[2,2])
  
  survival_rp_1 <- predict(rp, newdata = data.frame(eventtime = c(1,2,3,4,5), X = 0, X2 = 1), type = 'surv', se.fit = TRUE)
  survival_rp_0 <- predict(rp, newdata = data.frame(eventtime = c(1,2,3,4,5), X = 0, X2 = 0), type = 'surv', se.fit = TRUE)
  
  cumhaz_rp_1 <- predict(rp, newdata = data.frame(eventtime = c(1,2,3,4,5), X = 0, X2 = 1), type = 'cumhaz', se.fit = TRUE)
  cumhaz_rp_0 <- predict(rp, newdata = data.frame(eventtime = c(1,2,3,4,5), X = 0, X2 = 0), type = 'cumhaz', se.fit = TRUE)
  
  # GWR-Cox
  c = coxph(Surv(eventtime, status) ~ X + X2, data = data, weights = W)
  
  betac = coef(c)[1]
  sec = sqrt(c$var[1,1])
  
  survival_c_1 = predict(c, newdata = data.frame(eventtime = c(1,2,3,4,5), status = 1, X = 0, X2 = 1), type = 'surv', se.fit = TRUE)
  survival_c_0 = predict(c, newdata = data.frame(eventtime = c(1,2,3,4,5), status = 1, X = 0, X2 = 0), type = 'surv', se.fit = TRUE)
  cslower1 = survival_c_1$fit - 1.96 * survival_c_1$se.fit
  csupper1 = survival_c_1$fit + 1.96 * survival_c_1$se.fit
  cslower0 = survival_c_0$fit - 1.96 * survival_c_0$se.fit
  csupper0 = survival_c_0$fit + 1.96 * survival_c_0$se.fit
  
  
  cumhaz_c_1 <- predict(c, newdata = data.frame(eventtime = c(1,2,3,4,5), status = 1, X = 0, X2 = 1), type = "expected", se.fit = TRUE)
  cumhaz_c_0 <- predict(c, newdata = data.frame(eventtime = c(1,2,3,4,5), status = 1, X = 0, X2 = 0), type = "expected", se.fit = TRUE)
  chlower1 = cumhaz_c_1$fit - 1.96 * cumhaz_c_1$se.fit
  chupper1 = cumhaz_c_1$fit + 1.96 * cumhaz_c_1$se.fit
  chlower0 = cumhaz_c_0$fit - 1.96 * cumhaz_c_0$se.fit
  chupper0 = cumhaz_c_0$fit + 1.96 * cumhaz_c_0$se.fit
  
  # Global RP
  grp = stpm2(Surv(eventtime, status) ~ X + X2, tvc = list(X2 = tdf), data = data, df = df)
  
  betagrp = coef(grp)[2]
  segrp = sqrt(grp@vcov[2,2])
  
  survival_grp_1 <- predict(grp, newdata = data.frame(eventtime = c(1,2,3,4,5), X = 0, X2 = 1), type = 'surv', se.fit = TRUE)
  survival_grp_0 <- predict(grp, newdata = data.frame(eventtime = c(1,2,3,4,5), X = 0, X2 = 0), type = 'surv', se.fit = TRUE)
  
  cumhaz_grp_1 <- predict(grp, newdata = data.frame(eventtime = c(1,2,3,4,5), X = 0, X2 = 1), type = 'cumhaz', se.fit = TRUE)
  cumhaz_grp_0 <- predict(grp, newdata = data.frame(eventtime = c(1,2,3,4,5), X = 0, X2 = 0), type = 'cumhaz', se.fit = TRUE)
  
  # Global Cox
  gc = coxph(Surv(eventtime, status) ~ X + X2, data = data)
  
  betagc = coef(gc)[1]
  segc = sqrt(gc$var[1,1])
  
  survival_gc_1 = predict(c, newdata = data.frame(eventtime = c(1,2,3,4,5), status = 1, X = 0, X2 = 1), type = 'surv', se.fit = TRUE)
  survival_gc_0 = predict(c, newdata = data.frame(eventtime = c(1,2,3,4,5), status = 1, X = 0, X2 = 0), type = 'surv', se.fit = TRUE)
  gcslower1 = survival_gc_1$fit - 1.96 * survival_gc_1$se.fit
  gcsupper1 = survival_gc_1$fit + 1.96 * survival_gc_1$se.fit
  gcslower0 = survival_gc_0$fit - 1.96 * survival_gc_0$se.fit
  gcsupper0 = survival_gc_0$fit + 1.96 * survival_gc_0$se.fit
  
  
  cumhaz_gc_1 <- predict(c, newdata = data.frame(eventtime = c(1,2,3,4,5), status = 1, X = 0, X2 = 1), type = "expected", se.fit = TRUE)
  cumhaz_gc_0 <- predict(c, newdata = data.frame(eventtime = c(1,2,3,4,5), status = 1, X = 0, X2 = 0), type = "expected", se.fit = TRUE)
  gchlower1 = cumhaz_gc_1$fit - 1.96 * cumhaz_gc_1$se.fit
  gchupper1 = cumhaz_gc_1$fit + 1.96 * cumhaz_gc_1$se.fit
  gchlower0 = cumhaz_gc_0$fit - 1.96 * cumhaz_gc_0$se.fit
  gchupper0 = cumhaz_gc_0$fit + 1.96 * cumhaz_gc_0$se.fit
  
  iteration_df <- data.frame(
    iteration = i,
    # RP
    betarp = beta,
    serp = se,
    survival_rp_1 = survival_rp_1$Estimate,
    survival_rp_0 = survival_rp_0$Estimate,
    cumhaz_rp_1 = cumhaz_rp_1$Estimate,
    cumhaz_rp_0 = cumhaz_rp_0$Estimate,
    rpLs1 = survival_rp_1$lower,
    rpUs1 = survival_rp_1$upper,
    rpLs0 = survival_rp_0$lower,
    rpUs0 = survival_rp_0$upper,
    rpLh1 = cumhaz_rp_1$lower,
    rpUh1 = cumhaz_rp_1$upper,
    rpLh0 = cumhaz_rp_0$lower,
    rpUh0 = cumhaz_rp_0$upper,
    # cox
    betac = betac,
    sec = sec,
    survival_c_1 = survival_c_1$fit,
    survival_c_0 = survival_c_0$fit,
    cumhaz_c_1 = cumhaz_c_1$fit,
    cumhaz_c_0 = cumhaz_c_0$fit,
    cLs1 = cslower1,
    cUs1 = csupper1,
    cLs0 = cslower0,
    cUs0 = csupper1,
    cLh1 = chlower1,
    cUh1 = chupper1,
    cLh0 = chlower0,
    cUh0 = chupper0,
    #GRP
    betagrp = betagrp,
    serp = segrp,
    survival_grp_1 = survival_grp_1$Estimate,
    survival_grp_0 = survival_grp_0$Estimate,
    cumhaz_grp_1 = cumhaz_grp_1$Estimate,
    cumhaz_grp_0 = cumhaz_grp_0$Estimate,
    grpLs1 = survival_grp_1$lower,
    grpUs1 = survival_grp_1$upper,
    grpLs0 = survival_grp_0$lower,
    grpUs0 = survival_grp_0$upper,
    grpLh1 = cumhaz_grp_1$lower,
    grpUh1 = cumhaz_grp_1$upper,
    grpLh0 = cumhaz_grp_0$lower,
    grpUh0 = cumhaz_grp_0$upper,
    #GCox
    betagc = betagc,
    segc = segc,
    survival_gc_1 = survival_gc_1$fit,
    survival_gc_0 = survival_gc_0$fit,
    cumhaz_gc_1 = cumhaz_gc_1$fit,
    cumhaz_gc_0 = cumhaz_gc_0$fit,
    gcLs1 = gcslower1,
    gcUs1 = gcsupper1,
    gcLs0 = gcslower0,
    gcUs0 = gcsupper1,
    gcLh1 = gchlower1,
    gcUh1 = gchupper1,
    gcLh0 = gchlower0,
    gcUh0 = gchupper0
  )
  
  # Return the iteration's data frame (this will be combined in the main list)
  return(iteration_df)
}





# res = iteration_df

#GWR-RP
betas = matrix(0, nrow = nsim, ncol = 5)
se = matrix(0, nrow = nsim, ncol = 5)
surv1 = matrix(0, nrow = nsim, ncol = 5)
surv0 = matrix(0, nrow = nsim, ncol = 5)
cumhaz1 = matrix(0, nrow = nsim, ncol = 5)
cumhaz0 = matrix(0, nrow = nsim, ncol = 5)

rpLs1 = matrix(0, nrow = nsim, ncol = 5)
rpUs1 = matrix(0, nrow = nsim, ncol = 5)
rpLs0 = matrix(0, nrow = nsim, ncol = 5)
rpUs0 = matrix(0, nrow = nsim, ncol = 5)
rpLh1 = matrix(0, nrow = nsim, ncol = 5)
rpUh1 = matrix(0, nrow = nsim, ncol = 5)
rpLh0 = matrix(0, nrow = nsim, ncol = 5)
rpUh0 = matrix(0, nrow = nsim, ncol = 5)

#GWR-Cox
betasc = matrix(0, nrow = nsim, ncol = 5)
sec = matrix(0, nrow = nsim, ncol = 5)
surv1c = matrix(0, nrow = nsim, ncol = 5)
surv0c = matrix(0, nrow = nsim, ncol = 5)
cumhaz1c = matrix(0, nrow = nsim, ncol = 5)
cumhaz0c = matrix(0, nrow = nsim, ncol = 5)

cLs1 = matrix(0, nrow = nsim, ncol = 5)
cUs1 = matrix(0, nrow = nsim, ncol = 5)
cLs0 = matrix(0, nrow = nsim, ncol = 5)
cUs0 = matrix(0, nrow = nsim, ncol = 5)
cLh1 = matrix(0, nrow = nsim, ncol = 5)
cUh1 = matrix(0, nrow = nsim, ncol = 5)
cLh0 = matrix(0, nrow = nsim, ncol = 5)
cUh0 = matrix(0, nrow = nsim, ncol = 5)
# GRP
betasgrp = matrix(0, nrow = nsim, ncol = 5)
segrp = matrix(0, nrow = nsim, ncol = 5)
surv1grp = matrix(0, nrow = nsim, ncol = 5)
surv0grp = matrix(0, nrow = nsim, ncol = 5)
cumhaz1grp = matrix(0, nrow = nsim, ncol = 5)
cumhaz0grp = matrix(0, nrow = nsim, ncol = 5)

grpLs1 = matrix(0, nrow = nsim, ncol = 5)
grpUs1 = matrix(0, nrow = nsim, ncol = 5)
grpLs0 = matrix(0, nrow = nsim, ncol = 5)
grpUs0 = matrix(0, nrow = nsim, ncol = 5)
grpLh1 = matrix(0, nrow = nsim, ncol = 5)
grpUh1 = matrix(0, nrow = nsim, ncol = 5)
grpLh0 = matrix(0, nrow = nsim, ncol = 5)
grpUh0 = matrix(0, nrow = nsim, ncol = 5)
#GCox
betasgc = matrix(0, nrow = nsim, ncol = 5)
segc = matrix(0, nrow = nsim, ncol = 5)
surv1gc = matrix(0, nrow = nsim, ncol = 5)
surv0gc = matrix(0, nrow = nsim, ncol = 5)
cumhaz1gc = matrix(0, nrow = nsim, ncol = 5)
cumhaz0gc = matrix(0, nrow = nsim, ncol = 5)

gcLs1 = matrix(0, nrow = nsim, ncol = 5)
gcUs1 = matrix(0, nrow = nsim, ncol = 5)
gcLs0 = matrix(0, nrow = nsim, ncol = 5)
gcUs0 = matrix(0, nrow = nsim, ncol = 5)
gcLh1 = matrix(0, nrow = nsim, ncol = 5)
gcUh1 = matrix(0, nrow = nsim, ncol = 5)
gcLh0 = matrix(0, nrow = nsim, ncol = 5)
gcUh0 = matrix(0, nrow = nsim, ncol = 5)

for(i in 1:nsim){
  new = res[res$iteration == i,]
  #GWR-Rp
  betas[i,] <- new$betarp
  se[i,] = new$serp
  surv1[i,] = new$survival_rp_1
  surv0[i,] = new$survival_rp_0
  cumhaz1[i,] = new$cumhaz_rp_1
  cumhaz0[i,] = new$cumhaz_rp_0
  
  rpLs1[i,] = new$rpLs1
  rpUs1[i,] = new$rpUs1
  rpLs0[i,] = new$rpLs0
  rpUs0[i,] = new$rpUs0
  rpLh1[i,] = new$rpLh1
  rpUh1[i,] = new$rpUh1
  rpLh0[i,] = new$rpLh0
  rpUh0[i,] = new$rpUh0
  #GWR-Cox
  betasc[i,] = new$betac
  sec[i,] = new$sec
  surv1c[i,] = new$survival_c_1
  surv0c[i,] = new$survival_c_0
  cumhaz1c[i,] = new$cumhaz_c_1
  cumhaz0c[i,] = new$cumhaz_c_0
  
  cLs1[i,] = new$cLs1
  cUs1[i,] = new$cUs1
  cLs0[i,] = new$cLs0
  cUs0[i,] = new$cUs0
  cLh1[i,] = new$cLh1
  cUh1[i,] = new$cUh1
  cLh0[i,] = new$cLh0
  cUh0[i,] = new$cUh0
  #GRP
  betasgrp[i,] <- new$betagrp
  se[i,] = new$serp
  surv1grp[i,] = new$survival_grp_1
  surv0grp[i,] = new$survival_grp_0
  cumhaz1grp[i,] = new$cumhaz_grp_1
  cumhaz0grp[i,] = new$cumhaz_grp_0
  
  grpLs1[i,] = new$grpLs1
  grpUs1[i,] = new$grpUs1
  grpLs0[i,] = new$grpLh0
  grpUs0[i,] = new$grpUs0
  grpLh1[i,] = new$grpLh1
  grpUh1[i,] = new$grpUh1
  grpLh0[i,] = new$grpLh0
  grpUh0[i,] = new$grpLh0
  #GCox
  betasgc[i,] = new$betagc
  segc[i,] = new$segc
  surv1gc[i,] = new$survival_gc_1
  surv0gc[i,] = new$survival_gc_0
  cumhaz1gc[i,] = new$cumhaz_gc_1
  cumhaz0gc[i,] = new$cumhaz_gc_0
  
  gcLs1[i,] = new$gcLs1
  gcUs1[i,] = new$gcUs1
  gcLs0[i,] = new$gcLs0
  gcUs0[i,] = new$gcUs0
  gcLh1[i,] = new$gcLh1
  gcUh1[i,] = new$gcUh1
  gcLh0[i,] = new$gcLh0
  gcUh0[i,] = new$gcUh0
}

true = log(3)
true_hazard <- function(t, X, X2) {
  baseline_hazard <- exp(0.0130 * exp(0.2043 * t)) #exp(-18 + 7.3*t-11.5*t^0.5*log(t) + 9.5*t^0.5)
  time_dependent_effect <- beta1_t(t)
  hazard <- baseline_hazard * exp(log(3) * X + time_dependent_effect * X2)
  return(hazard)
}

cumulative_hazard <- function(t, X, X2) {
  integral <- integrate(function(u) true_hazard(u, X, X2), lower = 0, upper = t)$value
  return(integral)
}

specific_times <- c(1, 2, 3, 4, 5)
true_hazard_rates <- sapply(specific_times, function(t) cumulative_hazard(t,0,1))
true_hazard_rates0 <- sapply(specific_times, function(t) cumulative_hazard(t,0,0))

true_survival <- function(t, X, X2) {
  ch <- cumulative_hazard(t, X, X2)
  survival <- exp(-ch)
  return(survival)
}

true_survival_probabilities <- sapply(specific_times, function(t) true_survival(t,0, 1))
true_survival_probabilities0 <- sapply(specific_times, function(t) true_survival(t,0, 0))

coverage = function(true, L, U){
  coverage_count = 0
  for(i in 1:nsim){
    if(true >= L[i] && true <= U[i]){
      coverage_count = coverage_count + 1
    }
  }
  coverage_rate <- coverage_count / nsim
  return(coverage_rate)
}
# GWR-RP ----------------------------------------------

biasrp = mean(abs(betas[,1] - true))
covrp = mean(abs(betas[,1] - true) <= 1.96 * se)


s11 = abs(mean(surv1[,1]) - true_survival_probabilities[1])
s12 = abs(mean(surv1[,2]) - true_survival_probabilities[2])
s13 = abs(mean(surv1[,3]) - true_survival_probabilities[3])
s14 = abs(mean(surv1[,4]) - true_survival_probabilities[4])
s15 = abs(mean(surv1[,5]) - true_survival_probabilities[5])

rps1 = rbind(s11,s12,s13,s14,s15)
rownames(rps1) = c('1','2','3','4','5')
colnames(rps1) = 'GWRRP Survival X = 1'

s01 = abs(mean(surv0[,1]) - true_survival_probabilities0[1])
s02 = abs(mean(surv0[,2]) - true_survival_probabilities0[2])
s03 = abs(mean(surv0[,3]) - true_survival_probabilities0[3])
s04 = abs(mean(surv0[,4]) - true_survival_probabilities0[4])
s05 = abs(mean(surv0[,5]) - true_survival_probabilities0[5])

rps0 = rbind(s01,s02,s03,s04,s05)
rownames(rps0) = c('1','2','3','4','5')
colnames(rps0) = 'GWRRP Survival X = 0'

h11 = abs(mean(cumhaz1[,1]) - true_hazard_rates[1])
h12 = abs(mean(cumhaz1[,2]) - true_hazard_rates[2])
h13 = abs(mean(cumhaz1[,3]) - true_hazard_rates[3])
h14 = abs(mean(cumhaz1[,4]) - true_hazard_rates[4])
h15 = abs(mean(cumhaz1[,5]) - true_hazard_rates[5])

rph1 = rbind(h11,h12,h13,h14,h15)
rownames(rph1) = c('1','2','3','4','5')
colnames(rph1) = 'GWRRP Cumhaz X = 1'

h01 = abs(mean(cumhaz0[,1]) - true_hazard_rates0[1])
h02 = abs(mean(cumhaz0[,2]) - true_hazard_rates0[2])
h03 = abs(mean(cumhaz0[,3]) - true_hazard_rates0[3])
h04 = abs(mean(cumhaz0[,4]) - true_hazard_rates0[4])
h05 = abs(mean(cumhaz0[,5]) - true_hazard_rates0[5])

rph0 = rbind(h01,h02,h03,h04,h05)
rownames(rph0) = c('1','2','3','4','5')
colnames(rph0) = 'GWRRP Cumhaz X = 0'

cs11 = coverage(true_survival_probabilities[1], rpLs1[,1], rpUs1[,1])
cs12 = coverage(true_survival_probabilities[2], rpLs1[,2], rpUs1[,2])
cs13 = coverage(true_survival_probabilities[3], rpLs1[,3], rpUs1[,3])
cs14 = coverage(true_survival_probabilities[4], rpLs1[,4], rpUs1[,4])
cs15 = coverage(true_survival_probabilities[5], rpLs1[,5], rpUs1[,5])

rpcs1 = rbind(cs11,cs12,cs13,cs14,cs15)
rownames(rpcs1) <- c('1','2','3','4','5')
colnames(rpcs1) <- 'GWRRP Survival Coverage X = 1'

cs01 = coverage(true_survival_probabilities0[1], rpLs0[,1], rpUs0[,1])
cs02 = coverage(true_survival_probabilities0[2], rpLs0[,2], rpUs0[,2])
cs03 = coverage(true_survival_probabilities0[3], rpLs0[,3], rpUs0[,3])
cs04 = coverage(true_survival_probabilities0[4], rpLs0[,4], rpUs0[,4])
cs05 = coverage(true_survival_probabilities0[5], rpLs0[,5], rpUs0[,5])

rpcs0 = rbind(cs01,cs02,cs03,cs04,cs05)
rownames(rpcs0) <- c('1','2','3','4','5')
colnames(rpcs0) <- 'GWRRP Survival Coverage X = 0'

ch11 = coverage(true_hazard_rates[1], rpLh1[,1], rpUh1[,1])
ch12 = coverage(true_hazard_rates[2], rpLh1[,2], rpUh1[,2])
ch13 = coverage(true_hazard_rates[3], rpLh1[,3], rpUh1[,3])
ch14 = coverage(true_hazard_rates[4], rpLh1[,4], rpUh1[,4])
ch15 = coverage(true_hazard_rates[5], rpLh1[,5], rpUh1[,5])

rpch1 = rbind(ch11,ch12,ch13,ch14,ch15)
rownames(rpch1) <- c('1','2','3','4','5')
colnames(rpch1) <- 'GWRRP Cumhaz Coverage X = 1'

ch01 = coverage(true_hazard_rates0[1], rpLh0[,1], rpUh0[,1])
ch02 = coverage(true_hazard_rates0[2], rpLh0[,2], rpUh0[,2])
ch03 = coverage(true_hazard_rates0[3], rpLh0[,3], rpUh0[,3])
ch04 = coverage(true_hazard_rates0[4], rpLh0[,4], rpUh0[,4])
ch05 = coverage(true_hazard_rates0[5], rpLh0[,5], rpUh0[,5])

rpch0 = rbind(ch01,ch02,ch03,ch04,ch05)
rownames(rpch0) <- c('1','2','3','4','5')
colnames(rpch0) <- 'GWRRP Cumhaz Coverage X = 0'

# GWR-Cox -----------------------------------------------------------

biasc = mean(abs(betasc[,1] - true))
covc = mean(abs(betasc[,1] - true) <= 1.96 * sec)


s11c = abs(mean(surv1c[,1]) - true_survival_probabilities[1])
s12c = abs(mean(surv1c[,2]) - true_survival_probabilities[2])
s13c = abs(mean(surv1c[,3]) - true_survival_probabilities[3])
s14c = abs(mean(surv1c[,4]) - true_survival_probabilities[4])
s15c = abs(mean(surv1c[,5]) - true_survival_probabilities[5])

cs1 = rbind(s11c,s12c,s13c,s14c,s15c)
rownames(cs1) <- c('1','2','3','4','5')
colnames(cs1) <- 'GWC Survival X = 1'

s01c = abs(mean(surv0c[,1]) - true_survival_probabilities0[1])
s02c = abs(mean(surv0c[,2]) - true_survival_probabilities0[2])
s03c = abs(mean(surv0c[,3]) - true_survival_probabilities0[3])
s04c = abs(mean(surv0c[,4]) - true_survival_probabilities0[4])
s05c = abs(mean(surv0c[,5]) - true_survival_probabilities0[5])

cs0 = rbind(s01c,s02c,s03c,s04c,s05c)
rownames(cs0) <- c('1','2','3','4','5')
colnames(cs0) <- 'GWC Survival X = 0'

h11c = abs(mean(cumhaz1c[,1]) - true_hazard_rates[1])
h12c = abs(mean(cumhaz1c[,2]) - true_hazard_rates[2])
h13c = abs(mean(cumhaz1c[,3]) - true_hazard_rates[3])
h14c = abs(mean(cumhaz1c[,4]) - true_hazard_rates[4])
h15c = abs(mean(cumhaz1c[,5]) - true_hazard_rates[5])

ch1 = rbind(h11c,h12c,h13c,h14c,h15c)
rownames(ch1) <- c('1','2','3','4','5')
colnames(ch1) <- 'GWC Cumhaz X = 1'

h01c = abs(mean(cumhaz0c[,1]) - true_hazard_rates0[1])
h02c = abs(mean(cumhaz0c[,2]) - true_hazard_rates0[2])
h03c = abs(mean(cumhaz0c[,3]) - true_hazard_rates0[3])
h04c = abs(mean(cumhaz0c[,4]) - true_hazard_rates0[4])
h05c = abs(mean(cumhaz0c[,5]) - true_hazard_rates0[5])

ch0 = rbind(h01c,h02c,h03c,h04c,h05c)
rownames(ch0) <- c('1','2','3','4','5')
colnames(ch0) <- 'GWC Cumhaz X = 0'

cs11c = coverage(true_survival_probabilities[1], cLs1[,1], cUs1[,1])
cs12c = coverage(true_survival_probabilities[2], cLs1[,2], cUs1[,2])
cs13c = coverage(true_survival_probabilities[3], cLs1[,3], cUs1[,3])
cs14c = coverage(true_survival_probabilities[4], cLs1[,4], cUs1[,4])
cs15c = coverage(true_survival_probabilities[5], cLs1[,5], cUs1[,5])

ccs1 = rbind(cs11c,cs12c,cs13c,cs14c,cs15c)
rownames(ccs1) <- c('1','2','3','4','5')
colnames(ccs1) <- 'GWC Survival Coverage X = 1'

cs01c = coverage(true_survival_probabilities0[1], cLs0[,1], cUs0[,1])
cs02c = coverage(true_survival_probabilities0[2], cLs0[,2], cUs0[,2])
cs03c = coverage(true_survival_probabilities0[3], cLs0[,3], cUs0[,3])
cs04c = coverage(true_survival_probabilities0[4], cLs0[,4], cUs0[,4])
cs05c = coverage(true_survival_probabilities0[5], cLs0[,5], cUs0[,5])

ccs0 = rbind(cs01c,cs02c,cs03c,cs04c,cs05c)
rownames(ccs0) <- c('1','2','3','4','5')
colnames(ccs0) <- 'GWC Survival Coverage X = 0'

ch11c = coverage(true_hazard_rates[1], cLh1[,1], cUh1[,1])
ch12c = coverage(true_hazard_rates[2], cLh1[,2], cUh1[,2])
ch13c = coverage(true_hazard_rates[3], cLh1[,3], cUh1[,3])
ch14c = coverage(true_hazard_rates[4], cLh1[,4], cUh1[,4])
ch15c = coverage(true_hazard_rates[5], cLh1[,5], cUh1[,5])

cch1 = rbind(ch11c,ch12c,ch13c,ch14c,ch15c)
rownames(cch1) <- c('1','2','3','4','5')
colnames(cch1) <- 'GWC Cumhaz Coverage X = 1'

ch01c = coverage(true_hazard_rates0[1], cLh0[,1], cUh0[,1])
ch02c = coverage(true_hazard_rates0[2], cLh0[,2], cUh0[,2])
ch03c = coverage(true_hazard_rates0[3], cLh0[,3], cUh0[,3])
ch04c = coverage(true_hazard_rates0[4], cLh0[,4], cUh0[,4])
ch05c = coverage(true_hazard_rates0[5], cLh0[,5], cUh0[,5])

cch0 = rbind(ch01c,ch02c,ch03c,ch04c,ch05c)
rownames(cch0) <- c('1','2','3','4','5')
colnames(cch0) <- 'GWC Cumhaz Coverage X = 0'

# G RP -----------------------------------------

biasgrp = mean(abs(betasgrp[,1] - true))
covgrp = mean(abs(betasgrp[,1] - true) <= 1.96 * segrp)


s11grp = abs(mean(surv1grp[,1]) - true_survival_probabilities[1])
s12grp = abs(mean(surv1grp[,2]) - true_survival_probabilities[2])
s13grp = abs(mean(surv1grp[,3]) - true_survival_probabilities[3])
s14grp = abs(mean(surv1grp[,4]) - true_survival_probabilities[4])
s15grp = abs(mean(surv1grp[,5]) - true_survival_probabilities[5])

grps1 = rbind(s11grp,s12grp,s13grp,s14grp,s15grp)
rownames(grps1) = c('1','2','3','4','5')
colnames(grps1) = 'GRP Survival X = 1'

s01grp = abs(mean(surv0grp[,1]) - true_survival_probabilities0[1])
s02grp = abs(mean(surv0grp[,2]) - true_survival_probabilities0[2])
s03grp = abs(mean(surv0grp[,3]) - true_survival_probabilities0[3])
s04grp = abs(mean(surv0grp[,4]) - true_survival_probabilities0[4])
s05grp = abs(mean(surv0grp[,5]) - true_survival_probabilities0[5])

grps0 = rbind(s01grp,s02grp,s03grp,s04grp,s05grp)
rownames(grps0) = c('1','2','3','4','5')
colnames(grps0) = 'GRP Survival X = 0'

h11grp = abs(mean(cumhaz1grp[,1]) - true_hazard_rates[1])
h12grp = abs(mean(cumhaz1grp[,2]) - true_hazard_rates[2])
h13grp = abs(mean(cumhaz1grp[,3]) - true_hazard_rates[3])
h14grp = abs(mean(cumhaz1grp[,4]) - true_hazard_rates[4])
h15grp = abs(mean(cumhaz1grp[,5]) - true_hazard_rates[5])

grph1 = rbind(h11grp,h12grp,h13grp,h14grp,h15grp)
rownames(grph1) = c('1','2','3','4','5')
colnames(grph1) = 'GRP Cumhaz X = 1'

h01grp = abs(mean(cumhaz0grp[,1]) - true_hazard_rates0[1])
h02grp = abs(mean(cumhaz0grp[,2]) - true_hazard_rates0[2])
h03grp = abs(mean(cumhaz0grp[,3]) - true_hazard_rates0[3])
h04grp = abs(mean(cumhaz0grp[,4]) - true_hazard_rates0[4])
h05grp = abs(mean(cumhaz0grp[,5]) - true_hazard_rates0[5])

grph0 = rbind(h01grp,h02grp,h03grp,h04grp,h05grp)
rownames(grph0) = c('1','2','3','4','5')
colnames(grph0) = 'GRP Cumhaz X = 0'

cs11grp = coverage(true_survival_probabilities[1], grpLs1[,1], grpUs1[,1])
cs12grp = coverage(true_survival_probabilities[2], grpLs1[,2], grpUs1[,2])
cs13grp = coverage(true_survival_probabilities[3], grpLs1[,3], grpUs1[,3])
cs14grp = coverage(true_survival_probabilities[4], grpLs1[,4], grpUs1[,4])
cs15grp = coverage(true_survival_probabilities[5], grpLs1[,5], grpUs1[,5])

grpcs1 = rbind(cs11grp,cs12grp,cs13grp,cs14grp,cs15grp)
rownames(grpcs1) <- c('1','2','3','4','5')
colnames(grpcs1) <- 'GRP Survival Coverage X = 1'

cs01grp = coverage(true_survival_probabilities0[1], grpLs0[,1], grpUs0[,1])
cs02grp = coverage(true_survival_probabilities0[2], grpLs0[,2], grpUs0[,2])
cs03grp = coverage(true_survival_probabilities0[3], grpLs0[,3], grpUs0[,3])
cs04grp = coverage(true_survival_probabilities0[4], grpLs0[,4], grpUs0[,4])
cs05grp = coverage(true_survival_probabilities0[5], grpLs0[,5], grpUs0[,5])

grpcs0 = rbind(cs01grp,cs02grp,cs03grp,cs04grp,cs05grp)
rownames(grpcs0) <- c('1','2','3','4','5')
colnames(grpcs0) <- 'GRP Survival Coverage X = 0'

ch11grp = coverage(true_hazard_rates[1], grpLh1[,1], grpUh1[,1])
ch12grp = coverage(true_hazard_rates[2], grpLh1[,2], grpUh1[,2])
ch13grp = coverage(true_hazard_rates[3], grpLh1[,3], grpUh1[,3])
ch14grp = coverage(true_hazard_rates[4], grpLh1[,4], grpUh1[,4])
ch15grp = coverage(true_hazard_rates[5], grpLh1[,5], grpUh1[,5])

grpch1 = rbind(ch11grp,ch12grp,ch13grp,ch14grp,ch15grp)
rownames(grpch1) <- c('1','2','3','4','5')
colnames(grpch1) <- 'GRP Cumhaz Coverage X = 1'

ch01grp = coverage(true_hazard_rates0[1], grpLh0[,1], grpUh0[,1])
ch02grp = coverage(true_hazard_rates0[2], grpLh0[,2], grpUh0[,2])
ch03grp = coverage(true_hazard_rates0[3], grpLh0[,3], grpUh0[,3])
ch04grp = coverage(true_hazard_rates0[4], grpLh0[,4], grpUh0[,4])
ch05grp = coverage(true_hazard_rates0[5], grpLh0[,5], grpUh0[,5])

grpch0 = rbind(ch01grp,ch02grp,ch03grp,ch04grp,ch05grp)
rownames(grpch0) <- c('1','2','3','4','5')
colnames(grpch0) <- 'GRP Cumhaz Coverage X = 0'


# G Cox -----------------------------------

biasgc = mean(abs(betasgc[,1] - true))
covgc = mean(abs(betasgc[,1] - true) <= 1.96 * segc)


s11gc = abs(mean(surv1gc[,1]) - true_survival_probabilities[1])
s12gc = abs(mean(surv1gc[,2]) - true_survival_probabilities[2])
s13gc = abs(mean(surv1gc[,3]) - true_survival_probabilities[3])
s14gc = abs(mean(surv1gc[,4]) - true_survival_probabilities[4])
s15gc = abs(mean(surv1gc[,5]) - true_survival_probabilities[5])

gcs1 = rbind(s11gc,s12gc,s13gc,s14gc,s15gc)
rownames(gcs1) <- c('1','2','3','4','5')
colnames(gcs1) <- 'GC Survival X = 1'

s01gc = abs(mean(surv0gc[,1]) - true_survival_probabilities0[1])
s02gc = abs(mean(surv0gc[,2]) - true_survival_probabilities0[2])
s03gc = abs(mean(surv0gc[,3]) - true_survival_probabilities0[3])
s04gc = abs(mean(surv0gc[,4]) - true_survival_probabilities0[4])
s05gc = abs(mean(surv0gc[,5]) - true_survival_probabilities0[5])

gcs0 = rbind(s01gc,s02gc,s03gc,s04gc,s05gc)
rownames(gcs0) <- c('1','2','3','4','5')
colnames(gcs0) <- 'GC Survival X = 0'

h11gc = abs(mean(cumhaz1gc[,1]) - true_hazard_rates[1])
h12gc = abs(mean(cumhaz1gc[,2]) - true_hazard_rates[2])
h13gc = abs(mean(cumhaz1gc[,3]) - true_hazard_rates[3])
h14gc = abs(mean(cumhaz1gc[,4]) - true_hazard_rates[4])
h15gc = abs(mean(cumhaz1gc[,5]) - true_hazard_rates[5])

gch1 = rbind(h11gc,h12gc,h13gc,h14gc,h15gc)
rownames(gch1) <- c('1','2','3','4','5')
colnames(gch1) <- 'GC Cumhaz X = 1'

h01gc = abs(mean(cumhaz0gc[,1]) - true_hazard_rates0[1])
h02gc = abs(mean(cumhaz0gc[,2]) - true_hazard_rates0[2])
h03gc = abs(mean(cumhaz0gc[,3]) - true_hazard_rates0[3])
h04gc = abs(mean(cumhaz0gc[,4]) - true_hazard_rates0[4])
h05gc = abs(mean(cumhaz0gc[,5]) - true_hazard_rates0[5])

gch0 = rbind(h01gc,h02gc,h03gc,h04gc,h05gc)
rownames(gch0) <- c('1','2','3','4','5')
colnames(gch0) <- 'GC Cumhaz X = 0'

cs11gc = coverage(true_survival_probabilities[1], gcLs1[,1], gcUs1[,1])
cs12gc = coverage(true_survival_probabilities[2], gcLs1[,2], gcUs1[,2])
cs13gc = coverage(true_survival_probabilities[3], gcLs1[,3], gcUs1[,3])
cs14gc = coverage(true_survival_probabilities[4], gcLs1[,4], gcUs1[,4])
cs15gc = coverage(true_survival_probabilities[5], gcLs1[,5], gcUs1[,5])

gccs1 = rbind(cs11gc,cs12gc,cs13gc,cs14gc,cs15gc)
rownames(gccs1) <- c('1','2','3','4','5')
colnames(gccs1) <- 'GC Survival Coverage X = 1'

cs01gc = coverage(true_survival_probabilities0[1], gcLs0[,1], gcUs0[,1])
cs02gc = coverage(true_survival_probabilities0[2], gcLs0[,2], gcUs0[,2])
cs03gc = coverage(true_survival_probabilities0[3], gcLs0[,3], gcUs0[,3])
cs04gc = coverage(true_survival_probabilities0[4], gcLs0[,4], gcUs0[,4])
cs05gc = coverage(true_survival_probabilities0[5], gcLs0[,5], gcUs0[,5])

gccs0 = rbind(cs01gc,cs02gc,cs03gc,cs04gc,cs05gc)
rownames(gccs0) <- c('1','2','3','4','5')
colnames(gccs0) <- 'GC Survival Coverage X = 0'

ch11gc = coverage(true_hazard_rates[1], gcLh1[,1], gcUh1[,1])
ch12gc = coverage(true_hazard_rates[2], gcLh1[,2], gcUh1[,2])
ch13gc = coverage(true_hazard_rates[3], gcLh1[,3], gcUh1[,3])
ch14gc = coverage(true_hazard_rates[4], gcLh1[,4], gcUh1[,4])
ch15gc = coverage(true_hazard_rates[5], gcLh1[,5], gcUh1[,5])

gcch1 = rbind(ch11gc,ch12gc,ch13gc,ch14gc,ch15gc)
rownames(gcch1) <- c('1','2','3','4','5')
colnames(gcch1) <- 'GC Cumhaz Coverage X = 1'

ch01gc = coverage(true_hazard_rates0[1], gcLh0[,1], gcUh0[,1])
ch02gc = coverage(true_hazard_rates0[2], gcLh0[,2], gcUh0[,2])
ch03gc = coverage(true_hazard_rates0[3], gcLh0[,3], gcUh0[,3])
ch04gc = coverage(true_hazard_rates0[4], gcLh0[,4], gcUh0[,4])
ch05gc = coverage(true_hazard_rates0[5], gcLh0[,5], gcUh0[,5])

gcch0 = rbind(ch01gc,ch02gc,ch03gc,ch04gc,ch05gc)
rownames(gcch0) <- c('1','2','3','4','5')
colnames(gcch0) <- 'GC Cumhaz Coverage X = 0'

# --------------
betabias = cbind(biasrp, biasc, biasgrp, biasgc) 
betacov = cbind(covrp, covc, covgrp, covgc)

survival1 = cbind(rps1, cs1, grps1, gcs1)
survival0 = cbind(rps0, cs0, grps0, gcs0)
hazard1 = cbind(rph1, ch1, grph1, gch1)
hazard0 = cbind(rph0, ch0, grph0, gch0)
coveage_survival1 = cbind(rpcs1, ccs1, grpcs1, gccs1)
coverage_survival0 = cbind(rpcs0, ccs0, grpcs0, gccs0)
coverage_hazard1 = cbind(rpch1, cch1, grpch1, gcch1)
coverage_hazard0 = cbind(rpch0, cch0, grpch0, gcch0)

options(scipen=999)

stopCluster(cl)
