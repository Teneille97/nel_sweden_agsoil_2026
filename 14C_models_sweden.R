#SOM decomposition models for agricultural soils at four sites in Sweden

## load libraries
library(here)
library(SoilR)
library(FME)
library(ggplot2)
library(dplyr)

## import data
all_data=read.csv2(here("csv_files", "LTEnitrogen2.csv"))

all_data <- all_data %>%
  mutate(across(c("Year", "TONgkg", "Cgkg", "MolarCN", "F14C", "err", "d14C"), as.numeric)) %>%
  mutate(across(c("LTE", "Teatment", "Temperature...C."), as.factor)) 

## data subsets by sites
M2_df=subset(all_data,all_data$LTE=="M2-1957")
M4_df=subset(all_data,all_data$LTE=="M4-1957")
M6_df=subset(all_data,all_data$LTE=="M6-1957")
R94_df=subset(all_data,all_data$LTE=="R94-1966")

## C stocks
# Calculate C stocks (g m-2) = C content (g kg-1) x BD (g cm-3) x depth (0.2 m) x 1000
M2_df$C_stocks_gm2<-M2_df$Cgkg*1.51*0.2*1000 #BD = 1.51
M4_df$C_stocks_gm2<-M4_df$Cgkg*1.72*0.2*1000 #BD = 1.72
M6_df$C_stocks_gm2<-M6_df$Cgkg*1.44*0.2*1000 #BD = 1.44
R94_df$C_stocks_gm2<-R94_df$Cgkg*1.37*0.2*1000 #BD = 1.37

## Inputs
# Constant inputs (kg C m-2) from Bolinder et al. 2012, doi: 10.4141/cjss2012-036 Table 4
BolinderIn<-c(0.277,0.268,0.304,0.271,0.269,0.29, # Rotation A
              0.259, 0.291, 0.296, 0.288, 0.265, 0.256) # Rotation B

mean_C_inputs<-mean(BolinderIn)*1000 #transform to a single mean input in g m-2

## Atmospheric radiocarbon
Atm14C=Hua2021$NHZone1[,1:2]
fAtm14C=read.csv(here("csv_files", "NHZ1forecast.csv")) #forecast later than 2019
Atm14C<-rbind(Atm14C, data.frame(Year=fAtm14C$time, mean.Delta14C=Delta14C_from_AbsoluteFractionModern(fAtm14C$F14C)))

## Nested mod: First using 325 pool as 'bulk' and 400 pool as slow, 
## the mod "fast" pool is used in subsequent three-pool mod, as "obs" inter pool

## Models: Site M2

### Obs dfs for each SOM thermal pool
M2_325=subset(M2_df,M2_df$Temperature...C.=="325")
M2_Cobs_325 <- data.frame(Year = M2_325$Year, Ct_325 = M2_325$C_stocks_gm2)
M2_C14obs_325 <- data.frame(Year = M2_325$Year, C14t_325 = M2_325$d14C)

M2_400=subset(M2_df,M2_df$Temperature...C.=="400")
M2_Cobs_400 <- data.frame(Year = M2_400$Year, Ct_slow = M2_400$C_stocks_gm2)
M2_C14obs_400 <- data.frame(Year = M2_400$Year, C14t_slow = M2_400$d14C)

M2_Cobs_inter <- data.frame(Year = M2_400$Year, Ct_inter = M2_325$C_stocks_gm2-M2_400$C_stocks_gm2) #C inter = C 325 - C 400

### Initial C & Delta14C
C0_M2_slow<-M2_Cobs_400[1,2]
C0_M2_inter<-M2_Cobs_inter[1,2]
C0_M2_325<-M2_Cobs_325[1,2]

F0_M2_400<-M2_C14obs_400[1,2]
F0_M2_325<-M2_C14obs_325[1,2]


### STEP 1. Define function to run a two-pool series model, using 325 pool as 'bulk' and 400 pool as slow
#pars[1:6] = kf, ks, alpha sf, I

mf1=function(pars){
  md=TwopSeriesModel14(t=yr,ks=pars[1:2],C0=c(C0_M2_inter, C0_M2_slow), 
                       F0_Delta14C = c(-100, -140),  #set initial inter pool C based on preliminary investigation
                       In=pars[4], # mod inputs  
                       a21=pars[1]*pars[3], inputFc = Atm14C) #where pars[1] = kf and pars[3] = alpha sf
  Ct_pools = getC(md) # matrix with 2 columns for fast and slow pool C
  C14_pools = getF14(md) # matrix with 2 columns for fast and slow pool Delta 14C
  C14t = getF14C(md) #bulk 14C    
  return(data.frame(
    Year = yr,
    Ct_325 = rowSums(Ct_pools), # total C
    C14t_325 = C14t, # bulk Delta 14C     
    Ct_inter = Ct_pools[,1],    
    Ct_slow = Ct_pools[,2],
    C14t_inter = C14_pools[,1],
    C14t_slow = C14_pools[,2]
  ))
}

### Define function to add obs datasets to cost function

mc1=function(pars){ 
  out=mf1(pars) #out = df of Year, Ct and C14t from mf output 
  Cost1 = modCost(model = out, obs = M2_Cobs_325, # C 325
                  x = "Year")
  Cost2 = modCost(model = out, obs = M2_C14obs_325, # Delta14C 325
                  x = "Year", cost = Cost1)
  Cost3 = modCost(model = out, obs = M2_Cobs_400, # C slow pool
                  x = "Year", cost = Cost2)
  Cost4 = modCost(model = out, obs = M2_C14obs_400, # Delta14C slow pool
                  x = "Year", cost = Cost3)
  Cost5 = modCost(model = out, obs = M2_Cobs_inter, # C "fast" pool
                 x = "Year", cost = Cost4)
  return(Cost5)
} 

inipars1=c(0.05,0.001,0.01, 2) #pars[1:5] = kf, ks, alpha sf, I

yr <- as.numeric(seq(1957,2019, by = 1/12))

### Run model


## Uncomment the following to run again
 mFit1_M2=modFit(f=mc1,p=inipars1,method="Nelder-Mead",upper=c(1,0.5,1,mean_C_inputs),lower=c(0,0,0,0))
 bestpars1_M2=mFit1_M2$par
# save(bestpars1_M2, file="bestpars1_M2.RData")
## Otherwise load previous results
# load("bestpars1.RData")
bestModel1_M2<-TwopSeriesModel14(t=yr,ks=bestpars1_M2[1:2],C0=c(C0_M2_inter, C0_M2_slow), 
                  F0_Delta14C = c(-100, -140),
                  In=bestpars1_M2[4], 
                  a21=bestpars1_M2[1]*bestpars1_M2[3], inputFc = Atm14C) 


mod1_C14t_MF2_pools=data.frame(Year = yr, 
                          "C14t_inter"=getF14(bestModel1_M2)[,1], 
                          "C14t_slow" =getF14(bestModel1_M2)[,2]) # Delta14C for both pools
mod1_C14t_MF2_325=data.frame(Year = yr, 
                         "C14t_325" = getF14C(bestModel1_M2)) # Delta14C for bulk
mod1_Ct_MF2_pools=data.frame(Year = yr,
                       "Ct_inter" = getC(bestModel1_M2)[,1],
                       "Ct_slow" = getC(bestModel1_M2)[,2]) # C for both pools
mod1_Ct_MF2_325=data.frame(Year = yr, 
                      "Ct_325" = rowSums(getC(bestModel1_M2))) # C for bulk


# plot modeled vs observed
Delta14Clabel<-expression(Delta^14*C)
stocks_label <- expression(C ~ (g ~ m^{-2}))

## C stocks plot
y_min_C1 <- min(M2_Cobs_325$Ct, M2_Cobs_325$Ct_325, M2_Cobs_400$Ct_slow, 
               mod1_Ct_MF2_325$Ct_325, mod1_Ct_MF2_pools$Ct_inter, mod1_Ct_MF2_pools$Ct_slow, na.rm = TRUE)
y_max_C1 <- max(M2_Cobs_325$Ct, M2_Cobs_325$Ct_325, M2_Cobs_400$Ct_slow, 
               mod1_Ct_MF2_325$Ct_325, mod1_Ct_MF2_pools$Ct_inter, mod1_Ct_MF2_pools$Ct_slow, na.rm = TRUE)

par(mar = c(5, 5, 4, 2))  # bottom, left, top, right
plot(mod1_Ct_MF2_325$Year, mod1_Ct_MF2_325$Ct_325, #mod 325 C 
     type = "l", lwd = 2,
     xlab = "Year", ylab = stocks_label,
     main = "Modeled vs Observed C stocks",
     ylim = c(y_min_C1, y_max_C1))
lines(mod1_Ct_MF2_pools$Year, mod1_Ct_MF2_pools$Ct_inter, #mod fast pool C 
      col = "blue", lty = 1, lwd = 2)
lines(mod1_Ct_MF2_pools$Year, mod1_Ct_MF2_pools$Ct_slow, #mod slow pool C 
      col = "darkgreen", lty = 1, lwd = 2)
points(M2_Cobs_325$Year, M2_Cobs_325$Ct_325, #observed 325 C points
       pch = 16, col = "black")
points(M2_Cobs_400$Year, M2_Cobs_400$Ct_slow, #observed slow pool C points
       pch = 16, col = "darkgreen")
legend("topright", 
       legend = c("325 pool", "intermediate pool", "slow pool"),
       col = c("black", "blue", "darkgreen"),
       lty = c(1, 1, 1, 1), 
       pch = c(NA, NA, NA, NA))

##Delta14C plot
y_min_14C1 <- min(M2_C14obs_325$C14t, M2_C14obs_325$C14t_325, M2_C14obs_400$C14t_slow, 
                 mod1_C14t_MF2_325$C14t_325, mod1_C14t_MF2_pools$C14t_inter, mod1_C14t_MF2_pools$C14t_slow, na.rm = TRUE)
y_max_14C1 <- max(M2_C14obs_325$C14t, M2_C14obs_325$C14t_325, M2_C14obs_400$C14t_slow, 
                 mod1_C14t_MF2_325$C14t_325, mod1_C14t_MF2_pools$C14t_inter, mod1_C14t_MF2_pools$C14t_slow, na.rm = TRUE)

par(mar = c(5, 5, 4, 2))  # bottom, left, top, right
plot(mod1_C14t_MF2_325$Year, mod1_C14t_MF2_325$C14t_325, # mod bulk Delta14C 
     type = "l", lwd = 2,
     xlab = "Year", ylab = Delta14Clabel,
     main = "Modeled vs Observed Delta14C",
     ylim = c(y_min_14C1, y_max_14C1))
lines(mod1_C14t_MF2_pools$Year, mod1_C14t_MF2_pools$C14t_inter, # mod fast pool Delta14C
      col = "blue", lty = 1, lwd = 2)
lines(mod1_C14t_MF2_pools$Year, mod1_C14t_MF2_pools$C14t_slow, # mod slow pool Delta14C 
      col = "darkgreen", lty = 1, lwd = 2)
lines(Atm14C$Year, Atm14C$mean.Delta14C, #atm 14C 
      col = "purple", lty = 1, lwd = 2)
points(M2_C14obs_325$Year, M2_C14obs_325$C14t, #observed bulk Delta14C points
       pch = 16, col = "black")
points(M2_C14obs_400$Year, M2_C14obs_400$C14t_slow, #observed slow pool Delta14C points
       pch = 16, col = "darkgreen")
legend("topleft", 
       legend = c("325 pool", "intermediate pool", "slow pool", "Atm"),
       col = c("black", "blue", "darkgreen", "purple"),
       lty = c(1, 1, 1, 1), 
       pch = c(NA, NA, NA, NA))



### STEP 2. Define function to run a three-pool series model, using mod "fast" pool (from step 1) as "obs" inter pool

### Obs dfs for each SOM thermal pool
M2_bulk=subset(M2_df,M2_df$Temperature...C.=="Soil")
M2_Cobs_bulk <- data.frame(Year = M2_bulk$Year, Ct = M2_bulk$C_stocks_gm2)
M2_C14obs_bulk <- data.frame(Year = M2_bulk$Year, C14t = M2_bulk$d14C)

M2_C14obs_inter <- data.frame(Year = mod1_C14t_MF2_pools$Year, C14t_inter = mod1_C14t_MF2_pools$C14t_inter)

M2_Cobs_fast <- data.frame(Year = M2_400$Year, Ct_fast = M2_bulk$C_stocks_gm2-M2_400$C_stocks_gm2) # fast C = bulk C - slow C

### Initial C & Delta14C
C0_M2_fast<-M2_Cobs_fast[1,2]
C0_M2_slow<-M2_Cobs_400[1,2]
C0_M2_inter<-M2_Cobs_inter[1,2]

F0_M2_bulk<-M2_C14obs_bulk[1,2]
F0_M2_400<-M2_C14obs_400[1,2]


#pars[1:6] = kf, ki, ks, alpha 21, alpha 32, Ffb

mf2=function(pars){
  md=ThreepSeriesModel14(t=yr,ks=pars[1:3],C0=c(C0_M2_fast, C0_M2_inter, C0_M2_slow), 
                       F0_Delta14C = c(pars[6]*F0_M2_bulk, -100, -140),  #par[6] assigns proportion of bulk 14C to fast pool, initial inter & slow pool C based on preliminary investigation
                       In=mean_C_inputs, # mod inputs  
                       a21=pars[1]*pars[4], a32= pars[2]*pars[5], inputFc = Atm14C) 
  Ct_pools = getC(md) # matrix with 3 columns for C of pools
  C14_pools = getF14(md) # matrix with 3 columns for Delta 14C of pools
  C14t = getF14C(md) #bulk 14C    
  return(data.frame(
    Year = yr,
    Ct = rowSums(Ct_pools), # total C
    C14t = C14t, # bulk Delta 14C     
    Ct_fast = Ct_pools[,1], 
    Ct_inter = Ct_pools[,2],    
    Ct_slow = Ct_pools[,3],
    C14t_fast = C14_pools[,1],
    C14t_inter = C14_pools[,2],
    C14t_slow = C14_pools[,3]
  ))
}


### Define function to add obs datasets to cost function

mc2=function(pars){ 
  out=mf2(pars) #out = df of Year, Ct and C14t from mf output 
  Cost1 = modCost(model = out, obs = M2_Cobs_bulk, # C bulk
                  x = "Year")
  Cost2 = modCost(model = out, obs = M2_C14obs_bulk, # Delta14C bulk
                  x = "Year", cost = Cost1)
  Cost3 = modCost(model = out, obs = M2_Cobs_400, # C slow pool
                  x = "Year", cost = Cost2)
  Cost4 = modCost(model = out, obs = M2_C14obs_400, # Delta14C slow pool
                  x = "Year", cost = Cost3)
  Cost5 = modCost(model = out, obs = M2_Cobs_inter, # C inter pool
                  x = "Year", cost = Cost4)
  Cost6 = modCost(model = out, obs = M2_C14obs_inter, # Delta14C inter pool
                  x = "Year", cost = Cost5)
  Cost7 = modCost(model = out, obs = M2_Cobs_fast, # C fast pool
                  x = "Year", cost = Cost5)
  return(Cost7)
} 

inipars2=c(0.1, 0.05,0.001,0.01, 0.01, 0.9) #pars[1:6] = kf, ki, ks, alpha 21, alpha 32, Ffb

yr <- as.numeric(seq(1957,2019, by = 1/12))

### Run model
## Uncomment the following to run again
mFit2_M2=modFit(f=mc2,p=inipars2,method="Nelder-Mead",upper=c(1,0.5,0.1,0.9,0.9,0.9),lower=c(0,0,0,0,0,0))
bestpars2_M2=mFit2_M2$par
# save(bestpars2_M2, file="bestpars2_M2.RData")
## Otherwise load previous results
# load("bestpars2.RData")
bestModel2_M2<-ThreepSeriesModel14(t=yr,ks=bestpars2_M2[1:3],C0=c(C0_M2_fast, C0_M2_inter, C0_M2_slow), 
                                   F0_Delta14C = c(bestpars2_M2[6]*F0_M2_bulk, -100, -140), 
                                   In=mean_C_inputs, # mod inputs  
                                   a21=bestpars2_M2[1]*bestpars2_M2[4], a32= bestpars2_M2[2]*bestpars2_M2[5], inputFc = Atm14C) 

mod2_C14t_MF2_pools=data.frame(Year = yr, 
                               "C14t_fast"=getF14(bestModel2_M2)[,1], 
                               "C14t_inter" =getF14(bestModel2_M2)[,2],
                               "C14t_slow" =getF14(bestModel2_M2)[,2]) # Delta14C for pools
mod2_C14t_MF2_bulk=data.frame(Year = yr, 
                             "C14t_bulk" = getF14C(bestModel2_M2)) # Delta14C for bulk
mod2_Ct_MF2_pools=data.frame(Year = yr,
                             "Ct_fast" = getC(bestModel2_M2)[,1],
                             "Ct_inter" = getC(bestModel2_M2)[,2],
                             "Ct_slow" = getC(bestModel2_M2)[,3]) # C for pools
mod2_Ct_MF2_bulk=data.frame(Year = yr, 
                           "Ct_bulk" = rowSums(getC(bestModel2_M2))) # C for bulk

# plot modeled vs observed
Delta14Clabel<-expression(Delta^14*C)
stocks_label <- expression(C ~ (g ~ m^{-2}))

## C stocks plot
y_min_C2 <- min(M2_Cobs_bulk$Ct, M2_Cobs_fast$Ct_fast, M2_Cobs_inter$Ct_inter, M2_Cobs_400$Ct_slow, 
                mod2_Ct_MF2_bulk$Ct_bulk, 
                mod2_Ct_MF2_pools$Ct_fast, mod2_Ct_MF2_pools$Ct_inter, mod2_Ct_MF2_pools$Ct_slow, na.rm = TRUE)
y_max_C2 <- max(M2_Cobs_bulk$Ct, M2_Cobs_fast$Ct_fast, M2_Cobs_inter$Ct_inter, M2_Cobs_400$Ct_slow, 
                mod2_Ct_MF2_bulk$Ct_bulk, 
                mod2_Ct_MF2_pools$Ct_fast, mod2_Ct_MF2_pools$Ct_inter, mod2_Ct_MF2_pools$Ct_slow, na.rm = TRUE)

par(mar = c(5, 5, 4, 2))  # bottom, left, top, right
plot(mod2_Ct_MF2_325$Year, mod2_Ct_MF2_325$Ct_bulk, #mod bulk C 
     type = "l", lwd = 2,
     xlab = "Year", ylab = stocks_label,
     main = "Modeled vs Observed C stocks",
     ylim = c(y_min_C2, y_max_C2))
lines(mod2_Ct_MF2_pools$Year, mod2_Ct_MF2_pools$Ct_fast, #mod fast pool C 
      col = "blue", lty = 1, lwd = 2)
lines(mod2_Ct_MF2_pools$Year, mod2_Ct_MF2_pools$Ct_slow, #mod slow pool C 
      col = "darkgreen", lty = 1, lwd = 2)
lines(mod2_Ct_MF2_pools$Year, mod2_Ct_MF2_pools$Ct_inter, #mod inter pool C 
      col = "orange", lty = 1, lwd = 2)
points(M2_Cobs_bulk$Year, M2_Cobs_bulk$Ct, #observed bulk C points
       pch = 16, col = "black")
points(M2_Cobs_400$Year, M2_Cobs_400$Ct_slow, #observed slow pool C points
       pch = 16, col = "darkgreen")
points(M2_Cobs_inter$Year, M2_Cobs_inter$Ct_inter, #observed inter pool C points
       pch = 16, col = "orange")
points(M2_Cobs_fast$Year, M2_Cobs_fast$Ct_fast, #observed fast pool C points
       pch = 16, col = "blue")
legend("topright", 
       legend = c("bulk pool", "fast pool", "slow pool", "intermediate pool"),
       col = c("black", "blue", "darkgreen", "orange"),
       lty = c(1, 1, 1, 1), 
       pch = c(NA, NA, NA, NA))


##Delta14C plot
y_min_14C2 <- min(M2_C14obs_bulk$C14t, M2_C14obs_inter$C14t_inter, M2_C14obs_400$C14t_slow, 
                mod2_C14t_MF2_bulk$C14t_bulk, 
                mod2_C14t_MF2_pools$C14t_fast, mod2_C14t_MF2_pools$C14t_inter, mod2_C14t_MF2_pools$C14t_slow, na.rm = TRUE)
y_max_14C2 <- max(M2_C14obs_bulk$C14t, M2_C14obs_inter$C14t_inter, M2_C14obs_400$C14t_slow, 
                mod2_C14t_MF2_bulk$C14t_bulk, 
                mod2_C14t_MF2_pools$C14t_fast, mod2_C14t_MF2_pools$C14t_inter, mod2_C14t_MF2_pools$C14t_slow, na.rm = TRUE)


par(mar = c(5, 5, 4, 2))  # bottom, left, top, right
plot(mod2_C14t_MF2_bulk$Year, mod2_C14t_MF2_bulk$C14t_bulk, #mod bulk C 
     type = "l", lwd = 2,
     xlab = "Year", ylab = Delta14Clabel,
     main = "Modeled vs Observed Delta14C",
     ylim = c(y_min_14C2, y_max_14C2))
lines(mod2_C14t_MF2_pools$Year, mod2_C14t_MF2_pools$C14t_fast, #mod fast pool 14C 
      col = "blue", lty = 1, lwd = 2)
lines(mod2_C14t_MF2_pools$Year, mod2_C14t_MF2_pools$C14t_slow, #mod slow pool 14C 
      col = "darkgreen", lty = 1, lwd = 2)
lines(mod2_C14t_MF2_pools$Year, mod2_C14t_MF2_pools$C14t_inter, #mod inter pool 14C 
      col = "orange", lty = 1, lwd = 2)
lines(Atm14C$Year, Atm14C$mean.Delta14C, #atm 14C 
      col = "purple", lty = 1, lwd = 2)
points(M2_C14obs_bulk$Year, M2_C14obs_bulk$C14t, #observed bulk 14C points
       pch = 16, col = "black")
points(M2_C14obs_400$Year, M2_C14obs_400$C14t_slow, #observed slow pool 14C points
       pch = 16, col = "darkgreen")
legend("topright", 
       legend = c("bulk pool", "fast pool", "slow pool", "intermediate pool", "atm"),
       col = c("black", "blue", "darkgreen", "orange", "purple"),
       lty = c(1, 1, 1, 1, 1), 
       pch = c(NA, NA, NA, NA, NA))


### MCMC optimization
var0_M2 <- mFit_M2$var_ms_unweighted
cov0_M2 <- summary(mFit_M2)$cov.scaled # The covariance matrix can be used for the jump 
# uncomment to run again 
# MCMC_M2 <- modMCMC(f=mc, p = bestpars_M2, niter = 2500, jump = NULL, var0 = var0, wvar0 = 1) #tnel: jump default is 10% of par value; what is wvar and why is it zero?
# save(MCMC_M2, file="MCMC_M2.RData")
# alternatively load saved run
# load("MCMC_M2.RData")
parsMCMC_M2<-summary(MCMC_M2)
# uncomment to map uncertainties of predicted Ct and C14t
# sR_M2=sensRange(func=mf, parInput=MCMC_M2$par) 
# save(sR_M2, file="sR_M2.RData") 
# alternatively, load saved run
# load("sR_M2.RData")
summarysR_M2<-summary(sR_M2) # mean is the Ct (g m-2) predicted
nx_M2<-attributes(summarysR_M2)$nx
C14tR_M2<-as.data.frame(summarysR_M2[(nx_M2+1):(2*nx_M2),]) #mean is the Delta_C14t (per mille) predicted

pairs(MCMC,nsample=500) #check inter-dependence of parameters using 500 MCMC runs;