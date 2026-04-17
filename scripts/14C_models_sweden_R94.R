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
R94_df=subset(all_data,all_data$LTE=="R94-1957")
R94_df=subset(all_data,all_data$LTE=="R94-1957")
R94_df=subset(all_data,all_data$LTE=="R94-1957")
R94_df=subset(all_data,all_data$LTE=="R94-1966")

## C stocks
# Calculate C stocks (g m-2) = C content (g kg-1) x BD (g cm-3) x depth (0.2 m) x 1000
R94_df$C_stocks_gR94<-R94_df$Cgkg*1.51*0.2*1000 #BD = 1.51
R94_df$C_stocks_gR94<-R94_df$Cgkg*1.72*0.2*1000 #BD = 1.72
R94_df$C_stocks_gR94<-R94_df$Cgkg*1.44*0.2*1000 #BD = 1.44
R94_df$C_stocks_gR94<-R94_df$Cgkg*1.37*0.2*1000 #BD = 1.37

## Inputs
# Constant inputs (kg C m-2) from Bolinder et al. 2012, doi: 10.4141/cjss2012-036 Table 4
BolinderIn<-c(0.277,0.268,0.304,0.271,0.269,0.29, # Rotation A
              0.259, 0.291, 0.296, 0.288, 0.265, 0.256) # Rotation B

mean_C_inputs<-mean(BolinderIn)*1000 #transform to a single mean input in g m-2

## Atmospheric radiocarbon
Atm14C=Hua2021$NHZone1[,1:2]
fAtm14C=read.csv(here("csv_files", "NHZ1forecast.csv")) #forecast later than 2019
Atm14C<-rbind(Atm14C, data.frame(Year=fAtm14C$time, mean.Delta14C=Delta14C_from_AbsoluteFractionModern(fAtm14C$F14C)))


## Models: Site R94

### Obs dfs for each SOM thermal pool
R94_bulk=subset(R94_df,R94_df$Temperature...C.=="Soil")
R94_Cobs_bulk <- data.frame(Year = R94_bulk$Year, Ct = R94_bulk$C_stocks_gR94)
R94_C14obs_bulk <- data.frame(Year = R94_bulk$Year, C14t = R94_bulk$d14C)

R94_325=subset(R94_df,R94_df$Temperature...C.=="325")
R94_Cobs_325 <- data.frame(Year = R94_325$Year, Ct_fast = R94_325$C_stocks_gR94)
R94_C14obs_325 <- data.frame(Year = R94_325$Year, C14t_fast = R94_325$d14C)

R94_400=subset(R94_df,R94_df$Temperature...C.=="400")
R94_Cobs_400 <- data.frame(Year = R94_400$Year, Ct_slow = R94_400$C_stocks_gR94)
R94_C14obs_400 <- data.frame(Year = R94_400$Year, C14t_slow = R94_400$d14C)

### Initial C & Delta14C
C0_R94_bulk<-R94_Cobs_bulk[1,2]
F0_R94_bulk<-R94_C14obs_bulk[1,2]
F0_R94_400<--120

### Define function to run a two-pool series model
#pars[1:5] = kf, ks, alpha sf, C0fb, F0fb, I

mf=function(pars){
  md=TwopSeriesModel14(t=yr,ks=pars[1:2],C0=C0_R94_bulk*c(pars[4], 1-pars[4]), #par[4] allocates a portion of initial bulk C to fast pool
                       F0_Delta14C = c(F0_R94_bulk * pars[5], F0_R94_400), #par[5] allocates a portion of initial Delta14C of bulk to fast pool 
                       In=pars[6], #constant input scalar 
                       a21=pars[1]*pars[3], inputFc = Atm14C) #where pars[1] = kf and pars[3] = alpha sf
  Ct_pools = getC(md) # matrix with 2 columns for fast and slow pool C
  C14_pools = getF14(md) # matrix with 2 columns for fast and slow pool Delta 14C
  C14t = getF14C(md) #bulk 14C    
  return(data.frame(
    Year = yr,
    Ct = rowSums(Ct_pools), # Total C
    C14t = C14t, #Bulk Delta 14C     
    Ct_fast = Ct_pools[,1],    
    Ct_slow = Ct_pools[,2],
    C14t_fast = C14_pools[,1],
    C14t_slow = C14_pools[,2]
  ))
}

### Define function to add obs datasets to cost function

mc=function(pars){ 
  out=mf(pars) #out = df of Year, Ct and C14t from mf output 
  Cost1 = modCost(model = out, obs = R94_Cobs_bulk, # C bulk
                  x = "Year")
  Cost2 = modCost(model = out, obs = R94_C14obs_bulk, # Delta14C bulk
                  x = "Year", cost = Cost1)
  Cost3 = modCost(model = out, obs = R94_Cobs_325, # C fast pool
                  x = "Year", cost = Cost2)
  Cost4 = modCost(model = out, obs = R94_C14obs_325, # Delta14C fast pool 
                  x = "Year", cost = Cost3)
  Cost5 = modCost(model = out, obs = R94_Cobs_400, # C slow pool
                  x = "Year", cost = Cost4)
  return(modCost(model = out, obs = R94_C14obs_400, # Delta 14C slow pool
                 x = "Year", cost = Cost5))
} 

inipars=c(0.5,0.001,0.01, 0.4, 0.9, mean_C_inputs*0.1) #pars[1:6] = kf, ks, alpha sf, C0fb, F0fb,I

yr <- as.numeric(seq(1957,2019, by = 1/12))

### Run model

## Uncomment the following to run again
mFit_R94=modFit(f=mc,p=inipars,method="Nelder-Mead",upper=c(2,0.5,1,1,1, mean_C_inputs*1.2),lower=c(0,0,0,0,0,0))
bestpars_R94=mFit_R94$par
# save(bestpars_R94, file="bestpars_R94.RData")
## Otherwise load previously saved bestpars.RData
bestModel_R94<-TwopSeriesModel14(t=yr,ks=bestpars_R94[1:2],C0=C0_R94_bulk*c(bestpars_R94[4], 1-bestpars_R94[4]), 
                                 F0_Delta14C = c(F0_R94_bulk * bestpars_R94[5], F0_R94_400),
                                 In=bestpars_R94[6], 
                                 a21=bestpars_R94[1]*bestpars_R94[3], inputFc = Atm14C) 

mod_C14t_R94_pools=data.frame(Year = yr, 
                              "C14t_fast"=getF14(bestModel_R94)[,1], 
                              "C14t_slow" =getF14(bestModel_R94)[,2]) # Delta14C for both pools
mod_C14t_R94_bulk=data.frame(Year = yr, 
                             "C14t" = getF14C(bestModel_R94)) # Delta14C for bulk
mod_Ct_R94_pools=data.frame(Year = yr,
                            "Ct_fast" = getC(bestModel_R94)[,1],
                            "Ct_slow" = getC(bestModel_R94)[,2]) # C for both pools
mod_Ct_R94_bulk=data.frame(Year = yr, 
                           "Ct_bulk" = rowSums(getC(bestModel_R94))) # C for bulk


# plot modeled vs observed
Delta14Clabel<-expression(Delta^14*C)
stocks_label <- expression(C ~ (g ~ m^{-2}))

## C stocks plot
y_min_C <- min(R94_Cobs_bulk$Ct, R94_Cobs_325$Ct_fast, R94_Cobs_400$Ct_slow, 
               mod_Ct_R94_bulk$Ct_bulk, mod_Ct_R94_pools$Ct_fast, mod_Ct_R94_pools$Ct_slow, na.rm = TRUE)
y_max_C <- max(R94_Cobs_bulk$Ct, R94_Cobs_325$Ct_fast, R94_Cobs_400$Ct_slow, 
               mod_Ct_R94_bulk$Ct_bulk, mod_Ct_R94_pools$Ct_fast, mod_Ct_R94_pools$Ct_slow, na.rm = TRUE)

par(mar = c(5, 5, 4, 2))  # bottom, left, top, right
plot(mod_Ct_R94_bulk$Year, mod_Ct_R94_bulk$Ct, #mod bulk C 
     type = "l", lwd = 2,
     xlab = "Year", ylab = stocks_label,
     main = "Modeled vs Observed C stocks",
     ylim = c(y_min_C, y_max_C))
lines(mod_Ct_R94_pools$Year, mod_Ct_R94_pools$Ct_fast, #mod fast pool C 
      col = "blue", lty = 1, lwd = 2)
lines(mod_Ct_R94_pools$Year, mod_Ct_R94_pools$Ct_slow, #mod slow pool C 
      col = "darkgreen", lty = 1, lwd = 2)
points(R94_Cobs_bulk$Year, R94_Cobs_bulk$Ct, #observed bulk C points
       pch = 16, col = "black")
points(R94_Cobs_325$Year, R94_Cobs_325$Ct_fast, #observed fast pool C points
       pch = 16, col = "blue")
points(R94_Cobs_400$Year, R94_Cobs_400$Ct_slow, #observed slow pool C points
       pch = 16, col = "darkgreen")
legend("topright", 
       legend = c("Bulk", "Fast Pool", "Slow Pool"),
       col = c("black", "blue", "darkgreen"),
       lty = c(1, 1, 1, 1), 
       pch = c(NA, NA, NA, NA))

##Delta14C plot
y_min_14C <- min(R94_C14obs_bulk$C14t, R94_C14obs_325$C14t_fast, R94_C14obs_400$C14t_slow, 
                 mod_C14t_R94_bulk$C14t, mod_C14t_R94_pools$C14t_fast, mod_C14t_R94_pools$C14t_slow, na.rm = TRUE)
y_max_14C <- max(R94_C14obs_bulk$C14t, R94_C14obs_325$C14t_fast, R94_C14obs_400$C14t_slow, 
                 R94_C14obs_bulk$C14t, R94_C14obs_325$C14t_fast, R94_C14obs_400$C14t_slow, 
                 mod_C14t_R94_bulk$C14t, mod_C14t_R94_pools$C14t_fast, mod_C14t_R94_pools$C14t_slow, na.rm = TRUE)

par(mar = c(5, 5, 4, 2))  # bottom, left, top, right
plot(mod_C14t_R94_bulk$Year, mod_C14t_R94_bulk$C14t, # mod bulk Delta14C 
     type = "l", lwd = 2,
     xlab = "Year", ylab = Delta14Clabel,
     main = "Modeled vs Observed Delta14C",
     ylim = c(y_min_14C, y_max_14C))
lines(mod_C14t_R94_pools$Year, mod_C14t_R94_pools$C14t_fast, # mod fast pool Delta14C
      col = "blue", lty = 1, lwd = 2)
lines(mod_C14t_R94_pools$Year, mod_C14t_R94_pools$C14t_slow, # mod slow pool Delta14C 
      col = "darkgreen", lty = 1, lwd = 2)
lines(Atm14C$Year, Atm14C$mean.Delta14C, #atm 14C 
      col = "purple", lty = 1, lwd = 2)
points(R94_C14obs_bulk$Year, R94_C14obs_bulk$C14t, #observed bulk Delta14C points
       pch = 16, col = "black")
points(R94_C14obs_325$Year, R94_C14obs_325$C14t_fast, #observed fast pool Delta14C points
       pch = 16, col = "blue")
points(R94_C14obs_400$Year, R94_C14obs_400$C14t_slow, #observed slow pool Delta14C points
       pch = 16, col = "darkgreen")
legend("topright", 
       legend = c("Bulk", "Fast Pool", "Slow Pool", "Atm"),
       col = c("black", "blue", "darkgreen", "purple"),
       lty = c(1, 1, 1, 1), 
       pch = c(NA, NA, NA, NA))
