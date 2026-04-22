# SOM decomposition models for agricultural soils at multiple sites

## load libraries
library(here)
library(SoilR)
library(FME)
library(ggplot2)
library(dplyr)

## import data
all_data <- read.csv2(here("csv_files", "LTEnitrogen2.csv"))
load("~/nel_sweden_agsoil_2026/mod_runs/results_325asslow.Rdata")

## data cleaning
all_data <- all_data %>%
  mutate(across(c("Year", "TONgkg", "Cgkg", "MolarCN", "F14C", "err", "d14C"), as.numeric)) %>%
  mutate(across(c("LTE", "Teatment", "Temperature...C."), as.factor)) 

## site settings -------------------------------------------------------------

sites <- c("M2-1957", "M4-1957", "M6-1957", "R94-1966")

bulk_density <- mean(c(1.51, 1.72, 1.44, 1.37))
F0_site_400_vals <- mean(all_data$d14C)

## Inputs --------------------------------------------------------------------

BolinderIn <- c(
  0.277,0.268,0.304,0.271,0.269,0.29,
  0.259,0.291,0.296,0.288,0.265,0.256
)

mean_C_inputs <- mean(BolinderIn) * 1000  # g m-2

## Atmospheric radiocarbon ---------------------------------------------------

Atm14C <- Hua2021$NHZone1[,1:2]

fAtm14C <- read.csv(here("csv_files", "NHZ1forecast.csv"))

Atm14C <- rbind(
  Atm14C,
  data.frame(
    Year = fAtm14C$time,
    mean.Delta14C = Delta14C_from_AbsoluteFractionModern(fAtm14C$F14C)
  )
)

# add C stocks
all_data$C_stocks_gm2 <- all_data$Cgkg * bulk_density * 0.2 * 1000

summary_df <- all_data %>%
  group_by(Year, Temperature...C.) %>%
  summarise(
    C_stocks_mean = mean(C_stocks_gm2, na.rm = TRUE),
    C_stocks_sd = sd(C_stocks_gm2, na.rm = TRUE),
    d14C_mean     = mean(d14C, na.rm = TRUE),
    d14C_sd       = sd(d14C, na.rm = TRUE),
    .groups = "drop"
  )

bulk <- subset(summary_df, Temperature...C. == "Soil")
pool_325 <- subset(summary_df, Temperature...C. == "325")
pool_400<- subset(summary_df, Temperature...C. == "400")

# obs data for cost func
Cobs_bulk <- data.frame(Year = bulk$Year, Ct = bulk$C_stocks_mean, Ct_sd = bulk$C_stocks_sd)
C14obs_bulk <- data.frame(Year = bulk$Year, C14t = bulk$d14C_mean, C14t_sd = bulk$d14C_sd)
Cobs_325 <- data.frame(Year = pool_325$Year, Ct_325 = pool_325$C_stocks_mean, Ct_325_sd = pool_325$C_stocks_sd)
Cobs_slow <- data.frame(Year = pool_400$Year, Ct_slow = pool_400$C_stocks_mean, Ct_slow_sd = pool_400$C_stocks_sd)
C14obs_slow <- data.frame(Year = pool_400$Year, C14t_slow = pool_400$d14C_mean, C14t_slow_sd = pool_400$d14C_sd)
Cobs_fast <- data.frame(Year = Cobs_325$Year, Ct_fast = Cobs_bulk$Ct-Cobs_325$Ct_325, Ct_fast_sd = sqrt(bulk$C_stocks_sd^2 + pool_325$C_stocks_sd^2)) # C fast  = C bulk - C 325
Cobs_inter<-data.frame(Year = Cobs_325$Year, Ct_inter = Cobs_325$Ct_325-Cobs_slow$Ct_slow, Ct_inter_sd = sqrt(Cobs_325$Ct_325_sd^2 + Cobs_slow$Ct_slow_sd^2))

#remove NA sd values
Cobs_bulk$Ct_sd[is.na(Cobs_bulk$Ct_sd)] <- mean(Cobs_bulk$Ct_sd, na.rm = TRUE) 
C14obs_bulk$C14t_sd[is.na(C14obs_bulk$C14t_sd)] <- mean(C14obs_bulk$C14t_sd, na.rm = TRUE)
Cobs_fast$Ct_fast_sd[is.na(Cobs_fast$Ct_fast_sd)] <- mean(Cobs_fast$Ct_fast_sd, na.rm = TRUE)
C14obs_slow$C14t_slow_sd[is.na(C14obs_slow$C14t_slow_sd)] <- mean(C14obs_slow$C14t_slow_sd, na.rm = TRUE)
Cobs_inter$Ct_inter_sd[is.na(Cobs_inter$Ct_inter_sd)] <- mean(Cobs_inter$Ct_inter_sd, na.rm = TRUE)
Cobs_slow$Ct_slow_sd[is.na(Cobs_slow$Ct_slow_sd)] <- mean(Cobs_slow$Ct_slow_sd, na.rm = TRUE)

#initial values
yr <- seq(1957, 2019, by = 1/12)
C0_bulk <- mean(Cobs_bulk[Cobs_bulk$Year==1957,]$Ct)
#C0_fast <- mean(Cobs_fast[Cobs_bulk$Year==1957,]$Ct_fast)
C0_fast <- 5100
C0_400 <- mean(Cobs_slow[Cobs_slow$Year==1957,]$Ct_slow)
F0_400 <- mean(C14obs_slow[C14obs_slow$Year==1957,]$C14t_slow)
F0_bulk <- mean(C14obs_bulk[C14obs_bulk$Year==1957,]$C14t) 

# func to run mod
run_mod <- function(pars){ #kf, ki, ks, alpha 21, +F0b=Fi, I, alpha 32, alpha 41
  c14_atm <- BoundFc(Atm14C, format = "Delta14C")
  c14_initial <- ConstFc(
    values = c(F0_bulk + 10, F0_bulk + pars[5], F0_400, F0_bulk + 10), 
    format = "Delta14C"
  )
  A4 <- diag(-c(pars[1:3],0.000000001))
  A4[2,1] <- pars[1]*pars[4]
  A4[3,2] <- pars[2]*pars[7]
  A4[4,1] <- pars[1]*pars[8]
  mod<-GeneralModel_14(
    t = yr,
    A = A4,
    ivList = c(C0_fast, abs(C0_bulk-C0_fast-C0_400), C0_400, 0),
    initialValF = c14_initial,
    inputFluxes = c(pars[6],0,0,0),
    inputFc = c14_atm
  )
  
  Ct_pools <- getC(mod)
  C14_pools <- getF14(mod)
  C14t <- getF14C(mod)
  
  mod_results<-data.frame(
    Year = yr,
    Ct = rowSums(Ct_pools),
    C14t = C14t,
    Ct_fast = Ct_pools[,1],
    Ct_inter = Ct_pools[,2],
    Ct_slow = Ct_pools[,3],
    Ct_loss = Ct_pools[,4],
    C14t_fast = C14_pools[,1],
    C14t_inter = C14_pools[,2],
    C14t_slow = C14_pools[,3],
    C14t_loss = C14_pools[,4]
  )
  
  return(mod_results)
}

inipars <- c(0.1,0.05,0.001, 0.1, -50, mean_C_inputs*0.5, 0.1, 0.005) #kf, ki, ks, alpha 21, F0ib, I, alpha 32, alpha 41

# cost func
mc <- function(pars){
  out = run_mod(pars)
  Cost1 <- modCost(out, Cobs_bulk, x = "Year", err = "Ct_sd")
  Cost2 <- modCost(out, C14obs_bulk, x = "Year", cost = Cost1, err = "C14t_sd")
  Cost3 <- modCost(out, Cobs_fast, x = "Year", cost = Cost2, err = "Ct_fast_sd")
  Cost4 <- modCost(out, Cobs_inter, x = "Year", cost = Cost3, err = "Ct_inter_sd")
  modCost(out, Cobs_slow, x = "Year", cost = Cost4, err = "Ct_slow_sd")
}

#fit mod
mFit <- modFit(
  f = mc,
  p = inipars,
  method = "Nelder-Mead",
  upper = c(0.5,0.2,0.005, 0.5,0, mean_C_inputs*2,1, 0.01), #kf, ki, ks, alpha 21, F0ib, I, alpha 32, alpha 41
  lower = c(0.05,0.01,0.0001,0,-160,0,0, 0) 
)

bestpars <- mFit$par
out_best <- run_mod(bestpars)
Delta14Clabel <- expression(Delta^14*C)
stocks_label <- expression(C ~ (g ~ m^{-2}))

col_bulk <- "black"
col_fast <- "#1b9e77"
col_slow <- "#d95f02"
col_inter <- "#0000FF"
col_loss <- "#FFEA00"

## plots
plot_C <- ggplot() +
  
  ## --- ribbons FIRST (so they stay behind lines) ---
  geom_ribbon(data = Cobs_bulk,
              aes(Year, ymin = Ct - Ct_sd, ymax = Ct + Ct_sd, fill = "Bulk"),
              alpha = 0.15) +
  
  geom_ribbon(data = Cobs_fast,
              aes(Year, ymin = Ct_fast - Ct_fast_sd, ymax = Ct_fast + Ct_fast_sd, fill = "Fast"),
              alpha = 0.15) +
  
  geom_ribbon(data = Cobs_slow,
              aes(Year, ymin = Ct_slow - Ct_slow_sd, ymax = Ct_slow + Ct_slow_sd, fill = "Slow"),
              alpha = 0.15) +
  
  geom_ribbon(data = Cobs_inter,
              aes(Year, ymin = Ct_inter - Ct_inter_sd, ymax = Ct_inter + Ct_inter_sd, fill = "Inter"),
              alpha = 0.15) +
  
  
  geom_line(data = out_best, aes(Year, Ct, colour = "Bulk"), linewidth = 1) +
  geom_line(data = out_best, aes(Year, Ct_fast, colour = "Fast"), linetype = "dashed") +
  geom_line(data = out_best, aes(Year, Ct_slow, colour = "Slow"), linetype = "dotted") +
  geom_line(data = out_best, aes(Year, Ct_inter, colour = "Inter"), linetype = "dotted") +
  geom_line(data = out_best, aes(Year, Ct_loss, colour = "Loss"), linetype = "dotted") +
  
  geom_point(data = Cobs_bulk, aes(Year, Ct, colour = "Bulk")) +
  geom_point(data = Cobs_fast, aes(Year, Ct_fast, colour = "Fast")) +
  geom_point(data = Cobs_slow, aes(Year, Ct_slow, colour = "Slow")) +
  geom_point(data = Cobs_inter, aes(Year, Ct_inter, colour = "Inter")) +
  
  ## --- scales ---
  scale_colour_manual(
    name = "Pool",
    values = c("Bulk" = col_bulk,
               "Fast" = col_fast,
               "Slow" = col_slow,
               "Inter" = col_inter,
               "Loss" = col_loss)
  ) +
  
  scale_fill_manual(
    name = "Pool",
    values = c("Bulk" = col_bulk,
               "Fast" = col_fast,
               "Slow" = col_slow,
               "Inter" = col_inter)
  ) +
  
  ## --- labels ---
  ggtitle("C stocks") +
  ylab(stocks_label) +
  xlab("Year") +
  theme_minimal()

plot_C14 <- ggplot() +
  
  ## --- ribbons ---
  geom_ribbon(data = C14obs_bulk,
              aes(Year, ymin = C14t - C14t_sd, ymax = C14t + C14t_sd, fill = "Bulk"),
              alpha = 0.15) +
  
  geom_ribbon(data = C14obs_slow,
              aes(Year, ymin = C14t_slow - C14t_slow_sd, ymax = C14t_slow + C14t_slow_sd, fill = "Slow"),
              alpha = 0.15) +
  
  ## --- model lines ---
  geom_line(data = out_best, aes(Year, C14t, colour = "Bulk"), linewidth = 1) +
  geom_line(data = out_best, aes(Year, C14t_fast, colour = "Fast"), linetype = "dashed") +
  geom_line(data = out_best, aes(Year, C14t_slow, colour = "Slow"), linetype = "dotted") +
  geom_line(data = out_best, aes(Year, C14t_inter, colour = "Inter"), linetype = "dotted") +
  geom_line(data = out_best, aes(Year, C14t_loss, colour = "Loss"), linetype = "dotted") +
  
  geom_point(data = C14obs_bulk, aes(Year, C14t, colour = "Bulk")) +
  geom_point(data = C14obs_slow, aes(Year, C14t_slow, colour = "Slow")) +
  
  ## --- scales ---
  scale_colour_manual(
    name = "Pool",
    values = c("Bulk" = col_bulk,
               "Fast" = col_fast,
               "Slow" = col_slow,
               "Inter" = col_inter,
               "Loss" = col_loss)
  ) +
  
  scale_fill_manual(
    name = "Pool",
    values = c("Bulk" = col_bulk,
               "Fast" = col_fast,
               "Slow" = col_slow,
               "Inter" = col_inter)
  ) +
  
  ## --- labels ---
  ggtitle(expression(Delta^14*C)) +
  ylab(Delta14Clabel) +
  xlab("Year") +
  theme_minimal()

# save results and plots
output_dir <- "plots/allsites/4_pool"

file_C <- file.path(output_dir, "C.png")
file_C14 <- file.path(output_dir, "C14.png")

ggsave(filename = file_C, plot = plot_C,
       width = 7, height = 5, dpi = 300, bg = "white")
ggsave(filename = file_C14, plot = plot_C14,
       width = 7, height = 5, dpi = 300, bg = "white")

saveRDS(plot_C,
        file = file.path(output_dir, "C.rds"))
saveRDS(plot_C14,
        file = file.path(output_dir, "14C.rds"))

results_4p<-list(
  pars = bestpars,
  output = out_best,
  plot_C = plot_C,
  plot_C14 = plot_C14
)

save(results_4p, file = file.path("mod_runs", "results_allsites_4pool.Rdata"))
#load(here::here("mod_runs/results_allsites_4pool.Rdata"))

## MCMC

var0 <- mFit$var_ms_unweighted
cov0 <- summary(mFit)$cov.scaled #  cov matrix can be used for jump 
# uncomment to run again 
MCMC_M2 <- modMCMC(f=mc, p = bestpars_M2, niter = 2500, jump = cov0, var0 = var0, wvar0 = 1) 
save(MCMC_M2, file="MCMC_M2.RData")
# alternatively load saved run "MCMC_M2.RData"
parsMCMC_M2<-summary(MCMC_M2)
# uncomment to map uncertainties of predicted Ct and C14t
sR_M2=sensRange(func=mf, parInput=MCMC_M2$par) 
save(sR_M2, file="sR_M2.RData") 
# alternatively, load saved run
# load("sR_M2.RData")
summarysR_M2<-summary(sR_M2) 
nx2_M2<-attributes(summarysR_M2)$nx
CtR2<-as.data.frame(summarysR_M2[1:nx2_M2,]) # mean is the Ct (g m-2) predicted
C14tR2<-as.data.frame(summarysR_M2[(nx2_M2+1):(2*nx2_M2),])#mean is the Delta_C14t (per mille) predicted

pairs(MCMC_M2,nsample=500) #check inter-dependence of parameters using 500 MCMC runs

plot(CtR2[,1:2], type="l", xlab="Year", ylab="C stock", ylim=c(0,6000)) 
polygon(c(CtR2$x,rev(CtR2$x)), c(CtR2$Min, rev(CtR2$Max)) ,col=gray(0.8), border="NA")
polygon(c(CtR2$x,rev(CtR2$x)), c(CtR2$Mean+CtR2$Sd, rev(CtR2$Mean-CtR2$Sd)) ,col=gray(0.5), border="NA")
lines(CtR2[,1:2])
points(Cobs,pch=20)

plot(C14tR2[,1:2], type="l", xlim=c(1950,2020), ylim=c(-120,400), ylab=expression(paste(Delta^14,"C ","(\u2030)")),xlab="Years")
polygon(c(C14tR2$x,rev(C14tR2$x)), c(C14tR2$Min, rev(C14tR2$Max)) ,col=gray(0.8), border="NA")
polygon(c(C14tR2$x,rev(C14tR2$x)), c(C14tR2$Mean+C14tR2$Sd, rev(C14tR2$Mean-C14tR2$Sd)) ,col=gray(0.5), border="NA")
lines(C14tR2[,1:2])
points(C14obs, pch=20)
lines(Atm14C,col=4)

meanpars2<-as.numeric(parsMCMC2[1,1:5])
meanModel2<-TwopSeriesModel14(t=yr,ks=meanpars2[1:2],C0=C0*c(meanpars[4], 1-meanpars[4]),
                              F0_Delta14C = c(F0*meanpars[5], -200), 
                              #In=mean(BolinderIn)*1000, 
                              In=mean(avgYld[,2], na.rm = TRUE)*1.5,
                              a21=meanpars[1]*meanpars[3], inputFc = Atm14C)
meanF14C2=getF14(meanModel2)
meanC2=getC(meanModel2)



# Age and transit time
tau=seq(0,500)
A1=meanModel@mat@map(1970) #tnel: matrix created using arbitrary year, since it remains the same for all years
A2=meanModel2@mat@map(1970)

#tnel: testing scenarios of min and max turnover rates
A1min<- -1*diag(parsMCMC[3,1:2]) #tnel: constructs a diag matrix of decomp rates using min p1 and p2
A1min[2,1]<- abs(parsMCMC[3,3])*parsMCMC[3,2] #tnel: alpha21 was neg.ve - even though we provided restrictions in SoilR model, the MCMC optimization gives a mean from distribution
A1max<- -1*diag(parsMCMC[4,1:2]) 
A1max[2,1]<- parsMCMC[4,3]*parsMCMC[4,2] 
A2min<- -1*diag(parsMCMC2[3,1:2])
A2min[2,1]<- parsMCMC2[3,3]*parsMCMC2[3,2] 
A2max<- -1*diag(parsMCMC2[4,1:2])
A2max[2,1]<- parsMCMC2[4,3]*parsMCMC2[4,2] 


SA1=systemAge(A=A1,u=c(1,0),a=tau) #here u are normalized inputs i.e. simply divided among pools according to system e.g. 0.5, 0.5 could be parallel
TT1=transitTime(A=A1,u=c(1,0), a=tau, q=c(0.1,0.5,0.9))
TT1min<-transitTime(A=A1min,u=c(1,0),q=c(0.1,0.5,0.9))
TT1max<-transitTime(A=A1max,u=c(1,0),q=c(0.1,0.5,0.9))

SA2=systemAge(A=A2,u=c(1,0),a=tau)
TT2=transitTime(A=A2,u=c(1,0), a=tau, q=c(0.1,0.5,0.9))
TT2min<-transitTime(A=A2min,u=c(1,0),q=c(0.1,0.5,0.9))
TT2max<-transitTime(A=A2max,u=c(1,0),q=c(0.1,0.5,0.9))

In=data.frame(avgYld[,1],avgYld[,2]*bestpars[4]) #tnel: In makes a df of variable inputs, where bestpars[4] = gamma yield to inputs

pdf("C:/Users/tnel/OneDrive - Universiteit Antwerpen/Teneille UAntwerp/Literature/carbon persistence/spohn 2023/code spohn 2023/Figures/Fig2.pdf", encoding = 'WinAnsi.enc', width=7.5*sqrt(2), height=7.5)
par(mfrow=c(2,2), mar=c(4,4.5,1,1), cex.lab=1.2)
#plot(avgYld$Year, avgYld$Yield*2,type="o", xlab="Calendar year", ylab=expression(paste("Mean yields [g DW ", m^-2, " y", r^-1, "]" )), xlim=c(1937,2022))
plot(In, type="b", ylim=c(0,700), xlab="Calendar year", ylab=expression(paste("Carbon inputs [g C ", m^-2, " y", r^-1,"]")), col=rgb(0,0,1), pch=20)
#abline(h=mean(BolinderIn)*1000, col=2)
abline(h=mean(avgYld[,2], na.rm = TRUE)*1.5, col=rgb(0,1,0))
legend("bottomright", c("Variable inputs predicted by the model", "Constant inputs from Bolinder et al. (2012)"), lty=1, 
       col=c(rgb(0,0,1), rgb(0,1,0)), bty="n")
legend("topleft", "a", cex=1.2, text.font=2, bty="n")

par(mar=c(4,4.5,1,1))
plot(C14tR[,1:2],type="l",ylim=c(-200,400), xlab="Calendar year", ylab=expression(paste(Delta^14,"C [\U2030]")), xlim=c(1950,2022))
polygon(c(C14tR$x,rev(C14tR$x)), c(C14tR$Min, rev(C14tR$Max)) ,col=rgb(0,0,1, alpha=0.3), border="NA") #tnel: part of plotting sd
polygon(c(C14tR$x,rev(C14tR$x)), c(C14tR$Mean+C14tR$Sd, rev(C14tR$Mean-C14tR$Sd)) ,col=rgb(0,0,1,alpha=0.5), border="NA")
lines(C14tR[,1:2], col=rgb(0,0,1))
polygon(c(C14tR2$x,rev(C14tR2$x)), c(C14tR2$Min, rev(C14tR2$Max)) ,col=rgb(0,1,0, alpha=0.3), border="NA")
polygon(c(C14tR2$x,rev(C14tR2$x)), c(C14tR2$Mean+C14tR2$Sd, rev(C14tR2$Mean-C14tR2$Sd)) ,col=rgb(0,1,0, alpha=0.5), border="NA")
lines(C14tR2[,1:2], col=rgb(0,1,0))
points(C14obs, pch=20)
lines(Atm14C)
legend("topright", c("Measurements", "Model with constant inputs", "Model with variable inputs", "Atmosphere"), 
       lty=c(NA,1,1,1), pch=c(19, NA, NA, NA), col=c(1,rgb(0,1,0), rgb(0,0,1),1), bty="n")
legend("topleft", "b", cex=1.2, text.font=2, bty="n")

matplot(yr,C1,type="l",col=rgb(0,0,1),lty=2:3, lwd=2,ylim=c(0,10000), xlim=c(1966,2022), xlab="Calendar year",
        ylab=expression(paste("Carbon stock [g C ", m^-2, "]")))
polygon(c(CtR$x,rev(CtR$x)), c(CtR$Min, rev(CtR$Max)) ,col=rgb(0,0,1, alpha=0.3), border="NA")
polygon(c(CtR$x,rev(CtR$x)), c(CtR$Mean+CtR$Sd, rev(CtR$Mean-CtR$Sd)) ,col=rgb(0,0,1, alpha=0.5), border="NA")
lines(CtR[,1:2], col=rgb(0,0,1))
matlines(yr, C2, col=rgb(0,1,0), lty=2:3, lwd=2)
polygon(c(CtR2$x,rev(CtR2$x)), c(CtR2$Min, rev(CtR2$Max)) ,col=rgb(0,1,0, alpha=0.3), border="NA")
polygon(c(CtR2$x,rev(CtR2$x)), c(CtR2$Mean+CtR2$Sd, rev(CtR2$Mean-CtR2$Sd)) ,col=rgb(0,1,0, alpha=0.5), border="NA")
lines(CtR2[,1:2], col=rgb(0,1,0))
points(Cobs, pch=20)
legend("top", c("Measurements", "TOC, constant inputs", "Fast pool", "Slow pool"), lty=c(NA,1,2,3), 
       pch=c(19, NA, NA, NA), col=c(1,rep(rgb(0,1,0),3)), bty="n")
legend("topright", c("TOC, variable inputs", "Fast pool", "Slow pool"), lty=c(1,2,3), 
       col=c(rep(rgb(0,0,1),3)), bty="n")
legend("topleft", "c", cex=1.2, text.font=2, bty="n")

plot(tau,log(TT1$transitTimeDensity), type="l", xlim=c(0,200), xlab="Transit time (yr)", ylab="Log probability density", col=rgb(0,0,1))
abline(v=TT1$meanTransitTime, lty=2, col=rgb(0,0,1))
#abline(v=TT1$quantiles[2],lty=3, col=4)
lines(tau, log(TT2$transitTimeDensity), col=rgb(0,1,0))
abline(v=TT2$meanTransitTime, lty=2, col=rgb(0,1,0))
#abline(v=TT2$quantiles[2],lty=3, col=2)
legend(x=60,y=-4,c("Transit time distribution, variable inputs", paste("Mean transit time = ", round(TT1$meanTransitTime,1), " yr"))
       ,lty=1:3, col=rgb(0,0,1), bty="n")
legend(x=60,y=-6,c("Transit time distribution, constant inputs", paste("Mean transit time = ", round(TT2$meanTransitTime,1), " yr"))
       ,lty=1:3, col=rgb(0,1,0), bty="n")
legend(x=5,y=-1, "d", cex=1.2, text.font=2, bty="n")
par(mfrow=c(1,1))
dev.off()


data.frame(Parameter=c("kf", "ks", "alpha_sf", "gamma", "beta"), 
           Model1=paste(round(bestpars[1:5], 3), "+-", round(parsMCMC[2,1:5], 3)), #tnel: note here the best pars, not the mean pars, are reported with the MCMC uncertainty
           Model2=paste(c(round(bestpars2[1:3], 3), NA, round(bestpars2[4], 3)), "+-", c(round(parsMCMC2[2,1:3],3), NA, round(parsMCMC2[2,4], 3)))
)

paste("Mean TT (uncertainty): ",round(TT1$meanTransitTime, 2), "(", round(TT1max$meanTransitTime, 2), "--", round(TT1min$meanTransitTime, 2),")")
paste("Mean TT (uncertainty): ",round(TT2$meanTransitTime, 2), "(", round(TT2max$meanTransitTime, 2), "--", round(TT2min$meanTransitTime, 2),")")

paste("Median TT (uncertainty): ",round(TT1$quantiles[2], 2), "(", round(TT1max$quantiles[2], 2), "--", round(TT1min$quantiles[2], 2),")")
paste("Median TT (uncertainty): ",round(TT2$quantiles[2], 2), "(", round(TT2max$quantiles[2], 2), "--", round(TT2min$quantiles[2], 2),")")

round(TT1$quantiles, 2) #tnel: just making extraction easy for table
round(TT1min$meanTransitTime, 2)
round(TT1min$quantiles, 2)
round(TT1max$meanTransitTime, 2)
round(TT1max$quantiles, 2)

round(TT2$meanTransitTime, 2)
round(TT2$quantiles, 2)

