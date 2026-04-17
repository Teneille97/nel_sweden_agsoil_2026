# SOM decomposition models for agricultural soils at multiple sites

## load libraries
library(here)
library(SoilR)
library(FME)
library(ggplot2)
library(dplyr)

## import data
all_data <- read.csv2(here("csv_files", "LTEnitrogen2.csv"))

all_data <- all_data %>%
  mutate(across(c("Year", "TONgkg", "Cgkg", "MolarCN", "F14C", "err", "d14C"), as.numeric)) %>%
  mutate(across(c("LTE", "Teatment", "Temperature...C."), as.factor)) 

## site settings -------------------------------------------------------------

sites <- c("M2-1957", "M4-1957", "M6-1957", "R94-1966")

bulk_density <- c(
  "M2-1957" = 1.51,
  "M4-1957" = 1.72,
  "M6-1957" = 1.44,
  "R94-1966" = 1.37
)

C0_site_400_vals <- c(
  "M2-1957" = -150,
  "M4-1957" = -500,
  "M6-1957" = -250,
  "R94-1966" = -120
)

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
run_model_site <- function(site_name) {
  
  cat("Running site:", site_name, "\n")
  
  # subset data
  site_df <- subset(all_data, LTE == site_name)
  
  # add C stocks
  bd <- bulk_density[site_name]
  site_df$C_stocks_gm2 <- site_df$Cgkg * bd * 0.2 * 1000
  
  ## split thermal pools -----------------------------------------------------
  
  bulk <- subset(site_df, Temperature...C. == "Soil")
  pool_325 <- subset(site_df, Temperature...C. == "325")
  pool_400 <- subset(site_df, Temperature...C. == "400")
  
  ## observed data -----------------------------------------------------------
  
  Cobs_bulk <- data.frame(Year = bulk$Year, Ct = bulk$C_stocks_gm2)
  C14obs_bulk <- data.frame(Year = bulk$Year, C14t = bulk$d14C)
  
  Cobs_325 <- data.frame(Year = pool_325$Year, Ct_fast = pool_325$C_stocks_gm2)
  C14obs_325 <- data.frame(Year = pool_325$Year, C14t_fast = pool_325$d14C)
  
  Cobs_400 <- data.frame(Year = pool_400$Year, Ct_slow = pool_400$C_stocks_gm2)
  C14obs_400 <- data.frame(Year = pool_400$Year, C14t_slow = pool_400$d14C)
  
  ## time + initial parameters ----------------------------------------------
  
  yr <- seq(1957, 2019, by = 1/12)
  inipars <- c(0.5,0.001,0.01, 0.4, 0.9, mean_C_inputs*0.1)
  C0_bulk <- Cobs_bulk[1,2]
  F0_400 <- C0_site_400_vals[site_name]
  F0_bulk <- C14obs_bulk[1,2]
  
  ## model function ----------------------------------------------------------
  
  mf <- function(pars){
    
    md <- TwopSeriesModel14(
      t = yr,
      ks = pars[1:2],
      C0 = C0_bulk * c(pars[4], 1 - pars[4]),
      F0_Delta14C = c(F0_bulk * pars[5], F0_400),
      In = pars[6],
      a21 = pars[1] * pars[3],
      inputFc = Atm14C
    )
    
    Ct_pools <- getC(md)
    C14_pools <- getF14(md)
    C14t <- getF14C(md)
    
    data.frame(
      Year = yr,
      Ct = rowSums(Ct_pools),
      C14t = C14t,
      Ct_fast = Ct_pools[,1],
      Ct_slow = Ct_pools[,2],
      C14t_fast = C14_pools[,1],
      C14t_slow = C14_pools[,2]
    )
  }
  
  ## cost function -----------------------------------------------------------
  
  mc <- function(pars){
    
    out <- mf(pars)
    
    Cost1 <- modCost(out, Cobs_bulk, x = "Year")
    Cost2 <- modCost(out, C14obs_bulk, x = "Year", cost = Cost1)
    Cost3 <- modCost(out, Cobs_325, x = "Year", cost = Cost2)
    Cost4 <- modCost(out, C14obs_325, x = "Year", cost = Cost3)
    Cost5 <- modCost(out, Cobs_400, x = "Year", cost = Cost4)
    
    modCost(out, C14obs_400, x = "Year", cost = Cost5)
  }
  
  ## optimize ---------------------------------------------------------------
  
  mFit <- modFit(
    f = mc,
    p = inipars,
    method = "Nelder-Mead",
    upper = c(2,0.5,1,1,1, mean_C_inputs*1.2),
    lower = c(0,0,0,0,0,0)
  )
  
  bestpars <- mFit$par
  out_best <- mf(bestpars)
  
  ## ------------------------------------------------------------------------
  ## PLOTS
  ## ------------------------------------------------------------------------
  
  ## labels
  Delta14Clabel <- expression(Delta^14*C)
  stocks_label <- expression(C ~ (g ~ m^{-2}))
  
  ## define colours
  col_bulk <- "black"
  col_fast <- "#1b9e77"
  col_slow <- "#d95f02"
  
  ## Carbon stocks ------------------------------------------------------------
  
  plot_C <- ggplot() +
    
    # model lines
    geom_line(data = out_best, aes(Year, Ct, colour = "Bulk"), linewidth = 1) +
    geom_line(data = out_best, aes(Year, Ct_fast, colour = "Fast"), linetype = "dashed") +
    geom_line(data = out_best, aes(Year, Ct_slow, colour = "Slow"), linetype = "dotted") +
    
    # observations
    geom_point(data = Cobs_bulk, aes(Year, Ct, colour = "Bulk")) +
    geom_point(data = Cobs_325, aes(Year, Ct_fast, colour = "Fast")) +
    geom_point(data = Cobs_400, aes(Year, Ct_slow, colour = "Slow")) +
    
    scale_colour_manual(
      name = "Pool",
      values = c("Bulk" = col_bulk,
                 "Fast" = col_fast,
                 "Slow" = col_slow)
    ) +
    
    ggtitle(paste("C stocks -", site_name)) +
    ylab(stocks_label) +
    xlab("Year") +
    theme_minimal()
  
  
  ## Delta 14C ---------------------------------------------------------------
  
  plot_C14 <- ggplot() +
    
    # model lines
    geom_line(data = out_best, aes(Year, C14t, colour = "Bulk"), linewidth = 1) +
    geom_line(data = out_best, aes(Year, C14t_fast, colour = "Fast"), linetype = "dashed") +
    geom_line(data = out_best, aes(Year, C14t_slow, colour = "Slow"), linetype = "dotted") +
    
    # observations
    geom_point(data = C14obs_bulk, aes(Year, C14t, colour = "Bulk")) +
    geom_point(data = C14obs_325, aes(Year, C14t_fast, colour = "Fast")) +
    geom_point(data = C14obs_400, aes(Year, C14t_slow, colour = "Slow")) +
    
    scale_colour_manual(
      name = "Pool",
      values = c("Bulk" = col_bulk,
                 "Fast" = col_fast,
                 "Slow" = col_slow)
    ) +
    
    ggtitle(paste(expression(Delta^14*C), "-", site_name)) +
    ylab(Delta14Clabel) +
    xlab("Year") +
    theme_minimal()
  
  ## return ------------------------------------------------------------------
  
  return(list(
    site = site_name,
    pars = bestpars,
    output = out_best,
    plot_C = plot_C,
    plot_C14 = plot_C14
  ))
}
## --------------------------------------------------------------------------
## RUN ALL SITES
## --------------------------------------------------------------------------

results <- lapply(sites, run_model_site)
names(results) <- sites

## save all plots to local dir ----------------------------------------------

output_dir <- "plots/simple_twopool"
for (site in names(results)) {
  
  res <- results[[site]]
  
  # file names
  file_C <- file.path(output_dir, paste0(site, "_C.png"))
  file_C14 <- file.path(output_dir, paste0(site, "_C14.png"))
  
  # save plots
  ggsave(filename = file_C, plot = res$plot_C,
         width = 7, height = 5, dpi = 300, bg = "white")
  
  ggsave(filename = file_C14, plot = res$plot_C14,
         width = 7, height = 5, dpi = 300, bg = "white")
}

## save all plots as rds for loading to qmd ---------------------------------
for (name in names(results)) {
  saveRDS(results[[name]]$plot_C,
          file = file.path(output_dir, paste0(name, "_C.rds")))
  
  saveRDS(results[[name]]$plot_C14,
          file = file.path(output_dir, paste0(name, "_14C.rds")))
}

save(results, file   = "modruns/results_simple_twop.Rdata")

## --------------------------------------------------------------------------
## combine outputs
## --------------------------------------------------------------------------

# library(purrr)
# summary_df <- map_dfr(results, ~ data.frame(site = .x$site))

### MCMC optimization
var0_M2 <- mFit_M2$var_ms_unweighted
cov0_M2 <- summary(mFit_M2)$cov.scaled # The covariance matrix can be used for the jump 
# uncomment to run again 
 MCMC_M2 <- modMCMC(f=mc, p = bestpars_M2, niter = 2500, jump = NULL, var0 = var0_M2, wvar0 = 1) #tnel: jump default is 10% of par value; what is wvar and why is it zero?
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

#edit from here onwards

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
