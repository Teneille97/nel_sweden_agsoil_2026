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

#mean_C_inputs <- mean(BolinderIn) * 1000  # g m-2#
mean_C_inputs <- 80
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
C0_fast <- 5100
C0_400 <- mean(Cobs_slow[Cobs_slow$Year==1957,]$Ct_slow)
F0_400 <- mean(C14obs_slow[C14obs_slow$Year==1957,]$C14t_slow)
F0_bulk <- mean(C14obs_bulk[C14obs_bulk$Year==1957,]$C14t) 

# func to run mod
run_mod <- function(pars){ #kf, ki, ks, alpha 21, +F0b=Fi, alpha 32, alpha 41, alpha 42
  c14_atm <- BoundFc(Atm14C, format = "Delta14C")
  c14_initial <- ConstFc(
    values = c(F0_bulk + 10, F0_bulk + pars[5], F0_400, F0_bulk + 10), 
    format = "Delta14C"
  )
  A4 <- diag(-c(pars[1:3],0.000000001))
  A4[2,1] <- pars[1]*pars[4]
  A4[3,2] <- pars[2]*pars[6]
  A4[4,1] <- pars[1]*pars[7]
  A4[4,2] <- pars[2]*pars[8]
  mod<-GeneralModel_14(
    t = yr,
    A = A4,
    ivList = c(C0_fast, abs(C0_bulk-C0_fast-C0_400), C0_400, 0),
    initialValF = c14_initial,
    inputFluxes = c(mean_C_inputs,0,0,0),
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

inipars <- c(0.1,0.05,0.001, 0.1, -50, 0.1, 0.005, 0.005) #kf, ki, ks, alpha 21, F0ib, I, alpha 32, alpha 41

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
  upper = c(0.5,0.2,0.005, 0.5,0,1, 0.01, 0.01), #kf, ki, ks, alpha 21, F0ib, alpha 32, alpha 41, alpha 42
  lower = c(0.05,0.01,0.0001,0,-160,0, 0, 0) 
)

bestpars <- mFit$par
out_best <- run_mod(bestpars)
Delta14Clabel <- expression(Delta^14*C)
stocks_label <- expression(C ~ (g ~ m^{-2}))

col_bulk <- "black"
col_fast <- "#1b9e77"
col_slow <- "#BF40BF"
col_inter <- "#0000FF"
col_loss <- "#FF0000"

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

#save(results_4p, file = file.path("mod_runs", "results_allsites_4pool.Rdata"))
load(here::here("mod_runs/results_allsites_4pool.Rdata"))
out_best<-results_4p$output
bestpars<-results_4p$pars
## MCMC

var0 <- mFit$var_ms_unweighted
cov0 <- summary(mFit)$cov.scaled #  cov matrix can be used for jump 

# run 
MCMC <- modMCMC(f=mc, p = bestpars, niter = 50000, jump = cov0*0.001, var0 = var0, wvar0 = 1, updatecov = 1000, burninlength =  1000, 
                upper = c(0.1,0.2,0.01, 0.5,0,1, 0.01, 0.01), #kf, ki, ks, alpha 21, F0ib, alpha 32, alpha 41, alpha 42
                lower = c(0.05,0.005,0.0001,0,-20,0, 0, 0)) 
#save(MCMC, file = file.path("mod_runs", "MCMC_allsites_4pools.Rdata"))
load(here::here("mod_runs/MCMC_allsites_4pools.Rdata"))

# view distribution
class(MCMC) <- "modMCMC" #make sure FME library is loaded
parsMCMC<-summary(MCMC)

# plot convergence
convergence_plot<-plot(MCMC) 

# plot densitites
par(mfrow = c(3, 3)) 
for (p in colnames(MCMC$pars)) {
  plot(density(MCMC$pars[, p]), main = paste("Density of", p))
}

# check inter-dependence of parameters using 500 MCMC runs
pairs(MCMC,nsample=500) 

# performance
par(mfrow = c(1, 1))  
hist(MCMC$SS, breaks = 50)
mc(mFit$par)
percentage_accepted<-(100*MCMC$naccepted)/50000
MCMC_bestpars<-MCMC$bestpar
cost_modfit <- mc(bestpars)$model #should be lower than mcmc
cost_mcmc   <- mc(MCMC_bestpars)$model  

# plots showing 95% credible interval of model outputs given the posterior parameter distribution
# ie posterior predictive uncertainty conditional on sampled MCMC chains

set.seed(1)
pars_sub <- MCMC$pars[sample(1:nrow(MCMC$pars), 500), ]

runs <- lapply(1:nrow(pars_sub), function(i) run_mod(pars_sub[i, ]))
#save(runs, file = file.path("mod_runs", "runs.Rdata"))
load(here::here("mod_runs/runs.Rdata"))


# convert to array-like structure
extract_var <- function(var){
  sapply(runs, function(x) x[[var]])
}

vars <- c("Ct","Ct_fast","Ct_slow","Ct_inter", "Ct_loss",
          "C14t","C14t_fast","C14t_slow","C14t_inter", "C14t_loss")

unc_list <- lapply(vars, function(v){
  mat <- extract_var(v)
  
  data.frame(
    Year = runs[[1]]$Year,
    var = v,
    # Use na.rm = TRUE to skip the failed runs
    Mean = rowMeans(mat, na.rm = TRUE),
    Low  = apply(mat, 1, quantile, 0.025, na.rm = TRUE),
    High = apply(mat, 1, quantile, 0.975, na.rm = TRUE)
  )
})

unc_df <- do.call(rbind, unc_list)

# var names for plotting
unc_df$Pool <- dplyr::case_when(
  unc_df$var == "Ct" ~ "Bulk",
  unc_df$var == "Ct_fast" ~ "Fast",
  unc_df$var == "Ct_slow" ~ "Slow",
  unc_df$var == "Ct_inter" ~ "Inter",
  unc_df$var == "Ct_loss" ~ "Loss",
  unc_df$var == "C14t" ~ "Bulk",
  unc_df$var == "C14t_fast" ~ "Fast",
  unc_df$var == "C14t_slow" ~ "Slow",
  unc_df$var == "C14t_inter" ~ "Inter",
  unc_df$var == "C14t_loss" ~ "Loss"
)
unc_C    <- subset(unc_df, grepl("^Ct", var))
unc_C14  <- subset(unc_df, grepl("^C14t", var))

# C stocks plot final
plot_C_final <- ggplot() +
  
  ## --- MCMC ribbons (model uncertainty) ---
  geom_ribbon(data = unc_C,
              aes(x = Year, ymin = Low, ymax = High, fill = Pool),
              alpha = 0.2) +
  
  ## --- model lines ---
  geom_line(data = out_best, aes(Year, Ct, colour = "Bulk"), linewidth = 1) +
  geom_line(data = out_best, aes(Year, Ct_fast, colour = "Fast"), linetype = "dashed") +
  geom_line(data = out_best, aes(Year, Ct_slow, colour = "Slow"), linetype = "dotted") +
  geom_line(data = out_best, aes(Year, Ct_inter, colour = "Inter"), linetype = "dotdash") +
  geom_line(data = out_best, aes(Year, Ct_loss, colour = "Loss"), linetype = "dotted") +
  
  ## --- observations: points + error bars ---
  geom_point(data = Cobs_bulk, aes(Year, Ct, colour = "Bulk")) +
  geom_errorbar(data = Cobs_bulk,
                aes(Year, ymin = Ct - Ct_sd, ymax = Ct + Ct_sd, colour = "Bulk"),
                width = 0.5) +
  
  geom_point(data = Cobs_fast, aes(Year, Ct_fast, colour = "Fast")) +
  geom_errorbar(data = Cobs_fast,
                aes(Year, ymin = Ct_fast - Ct_fast_sd, ymax = Ct_fast + Ct_fast_sd, colour = "Fast"),
                width = 0.5) +
  
  geom_point(data = Cobs_slow, aes(Year, Ct_slow, colour = "Slow")) +
  geom_errorbar(data = Cobs_slow,
                aes(Year, ymin = Ct_slow - Ct_slow_sd, ymax = Ct_slow + Ct_slow_sd, colour = "Slow"),
                width = 0.5) +
  
  geom_point(data = Cobs_inter, aes(Year, Ct_inter, colour = "Inter")) +
  geom_errorbar(data = Cobs_inter,
                aes(Year, ymin = Ct_inter - Ct_inter_sd, ymax = Ct_inter + Ct_inter_sd, colour = "Inter"),
                width = 0.5) +
  
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
               "Inter" = col_inter,
               "Loss" = col_loss)
  ) +
  
  ggtitle("C stocks") +
  ylab(stocks_label) +
  xlab("Year") +
  theme_minimal()

# 14C plot final

plot_C14_final <- ggplot() +
  
  ## --- MCMC ribbons ---
  geom_ribbon(data = unc_C14,
              aes(x = Year, ymin = Low, ymax = High, fill = Pool),
              alpha = 0.2) +
  
  ## --- model lines ---
  geom_line(data = out_best, aes(Year, C14t, colour = "Bulk"), linewidth = 1) +
  geom_line(data = out_best, aes(Year, C14t_fast, colour = "Fast"), linetype = "dashed") +
  geom_line(data = out_best, aes(Year, C14t_slow, colour = "Slow"), linetype = "dotted") +
  geom_line(data = out_best, aes(Year, C14t_inter, colour = "Inter"), linetype = "dotdash") +
  geom_line(data = out_best, aes(Year, C14t_loss, colour = "Loss"), linetype = "dotted") +
  
  ## --- observations ---
  geom_point(data = C14obs_bulk, aes(Year, C14t, colour = "Bulk")) +
  geom_errorbar(data = C14obs_bulk,
                aes(Year, ymin = C14t - C14t_sd, ymax = C14t + C14t_sd, colour = "Bulk"),
                width = 0.5) +
  
  geom_point(data = C14obs_slow, aes(Year, C14t_slow, colour = "Slow")) +
  geom_errorbar(data = C14obs_slow,
                aes(Year, ymin = C14t_slow - C14t_slow_sd, ymax = C14t_slow + C14t_slow_sd, colour = "Slow"),
                width = 0.5) +
  
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
               "Inter" = col_inter,
               "Loss" = col_loss)
  ) +
  
  ggtitle(expression(Delta^14*C)) +
  ylab(Delta14Clabel) +
  xlab("Year") +
  theme_minimal()

# save plots
output_dir <- "plots/allsites/4_pool"

file_C_final <- file.path(output_dir, "C_final.png")
file_C14_final <- file.path(output_dir, "C14_final.png")

ggsave(filename = file_C_final, plot = plot_C_final,
       width = 7, height = 5, dpi = 300, bg = "white")
ggsave(filename = file_C14_final, plot = plot_C14_final,
       width = 7, height = 5, dpi = 300, bg = "white")

saveRDS(plot_C_final,
        file = file.path(output_dir, "C.rds"))
saveRDS(plot_C14_final,
        file = file.path(output_dir, "14C.rds"))


# Age and transit time
meanpars<-as.numeric(parsMCMC[1,1:8]) # mean 
minpars<-as.numeric(parsMCMC[3,1:8]) # min
maxpars<-as.numeric(parsMCMC[4,1:8]) # max

#A_mean <- diag(-c(meanpars[1:3],0.000000001)) #kf, ki, ks, alpha 21, F0ib, alpha 32, alpha 41, alpha 42
A_mean <- diag(-c(meanpars[1:3])) # the loss pool has infinite age because negligible turnover
A_mean[2,1] <- meanpars[1]*meanpars[4]
A_mean[3,2] <- meanpars[2]*meanpars[6]
#A_mean[4,1] <- meanpars[1]*meanpars[7]
#A_mean[4,2] <- meanpars[2]*meanpars[8]

tau=seq(0,500)

SA=systemAge(A=A_mean,u=c(1,0,0),a=tau) #here u are normalized inputs i.e. simply divided among pools
SA$meanSystemAge
par(mfrow = c(1, 1))  
plot(SA$systemAgeDensity)
plot(SA$poolAgeDensity[,1])
plot(SA$poolAgeDensity[,2])
plot(SA$poolAgeDensity[,3])

SA$meanPoolAge

TT=transitTime(A=A_mean,u=c(1,0,0), a=tau, q=c(0.1,0.5,0.9)) # system info only; we need to get these for pools?
plot(TT$transitTimeDensity)
TT$transitTimeDensity


## --- Summarize Age and Transit Time Distribution ---

# Helper function for formatting
sum_age_tt_fun <- function(vector) {
  data.frame(
    mean_age    = mean(vector),
    mean_age_sd = sd(vector),
    q05         = as.numeric(quantile(vector, 0.05)),
    median      = median(vector),
    q95         = as.numeric(quantile(vector, 0.95))
  )
}

# Bind the results from your previously run 'age_dist'
SA_TT_sum_final <- data.frame(
  rbind(
    cbind(Variable = "System",       sum_age_tt_fun(age_dist$system_age)),
    cbind(Variable = "Fast Pool",    sum_age_tt_fun(age_dist$fast_age)),
    cbind(Variable = "Intermediate", sum_age_tt_fun(age_dist$inter_age)),
    cbind(Variable = "Slow Pool",    sum_age_tt_fun(age_dist$slow_age)),
    cbind(Variable = "Transit Time", sum_age_tt_fun(age_dist$transit_time))
  )
)

print(SA_TT_sum_final)

## --- Plot Age Densities 

age_dens_long <- age_df %>%
  pivot_longer(c(fast, inter, slow), names_to = "Pool", values_to = "Density") %>%
  mutate(Pool = factor(Pool, levels = c("fast", "inter", "slow"), 
                       labels = c("Fast (kf)", "Intermediate (ki)", "Slow (ks)")))

plot_age_densities <- ggplot(age_dens_long, aes(x = age, y = Density, color = Pool)) +
  geom_line(linewidth = 1) +
  facet_wrap(~Pool, scales = "free") +
  theme_classic() +
  scale_color_manual(values = c("#33a02c", "#b2df8a", "#a6cee3")) +
  labs(x = "Age (years)", y = "Density Function", title = "Pool Age Distributions")

print(plot_age_densities)

## --- Calculate Transit Time Densities for Individual Pools ---

tau <- seq(0, 1000)

# Calculate transit time density for the system
TT_sys <- transitTime(A = A_med, u = c(1,0,0), a = tau)

# For individual pools
pool_transit_list <- lapply(1:3, function(i) {
  u_pool <- rep(0, 3)   # unit vector for the pool
  u_pool[i] <- 1
  TT_pool <- transitTime(A = A_med, u = u_pool, a = tau)
  return(TT_pool$transitTimeDensity)   # transit time starting from pool i
})

tt_density_df <- data.frame(
  age = tau,
  Fast = pool_transit_list[[1]],
  Intermediate = pool_transit_list[[2]],
  Slow = pool_transit_list[[3]],
  System = TT_sys$transitTimeDensity
)

## plot tt densities

tt_long <- tt_density_df %>%
  pivot_longer(-age, names_to = "Component", values_to = "Density")

ggplot(tt_long, aes(x = age, y = Density, color = Component)) +
  geom_line(linewidth = 1) +
  facet_wrap(~Component, scales = "free_y") +
  theme_classic() +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Transit Time Densities (Exit Time Distributions)",
    x = "Time in System (years)",
    y = "Probability Density",
    caption = "Note: System Transit Time is the weighted average of pool transit times based on inputs (u)."
  )

# other approach

propagate_age_tt_full <- function(MCMC_pars, n_iter = 5000) {
  
  results <- data.frame(
    system_age = numeric(n_iter),
    fast_age   = numeric(n_iter),
    inter_age  = numeric(n_iter),
    slow_age   = numeric(n_iter),
    # Transit times for the whole system vs starting from specific pools
    system_tt  = numeric(n_iter),
    fast_tt    = numeric(n_iter),
    inter_tt   = numeric(n_iter),
    slow_tt    = numeric(n_iter)
  )
  
  for(i in 1:n_iter) {
    pars <- MCMC_pars[sample(1:nrow(MCMC_pars), 1), ]
    A <- build_A(pars)
    
    # 1. System/Pool Ages
    SA <- systemAge(A = A, u = c(1,0,0))
    results$system_age[i] <- SA$meanSystemAge
    results$fast_age[i]   <- SA$meanPoolAge[1]
    results$inter_age[i]  <- SA$meanPoolAge[2]
    results$slow_age[i]   <- SA$meanPoolAge[3]
    
    # 2. Transit Times (Exit times starting from each pool)
    results$system_tt[i] <- transitTime(A = A, u = c(1,0,0))$meanTransitTime
    results$fast_tt[i]   <- transitTime(A = A, u = c(1,0,0))$meanTransitTime # Usually same as system if all input enters Fast
    results$inter_tt[i]  <- transitTime(A = A, u = c(0,1,0))$meanTransitTime
    results$slow_tt[i]   <- transitTime(A = A, u = c(0,0,1))$meanTransitTime
  }
  return(results)
}

# Run the updated propagation
age_tt_dist <- propagate_age_tt_full(MCMC$pars, n_iter = 5000)


# Calculate densities using median parameters
tau <- seq(0, 1500)
# Exit distribution starting from Fast (Pool 1)
TT_fast  <- transitTime(A = A_med, u = c(1,0,0), a = tau)
# Exit distribution starting from Inter (Pool 2)
TT_inter <- transitTime(A = A_med, u = c(0,1,0), a = tau)
# Exit distribution starting from Slow (Pool 3)
TT_slow  <- transitTime(A = A_med, u = c(0,0,1), a = tau)

tt_density_df <- data.frame(
  age = tau,
  Fast = TT_fast$transitTimeDensity,
  Intermediate = TT_inter$transitTimeDensity,
  Slow = TT_slow$transitTimeDensity
)

# Plotting
tt_density_df %>%
  pivot_longer(-age, names_to = "Starting_Pool", values_to = "Density") %>%
  ggplot(aes(x = age, y = Density, color = Starting_Pool)) +
  geom_line(linewidth = 1) +
  facet_wrap(~Starting_Pool, scales = "free_y") +
  theme_minimal() +
  labs(title = "Transit Time Densities by Starting Pool",
       subtitle = "Probability of exiting the system at time t, given starting pool",
       x = "Years", y = "Density")


## --- Calculate fluxes and stocks for all pools 

# 1. par data
pf_4ps <- MCMC$pars %>%
  as.data.frame() %>%
  mutate(
    C_fast  = mean(unc_C$Mean[unc_C$Pool == "Fast"]),
    C_inter = mean(unc_C$Mean[unc_C$Pool == "Inter"]),
    C_slow  = mean(unc_C$Mean[unc_C$Pool == "Slow"]),
    C_loss  = mean(unc_C$Mean[unc_C$Pool == "Loss"])
  ) %>%
  mutate(
    out_fast = (1 - MCMC$bestpar[4] - MCMC$bestpar[7]) * MCMC$bestpar[1] * C_fast, #respiration
    trans_fast_to_inter = MCMC$bestpar[4] * MCMC$bestpar[1] * C_fast,
    trans_fast_to_loss = MCMC$bestpar[7] * MCMC$bestpar[1] * C_fast,
    out_inter = (1 - MCMC$bestpar[6] - MCMC$bestpar[8]) * MCMC$bestpar[2] * C_inter,
    trans_inter_to_slow = MCMC$bestpar[6] * MCMC$bestpar[2] * C_inter,
    trans_inter_to_loss = MCMC$bestpar[8] * MCMC$bestpar[2] * C_inter,
    out_slow = MCMC$bestpar[3] * C_slow
  )

flux_summary <- pf_4ps %>%
  summarise(across(starts_with("out") | starts_with("trans"), list(
    median = ~round(median(.), 2),
    lci    = ~round(quantile(., 0.05), 2),
    uci    = ~round(quantile(., 0.95), 2)
  ))) %>%
  pivot_longer(everything(), names_to = "Metric", values_to = "Value")

print(flux_summary)