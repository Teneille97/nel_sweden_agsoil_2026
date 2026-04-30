# SOM decomposition models for agricultural soils at multiple sites

## load libraries
library(here)
library(SoilR)
library(FME)
library(ggplot2)
library(dplyr)
library(tidyr)

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
# --- Add N stocks ---
all_data$N_stocks_gm2 <- all_data$TONgkg * bulk_density * 0.2 * 1000

summary_df <- all_data %>%
  group_by(Year, Temperature...C.) %>%
  summarise(
    C_stocks_mean = mean(C_stocks_gm2, na.rm = TRUE),
    C_stocks_sd   = sd(C_stocks_gm2, na.rm = TRUE),
    N_stocks_mean = mean(N_stocks_gm2, na.rm = TRUE),
    N_stocks_sd   = sd(N_stocks_gm2, na.rm = TRUE),
    d14C_mean     = mean(d14C, na.rm = TRUE),
    d14C_sd       = sd(d14C, na.rm = TRUE),
    CN_mean       = mean(MolarCN, na.rm = TRUE),
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


# obs N stocks (not for models, only final age distributions)
Nobs_bulk <- data.frame(Year = bulk$Year, Nt = bulk$N_stocks_mean, Nt_sd = bulk$N_stocks_sd)
Nobs_325 <- data.frame(Year = pool_325$Year, Nt_325 = pool_325$N_stocks_mean, Nt_325_sd = pool_325$N_stocks_sd)
Nobs_slow <- data.frame(Year = pool_400$Year, Nt_slow = pool_400$N_stocks_mean, Nt_slow_sd = pool_400$N_stocks_sd)
Nobs_fast <- data.frame(Year = Nobs_325$Year, Nt_fast = Nobs_bulk$Nt-Nobs_325$Nt_325, Nt_fast_sd = sqrt(bulk$N_stocks_sd^2 + pool_325$N_stocks_sd^2)) # C fast  = C bulk - C 325
Nobs_inter<-data.frame(Year = Nobs_325$Year, Nt_inter = Nobs_325$Nt_325-Nobs_slow$Nt_slow, Nt_inter_sd = sqrt(Nobs_325$Nt_325_sd^2 + Nobs_slow$Nt_slow_sd^2))


#remove NA sd values
Cobs_bulk$Ct_sd[is.na(Cobs_bulk$Ct_sd)] <- mean(Cobs_bulk$Ct_sd, na.rm = TRUE) 
C14obs_bulk$C14t_sd[is.na(C14obs_bulk$C14t_sd)] <- mean(C14obs_bulk$C14t_sd, na.rm = TRUE)
Cobs_fast$Ct_fast_sd[is.na(Cobs_fast$Ct_fast_sd)] <- mean(Cobs_fast$Ct_fast_sd, na.rm = TRUE)
C14obs_slow$C14t_slow_sd[is.na(C14obs_slow$C14t_slow_sd)] <- mean(C14obs_slow$C14t_slow_sd, na.rm = TRUE)
Cobs_inter$Ct_inter_sd[is.na(Cobs_inter$Ct_inter_sd)] <- mean(Cobs_inter$Ct_inter_sd, na.rm = TRUE)
Cobs_slow$Ct_slow_sd[is.na(Cobs_slow$Ct_slow_sd)] <- mean(Cobs_slow$Ct_slow_sd, na.rm = TRUE)


Nobs_bulk$Nt_sd[is.na(Nobs_bulk$Nt_sd)] <- mean(Nobs_bulk$Nt_sd, na.rm = TRUE) 
Nobs_fast$Nt_fast_sd[is.na(Nobs_fast$Nt_fast_sd)] <- mean(Nobs_fast$Nt_fast_sd, na.rm = TRUE)
Nobs_inter$Nt_inter_sd[is.na(Nobs_inter$Nt_inter_sd)] <- mean(Nobs_inter$Nt_inter_sd, na.rm = TRUE)
Nobs_slow$Nt_slow_sd[is.na(Nobs_slow$Nt_slow_sd)] <- mean(Nobs_slow$Nt_slow_sd, na.rm = TRUE)


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
#mFit <- modFit(
 # f = mc,
#  p = inipars,
#  method = "Nelder-Mead",
#  upper = c(0.5,0.2,0.005, 0.5,0,1, 0.01, 0.01), #kf, ki, ks, alpha 21, F0ib, alpha 32, alpha 41, alpha 42
#  lower = c(0.05,0.01,0.0001,0,-160,0, 0, 0) 
#)


#bestpars <- mFit$par
#out_best <- run_mod(bestpars)
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
#MCMC <- modMCMC(f=mc, p = bestpars, niter = 50000, jump = cov0*0.001, var0 = var0, wvar0 = 1, updatecov = 1000, burninlength =  1000, 
#                upper = c(0.1,0.2,0.01, 0.5,0,1, 0.01, 0.01), #kf, ki, ks, alpha 21, F0ib, alpha 32, alpha 41, alpha 42
#                lower = c(0.05,0.005,0.0001,0,-20,0, 0, 0)) 
#save(MCMC, file = file.path("mod_runs", "MCMC_allsites_4pools.Rdata"))
load(here::here("mod_runs/MCMC_allsites_4pools.Rdata"))

# view distribution
class(MCMC) <- "modMCMC" #make sure FME library is loaded
parsMCMC<-summary(MCMC)

# plot convergence
convergence_plot<-plot(MCMC) 

output_dir <- "plots/allsites/4_pool"
file_convergence_plot <- file.path(output_dir, "convergence_plot.png")

png(file_convergence_plot, width = 10, height = 8, units = "in", res = 300)
plot(MCMC)
dev.off()

# plot densitites
mcmc_long <- as.data.frame(MCMC$pars) %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value")

density_plot <- ggplot(mcmc_long, aes(x = value)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  facet_wrap(~ parameter, scales = "free") +
  theme_minimal() +
  labs(title = "MCMC Parameter Densities", x = "Value", y = "Density")

file_density_plot <- file.path(output_dir, "density_plot.png")
ggsave(filename = file_density_plot, plot = density_plot,
       width = 10, height = 8, dpi = 300, bg = "white")
saveRDS(density_plot,
        file = file.path(output_dir, "density_plot.rds"))


# check inter-dependence of parameters using 500 MCMC runs
pairs_plot<-pairs(MCMC,nsample=500) 

file_pairs_plot <- file.path(output_dir, "pairs_plot.png")

png(file_pairs_plot, width = 10, height = 8, units = "in", res = 300)
pairs(MCMC, nsample = 500)
dev.off()
saveRDS(pairs_plot,
        file = file.path(output_dir, "pairs_plot.rds"))

# performance
par(mfrow = c(1, 1))  
hist(MCMC$SS, breaks = 50)
percentage_accepted<-(100*MCMC$naccepted)/50000
MCMC_bestpars<-MCMC$bestpar
cost_modfit <- mc(bestpars)$model 
cost_mcmc   <- mc(MCMC_bestpars)$model  

# plots showing 95% credible interval of model outputs given the posterior parameter distribution
# i.e. posterior predictive uncertainty conditional on sampled MCMC chains

set.seed(1)
pars_sub <- MCMC$pars[sample(1:nrow(MCMC$pars), 500), ]

#runs <- lapply(1:nrow(pars_sub), function(i) run_mod(pars_sub[i, ]))
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
#### run mFit using MCMC pars

# first need to redefine upper and lower as that used for MCMC 
#mFit <- modFit(
 # f = mc,
#  p = MCMC_bestpars,
#  method = "Nelder-Mead",
#  upper = c(0.1,0.2,0.01, 0.5,0,1, 0.01, 0.01), #kf, ki, ks, alpha 21, F0ib, alpha 32, alpha 41, alpha 42
#  lower = c(0.05,0.005,0.0001,0,-20,0, 0, 0) 
#)

# run mod
#out_best_MCMC<- run_mod(MCMC_bestpars)
#save(out_best_MCMC, file = file.path("mod_runs", "out_best_MCMC.Rdata"))
load(here::here("mod_runs/out_best_MCMC.Rdata"))

# plot using predictions made by mcmc
plot_C_final <- ggplot() +
  
  ## --- MCMC ribbons (model uncertainty) ---
  geom_ribbon(data = unc_C,
              aes(x = Year, ymin = Low, ymax = High, fill = Pool),
              alpha = 0.2) +
  
  ## --- model lines ---
  geom_line(data = out_best_MCMC, aes(Year, Ct, colour = "Bulk"), linewidth = 1) +
  geom_line(data = out_best_MCMC, aes(Year, Ct_fast, colour = "Fast"), linetype = "dashed") +
  geom_line(data = out_best_MCMC, aes(Year, Ct_slow, colour = "Slow"), linetype = "dotted") +
  geom_line(data = out_best_MCMC, aes(Year, Ct_inter, colour = "Inter"), linetype = "dotdash") +
  geom_line(data = out_best_MCMC, aes(Year, Ct_loss, colour = "Loss"), linetype = "dotted") +
  
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
  geom_line(data = out_best_MCMC, aes(Year, C14t, colour = "Bulk"), linewidth = 1) +
  geom_line(data = out_best_MCMC, aes(Year, C14t_fast, colour = "Fast"), linetype = "dashed") +
  geom_line(data = out_best_MCMC, aes(Year, C14t_slow, colour = "Slow"), linetype = "dotted") +
  geom_line(data = out_best_MCMC, aes(Year, C14t_inter, colour = "Inter"), linetype = "dotdash") +
  geom_line(data = out_best_MCMC, aes(Year, C14t_loss, colour = "Loss"), linetype = "dotted") +
  
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
        file = file.path(output_dir, "C_final.rds"))
saveRDS(plot_C14_final,
        file = file.path(output_dir, "14C_final.rds"))

### Ages and transit times ###

# sample mcmc pars
set.seed(1)

burnin <- 20000  # or more, see below
pars_full <- as.data.frame(MCMC$pars)
pars_post <- pars_full[-(1:burnin), ]
pars_sub <- pars_post[sample(1:nrow(pars_post), 700), ]

# sample ~700 parameter sets
n_samp <- 700
pars_sub <- pars_full[sample(1:nrow(pars_full), n_samp), ]

# func to build 3 pool system matrix
build_A_u <- function(pars){
  # kf, ki, ks, alpha_fi, x, alpha_is, alpha_fl, alpha_il
  kf <- pars[1]
  ki <- pars[2]
  ks <- pars[3]
  alpha_fi <- pars[4]
  alpha_is <- pars[6]
  A <- diag(-c(kf, ki, ks)) # 3 pool system ---
  A[2,1] <- kf * alpha_fi
  A[3,2] <- ki * alpha_is
  u <- matrix(c(mean_C_inputs, 0, 0), ncol = 1) # input vector (only to fast pool)
  return(list(A = A, u = u))
}

# compute densities for one par set
get_age_tt_dens <- function(pars, ages){
  AU <- build_A_u(pars)
  SA <- systemAge(A = AU$A, u = AU$u, a = ages)
  TT <- transitTime(A = AU$A, u = AU$u, a = ages)
  data.frame(
    age = ages,
    system_age = SA$systemAgeDensity,
    fast = SA$poolAgeDensity[,1],
    inter = SA$poolAgeDensity[,2],
    slow = SA$poolAgeDensity[,3],
    transit_time = TT$transitTimeDensity
  )
}

# run for pars 
MCMC_bestpars<-MCMC$bestpar
ages <- seq(0, 500, by = 1)
dens_best <- get_age_tt_dens(MCMC_bestpars, ages)

# run for sampled mcmc pars 
dens_list <- lapply(1:nrow(pars_sub), function(i){
  get_age_tt_dens(as.numeric(pars_sub[i, ]), ages)
})

# save(dens_list, file = file.path("mod_runs", "dens_list.Rdata"))
load(here::here("mod_runs/dens_list.Rdata"))

# uncertainty envelopes
extract_var <- function(var){
  sapply(dens_list, function(x) x[[var]])
}

vars <- c("system_age","fast","inter","slow","transit_time")

unc_dens <- lapply(vars, function(v){
  
  mat <- extract_var(v)
  
  data.frame(
    age = ages,
    variable = v,
    mean = rowMeans(mat, na.rm = TRUE),
    median = apply(mat, 1, median, na.rm = TRUE),
    low  = apply(mat, 1, quantile, 0.025, na.rm = TRUE),
    high = apply(mat, 1, quantile, 0.975, na.rm = TRUE)
  )
})

unc_dens_df <- do.call(rbind, unc_dens)

# best pars for plotting
dens_best_long <- dens_best %>%
  pivot_longer(-age, names_to = "variable", values_to = "value")

# summary stats func
sum_fun <- function(age, dens){
  dens <- dens / sum(dens)   # <-- THIS is the key
  mean_age <- sum(age * dens)
  cdf <- cumsum(dens)
  
  data.frame(
    mean   = mean_age,
    median = age[which.min(abs(cdf - 0.5))]
  )
}

# summary stats from mean of MCMC densities
summary_list <- lapply(unique(unc_dens_df$variable), function(v){
  df <- unc_dens_df %>% filter(variable == v)
  stats <- sum_fun(df$age, df$mean)
  cbind(variable = v, stats)
})

summary_table <- do.call(rbind, summary_list)

summary_table_fmt <- summary_table %>%
  mutate(across(-variable, ~round(., 1)))


save(summary_table_fmt, file = file.path("mod_runs", "summary_table_fmt.Rdata"))
load(here::here("mod_runs/summary_table_fmt.Rdata"))


# plot
plot_age_tt <- ggplot() +
  # --- ribbon ---
  geom_ribbon(data = unc_dens_df,
              aes(x = log(age), ymin = low, ymax = high),
              fill = "grey70",
              alpha = 0.5) +
  # --- mean density curve (optional) ---
  geom_line(data = unc_dens_df,
            aes(x = log(age), y = mean),
            colour = "black",
            linewidth = 0.8) +
  facet_wrap(~variable, scales = "free", ncol = 2) +
  
  # --- MEDIAN ---
  geom_vline(data = summary_table,
             aes(xintercept = log(median), colour = "Median"),
             linetype = "dotted",
             linewidth = 1) +
  # --- MEAN ---
  geom_vline(data = summary_table,
             aes(xintercept = log(mean), colour = "Mean"),
             linetype = "dashed",
             linewidth = 1) +
  # --- legend control ---
  scale_colour_manual(
    name = "Statistic",
    values = c("Median" = "red", "Mean" = "blue")
  ) +
  labs(x = "Log age (years)", y = "Density") +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )


output_dir <- "plots/allsites/4_pool"
file_plot_age_tt <- file.path(output_dir, "plot_age_tt.png")

png(file_plot_age_tt, width = 10, height = 8, units = "in", res = 300)
print(plot_age_tt)
dev.off()

saveRDS(plot_age_tt,
        file = file.path(output_dir, "plot_age_tt.rds"))

# log age plots
plot_log_age_tt <- ggplot() +
  # --- ribbon ---
  geom_ribbon(data = unc_dens_df,
              aes(x = age, ymin = low, ymax = high),
              fill = "grey70",
              alpha = 0.5) +
  # --- mean density curve (optional) ---
  geom_line(data = unc_dens_df,
            aes(x = age, y = mean),
            colour = "black",
            linewidth = 0.8) +
  facet_wrap(~variable, scales = "free", ncol = 2) +
  
  # --- MEDIAN ---
  geom_vline(data = summary_table,
             aes(xintercept = median, colour = "Median"),
             linetype = "dotted",
             linewidth = 1) +
  # --- MEAN ---
  geom_vline(data = summary_table,
             aes(xintercept = mean, colour = "Mean"),
             linetype = "dashed",
             linewidth = 1) +
  # --- legend control ---
  scale_colour_manual(
    name = "Statistic",
    values = c("Median" = "red", "Mean" = "blue")
  ) +
  scale_x_log10(
    limits = c(1, 500),
    breaks = c(1, 10, 50, 100, 200, 500)
  )+
  labs(x = "Age (years), log-scale", y = "Density") +
  theme_minimal(base_size = 14) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

output_dir <- "plots/allsites/4_pool"
file_plot_log_age_tt <- file.path(output_dir, "plot_log_age_tt.png")

png(file_plot_log_age_tt, width = 10, height = 8, units = "in", res = 300)
print(plot_log_age_tt)
dev.off()

saveRDS(plot_log_age_tt,
        file = file.path(output_dir, "plot_log_age_tt.rds"))

## --- Calculate fluxes and stocks for all pools 

flux_list <- lapply(1:length(runs), function(i){
  
  pars <- as.numeric(MCMC$pars[i, ])
  run  <- runs[[i]]
  
  # take pool sizes at final year (quasi steady state)
  C_fast  <- tail(run$Ct_fast, 1)
  C_inter <- tail(run$Ct_inter, 1)
  C_slow  <- tail(run$Ct_slow, 1)
  
  kf <- pars[1]
  ki <- pars[2]
  ks <- pars[3]
  alpha_fi <- pars[4]
  alpha_is <- pars[6]
  alpha_fl <- pars[7]
  alpha_il <- pars[8]
  
  data.frame(
    fast_to_loss  = alpha_fl * kf * C_fast,
    fast_to_inter = alpha_fi * kf * C_fast,
    
    inter_to_loss = alpha_il * ki * C_inter,
    inter_to_slow = alpha_is * ki * C_inter,
    
    slow_to_loss  = ks * C_slow
  )
})

flux_df <- do.call(rbind, flux_list)

flux_summary <- flux_df %>%
  summarise(across(everything(),
                   list(
                     median = ~median(.),
                     lci = ~quantile(., 0.05),
                     uci = ~quantile(., 0.95)
                   )
  )) %>%
  pivot_longer(everything(),
               names_to = "Metric",
               values_to = "Value")


# ============================================================
# Weighted C and N age density distributions (RECONSTRUCTED BULK)
# ============================================================

# 1. Final stocks for best parameter set
best_run <- out_best_MCMC
final_yr_idx <- nrow(best_run)

# Best-fit C pools
best_C_stocks <- c(
  fast  = best_run$Ct_fast[final_yr_idx],
  inter = best_run$Ct_inter[final_yr_idx],
  slow  = best_run$Ct_slow[final_yr_idx]
)

# observed final C and N stocks

final_Cobs <- c(
  fast  = tail(Cobs_fast$Ct_fast, 1),
  inter = tail(Cobs_inter$Ct_inter, 1),
  slow  = tail(Cobs_slow$Ct_slow, 1)
)

final_Nobs <- c(
  fast  = tail(Nobs_fast$Nt_fast, 1),
  inter = tail(Nobs_inter$Nt_inter, 1),
  slow  = tail(Nobs_slow$Nt_slow, 1)
)

# ============================================================
# Pool-specific CN ratios (from observed decomposed pools)
# ============================================================

CN_bulk_final <- tail(bulk$CN_mean, 1)
CN_fast_final <- final_Cobs[1]  / final_Nobs[1]
CN_inter_final <- final_Cobs[2]  / final_Nobs[2]
CN_slow_final <- final_Cobs[3]  / final_Nobs[3]

# ============================================================
# 2. BEST-FIT density-weighted system
# ============================================================

dens_best_weighted <- dens_best %>%
  mutate(
    
    # density x Carbon stocks 
    C_fast_stock  = fast  * best_C_stocks["fast"],
    C_inter_stock = inter * best_C_stocks["inter"],
    C_slow_stock  = slow  * best_C_stocks["slow"],
    C_system_stock = C_fast_stock + C_inter_stock + C_slow_stock
  ) %>%
  
  mutate(
    
    # density x Nitrogen stocks
    N_fast_stock  = fast * best_C_stocks["fast"]* 1/CN_fast_final,
    N_inter_stock = inter * best_C_stocks["inter"] * 1/CN_inter_final,
    N_slow_stock  = slow * best_C_stocks["slow"] * 1/CN_slow_final,
    N_system_stock = N_fast_stock + N_inter_stock + N_slow_stock
  )

# ============================================================
# 3. MCMC ENSEMBLE (mass-conserving N)
# ============================================================

stock_vars <- c(
  "C_fast_stock", "C_inter_stock", "C_slow_stock", "C_system_stock",
  "N_fast_stock", "N_inter_stock", "N_slow_stock", "N_system_stock"
)

mcmc_stock_list <- lapply(seq_along(runs), function(i) {
  
  r <- runs[[i]]
  d <- dens_list[[i]]
  f_idx <- nrow(r)
  
  # Carbon pools
  cf <- d$fast  * r$Ct_fast[f_idx]
  ci <- d$inter * r$Ct_inter[f_idx]
  cs <- d$slow  * r$Ct_slow[f_idx]
  
  C_sys <- cf + ci + cs
  
  # ============================================================
  # using final observed CN to estimate N densities
  # ============================================================
  
  Nf <- d$fast * r$Ct_fast[f_idx] * 1/CN_fast_final
  Ni <- d$inter * r$Ct_inter[f_idx] * 1/CN_inter_final
  Ns <- d$slow * r$Ct_slow[f_idx] * 1/CN_slow_final
  
  N_sys <- Nf + Ni + Ns
  
  data.frame(
    age = d$age,
    
    # Carbon
    C_fast_stock  = cf,
    C_inter_stock = ci,
    C_slow_stock  = cs,
    C_system_stock = C_sys,
    
    # Nitrogen 
    N_fast_stock  = Nf,
    N_inter_stock = Ni,
    N_slow_stock  = Ns,
    N_system_stock = N_sys
  )
})

# ============================================================
# 4. Uncertainty envelopes
# ============================================================

ages <- dens_list[[1]]$age

unc_stocks_list <- lapply(stock_vars, function(v){
  
  mat <- sapply(mcmc_stock_list, function(x) x[[v]])
  
  data.frame(
    age = ages,
    variable = v,
    low  = apply(mat, 1, quantile, 0.025, na.rm = TRUE),
    high = apply(mat, 1, quantile, 0.975, na.rm = TRUE)
  )
})

unc_stocks_df <- do.call(rbind, unc_stocks_list)

unc_stocks_df <- unc_stocks_df %>%
  mutate(type = ifelse(grepl("^C_", variable), "C", "N"))

# ============================================================
# 5. Merge best-fit + uncertainty
# ============================================================

plot_df_final <- dens_best_weighted %>%
  select(age, contains("_stock")) %>%
  pivot_longer(-age, names_to = "variable", values_to = "best_val") %>%
  mutate(type = ifelse(grepl("^C_", variable), "C", "N")) %>%
  left_join(unc_stocks_df, by = c("age", "variable", "type"))


plot_df_final$variable <- recode(plot_df_final$variable,
                                 "C_system_stock" = "C system (reconstructed)",
                                 "N_system_stock" = "N system (reconstructed)",
                                 "C_fast_stock"   = "C fast",
                                 "C_inter_stock"  = "C intermediate",
                                 "C_slow_stock"   = "C slow",
                                 "N_fast_stock"   = "N fast",
                                 "N_inter_stock"  = "N intermediate",
                                 "N_slow_stock"   = "N slow"
)

# ============================================================
# 6. Summary statistics
# ============================================================

fast_ref <- summary_table_fmt$median[
  summary_table_fmt$variable == "fast"
]

slow_ref <- summary_table_fmt$median[
  summary_table_fmt$variable == "slow"
]

final_stats <- dens_best_weighted %>%
  select(age, contains("_stock")) %>%
  pivot_longer(-age, names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    total_stock = sum(value),
    mean_val   = sum(age * value) / sum(value),
    median_val = age[which.min(abs(cumsum(value) / sum(value) - 0.5))],
    frac_younger_fast = sum(value[age <= fast_ref]) / sum(value),
    frac_older_slow   = sum(value[age >= slow_ref]) / sum(value)
  )

final_stats$variable <- recode(final_stats$variable,
                               "C_system_stock" = "C system (reconstructed)",
                               "N_system_stock" = "N system (reconstructed)",
                               "C_fast_stock"   = "C fast",
                               "C_inter_stock"  = "C intermediate",
                               "C_slow_stock"   = "C slow",
                               "N_fast_stock"   = "N fast",
                               "N_inter_stock"  = "N intermediate",
                               "N_slow_stock"   = "N slow")

print(final_stats)
save(final_stats, file = file.path("mod_runs", "final_stats_weighted.Rdata"))

# ============================================================
# 7. PLOTTING 
# ============================================================

col_fast  <- "#1b9e77"
col_slow  <- "#BF40BF"
col_inter <- "#0000FF"

pool_colors <- c(
  "C fast" = col_fast,
  "C intermediate" = col_inter,
  "C slow" = col_slow,
  "N fast" = col_fast,
  "N intermediate" = col_inter,
  "N slow" = col_slow,
  "C system (reconstructed)" = "black",
  "N system (reconstructed)" = "black",
  "Mean" = "blue",
  "Median" = "red"
)

create_stock_age_plot <- function(data_subset, title_text, y_label) {
  
  stats_subset <- final_stats %>%
    filter(variable %in% unique(data_subset$variable))
  
  ggplot(data_subset, aes(x = age)) +
    
    geom_ribbon(aes(ymin = low, ymax = high, fill = variable), alpha = 0.4) +
    geom_line(aes(y = best_val, colour = variable), linewidth = 0.8) +
    
    geom_vline(data = stats_subset,
               aes(xintercept = mean_val, colour = "Mean"),
               linetype = "dashed") +
    
    geom_vline(data = stats_subset,
               aes(xintercept = median_val, colour = "Median"),
               linetype = "dotted") +
    
    facet_wrap(~variable, scales = "free_y", ncol = 2) +
    
    scale_x_log10(
      limits = c(1, 500),
      breaks = c(1, 10, 50, 100, 200, 500)
    ) +
    scale_fill_manual(values = pool_colors, guide = "none") +
    
    scale_colour_manual(values = pool_colors,
                        breaks = c("Mean", "Median"),
                        name = "Statistics") +
    
    labs(x = "Age (years, log-scale)",
         y = y_label,
         title = title_text) +
    
    theme_minimal(base_size = 14)
}

# ============================================================
# 8. Generate plots
# ============================================================

plot_C_weighted_final <- create_stock_age_plot(
  subset(plot_df_final, type == "C"),
  "Reconstructed Carbon Stock Distribution by Age",
  expression(Carbon~stocks~(g~m^{-2}))
)

plot_N_weighted_final <- create_stock_age_plot(
  subset(plot_df_final, type == "N"),
  "Reconstructed Nitrogen Stock Distribution by Age (Mass-conserving)",
  expression(Nitrogen~stocks~(g~m^{-2}))
)

# ============================================================
# 9. Save outputs
# ============================================================

output_dir <- "plots/allsites/4_pool"

ggsave(file.path(output_dir, "plot_C_weighted_final.png"),
       plot_C_weighted_final, width = 7, height = 5, dpi = 300, bg = "white")

ggsave(file.path(output_dir, "plot_N_weighted_final.png"),
       plot_N_weighted_final, width = 7, height = 5, dpi = 300, bg = "white")

saveRDS(plot_C_weighted_final,
        file.path(output_dir, "plot_C_weighted_final.rds"))

saveRDS(plot_N_weighted_final,
        file.path(output_dir, "plot_N_weighted_final.rds"))