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

output_dir <- "plots/allsites/4_pool"
file_convergence_plot <- file.path(output_dir, "convergence_plot.png")

png(file_convergence_plot, width = 10, height = 8, units = "in", res = 300)
plot(MCMC)
dev.off()

saveRDS(convergence_plot,
        file = file.path(output_dir, "convergence_plot.rds"))

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
mc(mFit$par)
percentage_accepted<-(100*MCMC$naccepted)/50000
MCMC_bestpars<-MCMC$bestpar
cost_modfit <- mc(bestpars)$model 
cost_mcmc   <- mc(MCMC_bestpars)$model  

# plots showing 95% credible interval of model outputs given the posterior parameter distribution
# i.e. posterior predictive uncertainty conditional on sampled MCMC chains

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


# Ages and transit times


## --- Calculate fluxes and stocks for all pools 

pf_4ps <- MCMC$pars %>%
  as.data.frame()

colnames(pf_4ps) <- c("kf","ki","ks","alpha_fi","x","alpha_is","alpha_fl","alpha_il")

pf_4ps_flux <- pf_4ps %>%
  mutate(
    C_fast  = sample(unc_C$Mean[unc_C$Pool == "Fast"],  n(), replace = TRUE),
    C_inter = sample(unc_C$Mean[unc_C$Pool == "Inter"], n(), replace = TRUE),
    C_slow  = sample(unc_C$Mean[unc_C$Pool == "Slow"],  n(), replace = TRUE),
    
    # Fast pool
    out_fast = (1 - alpha_fi - alpha_fl) * kf * C_fast,
    trans_fast_to_inter = alpha_fi * kf * C_fast,
    trans_fast_to_loss  = alpha_fl * kf * C_fast,
    
    # Intermediate pool
    out_inter = (1 - alpha_is - alpha_il) * ki * C_inter,
    trans_inter_to_slow = alpha_is * ki * C_inter,
    trans_inter_to_loss = alpha_il * ki * C_inter,
    
    # Slow pool
    out_slow = ks * C_slow
  )

flux_summary <- pf_4ps_flux %>%
  summarise(across(starts_with("out") | starts_with("trans"),
                   list(
                     median = ~round(median(.), 2),
                     lci    = ~round(quantile(., 0.05), 2),
                     uci    = ~round(quantile(., 0.95), 2)
                   ))) %>%
  pivot_longer(everything(),
               names_to = "Metric",
               values_to = "Value")

print(flux_summary)

# ages and transit times
get_age_tt <- function(pars){
  
  # Build A matrix (same structure as in run_mod)
  A4 <- diag(-c(pars[1:3], 1e-9))
  A4[2,1] <- pars[1]*pars[4]
  A4[3,2] <- pars[2]*pars[6]
  
  # Build model (steady-state approximation)
  mod <- GeneralModel(
    t = yr,
    A = A4,
    ivList = c(C0_fast, abs(C0_bulk-C0_fast-C0_400), C0_400),
    inputFluxes = c(mean_C_inputs,0,0)
  )
  
  # Ages
  system_age  <- systemAge(mod)$meanSystemAge
  pool_age    <- systemAge(mod$meanPoolAge)
  
  # Transit time
  transit_time <- getTransitTime(mod)
  
  return(c(
    system_age = system_age,
    fast_age   = pool_age[1],
    inter_age  = pool_age[2],
    slow_age   = pool_age[3],
    transit_time = transit_time
  ))
}
set.seed(1)

pars_sub <- MCMC$pars[sample(1:nrow(MCMC$pars), 3000), ]

age_tt_mat <- t(apply(pars_sub, 1, get_age_tt))
age_tt_df  <- as.data.frame(age_tt_mat)

age_tt_long <- age_tt_df %>%
  pivot_longer(cols = everything(),
               names_to = "metric",
               values_to = "value")

density_age_tt <- ggplot(age_tt_long, aes(x = value)) +
  geom_density(fill = "orange", alpha = 0.5) +
  facet_wrap(~ metric, scales = "free") +
  theme_minimal() +
  labs(title = "System Age, Pool Ages, and Transit Time Distributions",
       x = "Years",
       y = "Density")

ggsave(file.path(output_dir, "age_tt_density.png"),
       density_age_tt,
       width = 10, height = 8, dpi = 300)

saveRDS(density_age_tt,
        file = file.path(output_dir, "age_tt_density.rds"))

sum_age_tt <- function(x){
  c(
    mean   = mean(x),
    sd     = sd(x),
    median = median(x),
    q05    = quantile(x, 0.05),
    q25    = quantile(x, 0.25),
    q75    = quantile(x, 0.75),
    q95    = quantile(x, 0.95)
  )
}

summary_table <- as.data.frame(t(sapply(age_tt_df, sum_age_tt)))
summary_table$metric <- rownames(summary_table)

summary_table <- summary_table %>%
  select(metric, everything())

write.csv(summary_table,
          file = file.path(output_dir, "age_tt_summary.csv"),
          row.names = FALSE)

summary_table