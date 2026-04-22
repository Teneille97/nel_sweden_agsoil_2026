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

# import modeled fast pool
load("~/nel_sweden_agsoil_2026/mod_runs/results_allsites_325asslow.Rdata")

C14obs_fast <- data.frame(
  Year = results_325asslow$output$Year,
  C14t_fast = results_325asslow$output$C14t_fast
)

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
C0_inter <- mean(Cobs_inter[Cobs_inter$Year==1957,]$Ct_inter)
F0_400 <- mean(C14obs_slow[C14obs_slow$Year==1957,]$C14t_slow)
F0_bulk <- mean(C14obs_bulk[C14obs_bulk$Year==1957,]$C14t) 

# func to run mod
run_mod <- function(pars){ #kf, ki, ks, alpha 21, +F0b=Fi, I, alpha 32, alpha 41, alpha 42
  c14_atm <- BoundFc(Atm14C, format = "Delta14C")
  c14_initial <- ConstFc(
    values = c(F0_bulk + 10, F0_bulk + pars[5], F0_400, F0_bulk), 
    format = "Delta14C"
  )
  A4 <- diag(-c(pars[1:3],0.000000001))
  A4[2,1] <- pars[1]*pars[4]
  A4[3,2] <- pars[2]*pars[7]
  A4[4,1] <- pars[1]*pars[8]
  A4[4,2] <- pars[2]*(pars[9])
  mod<-GeneralModel_14(
    t = yr,
    A = A4,
    ivList = c(C0_fast, C0_inter, C0_400, 0),
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

inipars <- c(0.1,0.05,0.001, 0.05, -50, mean_C_inputs*0.5, 0.1, 0.005, 0.005) #kf, ki, ks, alpha 21, F0ib, I, alpha 32, alpha 41

# cost func
mc <- function(pars){
  out = run_mod(pars)
  Cost1 <- modCost(out, Cobs_bulk, x = "Year", err = "Ct_sd")
  Cost2 <- modCost(out, C14obs_bulk, x = "Year", cost = Cost1, err = "C14t_sd")
  Cost3 <- modCost(out, Cobs_fast, x = "Year", cost = Cost2, err = "Ct_fast_sd")
  Cost4 <- modCost(out, Cobs_inter, x = "Year", cost = Cost3, err = "Ct_inter_sd")
  Cost5 <- modCost(out, Cobs_slow, x = "Year", cost = Cost4, err = "Ct_slow_sd")
  Cost6 <- modCost(out, C14obs_fast, x = "Year", cost = Cost5) #maybe add uncertainty from MCMC of 325asslow mod?
  modCost(out, C14obs_slow, x = "Year", cost = Cost6, err = "C14t_slow_sd")
}

#fit mod
mFit <- modFit(
  f = mc,
  p = inipars,
  method = "Nelder-Mead",
  upper = c(0.5,0.2,0.005, 0.1,0, mean_C_inputs*2,1, 0.01, 0.01), #kf, ki, ks, alpha 21, F0ib, I, alpha 32, alpha 41, alpha 42
  lower = c(0.05,0.01,0.0001,0,-160,0,0, 0,0) 
)

bestpars <- mFit$par
out_best <- run_mod(bestpars)
Delta14Clabel <- expression(Delta^14*C)
stocks_label <- expression(C ~ (g ~ m^{-2}))

col_bulk <- "black"
col_fast <- "#1b9e77"
col_slow <- "#d95f02"
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
  geom_point(data = C14obs_fast, aes(Year, C14t_fast, colour = "Fast")) +
  
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
output_dir <- "plots/allsites/325asslow/nested"

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

results_nested<-list(
  pars = bestpars,
  output = out_best,
  plot_C = plot_C,
  plot_C14 = plot_C14
)

save(results_nested, file = file.path("mod_runs", "results_allsites_nested.Rdata"))
#load(here::here("mod_runs/results_allsites_nested.Rdata"))

