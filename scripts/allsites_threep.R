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
  

Cobs_bulk <- data.frame(Year = bulk$Year, Ct = bulk$C_stocks_mean, Ct_sd = bulk$C_stocks_sd)
C14obs_bulk <- data.frame(Year = bulk$Year, C14t = bulk$d14C_mean, C14t_sd = bulk$d14C_sd)
Cobs_325 <- data.frame(Year = pool_325$Year, Ct_325 = pool_325$C_stocks_mean, Ct_325_sd = pool_325$C_stocks_sd)
Cobs_slow <- data.frame(Year = pool_400$Year, Ct_slow = pool_400$C_stocks_mean, Ct_slow_sd = pool_400$C_stocks_sd)
C14obs_slow <- data.frame(Year = pool_400$Year, C14t_slow = pool_400$d14C_mean, C14t_slow_sd = pool_400$d14C_sd)
Cobs_fast <- data.frame(Year = Cobs_325$Year, Ct_fast = Cobs_bulk$Ct-Cobs_325$Ct_325, Ct_fast_sd = sqrt(bulk$C_stocks_sd^2 + pool_325$C_stocks_sd^2)) # C fast  = C bulk - C 325

Cobs_bulk$Ct_sd[is.na(Cobs_bulk$Ct_sd)] <- mean(Cobs_bulk$Ct_sd, na.rm = TRUE)
C14obs_bulk$C14t_sd[is.na(C14obs_bulk$C14t_sd)] <- mean(C14obs_bulk$C14t_sd, na.rm = TRUE)
Cobs_fast$Ct_fast_sd[is.na(Cobs_fast$Ct_fast_sd)] <- mean(Cobs_fast$Ct_fast_sd, na.rm = TRUE)
C14obs_slow$C14t_slow_sd[is.na(C14obs_slow$C14t_slow_sd)] <- mean(C14obs_slow$C14t_slow_sd, na.rm = TRUE)

yr <- seq(1957, 2019, by = 1/12)
C0_bulk <- mean(Cobs_bulk[Cobs_bulk$Year==1957,]$Ct)
C0_fast <- mean(Cobs_fast[Cobs_fast$Year==1957,]$Ct_fast)
C0_400 <- mean(Cobs_slow[Cobs_slow$Year==1957,]$Ct_slow)
F0_400 <- mean(C14obs_slow[C14obs_slow$Year==1957,]$C14t_slow)
F0_bulk <- mean(C14obs_bulk[C14obs_bulk$Year==1957,]$C14t) 

run_mod <- function(pars){
  mod<-TwopSeriesModel14(
    t = yr,
    ks = pars[1:2],
    C0 = c(C0_fast, C0_400),
    F0_Delta14C = c(F0_bulk * pars[4], F0_400),
    In = pars[5],
    a21 = pars[1] * pars[3],
    inputFc = Atm14C
  )

  Ct_pools <- getC(mod)
  C14_pools <- getF14(mod)
  C14t <- getF14C(mod)
    
  mod_results<-data.frame(
    Year = yr,
    Ct = rowSums(Ct_pools),
    C14t = C14t,
    Ct_fast = Ct_pools[,1],
    Ct_slow = Ct_pools[,2],
    C14t_fast = C14_pools[,1],
    C14t_slow = C14_pools[,2]
  )
  
  return(mod_results)
}

inipars <- c(0.5, 0.05, 0.1, 0.9, mean_C_inputs*0.5) #kf, ks, alpha 21, F0fb, I
  
mc <- function(pars){
  out <- run_mod(pars)
  Cost1 <- modCost(model = out, obs = Cobs_bulk, x = "Year", err = "Ct_sd", weight = "std")
  Cost2 <- modCost(model = out, obs = C14obs_bulk, x = "Year", cost = Cost1, err = "C14t_sd", weight = "std")
  Cost3 <- modCost(model = out, obs = Cobs_fast, x = "Year", cost = Cost2, err = "Ct_fast_sd", weight = "std")
  return(modCost(model = out, obs = C14obs_slow, x = "Year", cost = Cost3, err = "C14t_slow_sd", weight = "std"))
}

mFit <- modFit(
  f = mc,
  p = inipars,
  method = "Nelder-Mead",
  upper = c(2,0.2, 1,3, mean_C_inputs*2), #kf, ki, alpha 21, F0fb, I
  lower = c(0.1,0.01,0,0,0)
)
  
bestpars <- mFit$par
out_best <- run_mod(bestpars)
Delta14Clabel <- expression(Delta^14*C)
stocks_label <- expression(C ~ (g ~ m^{-2}))
  
col_bulk <- "black"
col_fast <- "#1b9e77"
col_slow <- "#d95f02"

plot_C <- ggplot() +
  geom_line(data = out_best, aes(Year, Ct, colour = "Bulk"), linewidth = 1) +
  geom_line(data = out_best, aes(Year, Ct_fast, colour = "Fast"), linetype = "dashed") +
  geom_line(data = out_best, aes(Year, Ct_slow, colour = "Slow"), linetype = "dotted") +

  geom_point(data = Cobs_bulk, aes(Year, Ct, colour = "Bulk")) +
  geom_point(data = Cobs_fast, aes(Year, Ct_fast, colour = "Fast")) +
  geom_point(data = Cobs_slow, aes(Year, Ct_slow, colour = "Slow")) +
    
  scale_colour_manual(
    name = "Pool",
    values = c("Bulk" = col_bulk,
                 "Fast" = col_fast,
                 "Slow" = col_slow)
    ) +
    
    ggtitle("C stocks") +
    ylab(stocks_label) +
    xlab("Year") +
    theme_minimal()
  
  
plot_C14 <- ggplot() +
    
    geom_line(data = out_best, aes(Year, C14t, colour = "Bulk"), linewidth = 1) +
    geom_line(data = out_best, aes(Year, C14t_fast, colour = "Fast"), linetype = "dashed") +
    geom_line(data = out_best, aes(Year, C14t_slow, colour = "Slow"), linetype = "dotted") +

    geom_point(data = C14obs_bulk, aes(Year, C14t, colour = "Bulk")) +
    geom_point(data = C14obs_slow, aes(Year, C14t_slow, colour = "Slow")) +
    
    scale_colour_manual(
      name = "Pool",
      values = c("Bulk" = col_bulk,
                 "Fast" = col_fast,
                 "Slow" = col_slow)
    ) +
    
    ggtitle(expression(Delta^14*C)) +
    ylab(Delta14Clabel) +
    xlab("Year") +
    theme_minimal()


output_dir <- "plots/allsites/3_pool"
for (site in names(results)) {
  
  res <- results_3p[[site]]
  
  # file names
  file_C <- file.path(output_dir, paste0(site, "_C.png"))
  file_C14 <- file.path(output_dir, paste0(site, "_C14.png"))
  
  # save plots
  ggsave(filename = file_C, plot = res$plot_C,
         width = 7, height = 5, dpi = 300, bg = "white")
  
  ggsave(filename = file_C14, plot = res$plot_C14,
         width = 7, height = 5, dpi = 300, bg = "white")
}


for (name in names(results_3p)) {
  saveRDS(results_3p[[name]]$plot_C,
          file = file.path(output_dir, paste0(name, "_C.rds")))
  
  saveRDS(results_3p[[name]]$plot_C14,
          file = file.path(output_dir, paste0(name, "_14C.rds")))
}

save(results_3p, file = file.path("mod_runs", "results_allsites_3pool.Rdata"))
