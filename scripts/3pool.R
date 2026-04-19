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

bulk_density <- c(
  "M2-1957" = 1.51,
  "M4-1957" = 1.72,
  "M6-1957" = 1.44,
  "R94-1966" = 1.37
)

F0_site_400_vals <- c(
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
  pool_400<- subset(site_df, Temperature...C. == "400")
  
  ## observed data -----------------------------------------------------------
  
  Cobs_bulk <- data.frame(Year = bulk$Year, Ct = bulk$C_stocks_gm2)
  C14obs_bulk <- data.frame(Year = bulk$Year, C14t = bulk$d14C)
  
  Cobs_325 <- data.frame(Year = pool_325$Year, Ct_325 = pool_325$C_stocks_gm2)
  
  Cobs_slow <- data.frame(Year = pool_400$Year, Ct_slow = pool_400$C_stocks_gm2)
  C14obs_slow <- data.frame(Year = pool_400$Year, C14t_slow = pool_400$d14C)
  
  Cobs_fast <- data.frame(Year = pool_400$Year, Ct_fast = Cobs_bulk$Ct-Cobs_325$Ct_325) # C fast  =C bulk - C 325
  C14obs_fast <- data.frame(
    Year = results[[site_name]]$output$Year,
    C14t_fast = results[[site_name]]$output$C14t_fast
  )
  
  ## time + initial parameters ----------------------------------------------
  
  yr <- seq(1957, 2019, by = 1/12)
  inipars <- c(0.5,0.05,0.001, 0.1, 1.5, mean_C_inputs*0.5, 0.1) #kf, ki, ks, alpha 21, F0ib, I, alpha 32
  C0_bulk <- Cobs_bulk[1,2]
  C0_fast <- Cobs_fast[1,2]
  C0_400 <- Cobs_slow[1,2]
  F0_400 <- F0_site_400_vals[site_name]
  F0_bulk <- C14obs_bulk[1,2]
  F0_fast <- C14obs_fast[1,2]
  
  ## model function ----------------------------------------------------------
  
  mf <- function(pars){
    
    md <- ThreepSeriesModel14(
      t = yr,
      ks = pars[1:3],
      C0 = c(C0_fast, abs(C0_bulk-C0_fast-C0_400), C0_400),
      F0_Delta14C = c(F0_fast, F0_bulk * pars[5], F0_400),
      In = pars[6],
      a21 = pars[1] * pars[4],
      a32 = pars[2] * pars[7],
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
      Ct_inter = Ct_pools[,2],
      Ct_slow = Ct_pools[,3],
      C14t_fast = C14_pools[,1],
      C14t_inter = C14_pools[,2],
      C14t_slow = C14_pools[,3]
    )
  }
  
  ## cost function -----------------------------------------------------------
  
  mc <- function(pars){
    
    out <- mf(pars)
    
    Cost1 <- modCost(out, Cobs_bulk, x = "Year")
    Cost2 <- modCost(out, C14obs_bulk, x = "Year", cost = Cost1)
    Cost3 <- modCost(out, Cobs_fast, x = "Year", cost = Cost2)
    Cost4 <- modCost(out, C14obs_fast, x = "Year", cost = Cost3)
    Cost5 <- modCost(out, Cobs_slow, x = "Year", cost = Cost4)
    modCost(out, C14obs_slow, x = "Year", cost = Cost5)
  }
  
  ## optimize ---------------------------------------------------------------
  
  mFit <- modFit(
    f = mc,
    p = inipars,
    method = "Nelder-Mead",
    upper = c(2,0.2,0.005, 1,3, mean_C_inputs*1.5,1), #kf, ki, ks, alpha 21, F0ib, I, alpha 32
    lower = c(0.1,0.01,0.0001,0,0,0,0)
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
  col_inter <- "#0000FF"
  
  ## Carbon stocks ------------------------------------------------------------
  
  plot_C <- ggplot() +
    
    # model lines
    geom_line(data = out_best, aes(Year, Ct, colour = "Bulk"), linewidth = 1) +
    geom_line(data = out_best, aes(Year, Ct_fast, colour = "Fast"), linetype = "dashed") +
    geom_line(data = out_best, aes(Year, Ct_slow, colour = "Slow"), linetype = "dotted") +
    geom_line(data = out_best, aes(Year, Ct_inter, colour = "Inter"), linetype = "dotted") +
    
    # observations
    geom_point(data = Cobs_bulk, aes(Year, Ct, colour = "Bulk")) +
    geom_point(data = Cobs_fast, aes(Year, Ct_fast, colour = "Fast")) +
    geom_point(data = Cobs_slow, aes(Year, Ct_slow, colour = "Slow")) +

    scale_colour_manual(
      name = "Pool",
      values = c("Bulk" = col_bulk,
                 "Fast" = col_fast,
                 "Slow" = col_slow,
                 "Inter" = col_inter)
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
    geom_line(data = out_best, aes(Year, C14t_inter, colour = "Inter"), linetype = "dotted") +
    # observations
    geom_point(data = C14obs_bulk, aes(Year, C14t, colour = "Bulk")) +
    geom_point(data = C14obs_slow, aes(Year, C14t_slow, colour = "Slow")) +

    scale_colour_manual(
      name = "Pool",
      values = c("Bulk" = col_bulk,
                 "Fast" = col_fast,
                 "Slow" = col_slow,
                 "Inter" = col_inter)
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

results_3p <- lapply(sites, run_model_site)
names(results_3p) <- sites

## save all plots to local dir ----------------------------------------------

output_dir <- "plots/325asslow/3_pool"
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

## save all plots as rds for loading to qmd ---------------------------------
for (name in names(results_3p)) {
  saveRDS(results_3p[[name]]$plot_C,
          file = file.path(output_dir, paste0(name, "_C.rds")))
  
  saveRDS(results_3p[[name]]$plot_C14,
          file = file.path(output_dir, paste0(name, "_14C.rds")))
}

save(results_3p, file = file.path("mod_runs", "results_3pool.Rdata"))
