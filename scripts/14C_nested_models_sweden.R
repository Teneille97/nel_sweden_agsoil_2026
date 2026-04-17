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
  
  bulk <- subset(site_df, Temperature...C. == "325")
  pool_400 <- subset(site_df, Temperature...C. == "400")
  
  ## observed data -----------------------------------------------------------
  
  Cobs_bulk <- data.frame(Year = bulk$Year, Ct = bulk$C_stocks_gm2)
  C14obs_bulk <- data.frame(Year = bulk$Year, C14t = bulk$d14C)
  
  Cobs_400 <- data.frame(Year = pool_400$Year, Ct_slow = pool_400$C_stocks_gm2)
  C14obs_400 <- data.frame(Year = pool_400$Year, C14t_slow = pool_400$d14C)
  
  Cobs_fast <- data.frame(Year = pool_400$Year, Ct_fast = Cobs_bulk$Ct-Cobs_400$Ct_slow) #C fast = C 325 - C 400
  
  
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
    Cost3 <- modCost(out, Cobs_fast, x = "Year", cost = Cost2)
    Cost4 <- modCost(out, Cobs_400, x = "Year", cost = Cost3)
    
    modCost(out, C14obs_400, x = "Year", cost = Cost4)
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
    geom_point(data = Cobs_fast, aes(Year, Ct_fast, colour = "Fast")) +
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

output_dir <- "plots/325asbulk"
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
