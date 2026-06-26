#tnel edits
##LTE

library(here)

Long=read.csv2(here("csv_files", "LTEnitrogen2.csv"))
head(Long)
unique(Long$LTE)

M2=subset(Long,Long$LTE=="M2-1957")
M4=subset(Long,Long$LTE=="M4-1957")
M6=subset(Long,Long$LTE=="M6-1957")
R94=subset(Long,Long$LTE=="R94-1966")

M2_Soil=subset(M2,M2$Temperature...C.=="Soil")
M2_325=subset(M2,M2$Temperature...C.=="325")
M2_400=subset(M2,M2$Temperature...C.=="400")

M4_Soil=subset(M4,M4$Temperature...C.=="Soil")
M4_325=subset(M4,M4$Temperature...C.=="325")
M4_400=subset(M4,M4$Temperature...C.=="400")

M6_Soil=subset(M6,M6$Temperature...C.=="Soil")
M6_325=subset(M6,M6$Temperature...C.=="325")
M6_400=subset(M6,M6$Temperature...C.=="400")

R94_Soil=subset(R94,R94$Temperature...C.=="Soil")
R94_175=subset(R94,R94$Temperature...C.=="175")
R94_250=subset(R94,R94$Temperature...C.=="250")
R94_325=subset(R94,R94$Temperature...C.=="325")
R94_400=subset(R94,R94$Temperature...C.=="400")
R94_445=subset(R94,R94$Temperature...C.=="445")
R94_500=subset(R94,R94$Temperature...C.=="500")

# Function to calculate summary statistics
summarise_d14C <- function(df, site_name, fraction_name) {
  data.frame(
    Site = site_name,
    Fraction = fraction_name,
    n = nrow(df),
    Mean_d14C = mean(as.numeric(df$d14C), na.rm = TRUE),
    Median_d14C = median(as.numeric(df$d14C), na.rm = TRUE)
  )
}

# Create summary table
d14C_summary <- rbind(
  summarise_d14C(M2_Soil, "M2", "Soil"),
  summarise_d14C(M2_325,  "M2", "325"),
  summarise_d14C(M2_400,  "M2", "400"),
  
  summarise_d14C(M4_Soil, "M4", "Soil"),
  summarise_d14C(M4_325,  "M4", "325"),
  summarise_d14C(M4_400,  "M4", "400"),
  
  summarise_d14C(M6_Soil, "M6", "Soil"),
  summarise_d14C(M6_325,  "M6", "325"),
  summarise_d14C(M6_400,  "M6", "400"),
  
  summarise_d14C(R94_Soil, "R94", "Soil"),
  summarise_d14C(R94_325,  "R94", "325"),
  summarise_d14C(R94_400,  "R94", "400")
)

# Round for presentation
d14C_summary$Mean_d14C <- round(d14C_summary$Mean_d14C, 1)
d14C_summary$Median_d14C <- round(d14C_summary$Median_d14C, 1)

d14C_summary

#compact summary
library(tidyr)

d14C_summary %>%
  select(Site, Fraction, Mean_d14C, Median_d14C) %>%
  pivot_wider(
    names_from = Fraction,
    values_from = c(Mean_d14C, Median_d14C)
  )