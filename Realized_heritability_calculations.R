### Calculate the realized heritability according to the Breeder's equation ----
rm(list = ls())
windowsFonts(my = windowsFont('Calibri Light'))
library(ggplot2)
library(stringr)
library(data.table)

setwd("/path/to/your/own/working/directory/")
data <- read.table("data_all_years.txt")
data <-  as.data.table(data)

Short_1_2016 <- data[Population == "Short_plants_1" & year == "2016",]
Short_1_2017 <- data[Population == "Short_plants_1" & year == "2017",]

calculated_realized_herit <- function(Data_Gen1,
                                      Data_Gen2){
  mean_Gen_1 <- mean(Data_Gen1$PlantHeight, na.rm = TRUE)
  mean_Gen_2 <- mean(Data_Gen2$PlantHeight, na.rm = TRUE)
  Selected_prop <- quantile(Data_Gen1$PlantHeight, 0.05, na.rm = TRUE)
  Selection_Differential <- Selected_prop - mean_Gen_1
  Response_to_Selection <- mean_Gen_2 - mean_Gen_1
  herit_BS <- as.numeric(Response_to_Selection)/as.numeric(Selection_Differential)
  return(herit_BS)
}
calculated_realized_herit(Data_Gen1 = Short_1_2016,
                          Data_Gen2 = Short_1_2017)
