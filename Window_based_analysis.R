# FDRfS significance thresholds based on Wstat based on the FST OD and SD ------
## Installing the packages -----------------------------------------------------
cat("Package installation starts.","\n")
install.packages("data.table",
                 repos='http://cran.rstudio.com/')
install.packages("parallel",
                 repos='http://cran.rstudio.com/')
install.packages("stringr",
                 repos='http://cran.rstudio.com/')
install.packages("doMC",
                 repos='http://cran.rstudio.com/')
install.packages("tidyverse",
                 repos='http://cran.rstudio.com/')
install.packages("gridExtra",
                 repos='http://cran.rstudio.com/')
### Loading of the packages ----------------------------------------------------
library(data.table)
library(stringr)
library(foreach)
library(parallel)
library(tidyverse)
library(doMC)
library(gridExtra)
cat("Packages are loaded.","\n")
cat("Core registration starts.","\n")
cores <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK'))
cores <- detectCores(all.tests = FALSE, logical = TRUE)
registerDoMC(cores)
### Load the data --------------------------------------------------------------
result_directory <- "/path/to/your/own/working/directory/"
input_directory <- "/path/to/your/own/working/directory/"
setwd(input_directory)
FST_values_od_cor <- fread("2022_02_16_GB1005_corrected_dt_Fst_values.txt",
                           skip = 1)
colnames_for_dt <- read.table("2022_02_16_GB1005_corrected_dt_Fst_values.txt",
                              nrows = 1, header = TRUE)
FST_values_od_cor <- FST_values_od_cor[,2:12]
colnames(FST_values_od_cor) <- colnames(colnames_for_dt)
FST_values_od_cor <- as.data.frame(FST_values_od_cor)
source("/usr/users/mtost/Shoepeg_resubmission_new_analysis/new_functions/2022_Wstat_calc_changed_by_Mila.R")
#### Prepare the data set ------------------------------------------------------
FST_values_od_cor <- as.data.table(FST_values_od_cor)
FST_values_od_cor[, FstSum_OD := FST_value_opposite_dir1 + FST_value_opposite_dir2]
FST_values_od_cor[, FstSum_SD := FST_value_same_dir1 + FST_value_same_dir2]
#### Divide the data set into the different chromosomes ------------------------
Fst_Chr_1 <-  FST_values_od_cor[which(data$Chromosome=="chr1"),]
Fst_Chr_2 <-  FST_values_od_cor[which(data$Chromosome=="chr2"),]
Fst_Chr_3 <-  FST_values_od_cor[which(data$Chromosome=="chr3"),]
Fst_Chr_4 <-  FST_values_od_cor[which(data$Chromosome=="chr4"),]
Fst_Chr_5 <-  FST_values_od_cor[which(data$Chromosome=="chr5"),]
Fst_Chr_6 <-  FST_values_od_cor[which(data$Chromosome=="chr6"),]
Fst_Chr_7 <-  FST_values_od_cor[which(data$Chromosome=="chr7"),]
Fst_Chr_8 <-  FST_values_od_cor[which(data$Chromosome=="chr8"),]
Fst_Chr_9 <-  FST_values_od_cor[which(data$Chromosome=="chr9"),]
Fst_Chr_10 <-  FST_values_od_cor[which(data$Chromosome=="chr10"),]
# Smoothness = 2000, method = 3 ------------------------------------------------
#### Load the parameters -------------------------------------------------------
ncol_FstSum_OD = which(colnames(FST_values_od_cor) == "FstSum_OD")
ncol_FstSum_SD = which(colnames(FST_values_od_cor) == "FstSum_SD")
ncol_Position = 2
smoothness = 2000
method_GenWin = 3
### Define window boundaries per chromosome ------------------------------------
run_Spline_Analyze_per_chr <- function(data,
                                       Chromosome,
                                       smoothness){
  cat("The window boundaries for", Chromosome, "are getting calculated!","\n")
  data <- as.data.frame(data)
  # Define window boundaries per chromosome
  win_FST_SD <- splineAnalyze_Mila(data[ ,ncol_FstSum_OD], data[ ,ncol_Position],
                                    smoothness = smoothness)
  win_FST_OD <- splineAnalyze_Mila(data[ ,ncol_FstSum_SD], data[ ,ncol_Position],
                                    smoothness = smoothness)
  # Prepare the data table
  win_FST_OD_dt <- win_FST_OD$windowData
  win_FST_SD_dt <- win_FST_SD$windowData
  win_FST_SD_dt_new <- cbind(rep(Chromosome, nrow(win_FST_SD_dt)),win_FST_SD_dt)
  win_FST_OD_dt_new <- cbind(rep(Chromosome, nrow(win_FST_OD_dt)),win_FST_OD_dt)
  colnames(win_FST_SD_dt_new) <- c("Chromosome", colnames(win_FST_SD_dt))
  colnames(win_FST_OD_dt_new) <- c("Chromosome", colnames(win_FST_OD_dt))
  return(list("Window_boundaries_SD" = win_FST_SD_dt_new,
              "Window_boundaries_OD" = win_FST_OD_dt_new))
}
#### Run window boundaries analysis --------------------------------------------
win_FST_allCalc_Chr1 <- run_Spline_Analyze_per_chr(data = Fst_Chr_1,
                                                   Chromosome = "1",
                                                   smoothness = smoothness)
win_FST_allCalc_Chr2 <- run_Spline_Analyze_per_chr(data = Fst_Chr_2,
                                                   Chromosome = "2",
                                                   smoothness = smoothness)
win_FST_allCalc_Chr3 <- run_Spline_Analyze_per_chr(data = Fst_Chr_3,
                                                   Chromosome = "3",
                                                   smoothness = smoothness)
win_FST_allCalc_Chr4 <- run_Spline_Analyze_per_chr(data = Fst_Chr_4,
                                                   Chromosome = "4",
                                                   smoothness = smoothness)
win_FST_allCalc_Chr5 <- run_Spline_Analyze_per_chr(data = Fst_Chr_5,
                                                   Chromosome = "5",
                                                   smoothness = smoothness)
win_FST_allCalc_Chr6 <- run_Spline_Analyze_per_chr(data = Fst_Chr_6,
                                                   Chromosome = "6",
                                                   smoothness = smoothness)
win_FST_allCalc_Chr7 <- run_Spline_Analyze_per_chr(data = Fst_Chr_7,
                                                   Chromosome = "7",
                                                   smoothness = smoothness)
win_FST_allCalc_Chr8 <- run_Spline_Analyze_per_chr(data = Fst_Chr_8,
                                                   Chromosome = "8",
                                                   smoothness = smoothness)
win_FST_allCalc_Chr9 <- run_Spline_Analyze_per_chr(data = Fst_Chr_9,
                                                   Chromosome = "9",
                                                   smoothness = smoothness)
win_FST_allCalc_Chr10 <- run_Spline_Analyze_per_chr(data = Fst_Chr_10,
                                                    Chromosome = "10",
                                                    smoothness = smoothness)
cat("The window boundaries were calculated for all chromosomes","\n","\n")
#### Combine data from all chromosomes -----------------------------------------
win_FST_SD_allChr <- rbind(win_FST_allCalc_Chr1$Window_boundaries_SD,
                            win_FST_allCalc_Chr2$Window_boundaries_SD,
                            win_FST_allCalc_Chr3$Window_boundaries_SD,
                            win_FST_allCalc_Chr4$Window_boundaries_SD,
                            win_FST_allCalc_Chr5$Window_boundaries_SD,
                            win_FST_allCalc_Chr6$Window_boundaries_SD,
                            win_FST_allCalc_Chr7$Window_boundaries_SD,
                            win_FST_allCalc_Chr8$Window_boundaries_SD,
                            win_FST_allCalc_Chr9$Window_boundaries_SD,
                            win_FST_allCalc_Chr10$Window_boundaries_SD)
win_FST_OD_allChr <- rbind(win_FST_allCalc_Chr1$Window_boundaries_OD,
                            win_FST_allCalc_Chr2$Window_boundaries_OD,
                            win_FST_allCalc_Chr3$Window_boundaries_OD,
                            win_FST_allCalc_Chr4$Window_boundaries_OD,
                            win_FST_allCalc_Chr5$Window_boundaries_OD,
                            win_FST_allCalc_Chr6$Window_boundaries_OD,
                            win_FST_allCalc_Chr7$Window_boundaries_OD,
                            win_FST_allCalc_Chr8$Window_boundaries_OD,
                            win_FST_allCalc_Chr9$Window_boundaries_OD,
                            win_FST_allCalc_Chr10$Window_boundaries_OD)
### Save the data sets ---------------------------------------------------------
write.table(win_FST_SD_allChr, paste0(result_directory, Sys.Date(),"_Method_",method_GenWin,"_Smoothness_",smoothness,"_window_boundaries_based_on_FstSum_SD.txt"),
            sep = "\t",
            row.names = FALSE,
            quote = TRUE)
write.table(win_FST_OD_allChr, paste0(result_directory, Sys.Date(),"_Method_",method_GenWin,"_Smoothness_",smoothness,"_window_boundaries_based_on_FstSum_OD.txt"),
            sep = "\t",
            row.names = FALSE,
            quote = TRUE)
cat("The data sets for all chromsomes were created","\n","\n")
q()
