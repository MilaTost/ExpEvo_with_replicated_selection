######################## Selection_signature_mapping ###########################
# load required packages
library(stringr)
library(data.table)
library(doMC)
library(vcfR)
library(parallel)

cores <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK'))
cores <- detectCores(all.tests = FALSE, logical = TRUE)
registerDoMC(cores)
time_0 <- Sys.time()
cat("The data set loading starts")
# vcf_S4 <- read.vcfR("YOUR_DATA_after_filtering.vcf.gz", verbose = FALSE)
#setwd("/path/to/your/own/working/directory/")
vcf_S4 <- read.vcfR("/usr/users/mtost/wd_GBeasy_rerun/Results/2022_GB10_Shoepeg_after_filtering.vcf.gz", verbose = FALSE)
setwd("/usr/users/mtost/GB_easy_results_analysis_wd/final_data_analysis_paper/")
cat("The data set is loaded.","\n")
time_1 <- Sys.time()
print(time_1 - time_0)

########### Selection signature mapping analysis -------------------------------
### mean Fst value calculation -------------------------------------------------
### This is a function that will calculate Fst over two or three             ###
### populations.  The inputs are vectors of allele frequencies at each locus ###
### for each population.  The function is explained in the documentation     ###
### below.                                                                   ###
### Adapted from Beissinger (2012)

### The function below will calculate Fst between population which got selected
### in the same direction and populations which got selected in opposite direct-
### ions.
### Formula:   Should the simplest formula for Fst be used?  If T, the
###            formula used is Fst = s^2 / (mean(p) * (1- mean(p)) .
###            if false, a better formula which corrects for small number
###            of populations is used: Weir and Cockerham, 1984

###### Fst value calculation ---------------------------------------------------
### References: Weir and Cockerham, 1984, Estimating F-statistics for the 
### analysis of population structure. Evoluation 38(6): 1358-1390.

### Calculate Fst value between the different populations all togehter
### and separately. Furthermore this function calculates the FST Sum
### statistic between the subpopulations selected in opposite directions

# The data table with the FST values is called "FST_values_od_cor".
# This function calculates the FST between subpopulations selected 
# in the same and opposite directions, which are needed to calculate
# the significance threshold based on the FDR for selection. 
# Because this function is calculating additionally the absolute 
# FST observed between the subpopulations selected in opposite directions.
# This  adjustment can be used to identify truly selected sites.
### Fst value calculation -------------------------------------------------
calculate_FST_value <- function(data,
                                pop_low_phenotype_sel_1,
                                pop_low_phenotype_sel_2,
                                pop_high_phenotype_sel_1,
                                pop_high_phenotype_sel_2){
  time_0 <- Sys.time()
  cat("Calculation of the FST value started.","\n")
  gt <- extract.gt(data, element = "GT", as.numeric = FALSE)
  gt_file <- gt
  dt_pop_low_phenotype_1 <- gt_file[,which(str_detect(colnames(gt_file),pop_low_phenotype_sel_1))]
  dt_pop_low_phenotype_2 <- gt_file[,which(str_detect(colnames(gt_file),pop_low_phenotype_sel_2))]
  dt_pop_high_phenotype_1 <- gt_file[,which(str_detect(colnames(gt_file),pop_high_phenotype_sel_1))]
  dt_pop_high_phenotype_2 <- gt_file[,which(str_detect(colnames(gt_file),pop_high_phenotype_sel_2))]
  dp <- extract.gt(data, element = "DP")
  SNP_IDs <- rownames(dp)
  get_allele_freq_per_marker_per_pop <- function(data, population_name){
    get_allele_freq_per_marker <- function(i){
      N_alleles_total <- 2*length(which(!is.na(data[i,])))
      N_ALT_alleles <- 2*length(which(data[i,] == "1|1")) + 2*length(which(data[i,] == "1/1"))
      N_REF_alleles <- 2*length(which(data[i,] == "0|0")) + 2*length(which(data[i,] == "0/0"))
      N_HET_alleles <- length(which(data[i,] == "1|0")) + length(which(data[i,] == "1/0")) +
        length(which(data[i,] == "0|1")) + length(which(data[i,] == "0/1"))
      N_GT_tot <- length(which(!is.na(data[i,])))
      N_P_GT <- length(which(data[i,] == "0|0")) + length(which(data[i,] == "0/0"))
      N_Q_GT <- length(which(data[i,] == "1|1")) + length(which(data[i,] == "1/1"))
      freq_REF <- (N_REF_alleles + N_HET_alleles)/N_alleles_total
      freq_ALT <- (N_ALT_alleles + N_HET_alleles)/N_alleles_total
      freq_HET <- N_HET_alleles/N_GT_tot
      freq_P <- N_P_GT/N_GT_tot
      freq_Q <- N_Q_GT/N_GT_tot
      Exp_Het <- 2*N_P_GT/N_GT_tot*N_Q_GT/N_GT_tot
      marker_info <- cbind(freq_REF,freq_ALT,freq_HET,freq_P,freq_Q,N_GT_tot,Exp_Het)
      return(marker_info)
    }
    af_per_Marker_pop <- unlist(mclapply(1:nrow(data),get_allele_freq_per_marker))
    marker_tab_pop <- matrix(af_per_Marker_pop, ncol = 7, byrow = TRUE)
    marker_tab_pop <- as.data.frame(marker_tab_pop)
    colnames(marker_tab_pop) <- c(paste("freq_REF",population_name,sep="_"),paste("freq_ALT",population_name,sep="_"),
                                  paste("freq_HET",population_name,sep="_"),paste("freq_P",population_name,sep="_"),
                                  paste("freq_Q",population_name,sep="_"),paste("N_GT_tot",population_name,sep="_"),
                                  paste("Exp_HET",population_name,sep="_"))
    rownames(marker_tab_pop) <- SNP_IDs
    return(marker_tab_pop)
  }
  dt_pop_low_phenotype_1_dt <- get_allele_freq_per_marker_per_pop(data=dt_pop_low_phenotype_1, population_name = pop_low_phenotype_sel_1)
  dt_pop_low_phenotype_2_dt <- get_allele_freq_per_marker_per_pop(data=dt_pop_low_phenotype_2, population_name = pop_low_phenotype_sel_2)
  dt_pop_high_phenotype_1_dt <- get_allele_freq_per_marker_per_pop(data=dt_pop_high_phenotype_1, population_name = pop_high_phenotype_sel_1)
  dt_pop_high_phenotype_2_dt <- get_allele_freq_per_marker_per_pop(data=dt_pop_high_phenotype_2, population_name = pop_high_phenotype_sel_2)
  marker_table_per_pop <- cbind(dt_pop_low_phenotype_1_dt,dt_pop_low_phenotype_2_dt,
                                dt_pop_high_phenotype_1_dt,dt_pop_high_phenotype_2_dt)
  marker_table_per_pop <- as.data.frame(marker_table_per_pop)
  fst_value_calc <- function(population_1,population_2){
    mean <- (population_1+population_2) / 2
    var <- (population_1-mean)^2 + (population_2-mean)^2
    fst_pop <- var / (mean*(1-mean)+(var/2))
    return(fst_pop)
  }
  Fst_value_sd_1 <- fst_value_calc(population_1 = marker_table_per_pop[,1], 
                                   population_2 = marker_table_per_pop[,22])
  Fst_value_sd_2 <- fst_value_calc(population_1 = marker_table_per_pop[,15], 
                                   population_2 = marker_table_per_pop[,8])
  Fst_value_od_1 <- fst_value_calc(population_1 = marker_table_per_pop[,1], 
                                   population_2 = marker_table_per_pop[,15])
  Fst_value_od_2 <- fst_value_calc(population_1 = marker_table_per_pop[,22], 
                                   population_2 = marker_table_per_pop[,8])
  FSTSum_OD <- Fst_value_od_1 + Fst_value_od_2
  FSTSum_SD <- Fst_value_sd_1 + Fst_value_sd_2
  FST_value_dt <- cbind(data@fix[,c(1,2)],SNP_IDs,
                        Fst_value_sd_1,
                        Fst_value_sd_2,
                        Fst_value_od_1,
                        Fst_value_od_2,
                        FSTSum_OD,
                        FSTSum_SD)
  FST_value_dt <- as.data.frame(FST_value_dt)
  colnames(FST_value_dt) <- c("Chromosome","Position","SNP_ID",
                              "FST_value_opposite_dir1",
                              "FST_value_opposite_dir2",
                              "FST_value_same_dir1",
                              "FST_value_same_dir2",
                              "FST_sum_OD",
                              "FST_sum_SD")
  FST_value_dt <- cbind(marker_table_per_pop,FST_value_dt)
  cat("Calculation of the FST value is done.","\n")
  time_1 <- Sys.time()
  cat("This was done within",print(time_1 - time_0),"\n")
  return(FST_value_dt)
}
FST_values_od_cor <- calculate_FST_value(data = vcf_S4, 
                                         pop_low_phenotype_sel_1 = "Shoepag_1",
                                         pop_low_phenotype_sel_2 = "Shoepag_4",
                                         pop_high_phenotype_sel_1 = "Shoepag_3",
                                         pop_high_phenotype_sel_2 = "Shoepag_2")
write.table(FST_values_od_cor, "2022_03_24_GB1005_Fst_values.txt", row.names = TRUE, sep = "  ",
            quote = FALSE)