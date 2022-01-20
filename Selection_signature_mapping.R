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
vcf_S4 <- read.vcfR("YOUR_DATA_after_filtering.vcf.gz", verbose = FALSE)
setwd("/path/to/your/own/working/directory/")

########### Allele frequency calculation ---------------------------------------
### Get the allele frequencies per population ----------------------------------
calculate_the_allele_freq <- function(data,
                                      population_name_1,
                                      population_name_2,
                                      population_name_3,
                                      population_name_4){
  time_0 <- Sys.time()
  cat("Calculation of the allele frequency spectrum started.","\n")
  gt <- extract.gt(data, element = "GT", as.numeric = FALSE)
  gt_file <- gt
  dt_pop_1 <- gt_file[,which(str_detect(colnames(gt_file),population_name_1))]
  dt_pop_2 <- gt_file[,which(str_detect(colnames(gt_file),population_name_2))]
  dt_pop_3 <- gt_file[,which(str_detect(colnames(gt_file),population_name_3))]
  dt_pop_4 <- gt_file[,which(str_detect(colnames(gt_file),population_name_4))]
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
  dt_pop_1 <- get_allele_freq_per_marker_per_pop(data=dt_pop_1, population_name = population_name_1)
  dt_pop_2 <- get_allele_freq_per_marker_per_pop(data=dt_pop_2, population_name = population_name_2)
  dt_pop_3 <- get_allele_freq_per_marker_per_pop(data=dt_pop_3, population_name = population_name_3)
  dt_pop_4 <- get_allele_freq_per_marker_per_pop(data=dt_pop_4, population_name = population_name_4)
  marker_table_per_pop <- cbind(dt_pop_1,dt_pop_2,
                                dt_pop_3,dt_pop_4)
  marker_table_per_pop <- as.data.frame(marker_table_per_pop)
  cat("Pop 1: For",length(which(marker_table_per_pop[,3] < marker_table_per_pop[,7])),
      "markers the observed heterozygosity was higher than the expected heterozygosity!","\n")
  cat("Pop 2: For",length(which(marker_table_per_pop[,10] < marker_table_per_pop[,14])),
      "markers the observed heterozygosity was higher than the expected heterozygosity!","\n")
  cat("Pop 3: For",length(which(marker_table_per_pop[,17] < marker_table_per_pop[,22])),
      "markers the observed heterozygosity was higher than the expected heterozygosity!","\n")
  cat("Pop 4: For",length(which(marker_table_per_pop[,24] < marker_table_per_pop[,28])),
      "markers the observed heterozygosity was higher than the expected heterozygosity!","\n")
  cat("Pop 1:",length(which(marker_table_per_pop[,1] < marker_table_per_pop[,2])),"times the REF allele was the minor allele","\n")
  cat("Pop 1:",length(which(marker_table_per_pop[,1] > marker_table_per_pop[,2])),"times the REF allele was the major allele","\n")
  cat("Pop 2:",length(which(marker_table_per_pop[,8] < marker_table_per_pop[,9])),"times the REF allele was the minor allele","\n")
  cat("Pop 2:",length(which(marker_table_per_pop[,8] > marker_table_per_pop[,9])),"times the REF allele was the major allele","\n")
  cat("Pop 3:",length(which(marker_table_per_pop[,15] < marker_table_per_pop[,16])),"times the REF allele was the minor allele","\n")
  cat("Pop 3:",length(which(marker_table_per_pop[,15] > marker_table_per_pop[,16])),"times the REF allele was the major allele","\n")
  cat("Pop 4:",length(which(marker_table_per_pop[,22] < marker_table_per_pop[,22])),"times the REF allele was the minor allele","\n")
  cat("Pop 4:",length(which(marker_table_per_pop[,22] > marker_table_per_pop[,22])),"times the REF allele was the major allele","\n")
  cat("Calculation of the allele frequency spectrum is done.","\n")
  time_1 <- Sys.time()
  cat("This was done within",print(time_1 - time_0),"\n")
  return(marker_table_per_pop)
}
marker_table_per_pop <- calculate_the_allele_freq_diff(data = vcf_S4, 
                                                       population_name_1 = "Shoepag_1",
                                                       population_name_2 = "Shoepag_2",
                                                       population_name_3 = "Shoepag_3",
                                                       population_name_4 = "Shoepag_4")
write.table(marker_table_per_pop,"2022_01_17_GB1004_Allele_and_GT_freq_per_pop.txt", sep = "  ", row.names = TRUE,
            quote = FALSE)
########### Selection signature mapping analysis -------------------------------
### Fst_value calculation ------------------------------------------------------
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
### and separately. Furthermore this function calculates the absolute FST
### statistic between the subpopulations selected in opposite directions
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
  fst_value_calc <- function(population_1,population_2,population_3,population_4){
    mean <- (population_1+population_2+population_3+population_4) / 4
    var <- (population_1-mean)^2 + (population_2-mean)^2 + (population_3-mean)^2 + (population_4-mean)^2
    fst_pop <- var / (mean*(1-mean)+(var/4))
    return(fst_pop)
  }
  FST_value <- fst_value_calc(population_1 = marker_table_per_pop[,1], 
                              population_2 = marker_table_per_pop[,8],
                              population_3 = marker_table_per_pop[,15],
                              population_4 = marker_table_per_pop[,22])
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
  od_stat <- abs(Fst_value_sd_1 + Fst_value_sd_2)
  FST_value_dt <- cbind(data@fix[,c(1,2)],SNP_IDs,
                        FST_value,
                        Fst_value_sd_1,
                        Fst_value_sd_2,
                        Fst_value_od_1,
                        Fst_value_od_2,
                        od_stat)
  FST_value_dt <- as.data.frame(FST_value_dt)
  colnames(FST_value_dt) <- c("Chromosome","Position","SNP_ID",
                              "FST_value","FST_value_opposite_dir1",
                              "FST_value_opposite_dir2",
                              "FST_value_same_dir1",
                              "FST_value_same_dir2",
                              "Abs_sum_FST")
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
# The data table with the FST values is called "FST_values_od_cor".
# This function calculates the FST between subpopulations selected 
# in the same and opposite directions, which are needed to calculate
# the significance threshold based on the FDR for selection. 
# Because this function is calculating additionally the absolute 
# FST observed between the subpopulations selected in opposite directions.
# This  adjustment can be used to identify truly selected sites.
write.table(FST_values_od_cor, "GB1005_Fst_values.txt", row.names = TRUE, sep = "  ",
            quote = FALSE)
###### Allele frequency differences ---------------------------------------------------
## Calculate allele frequeny differences between the different 
## subpopulations and the absolute allele frequency difference 
## between all subpopulations together according to 
## Turner et al. 2012 :
calculate_the_allele_freq_diff <- function(data,
                                           pop_low_phenotype_sel_1,
                                           pop_low_phenotype_sel_2,
                                           pop_high_phenotype_sel_1,
                                           pop_high_phenotype_sel_2){
  time_0 <- Sys.time()
  cat("Calculation of the allele frequency differences started.","\n")
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
  selected_opposite_dir_1 <- as.numeric(marker_table_per_pop[,1] - marker_table_per_pop[,15])
  selected_opposite_dir_2 <- as.numeric(marker_table_per_pop[,8] - marker_table_per_pop[,22])
  selected_same_dir_1 <- as.numeric(marker_table_per_pop[,1] - marker_table_per_pop[,8])
  selected_same_dir_2 <- as.numeric(marker_table_per_pop[,15] - marker_table_per_pop[,22])
  od_stat <- abs(selected_opposite_dir_1 + selected_opposite_dir_2)
  allele_freq_diff_REF <- cbind(data@fix[,c(1,2)],SNP_IDs,
                                selected_opposite_dir_1,
                                selected_opposite_dir_2,
                                selected_same_dir_1,
                                selected_same_dir_2,
                                od_stat)
  allele_freq_diff_REF <- as.data.frame(allele_freq_diff_REF)
  colnames(allele_freq_diff_REF) <- c("Chromosome","Position","SNP_ID",
                                      "Low1_vs_Low2","High1_vs_High2",
                                      "Low1_vs_High1","Low2_vs_High2",
                                      "Abs_allele_freq_diff")
  cat("Calculation of the allele frequency differences is done.","\n")
  time_1 <- Sys.time()
  cat("This was done within",print(time_1 - time_0),"\n")
  return(allele_freq_diff_REF)
}
allele_freq_diff <- calculate_the_allele_freq_diff(data = vcf_S4, 
                                                   pop_low_phenotype_sel_1 = "Shoepag_1",
                                                   pop_low_phenotype_sel_2 = "Shoepag_4",
                                                   pop_high_phenotype_sel_1 = "Shoepag_3",
                                                   pop_high_phenotype_sel_2 = "Shoepag_2")
write.table(allele_freq_diff,"GB1006_Allele_freq_diff.txt", sep = "  ", row.names = TRUE,
            quote = FALSE)

q()