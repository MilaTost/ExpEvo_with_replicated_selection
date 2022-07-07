#### DRIFT SIMULATOR CHANGED BY MILA TOST ######################################
rm(list = ls())
library(plyr)
library(foreach)
library(parallel)
library(data.table)
library(doMC)
library(stringr)

cores <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK'))
cores <- detectCores(all.tests = FALSE, logical = TRUE)
registerDoMC(cores)
setwd("/usr/users/mtost/new_drift_simulations/")
# use the real allele frequencies from generation 0
marker_data_tab <- read.table("/usr/users/mtost/new_drift_simulations/2022_06_29_Allele_freq_Gen0.txt",
                              row.names = 1)
### SHORTCUT: read in directly the allele frequencies --------------------------
initial_Short1 <- marker_data_tab$freq_P_Gen1_Shoe0
initial_Short2 <- marker_data_tab$freq_P_Gen2_Shoe0
initial_Tall1 <- marker_data_tab$freq_P_Rol1_Shoe0
initial_Tall2 <- marker_data_tab$freq_P_Rol2_Shoe0

total_pop <- 5000
males_pop <- 5000
females_pop <- 250
### 1 M simulations ----------------------------------------------------------
sim <- 1000000
cat(sim, "simulations are conducted.","\n")
cycles <- 3
kernels_per_plant <- 500

initial_Short1_new <- sample(initial_Short1, sim, replace = TRUE)
initial_Short2_new <- sample(initial_Short2, sim, replace = TRUE)
initial_Tall1_new <- sample(initial_Tall1, sim, replace = TRUE)
initial_Tall2_new <- sample(initial_Tall2, sim, replace = TRUE)

drift_simulations <- function(total_pop_size,
                              no_males,
                              no_females,
                              simulations,
                              cycles,
                              kernels_per_plant,
                              initial_freq){
  drift_many_loci <- function(a){
    maleprog <- c()
    femaleprog <- c()
    frequency <- c()
    initial_frequency <- as.numeric(initial_freq[a])
    drift_cycles <- function(i){
      prog <- c(rep("A",initial_frequency*total_pop_size*2),rep("B",(1-initial_frequency)*total_pop_size*2))
      maleprog <- sample(prog, size = no_males, replace=T)  #sample male alleles from total population
      femaleprog <- sample(prog,size = no_females,replace=F) #sample females from total population
      kernels_from_female_mother_plants <- rep(femaleprog,kernels_per_plant)
      pollen_father_plants <- sample(maleprog, length(kernels_from_female_mother_plants),replace = TRUE)
      progenies <- c(kernels_from_female_mother_plants, pollen_father_plants)
      prog_1 <- sample(progenies, total_pop_size, replace = FALSE) #from all harvested ears, only 5000 are sown
      frequency <- length(which(prog_1=="A"))/(2*total_pop_size)
      return(frequency)
    }
    drifted_freq_over_cycles <- mclapply(1:cycles, drift_cycles)
    drifted_freq_over_cycles <- unlist(drifted_freq_over_cycles)
    frequency_0 <- initial_frequency
    drifted_freq <- c(frequency_0, drifted_freq_over_cycles)
    return(drifted_freq)
  }
  drift_many_loci_freq <- mclapply(1:sim, drift_many_loci)
  data <- unlist(drift_many_loci_freq)
  dt <- matrix(data = data, ncol={cycles+1},nrow=sim, byrow = TRUE)
  colnames(dt) <- paste("Generation",0:cycles, sep = "_")
  return(dt)
}
### RUN ------------------------------------------------------------------------
simulated_drift_Short1 <- drift_simulations(total_pop_size = total_pop,
                                            no_males = males_pop,
                                            no_females = females_pop,
                                            simulations = sim,
                                            cycles = cycles,
                                            kernels_per_plant = kernels_per_plant,
                                            initial_freq = initial_Short1_new)
simulated_drift_Short2 <- drift_simulations(total_pop_size = total_pop,
                                            no_males = males_pop,
                                            no_females = females_pop,
                                            simulations = sim,
                                            cycles = cycles,
                                            kernels_per_plant = kernels_per_plant,
                                            initial_freq = initial_Short2_new)
simulated_drift_Tall1 <- drift_simulations(total_pop_size = total_pop,
                                           no_males = males_pop,
                                           no_females = females_pop,
                                           simulations = sim,
                                           cycles = cycles,
                                           kernels_per_plant = kernels_per_plant,
                                           initial_freq = initial_Tall1_new)
simulated_drift_Tall2 <- drift_simulations(total_pop_size = total_pop,
                                           no_males = males_pop,
                                           no_females = females_pop,
                                           simulations = sim,
                                           cycles = cycles,
                                           kernels_per_plant = kernels_per_plant,
                                           initial_freq = initial_Tall2_new)

simulated_drift_Short1 <- as.data.frame(simulated_drift_Short1)
simulated_drift_Short2 <- as.data.frame(simulated_drift_Short2)
simulated_drift_Tall1 <- as.data.frame(simulated_drift_Tall1)
simulated_drift_Tall2 <- as.data.frame(simulated_drift_Tall2)

dt_drift_sim <- rbind(simulated_drift_Short1,
                      simulated_drift_Short2,
                      simulated_drift_Tall1,
                      simulated_drift_Tall2)
write.table(dt_drift_sim, "2022_06_29_Simulated_drift_Allele_freq_distribution_V4R1.txt", 
            row.names = TRUE, sep = "  ",
            quote = FALSE)
### Calculate FST values -------------------------------------------------------
fst_value_calc <- function(population_1,population_2){
  population_1 <- as.numeric(population_1)
  population_2 <- as.numeric(population_2)
  mean <- (population_1+population_2) / 2
  var <- (population_1-mean)^2 + (population_2-mean)^2
  fst_pop <- var / (mean*(1-mean)+(var/2))
  return(fst_pop)
}
Fst_value_od_1 <- fst_value_calc(population_1 = simulated_drift_Short1$Generation_3, 
                                 population_2 = simulated_drift_Tall1$Generation_3)
Fst_value_od_2 <- fst_value_calc(population_1 = simulated_drift_Short2$Generation_3, 
                                 population_2 = simulated_drift_Tall2$Generation_3)
all_FST_values <- cbind(Fst_value_od_1,
                        Fst_value_od_2)
### Calculate the FST sum value ------------------------------------------------
calculate_FST_sum_for_drift <- function(FST_values){
  calculate_FST_sum_for_drift_per_sim <- function(i){
    Fst_sum <- sum(FST_values[i,1],FST_values[i,2])
    all_fst_sum_values_per_sim <- Fst_sum
    return(all_fst_sum_values_per_sim)
  }
  FST_sum_values <- mclapply(1:nrow(FST_values),calculate_FST_sum_for_drift_per_sim)
  unlisted_FST_sum_values <- unlist(FST_sum_values)
  matrix_FST_sum_values <- matrix(data = unlisted_FST_sum_values, nrow = nrow(FST_values),
                                  ncol = 1, byrow = TRUE)
  dt_FST_sum_values <- as.data.frame(matrix_FST_sum_values)
  return(dt_FST_sum_values)
}
all_FST_sum_values <- calculate_FST_sum_for_drift(FST_values = all_FST_values)
# Significance threshold based on the 99.9999th percentile 
threshold_FSTSum <- quantile(all_FST_sum_values, probs = 0.999999, na.rm = TRUE)
cat("FstSum: The significance threshold based on drift simulations","\n",
    "and on the 99.9999th percentile is:",threshold_FSTSum,"\n")
# Significance threshold based on the 99.999th percentile 
threshold_FSTSum <- quantile(all_FST_sum_values, probs = 0.99999, na.rm = TRUE)
cat("FstSum: The significance threshold based on drift simulations","\n",
    "and on the 99.999th percentile is:",threshold_FSTSum,"\n")
# Significance threshold based on the 99.99th percentile 
threshold_FSTSum <- quantile(all_FST_sum_values, probs = 0.9999, na.rm = TRUE)
cat("FstSum: The significance threshold based on drift simulations","\n",
    "and on the 99.99th percentile is:",threshold_FSTSum,"\n")
# Significance threshold based on the 99.9th percentile 
threshold_FSTSum <- quantile(all_FST_sum_values, probs = 0.999, na.rm = TRUE)
cat("FstSum: The significance threshold based on drift simulations","\n",
    "and on the 99.9th percentile is:",threshold_FSTSum,"\n")
# Significance threshold based on the 99th percentile 
threshold_FSTSum <- quantile(all_FST_sum_values, probs = 0.99, na.rm = TRUE)
cat("FstSum: The significance threshold based on drift simulations","\n",
    "and on the 99th percentile is:",threshold_FSTSum,"\n")
# Create the table with all FSTSum values --------------------------------------
write.table(all_FST_sum_values, "2022_06_29_Simulated_drift_FST_sum_values_V4R2.txt", 
            row.names = TRUE, sep = "  ",
            quote = FALSE)