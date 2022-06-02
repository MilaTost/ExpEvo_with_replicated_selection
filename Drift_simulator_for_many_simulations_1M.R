#### DRIFT SIMULATOR CHANGED BY MILA TOST ######################################
rm(list = ls())
library(plyr)
library(foreach)
library(parallel)
library(data.table)
library(doMC)

cores <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK'))
cores <- detectCores(all.tests = FALSE, logical = TRUE)
registerDoMC(cores)
setwd("/usr/users/mtost/new_drift_simulations/")
# use the real allele frequencies from generation 0
## Local machine
#real_data <- read.table("C:/Users/mtost/Documents/Masterarbeit/Data_analysis/Simulations/Simulations_new/test_af_2016_for_simulations.txt")
#initial_Short1 <- real_data$freq_P

# use the real allele frequencies from generation 0
real_data <- read.table("/usr/users/mtost/wd_test_GB_easy_rerun/R_analysis_20162020/GB1307_Allele_and_GT_freq_per_pop_after_MAF005_filt.txt",
                        header = TRUE)
initial_Short1 <- real_data$freq_P_Short1Shoe0
initial_Short2 <- real_data$freq_P_Short2Shoe0
initial_Tall1 <- real_data$freq_P_Tall1Shoe0
initial_Tall2 <- real_data$freq_P_Tall2Shoe0

total_pop <- 5000
males_pop <- 5000
females_pop <- 250
sim <- 1000000
cycles <- 3
kernels_per_plant <- 500

initial_Short1_new <- sample(initial_Short1, sim, replace = TRUE)
initial_Short2_new <- sample(initial_Short2, sim, replace = TRUE)
initial_Tall1_new <- sample(initial_Short1, sim, replace = TRUE)
initial_Tall2_new <- sample(initial_Short1, sim, replace = TRUE)

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
  colnames(dt) <- 0:cycles
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

dt_drift_sim <- rbind(simulated_drift_Short1,
                      simulated_drift_Short2,
                      simulated_drift_Tall1,
                      simulated_drift_Tall2)
write.table(dt_drift_sim, "2022_05_27_Simulated_drift_Allele_freq_distribution_V4R1.txt", 
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
Fst_value_sd_1 <- fst_value_calc(population_1 = simulated_drift_Short1, 
                                 population_2 = simulated_drift_Short2)
Fst_value_sd_2 <- fst_value_calc(population_1 = simulated_drift_Tall1, 
                                 population_2 = simulated_drift_Tall2)
Fst_value_od_1 <- fst_value_calc(population_1 = simulated_drift_Short1, 
                                 population_2 = simulated_drift_Tall1)
Fst_value_od_2 <- fst_value_calc(population_1 = simulated_drift_Short2, 
                                 population_2 = simulated_drift_Tall2)
all_FST_values <- rbind(Fst_value_sd_1,
                        Fst_value_sd_2,
                        Fst_value_od_1,
                        Fst_value_od_2,
                        sum(Fst_value_sd_1,Fst_value_sd_2),
                        sum(Fst_value_od_1,Fst_value_od_2),
                        sum(Fst_value_od_1,Fst_value_sd_2),
                        sum(Fst_value_sd_1,Fst_value_od_2))
threshold <- quantile(all_FST_values, probs = 0.999999, na.rm = TRUE)
write.table(all_FST_values, "2022_05_27_Simulated_drift_FST_values_V4R1.txt", 
            row.names = TRUE, sep = "  ",
            quote = FALSE)
cat("The significance threshold based on drift simulations is:",threshold)

### Add sequencing -------------------------------------------------------------
drift_simulations_with_sequencing <- function(total_pop_size,
                                              no_males,
                                              no_females,
                                              simulations,
                                              cycles,
                                              kernels_per_plant,
                                              initial_freq,
                                              minobs,
                                              seqind){
  seed <- sample(1000,1)
  set.seed(seed)
  drift_many_loci <- function(a){
    seed <- sample(1000,1)
    set.seed(seed)
    maleprog <- c()
    femaleprog <- c()
    frequency <- c()
    initial_frequency <- initial_freq[a]
    drift_cycles <- function(i){
      seed <- sample(1000,1)
      set.seed(seed)
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
  coverage_sampling <- function(i){
    seed <- sample(1000,1)
    set.seed(seed)
    allele_freq <- dt[i,cycles+1]
    alleles <- c(rep("A",allele_freq*2*total_pop_size),rep("B",2*total_pop_size-(allele_freq*2*total_pop_size)))
    possible_marker_coverages <- seq(minobs, seqind)
    marker_coverage <- sample(possible_marker_coverages, 1, replace = TRUE)
    seq_sample <- sample(alleles, marker_coverage, replace = FALSE)
    allele_freq_after_seq <- length(which(seq_sample=="A"))/length(seq_sample)
    return(allele_freq_after_seq)
  }
  sequenced_frequencies <- mclapply(1:sim, coverage_sampling)
  dt_all <- cbind(dt,sequenced_frequencies)
  return(dt_all)
}

### RUN ------------------------------------------------------------------------
simulated_drift_Short1 <- drift_simulations_with_sequencing(total_pop_size = total_pop,
                                                 no_males = males_pop,
                                                 no_females = females_pop,
                                                 simulations = sim,
                                                 cycles = cycles,
                                                 kernels_per_plant = kernels_per_plant,
                                                 initial_freq = initial_Short1_new,
                                                 minobs = 40,
                                                 seqind = 96)
simulated_drift_Short2 <- drift_simulations_with_sequencing(total_pop_size = total_pop,
                                                 no_males = males_pop,
                                                 no_females = females_pop,
                                                 simulations = sim,
                                                 cycles = cycles,
                                                 kernels_per_plant = kernels_per_plant,
                                                 initial_freq = initial_Short2_new,
                                                 minobs = 40,
                                                 seqind = 96)
simulated_drift_Tall1 <- drift_simulations_with_sequencing(total_pop_size = total_pop,
                                                no_males = males_pop,
                                                no_females = females_pop,
                                                simulations = sim,
                                                cycles = cycles,
                                                kernels_per_plant = kernels_per_plant,
                                                initial_freq = initial_Tall1_new,
                                                minobs = 40,
                                                seqind = 96)
simulated_drift_Tall2 <- drift_simulations_with_sequencing(total_pop_size = total_pop,
                                                no_males = males_pop,
                                                no_females = females_pop,
                                                simulations = sim,
                                                cycles = cycles,
                                                kernels_per_plant = kernels_per_plant,
                                                initial_freq = initial_Tall2_new,
                                                minobs = 40,
                                                seqind = 96)
### Combine the different simulations to one data set ----------------------
dt_drift_sim <- rbind(simulated_drift_Short1,
                      simulated_drift_Short2,
                      simulated_drift_Tall1,
                      simulated_drift_Tall2)
write.table(dt_drift_sim, "2022_02_16_new_simulated_drift_with_sequencing_af_dist_V3R1.txt", 
            row.names = TRUE, sep = "  ",
            quote = FALSE)

fst_value_calc <- function(population_1,population_2){
  mean <- (population_1+population_2) / 2
  var <- (population_1-mean)^2 + (population_2-mean)^2
  fst_pop <- var / (mean*(1-mean)+(var/2))
  return(fst_pop)
}
Fst_value_sd_1 <- fst_value_calc(population_1 = simulated_drift_Short1, 
                                 population_2 = simulated_drift_Short2)
Fst_value_sd_2 <- fst_value_calc(population_1 = simulated_drift_Tall1, 
                                 population_2 = simulated_drift_Tall2)
Fst_value_od_1 <- fst_value_calc(population_1 = simulated_drift_Short1, 
                                 population_2 = simulated_drift_Tall1)
Fst_value_od_2 <- fst_value_calc(population_1 = simulated_drift_Short2, 
                                 population_2 = simulated_drift_Tall2)
all_FST_values <- rbind(Fst_value_sd_1,
                        Fst_value_sd_2,
                        Fst_value_od_1,
                        Fst_value_od_2,
                        sum(Fst_value_sd_1,Fst_value_sd_2),
                        sum(Fst_value_od_1,Fst_value_od_2),
                        sum(Fst_value_od_1,Fst_value_sd_2),
                        sum(Fst_value_sd_1,Fst_value_od_2))
threshold <- quantile(all_FST_values, probs = 0.999999, na.rm = TRUE)
cat("The significance threshold based on drift simulations with sequenicng is:",threshold)
write.table(all_FST_values, "2022_05_27_FST_values_V4.txt", 
            row.names = TRUE, sep = "  ",
            quote = FALSE)
### Note: This function is free to use and distribute.  Please provide appropriate credit when doing so.