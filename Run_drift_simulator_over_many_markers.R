library(plyr)
library(foreach)
library(parallel)
library(data.table)
library(doMC)

#provide the path to this script; this is the work-horse in which neutral drift is simulated
setwd("/usr/users/mtost/Drift_simulations/")
source("Sim10_DriftSimulator.R")

cores <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK'))
cores <- detectCores(all.tests = FALSE, logical = TRUE)
registerDoMC(cores)

digits=2 ##number of digits to round to 
step=0.01 ##interval between allele frequencies
sim=10##number of simulations to run
cycles=3 #the number of cycles
markers=379485
minobs=40 #the number of minimal observations at a marker.
#This value is used as a minimal value to sample the marker coverage/missingness
# at a marker
seqind=96 #the number of individuals which are sequenced

### Population structure
total_pop_1=2500 ## total individuals/generation in population 1
## Population 1 only had 2500 individuals
total_pop_2=5000 ## total individuals/generation in population 2
total_pop_3=5000 ## total individuals/generation in population 3
total_pop_4=5000 ## total individuals/generation in population 4
# In each population ~250 individuals were selected.
# The selected indivduals are the mothers of the next generation.
# All 5000 plants could have contributed to the next generation
# as father plants, but due to the drift simulator function,
# the males + females have to equal the population size. 
# An ear can be pollinated by several fathers, and an ear consists
# out of many seeds.
males_pop_1=2250 ##total males/generation in population 1
females_pop_1=250 ##total females/generation in population 1
## In population 1 ~250 individuals were selected
males_pop_2=4750 ##total males/generation in population 1
females_pop_2=250 ##total females/generation in population 1
## In population 1 ~250 individuals were selected
males_pop_3=4750 ##total males/generation in population 1
females_pop_3=250 ##total females/generation in population 1
## In population 1 ~250 individuals were selected
males_pop_4=4750 ##total males/generation in population 1
females_pop_4=250 ##total females/generation in population 1
## In population 1 ~250 individuals were selected

fs<-round(seq(step,1-step,by=step),digits) ##create a vector of allele frequencies with step; exclude 0 and 1.

### Sample the possible range of allele frequencies
sample_af <- function(i){
  i = 1
  sample(fs,i,replace = TRUE)
}
fs_sampled <- unlist(mclapply(1:markers,sample_af))

sample_initial_af_base_pop <- unlist(mclapply(1:length(fs_sampled),sample_initial_allele_freq))
cat("Sampling of initial allele frequencies is done","\n")

run_simulations_pop1 <- function(i){
  unlist(mclapply(1:length(sample_initial_af_base_pop),seq_drifted_ind_pop1))
}
sampled_seq_ind_pop1 <- unlist(mclapply(1:sim,run_simulations_pop1))
cat("Simulation pop 1 is done","\n")
run_simulations_pop2 <- function(i){
  unlist(mclapply(1:length(sample_initial_af_base_pop),seq_drifted_ind_pop234))
}
sampled_seq_ind_pop2 <- unlist(mclapply(1:sim,run_simulations_pop2))
cat("Simulation pop 2 is done","\n")
run_simulations_pop3 <- function(i){
  unlist(mclapply(1:length(sample_initial_af_base_pop),seq_drifted_ind_pop234))
}
sampled_seq_ind_pop3 <- unlist(mclapply(1:sim,run_simulations_pop3))
cat("Simulation pop 3 is done","\n")
run_simulations_pop4 <- function(i){
  unlist(mclapply(1:length(sample_initial_af_base_pop),seq_drifted_ind_pop234))
}
sampled_seq_ind_pop4 <- unlist(mclapply(1:sim,run_simulations_pop4))
cat("Simulation pop 4 is done","\n")
cat("Simulation of drift is done","\n")
################################################################################
### Calculate quantiles for allele frequency change ----------------------------
################################################################################
# Quantiles from Kumar et al. 2020: 0.0000001  0.9999999                      ##
### AFC ########################################################################
vector_simulations <- c(sampled_seq_ind_pop1,sampled_seq_ind_pop2,sampled_seq_ind_pop3,sampled_seq_ind_pop4)
cat("The quantile of simulations is:","\n",
    quantile(vector_simulations, probs = c(0.0000001,0.999999), na.rm = TRUE),"\n")

### Fst ########################################################################
fst_value_calc <- function(population_1,population_2,population_3,population_4){
  # Fst value between two populations
  mean <- (population_1+population_2+population_3+population_4) / 4
  var <- (population_1-mean)^2 + (population_2-mean)^2 + (population_3-mean)^2 + (population_4-mean)^2
  fst_pop <- var / (mean*(1-mean)+(var/4))
  return(fst_pop)
}

Fst_value_all_pop <- fst_value_calc(population_1 = sampled_seq_ind_pop1,
                                    population_2 = sampled_seq_ind_pop2,
                                    population_3 = sampled_seq_ind_pop3,
                                    population_4 = sampled_seq_ind_pop4)
length(Fst_value_all_pop)
length(sampled_seq_ind_pop1)

Fst_value_all_pop <- matrix(data = Fst_value_all_pop, ncol = sim, nrow = markers, byrow = TRUE)
colnames(Fst_value_all_pop) <- paste("Fst_value_simulation", seq(1,sim,1), sep = "_")
quantiles_fst_whole_pop <- quantile(Fst_value_all_pop, probs = c(0.0000001,0.999999), na.rm = TRUE)

cat("The quantile of simulations for the Fst value calculation is:","\n",
    quantiles_fst_whole_pop,"\n")

quantile_per_marker_Fst <- function(i){
  quantile(Fst_value_all_pop[i,], probs = 0.999999, na.rm = TRUE)
}
quantile_marker_Fst <- unlist(mclapply(1:markers,quantile_per_marker_Fst))
quantile_marker_Fst <- matrix(data = quantile_marker_Fst, nrow = markers)
colnames(quantile_marker_Fst) <- "Fst_value"
quantile_marker_Fst <- as.data.table(quantile_marker_Fst)
dist_Fst <- quantile_marker_Fst[,.N, by = round(Fst_value,2)]
write.table(dist_Fst, "/usr/users/mtost/Drift_simulations/Sim14GB10_Fst_value_distribution.txt", 
            row.names = TRUE, sep = "  ",
            quote = FALSE)
  
################################################################################
### Calculate allele frequency differences between the populations -------------
################################################################################
##change list to a matrix to write.

sample_seq_ind_pop1 <- matrix(data = sampled_seq_ind_pop1, ncol = markers)
sample_seq_ind_pop2 <- matrix(data = sampled_seq_ind_pop2, ncol = markers)
sample_seq_ind_pop3 <- matrix(data = sampled_seq_ind_pop3, ncol = markers)
sample_seq_ind_pop4 <- matrix(data = sampled_seq_ind_pop4, ncol = markers)

colnames(sample_seq_ind_pop1) <- paste0("SNP",seq(1,markers,1))
colnames(sample_seq_ind_pop2) <- paste0("SNP",seq(1,markers,1))
colnames(sample_seq_ind_pop3) <- paste0("SNP",seq(1,markers,1))
colnames(sample_seq_ind_pop4) <- paste0("SNP",seq(1,markers,1))

sample_seq_ind_pop1 <- t(sample_seq_ind_pop1)
sample_seq_ind_pop2 <- t(sample_seq_ind_pop2)
sample_seq_ind_pop3 <- t(sample_seq_ind_pop3)
sample_seq_ind_pop4 <- t(sample_seq_ind_pop4)

colnames(sample_seq_ind_pop1) <- paste("pop1","sim",seq(1,sim,1), sep = "_")
colnames(sample_seq_ind_pop2) <- paste("pop2","sim",seq(1,sim,1), sep = "_")
colnames(sample_seq_ind_pop3) <- paste("pop3","sim",seq(1,sim,1), sep = "_")
colnames(sample_seq_ind_pop4) <- paste("pop4","sim",seq(1,sim,1), sep = "_")

dim(sample_seq_ind_pop1)

Short1 <- sample_seq_ind_pop1
Short2 <- sample_seq_ind_pop4
Tall1 <- sample_seq_ind_pop3
Tall2 <- sample_seq_ind_pop2

Short1_vs_Short2 <- Short1 -Short2
Tall1_vs_Tall2 <- Tall1 - Tall2
Short1_vs_Tall1 <- Short1 - Short2
Short2_vs_Tall2 <- Short1 - Short2

colnames(Short1_vs_Short2) <- paste("Short1_vs_Short2_simulation", seq(1,sim,1), sep = "_")
colnames(Tall1_vs_Tall2) <- paste("Tall1_vs_Tall2_simulation", seq(1,sim,1), sep = "_")
colnames(Short1_vs_Tall1) <- paste("Short1_vs_Tall1_simulation", seq(1,sim,1), sep = "_")
colnames(Short2_vs_Tall2) <- paste("Short2_vs_Tall2_simulation", seq(1,sim,1), sep = "_")

allele_freq_diff <- cbind(Short1_vs_Short2, Tall1_vs_Tall2,Short1_vs_Tall1,Short2_vs_Tall2)

cat("The quantile of simulations for allele frequency difference, first approach:","\n",
    quantile(allele_freq_diff, probs = c(0.0000001,0.999999), na.rm = TRUE),"\n")

quantile_per_marker_diffStat <- function(i){
  quantile(allele_freq_diff[i,], probs = 0.999999, na.rm = TRUE)
}
quantile_marker_diffStat <- unlist(mclapply(1:markers,quantile_per_marker_diffStat))
quantile_marker_diffStat <- matrix(data = quantile_marker_diffStat, nrow = markers)
quantile_marker_diffStat <- as.data.table(quantile_marker_diffStat)
colnames(quantile_marker_diffStat) <- "Allele_freq_diff"
dist_diffStat <- quantile_marker_diffStat[,.N, by = round(Allele_freq_diff,2)]

write.table(dist_diffStat, "/usr/users/mtost/Drift_simulations/Sim14GB10_DiffStat_value_distribution.txt", 
            row.names = TRUE, sep = "  ",
            quote = FALSE)
##### Different approach to get the quantile for allele freq difference --------
write.table(allele_freq_diff, "/usr/users/mtost/Drift_simulations/Sim14GB10_Allele_freq_diff_simulated_drift_diff_approach.txt", 
            row.names = TRUE, sep = "  ",
            quote = FALSE)
cat("The quantile of simulations for allele frequency difference, second approach:","\n",
    quantile(allele_freq_diff, probs = c(0.0000001,0.999999), na.rm = TRUE),"\n")



q()
n

