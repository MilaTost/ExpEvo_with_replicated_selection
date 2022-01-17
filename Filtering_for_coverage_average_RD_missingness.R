################ Filtering_for_coverage_average_RD_missingness #################
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
vcf_S4 <- read.vcfR("/path/to/your/own/VCF_file_directory/YOUR_DATA.vcf", verbose = FALSE)

setwd("/path/to/your/own/working/directory/")
# Here you need to insert the path to your own working directory
cat("The data set is loaded.","\n")
time_1 <- Sys.time()
print(time_1 - time_0)

#### Overview over VCF file with functions from the vcfR package ---------------
chrom <- create.chromR(name="Supercontig", vcf=vcf_S4, verbose=FALSE)
png(file='GB10P1_Quality_plot_distribution_before_filtering_vcf.png',width = 480, height = 480, units = "px", pointsize=12)
par(
  family='sans',
  cex.axis=1,
  bg='white')
chromoqc(chrom, dp.alpha = 66)
dev.off()

cat("VCF file before anything was done:","\n",
    "the vcf file S4 format contains:",nrow(vcf_S4@fix),"markers.","\n")
rm(chrom)

####### Filter for individuals with low coverage -------------------------------
calc_coverage_per_individual <- function(data){
  cat("The calculation coverage of individuals started.","\n")
  time_0 <- Sys.time()
  gt <- extract.gt(data, element = "GT", as.numeric=TRUE)
  coverage_ind <- function(i){
    length(which(!is.na(gt[,i])))
  }
  coverage_ind_tot <- mapply(coverage_ind,1:ncol(gt))
  coverage_ind_rel <- coverage_ind_tot/NROW(gt)
  NAs_per_ind <- cbind(coverage_ind_tot,coverage_ind_rel)
  rownames(NAs_per_ind) <- colnames(gt)
  NAs_per_ind <- as.data.frame(NAs_per_ind)
  cat("The calculation coverage of individuals is finished.","\n")
  time_1 <- Sys.time()
  cat("This was done within",print(time_1 - time_0),"\n")
  return(NAs_per_ind)
}
NAs_per_ind <- calc_coverage_per_individual(data = vcf_S4)

### Create a file for plotting  ------------------------------------------------
write.table(NAs_per_ind,"GB1001_Coverage_per_individual.txt", sep = "  ", row.names = TRUE,
            quote = FALSE)
### Start filtering-------------------------------------------------------------
filter_for_ind_coverage <- function(data, threshold){
  time_0 <- Sys.time()
  cat("Filtering for low coverage individuals started.","\n")
  cat("Before filtering the vcf_S4 file contained:",ncol(data@gt)-1,"individuals.","\n")
  gt <- extract.gt(data, element = "GT", as.numeric=TRUE)
  coverage_ind <- function(i){
    length(which(!is.na(gt[,i])))
  }
  coverage_ind_tot <- mapply(coverage_ind,1:ncol(gt))
  coverage_ind_rel <- coverage_ind_tot/NROW(gt)
  NAs_per_ind <- cbind(coverage_ind_tot,coverage_ind_rel)
  rownames(NAs_per_ind) <- colnames(gt)
  NAs_per_ind <- as.data.frame(NAs_per_ind)
  index_ind <- rownames(NAs_per_ind[which(NAs_per_ind$coverage_ind_rel < threshold),])
  look_up_index <- function(i){
    which(colnames(gt)==index_ind[i])
  }
  index_vcf_S4 <- unlist(mclapply(1:length(index_ind),look_up_index))
  data@gt <- data@gt[,-index_vcf_S4]
  cat("The vcf_S4 file still contains:",ncol(data@gt)-1,"individuals.","\n")
  cat("Filtering for low coverage individuals is done.","\n")
  time_1 <- Sys.time()
  cat("This was done within",print(time_1 - time_0),"\n")
  return(data)
}
vcf_S4 <- filter_for_ind_coverage(data = vcf_S4,
                                  threshold = 0.1)
#### Get average sample read depth with dp function from the vcfR package ######
calulate_average_read_depth <- function(data){
  cat("Calculation of average read depth started.","\n")
  time_0 <- Sys.time()
  dp <- extract.info(data, element = "DP", as.numeric=TRUE)
  an <- extract.info(data, element = "AN", as.numeric=TRUE)
  obs <- (an/2)
  average_depth <- dp/obs
  average_depth <- as.data.frame(average_depth)
  rownames(average_depth) <- paste(data@fix[,1],data@fix[,2], sep = "_")
  average_depth$SNP_ID <- paste(data@fix[,1],data@fix[,2], sep = "_")
  average_depth$tot_DP <- dp
  average_depth$observations <- obs
  average_depth$observations <- as.numeric(average_depth$observations)
  average_depth$tot_DP <- as.numeric(average_depth$tot_DP)
  average_depth$average_depth <- as.numeric(average_depth$average_depth)
  cat("Calculation of average read depth is finished.","\n")
  time_1 <- Sys.time()
  cat("This was done within",print(time_1 - time_0),"\n")
  return(average_depth)
}
average_depth <- calulate_average_read_depth(data = vcf_S4)

write.table(average_depth,"GB1002_Average_Read_Depth.txt", sep = "  ", row.names = TRUE,
            quote = FALSE)
### Start filtering ------------------------------------------------------------
filter_for_average_read_depth <- function(data, min_av_depth, max_av_depth){
  cat("Filtering for average read depth started.","\n")
  cat("VCF file before filtering contains:",nrow(data@fix),"markers.","\n")
  time_0 <- Sys.time()
  dp <- extract.info(data, element = "DP", as.numeric=TRUE)
  an <- extract.info(data, element = "AN", as.numeric=TRUE)
  obs <- (an/2)
  average_depth <- dp/obs
  average_depth <- as.data.frame(average_depth)
  rownames(average_depth) <- paste(data@fix[,1],data@fix[,2], sep = "_")
  average_depth$SNP_ID <- paste(data@fix[,1],data@fix[,2], sep = "_")
  average_depth$tot_DP <- dp
  average_depth$observations <- obs
  average_depth$observations <- as.numeric(average_depth$observations)
  average_depth$tot_DP <- as.numeric(average_depth$tot_DP)
  average_depth$average_depth <- as.numeric(average_depth$average_depth)
  index_depth <- which(average_depth$average_depth > min_av_depth & average_depth$average_depth < max_av_depth)
  data@gt <- data@gt[index_depth,]
  data@fix <- data@fix[index_depth,]
  cat(length(index_depth),"markers passed the threshold.","\n")
  cat("VCF file before filtering contains:",nrow(data@fix),"markers.","\n")
  cat("Filtering for average read depth is finished.","\n")
  time_1 <- Sys.time()
  cat("This was done within",print(time_1 - time_0),"\n")
  return(data)
}
vcf_S4 <- filter_for_average_read_depth(data = vcf_S4, 
                                        min_av_depth = 1, 
                                        max_av_depth = 10)
## Create again a quality control plot -----------------------------------------
chrom <- create.chromR(name="Supercontig", vcf=vcf_S4, verbose=FALSE)
png(file='GB10P2_Quality_plot_distribution_after_filtering_for_depth_vcf.png',width = 480, height = 480, units = "px", pointsize=12)
par(
  family='sans',
  cex.axis=1,
  bg='white')
chromoqc(chrom, dp.alpha = 66)
dev.off()
rm(chrom)
## Calculate amount of markers with too many missing observations --------------
calculate_missingness_per_pop <- function(data, 
                                          population_1,
                                          population_2,
                                          population_3,
                                          population_4){
  cat("The calculation of missingness started.","\n")
  time_0 <- Sys.time()
  gt_pop <- extract.gt(vcf_S4, element = "GT", as.numeric=TRUE)
  dp <- extract.info(data, element = "DP", as.numeric=TRUE)
  pop_1 <- gt_pop[,which(str_detect(colnames(gt_pop),population_1))]
  pop_2 <- gt_pop[,which(str_detect(colnames(gt_pop),population_2))]
  pop_3 <- gt_pop[,which(str_detect(colnames(gt_pop),population_3))]
  pop_4 <- gt_pop[,which(str_detect(colnames(gt_pop),population_4))]
  get_missing_obs_per_pop <- function(pop){
    get_missing_obs <- function(i){
      length(which(is.na(pop[i,])))
    }
    mis_obs_pop <- unlist(mclapply(1:nrow(pop), get_missing_obs))
    return(mis_obs_pop)
  }
  mis_obs_pop1 <- get_missing_obs_per_pop(pop = pop_1)
  mis_obs_pop2 <- get_missing_obs_per_pop(pop = pop_2)
  mis_obs_pop3 <- get_missing_obs_per_pop(pop = pop_3)
  mis_obs_pop4 <- get_missing_obs_per_pop(pop = pop_4)
  NAs_per_marker_per_pop <- cbind(mis_obs_pop1,mis_obs_pop2,mis_obs_pop3,mis_obs_pop4)
  NAs_per_marker_per_pop <- as.data.frame(NAs_per_marker_per_pop)
  NAs_per_marker_per_pop <- cbind(rownames(gt_pop),NAs_per_marker_per_pop)
  colnames(NAs_per_marker_per_pop) <- c("SNP_ID","mis_obs_pop1","mis_obs_pop2","mis_obs_pop3","mis_obs_pop4")
  NAs_per_marker_per_pop <- as.data.frame(NAs_per_marker_per_pop)
  get_range_of_missingness <- function(i){
    abs(range(NAs_per_marker_per_pop[i,2:5])[1]-range(NAs_per_marker_per_pop[i,2:5])[2])
  }
  markers_diff <- unlist(mclapply(1:nrow(NAs_per_marker_per_pop),get_range_of_missingness))
  NAs_per_marker_per_pop$range_of_missingness <- markers_diff
  cat("The calculation of missingness is finished.","\n")
  time_1 <- Sys.time()
  cat("This was done within",print(time_1 - time_0),"\n")
  return(NAs_per_marker_per_pop)
}
NAs_per_marker_per_pop <- calculate_missingness_per_pop(data = vcf_S4,
                                                        population_1 = "Shoepag_1",
                                                        population_2 = "Shoepag_2",
                                                        population_3 = "Shoepag_3",
                                                        population_4 = "Shoepag_4")
write.table(NAs_per_marker_per_pop,"GB1003_Missingness.txt", sep = "  ", row.names = TRUE,
            quote = FALSE)

min_40_obs <- length(rownames(NAs_per_marker_per_pop)[which(96-NAs_per_marker_per_pop[,2] > 40 &
                                                              96-NAs_per_marker_per_pop[,3] > 40 &
                                                              96-NAs_per_marker_per_pop[,4] > 40 &
                                                              96-NAs_per_marker_per_pop[,5] > 40)])
min_20_obs <- length(rownames(NAs_per_marker_per_pop)[which(96-NAs_per_marker_per_pop[,2] > 20 &
                                                              96-NAs_per_marker_per_pop[,3] > 20 &
                                                              96-NAs_per_marker_per_pop[,4] > 20 &
                                                              96-NAs_per_marker_per_pop[,5] > 20)])
cat("When you would require, SNPs with at least 20 observations","\n",
    "you would keep", min_20_obs,"SNPs in the data set.","\n",
    "When you would require, SNPs with at least 40 observations","\n",
    "you would keep", min_40_obs,"SNPs in the data set","\n")
### Start filtering for missingness --------------------------------------------
filter_missingness_per_pop <- function(data, 
                                       population_1,
                                       population_2,
                                       population_3,
                                       population_4,
                                       max_missing_obs){
  cat("Filtering for missing observations per markers started.","\n")
  time_0 <- Sys.time()
  cat("VCF file contains before filtering for missingness:",nrow(vcf_S4@fix),"markers.","\n")
  gt_pop <- extract.gt(vcf_S4, element = "GT", as.numeric=TRUE)
  dp <- extract.info(data, element = "DP", as.numeric=TRUE)
  pop_1 <- gt_pop[,which(str_detect(colnames(gt_pop),population_1))]
  pop_2 <- gt_pop[,which(str_detect(colnames(gt_pop),population_2))]
  pop_3 <- gt_pop[,which(str_detect(colnames(gt_pop),population_3))]
  pop_4 <- gt_pop[,which(str_detect(colnames(gt_pop),population_4))]
  get_missing_obs_per_pop <- function(pop){
    get_missing_obs <- function(i){
      length(which(is.na(pop[i,])))
    }
    mis_obs_pop <- unlist(mclapply(1:nrow(pop), get_missing_obs))
    return(mis_obs_pop)
  }
  mis_obs_pop1 <- get_missing_obs_per_pop(pop = pop_1)
  mis_obs_pop2 <- get_missing_obs_per_pop(pop = pop_2)
  mis_obs_pop3 <- get_missing_obs_per_pop(pop = pop_3)
  mis_obs_pop4 <- get_missing_obs_per_pop(pop = pop_4)
  NAs_per_marker_per_pop <- cbind(mis_obs_pop1,mis_obs_pop2,mis_obs_pop3,mis_obs_pop4)
  NAs_per_marker_per_pop <- as.data.frame(NAs_per_marker_per_pop)
  NAs_per_marker_per_pop <- cbind(rownames(gt_pop),NAs_per_marker_per_pop)
  colnames(NAs_per_marker_per_pop) <- c("SNP_ID","mis_obs_pop1","mis_obs_pop2","mis_obs_pop3","mis_obs_pop4")
  NAs_per_marker_per_pop <- as.data.frame(NAs_per_marker_per_pop)
  get_range_of_missingness <- function(i){
    abs(range(NAs_per_marker_per_pop[i,2:5])[1]-range(NAs_per_marker_per_pop[i,2:5])[2])
  }
  markers_diff <- unlist(mclapply(1:nrow(NAs_per_marker_per_pop),get_range_of_missingness))
  NAs_per_marker_per_pop$range_of_missingness <- markers_diff
  index_NA_markers_dp <- subset(rownames(dp),96-as.numeric(NAs_per_marker_per_pop[,2]) > max_missing_obs &
                                  96-as.numeric(NAs_per_marker_per_pop[,3]) > max_missing_obs &
                                  96-as.numeric(NAs_per_marker_per_pop[,4]) > max_missing_obs &
                                  96-as.numeric(NAs_per_marker_per_pop[,5]) > max_missing_obs)
  index_NA_markers_gt <- subset(rownames(gt_pop),96-as.numeric(NAs_per_marker_per_pop[,2]) > max_missing_obs &
                                  96-as.numeric(NAs_per_marker_per_pop[,3]) > max_missing_obs &
                                  96-as.numeric(NAs_per_marker_per_pop[,4]) > max_missing_obs &
                                  96-as.numeric(NAs_per_marker_per_pop[,5]) > max_missing_obs)
  cat(length(index_NA_markers_gt),"markers passed your threshold for missingness.","\n")
  data@gt <- data@gt[match(index_NA_markers_gt, rownames(gt_pop)),]
  data@fix <- data@fix[match(index_NA_markers_dp, rownames(dp)),]
  cat("VCF file contains after filtering for missingness:",nrow(data@gt),"markers.","\n")
  cat("Filtering for missingness is finished.","\n")
  time_1 <- Sys.time()
  cat("This was done within",print(time_1 - time_0),"\n")
  return(data)
}
vcf_S4 <- filter_missingness_per_pop(data = vcf_S4, 
                                     population_1 = "Shoepag_1",
                                     population_2 = "Shoepag_2",
                                     population_3 = "Shoepag_3",
                                     population_4 = "Shoepag_4",
                                     max_missing_obs = 40)
###### Filter for monomorphic and triallelic markers ---------------------------
cat("Filtering for diallelic markers started.","\n")
cat("VCF file before filtering for diallelic markers contains:,","\n",
    nrow(vcf_S4@fix),"markers.","\n")
filter_for_non_diallelic_non_polymorphic <- function(data){
  cat("Filtering for diallelic markers started.","\n")
  time_0 <- Sys.time()
  dp <- extract.gt(data, element = "DP", as.numeric=TRUE)
  gt <- extract.gt(data, element = "GT", as.numeric=TRUE)
  index_not_diallelic <- which(str_count(data@fix[,5])!=1)
  cat(length(index_not_diallelic),"markers were not diallelic.","\n")
  data@gt <- data@gt[-index_not_diallelic,]
  data@fix <- data@fix[-index_not_diallelic,]
  index_polymorphic_markers <- which(is.polymorphic(data, na.omit = TRUE))
  data@gt <- data@gt[index_polymorphic_markers,]
  data@fix <- data@fix[index_polymorphic_markers,]
  cat(length(index_polymorphic_markers),"markers were polymorphic.","\n")
  cat("The VCF file still contains:",nrow(vcf_S4@gt),"markers.","\n")
  time_1 <- Sys.time()
  cat("This was done within",print(time_1 - time_0),"\n")
  return(data)
}
vcf_S4 <- filter_for_non_diallelic_non_polymorphic(data = vcf_S4)

#### Write new file ------------------------------------------------------------
write.vcf(vcf_S4, file="YOUR_DATA_after_filtering.vcf.gz", mask=FALSE)
rm(dp,gt)

chrom <- create.chromR(name="Supercontig", vcf=vcf_S4, verbose=FALSE)
png(file='GB10P3_Quality_plot_distribution_after_filtering_vcf.png',width = 480, height = 480, units = "px", pointsize=12)
par(
  family='sans',
  cex.axis=1,
  bg='white')
chromoqc(chrom, dp.alpha = 66)
dev.off()
rm(chrom)
q()
