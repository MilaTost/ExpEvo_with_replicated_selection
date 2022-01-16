# Selection signature mapping scripts with replicated selection
## Table of contents
[0 Introduction](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#introduction) <br />
1 Pipeline for the analysis of GBS data adapted from [Wickland et al. 2013](https://github.com/dpwickland/GB-eaSy) <br />
[2 Filtering](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#filtering) <br />
&emsp;[2.1 Filtering for individual samples with low coverage](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#21-filtering-for-individual-samples-with-low-coverage)
&emsp;[2.2 Filtering for read depth per sample](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#22-filtering-for-read-depth-per-sample)
&emsp;[2.3 Filtering for missingness](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#23-filtering-for-missingness)
&emsp;[2.4 Removal of non-diallelic and non-polymorphic markers](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#24-removal-of-non-diallelic-and-non-polymorphic-markers)
[3 Selection signature mapping](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#3-selection-signature-mapping) <br />
&emsp; [3.1 <img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> leveraging replicated selection](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#31--leveraging-replicated-selection) <br />
&emsp; [3.2 Allele frequency differences](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#32-allele-frequency-differences) <br />
[4 Significance thresholds](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#4-significance-thresholds) <br />
&emsp; [4.1 Based on the empiric distribution](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#41-based-on-the-empiric-distribution) <br />
&emsp; [4.2 Based on drift simulations](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#42-based-on-drift-simulations) <br />
&emsp; [4.3 Simulation of drift](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#43-simulation-of-drift) <br />
&emsp; [4.4 Based on the FDR for selection](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#44-based-on-the-fdr-for-selection)<br />
[5 Sauron plot](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#5-sauron-plot) <br />
[6 Other plotting scripts](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#6-other-plotting-scripts) <br />


## 0 Introduction
This repository contains scripts for selection signature mapping with replicated selection.
Usually in experimental evolution studies a random-mating population is selected in divergent directions for one
trait of interest. The highly divergent subpopulations are then scanned for signatures of selection.
A previous problem in these studies was the determination of significance thresholds. 
One solution to this problem is the usage of **replicated selection**. 
This was introduced by Turner and Miller (2012), but only tested on *Drospohila melanogaster*. <br /> <br />
In this approach two subpopulations are selected in the same direction. 
By that we differentiate between changes caused by selection and changes caused by other factors like drift.
The **replicated selection** is used to calculate a *FDR for selection*. 
[Turner and Miller (2012)](http://www.genetics.org/content/suppl/2012/03/30/genetics.112.139337.DC1) applied this significance threshold only on the
statistic based on the allele frequency differences. 
Whereas we applied the this new significance threshold to the commonly used 
          <img src="https://render.githubusercontent.com/render/math?math=F_{ST}">
statistic.

## 1 Pipeline for the analysis of GBS data adapted from [Wickland et al. 2013](https://github.com/dpwickland/GB-eaSy)
For the analysis of our raw reads from paired-end genotyping-by-sequencing (GBS) with ApeKI according to [Elshire et al. 2011](https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0019379&type=printable), we used the GB-eaSy pipeline from [Wickland et al. 2013](https://github.com/dpwickland/GB-eaSy). <br />
The pipeline consists out of several steps, which comprises:
- [Demultiplexing and trimming of the adapter sequence](https://github.com/dpwickland/GB-eaSy#step-2-demultiplex-raw-reads)
- [Alignement to the reference genome](https://github.com/dpwickland/GB-eaSy#step-3-align-to-reference)
- [Create a list of sorted bam files](https://github.com/dpwickland/GB-eaSy#step-4-create-list-of-sorted-bam-files)
- [Generate pileup and SNP calling](https://github.com/dpwickland/GB-eaSy#step-5-and-6-generate-pileup-and-call-snps)
- Filtering for quality parameters <br /> <br />

This bash-script was run on every sequenced plate separately, since every sequenced plate had it's own adapters and the same set of barcodes was used for all sequenced plates in our case. The *Demultiplexing* and *Alignement to the reference genome* steps were conducted by the `Demultiplexing_and_alignement_to_the_reference_genome.bash` script. <br /> <br />
The *Generate Pileup*,*SNP calling* and *Filtering for quality parameters* were run for all sequencing runs together. The *Filtering for quality parameters* was changed a lot and differed from the GB-eaSy pipline introduced by [Wickland et al. 2013](https://github.com/dpwickland/GB-eaSy). The  *Generate Pileup*,*SNP calling* and *Filtering for quality parameters* steps were conducted by the `Pileup_SNP_calling_Filtering.bash` script.
          
## 2 Filtering
After the VCF file was created with the GB-eaSy pipeline from [Wickland et al. 2013](https://github.com/dpwickland/GB-eaSy) and filtered for some quality parameters we did some additional filtering in R. All  listed *Filterung* steps were conducted by the `Filtering_for_coverage_average_RD_missingness.R` script. The Rscript contains the *Filtering* functions for the following filtering steps.  
### The filtering steps:
- Filtering for individual samples with low coverage
- Filtering for read depth per sample
- Filtering for missingness
- Removal of non-diallelic and non-polymorphic markers 
<br /> <br />
The Rscript generates a new filter VCF file, when all filtering steps are passed. The script additionally contains functions which calculate the  filtering parameters per marker, so that these can be plotted or checked, if needed. Furthermore the script also creates plots with the quality parameters **total read depth**, **mapping quality**, and **phred-scaled quality** before and after filtering. This plot is created with a plotting function from the `vcfR` package from [Knaus and Grünwald 2018]

### Loading of required packages and the VCF file:
Most of the function come from the `vcfR` package from [Knaus and Grünwald 2018](https://github.com/knausb/vcfR#:~:text=VcfR%20is%20an%20R%20package%20intended%20to%20allow,rapidly%20read%20from%20and%20write%20to%20VCF%20files.). The <vcfR> package from [Knaus and Grünwald 2018] works with S4 objects, therefore the syntax is quite different from "normal" R jargon. 
Most calculations and filtering steps over the entire set of markers were run by using the `mclapply()` function from the `parallel` package. We observed that the `mclapply()` function was extremly fast and efficient. 
```{r}
library(stringr)
library(data.table)
library(doMC)
library(vcfR)
library(parallel)
cores <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK'))
cores <- detectCores(all.tests = FALSE, logical = TRUE)
registerDoMC(cores)
setwd("/your/personal/working/directory/")
vcf_S4 <- read.vcfR("DATA.vcf", verbose = FALSE)
cat("VCF contains before anything was done:",nrow(vcf_S4@fix),"markers.","\n")
```

#### 2.1 Filtering for individual samples with low coverage
In our case, we choosed a minimum coverage of 10% per individual sample.
```{r}
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
```
### 2.2 Filtering for read depth per sample
In our case, we choosed a minimum average read depth of 1 and a maximum average read depth of 10 per marker. We choosed these thresholds based on the distribution of average read depth across all markers. 
```{r}
filter_for_average_read_depth <- function(data, min_av_depth, max_av_depth){
  cat("Filtering for average read depth started.","\n")
  cat("VCF file before filtering contains:",nrow(data@fix),"markers.","\n")
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
  return(data)
}
vcf_S4 <- filter_for_average_read_depth(data = vcf_S4, 
                                        min_av_depth = 1, 
                                        max_av_depth = 10)  
``` 
**Distribution of average read depth across all markers**
![GB1002 1_Average_read_depth](https://user-images.githubusercontent.com/63467079/149306207-62129755-9db6-473d-a1a9-dc2795ec84c6.png)

### 2.3 Filtering for missingness
When we filter for missingness, we filter in every subpopulation for at least 40 observations. So that only markers pass the threshold, which have 40 observations in every subpopulation.          
```{r}
filter_missingness_per_pop <- function(data, 
                                       population_1,
                                       population_2,
                                       population_3,
                                       population_4,
                                       max_missing_obs){
  cat("Filtering for missing observations per markers started.","\n")
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
  return(data)
}
vcf_S4_new <- filter_missingness_per_pop(data = vcf_S4, 
                           population_1 = "Shoepag_1",
                           population_2 = "Shoepag_2",
                           population_3 = "Shoepag_3",
                           population_4 = "Shoepag_4",
                           max_missing_obs = 40)
```
In these steps the average read depth per sample is calculated. This function needs only be applied, when you want to **save the information about individual coverage** e.g. for plotting. 
```{r}
calculate_missingness_per_pop <- function(data, 
                                          population_1,
                                          population_2,
                                          population_3,
                                          population_4){
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
  return(NAs_per_marker_per_pop)
}
NAs_per_marker_per_pop <- calculate_missingness_per_pop(data = vcf_S4,
                              population_1 = "Shoepag_1",
                              population_2 = "Shoepag_2",
                              population_3 = "Shoepag_3",
                              population_4 = "Shoepag_4")
```
With this command we also calculate the range of missingness between the populations. The range of missingness can be evaluated to look for markers which are more abundant in one population than in the others, which could skew the results.      
![GB1003_plot_missingness_02](https://user-images.githubusercontent.com/63467079/149474916-930cb3f0-bc73-4846-b334-f803ee950eb0.png)

### 2.4 Removal of non-diallelic and non-polymorphic markers
```{r}
filter_for_non_diallelic_non_polymorphic <- function(data){
  cat("Filtering for diallelic markers started.","\n")
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
  return(data)
}
vcf_S4 <- filter_for_non_diallelic_non_polymorphic(data = vcf_S4)
```
### Save the filtered VCF file
This function also comes from the `vcfR` package from [Knaus and Grünwald 2018](https://github.com/knausb/vcfR#:~:text=VcfR%20is%20an%20R%20package%20intended%20to%20allow,rapidly%20read%20from%20and%20write%20to%20VCF%20files.). 
```{r}
write.vcf(vcf_S4, file="path/to/your/working/directory/filtered_DATA.vcf.gz", mask=FALSE)          
```          
## 3 Selection signature mapping
The function for the calculation of the **<img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> leveraging replicated selection** and the **allele frequency differences** are contained in the `selection_signature_mapping.R script`. 
### 3.1 <img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> leveraging replicated selection
The function below will calculate the <img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> between all four subpopulations as: <br /> <br />
<img src="https://render.githubusercontent.com/render/math?math=F_{ST}=\frac{s^2}{\mu(p)*(1-\mu(p))%2B(\frac{s^2}{4})}"> 
<br /> <br /> according to [Weir and Cockerham, 1984](https://doi.org/10.1111/j.1558-5646.1984.tb05657.x). <br /> <br />
          
Additionally this function is calculating the absolute allele frequency difference between the subpopulations selected in opposite directions. The absolute allele frequency difference can be used to exclude markers with significant <img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> values, but which did diverge between the subpopulations selected in opposite directions. Significant <img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> values can be observed when one subpopulations was selected for specific environmental conditions. Therefore, we use this adjustment to identify truly selected sites.   
```{r}
calculate_FST_value <- function(data,
                                pop_low_phenotype_sel_1,
                                pop_low_phenotype_sel_2,
                                pop_high_phenotype_sel_1,
                                pop_high_phenotype_sel_2){
  time_0 <- Sys.time()
  cat("Calculation of the FST value started.","\n")
  gt <- extract.gt(data, element = "GT", as.numeric = FALSE)
  gt_file <- gt
  pop_1 <- gt_file[,which(str_detect(colnames(gt_file),pop_low_phenotype_sel_1))]
  pop_2 <- gt_file[,which(str_detect(colnames(gt_file),pop_low_phenotype_sel_2))]
  pop_3 <- gt_file[,which(str_detect(colnames(gt_file),pop_high_phenotype_sel_1))]
  pop_4 <- gt_file[,which(str_detect(colnames(gt_file),pop_high_phenotype_sel_2))]
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
  pop_1_dt <- get_allele_freq_per_marker_per_pop(data=pop_1, population_name = pop_low_phenotype_sel_1)
  pop_2_dt <- get_allele_freq_per_marker_per_pop(data=pop_2, population_name = pop_low_phenotype_sel_2)
  pop_3_dt <- get_allele_freq_per_marker_per_pop(data=pop_3, population_name = pop_high_phenotype_sel_1)
  pop_4_dt <- get_allele_freq_per_marker_per_pop(data=pop_4, population_name = pop_high_phenotype_sel_2)
  marker_table_per_pop <- cbind(pop_1_dt,pop_2_dt,
                                pop_3_dt,pop_4_dt)
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
                              population_4 = marker_table_per_pop[,21])
  selected_opposite_dir_1 <- as.numeric(marker_table_per_pop[,1] - marker_table_per_pop[,15])
  selected_opposite_dir_2 <- as.numeric(marker_table_per_pop[,8] - marker_table_per_pop[,21])
  od_stat <- abs(selected_opposite_dir_1 + selected_opposite_dir_2)
  FST_value_dt <- cbind(data@fix[,c(1,2)],SNP_IDs,
                        FST_value,od_stat)
  FST_value_dt <- as.data.frame(FST_value_dt)
  colnames(FST_value_dt) <- c("Chromosome","Position","SNP_ID",
                              "FST_value","Abs_allele_freq_diff")
  cat("Calculation of the FST value is done.","\n")
  time_1 <- Sys.time()
  cat("This was done within",print(time_1 - time_0),"\n")
  return(FST_value_dt)
}
FST_value_all_pop <- calculate_FST_value(data = vcf_S4, 
                                         pop_low_phenotype_sel_1 = "Shoepag_1",
                                         pop_low_phenotype_sel_2 = "Shoepag_4",
                                         pop_high_phenotype_sel_1 = "Shoepag_3",
                                         pop_high_phenotype_sel_2 = "Shoepag_2")
```

### 3.2 Allele frequency differences
The function below will calculate the allele frequency differences between the subpopulations selected in the same and opposite directions and the absolute allele frequency difference. The absolute allele frequency difference is calculated as: <br /> <br />
<img src="https://render.githubusercontent.com/render/math?math=%7C(p_{Low1}%2Dp_{High1})%2B(p_{Low2}%2Dp_{High2})%7C">
<br /> 
whereas: <br /> <img src="https://render.githubusercontent.com/render/math?math=Low1"> = Population 1 selected for the low phenotype <br />
        <img src="https://render.githubusercontent.com/render/math?math=Low1"> = Population 2 selected for the low phenotype <br />
        <img src="https://render.githubusercontent.com/render/math?math=High1"> = Population 1 selected for the high phenotype <br />
        <img src="https://render.githubusercontent.com/render/math?math=High2"> = Population 2 selected for the high phenotype <br />
according to [Turner and Miller (2012)](http://www.genetics.org/content/suppl/2012/03/30/genetics.112.139337.DC1). <br />   
```{r}
calculate_the_allele_freq_diff <- function(data,
                                           pop_low_phenotype_sel_1,
                                           pop_low_phenotype_sel_2,
                                           pop_high_phenotype_sel_1,
                                           pop_high_phenotype_sel_2){
  time_0 <- Sys.time()
  cat("Calculation of the allele frequency differences started.","\n")
  gt <- extract.gt(data, element = "GT", as.numeric = FALSE)
  gt_file <- gt
  pop_1 <- gt_file[,which(str_detect(colnames(gt_file),pop_low_phenotype_sel_1))]
  pop_2 <- gt_file[,which(str_detect(colnames(gt_file),pop_low_phenotype_sel_2))]
  pop_3 <- gt_file[,which(str_detect(colnames(gt_file),pop_high_phenotype_sel_1))]
  pop_4 <- gt_file[,which(str_detect(colnames(gt_file),pop_high_phenotype_sel_2))]
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
  pop_1_dt <- get_allele_freq_per_marker_per_pop(data=pop_1, population_name = pop_low_phenotype_sel_1)
  pop_2_dt <- get_allele_freq_per_marker_per_pop(data=pop_2, population_name = pop_low_phenotype_sel_2)
  pop_3_dt <- get_allele_freq_per_marker_per_pop(data=pop_3, population_name = pop_high_phenotype_sel_1)
  pop_4_dt <- get_allele_freq_per_marker_per_pop(data=pop_4, population_name = pop_high_phenotype_sel_2)
  marker_table_per_pop <- cbind(pop_1_dt,pop_2_dt,
                                pop_3_dt,pop_4_dt)
  marker_table_per_pop <- as.data.frame(marker_table_per_pop)
  selected_opposite_dir_1 <- as.numeric(marker_table_per_pop[,1] - marker_table_per_pop[,15])
  selected_opposite_dir_2 <- as.numeric(marker_table_per_pop[,8] - marker_table_per_pop[,21])
  selected_same_dir_1 <- as.numeric(marker_table_per_pop[,1] - marker_table_per_pop[,8])
  selected_same_dir_2 <- as.numeric(marker_table_per_pop[,15] - marker_table_per_pop[,21])
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
```
## 4 Significance thresholds
### 4.1 Based on the empiric distribution


### 4.2 Based on drift simulations 

### 4.3 Simulation of Drift
          
### 4.4 Based on the FDR for selection





## 5 Sauron plot
          
![GB1006_Sauron_plots_combined_values_2021_12_17](https://user-images.githubusercontent.com/63467079/149146525-ce94e222-dff8-4ad4-8dcb-ad14f7530032.png)

## 6 Other plotting scripts

