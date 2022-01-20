# Selection signature mapping scripts with replicated selection
## Table of contents
[0 Introduction](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#introduction) <br />
[1 Phenotypic data analysis](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#1-phenotypic-data-analysis) <br />
2 Pipeline for the analysis of GBS data adapted from [Wickland et al. 2013](https://github.com/dpwickland/GB-eaSy) <br />
[3 Filtering](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#filtering) <br />
&emsp;[3.1 Filtering for individual samples with low coverage](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#21-filtering-for-individual-samples-with-low-coverage) <br />
&emsp;[3.2 Filtering for read depth per sample](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#22-filtering-for-read-depth-per-sample) <br />
&emsp;[3.3 Filtering for missingness](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#23-filtering-for-missingness) <br />
&emsp;[3.4 Removal of non-diallelic and non-polymorphic markers](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#24-removal-of-non-diallelic-and-non-polymorphic-markers) <br />
[4 Selection signature mapping](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#3-selection-signature-mapping) <br />
&emsp; [4.1 FST leveraging replicated selection](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#31--leveraging-replicated-selection) <br />
&emsp; [4.2 Allele frequency differences](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#32-allele-frequency-differences) <br />
[5 Significance thresholds](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#4-significance-thresholds) <br />
&emsp; [5.1 Based on the empirical distribution](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#41-based-on-the-empiric-distribution) <br />
&emsp; [5.2 Based on drift simulations](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#42-based-on-drift-simulations) <br />
&emsp; [5.3 Simulation of drift](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#43-simulation-of-drift) <br />
&emsp; [5.4 Based on the FDR for selection](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#44-based-on-the-fdr-for-selection)<br />
[6 Sauron plot](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#5-sauron-plot) <br />
[7 Manhatten plots](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#6-other-plotting-scripts) <br />


## 0 Introduction
This repository contains scripts for selection signature mapping with replicated selection.
Usually in experimental evolution studies a random-mating population is selected in divergent directions for one
trait of interest. The highly divergent subpopulations are then scanned for signatures of selection.
A previous problem in these studies was the determination of significance thresholds. 
One solution to this problem is the usage of **replicated selection**. 
This was introduced by [Turner and Miller (2012)](http://www.genetics.org/content/suppl/2012/03/30/genetics.112.139337.DC1), but only tested on *Drospohila melanogaster*. <br /> <br />
In this approach two subpopulations are selected in the same direction. 
By that we differentiate between changes caused by selection and changes caused by other factors like drift.
The **replicated selection** is used to calculate a *FDR for selection*. 
[Turner and Miller (2012)](http://www.genetics.org/content/suppl/2012/03/30/genetics.112.139337.DC1) applied this significance threshold only on the
statistic based on the allele frequency differences. 
Whereas we applied this new significance threshold to the commonly used 
          <img src="https://render.githubusercontent.com/render/math?math=F_{ST}">
statistic. <br /> <br /> 
**Experimental design of an experimental evolution study with replicated selection** <br />
<img src="https://user-images.githubusercontent.com/63467079/150107509-3f984ca3-3a61-4338-87d5-d920a5673727.png" width="500" height="555.7324841">

## 1 Phenotypic data analysis
The phenotypic data analysis was conducted to evaluate the effect of selection on the phenotype. Therefore, the trait measurements in the subpopulations selected in opposite directions can be compared. We conducted a t-test between the subpopulations selected in opposite directions to test for significance.  <br /> 
```{r}

```
Additionally we also used the trait measurments from the base population and all generations of selection to show the decrease and increase in the selected trait. In our case the selected trait was plant height. <br /> <br /> 
**Measured plant height in all years in the subpopulations selected for short plant (green) and tall plant height (purple).**
<img src="https://user-images.githubusercontent.com/63467079/150105713-a27b5365-4822-483e-abe4-6fff10332bc7.png" width="600" height="360">

## 2 Pipeline for the analysis of GBS data adapted from [Wickland et al. 2013](https://github.com/dpwickland/GB-eaSy)
For the analysis of our raw reads from paired-end genotyping-by-sequencing (GBS) with ApeKI according to [Elshire et al. 2011](https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0019379&type=printable), we used the GB-eaSy pipeline from [Wickland et al. 2013](https://github.com/dpwickland/GB-eaSy). <br />
The pipeline consists out of several steps, which comprises:
- [Demultiplexing and trimming of the adapter sequence](https://github.com/dpwickland/GB-eaSy#step-2-demultiplex-raw-reads)
- [Alignement to the reference genome](https://github.com/dpwickland/GB-eaSy#step-3-align-to-reference)
- [Create a list of sorted bam files](https://github.com/dpwickland/GB-eaSy#step-4-create-list-of-sorted-bam-files)
- [Generate pileup and SNP calling](https://github.com/dpwickland/GB-eaSy#step-5-and-6-generate-pileup-and-call-snps)
- Filtering for quality parameters <br /> <br />

This bash-script was run on every sequenced plate separately, since every sequenced plate had it's own adapters and the same set of barcodes was used for all sequenced plates in our case. The *Demultiplexing* and *Alignement to the reference genome* steps were conducted by the `GB_eaSy_demultiplexing_and_alignemnt_S1.bash`,`GB_eaSy_demultiplexing_and_alignemnt_S2`,`GB_eaSy_demultiplexing_and_alignemnt_S3` and `GB_eaSy_demultiplexing_and_alignemnt_S4` scripts. <br /> <br />
The *Generate Pileup*,*SNP calling* and *Filtering for quality parameters* were run for all sequencing runs together. The *Filtering for quality parameters* was changed a lot and differed from the GB-eaSy pipline introduced by [Wickland et al. 2013](https://github.com/dpwickland/GB-eaSy). The  *Generate Pileup*,*SNP calling* and *Filtering for quality parameters* steps were conducted by the `GB_eaSy_pileup_SNP_calling_filtering.bash` script.
          
## 3 Filtering
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

### 3.1 Filtering for individual samples with low coverage
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
### 3.2 Filtering for read depth per sample
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
#### Distribution of average read depth across all markers
<img src="https://user-images.githubusercontent.com/63467079/149306207-62129755-9db6-473d-a1a9-dc2795ec84c6.png" width="500" height="300">
          
### 3.3 Filtering for missingness
When we filter for missingness, we filter in every subpopulation for at least 40 observations. So that only markers pass the threshold, which have 40 observations in every subpopulation.          
```{r}
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
<img src="https://user-images.githubusercontent.com/63467079/149474916-930cb3f0-bc73-4846-b334-f803ee950eb0.png" width="800" height="320">

### 3.4 Removal of non-diallelic and non-polymorphic markers
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
## 4 Selection signature mapping
The function for the calculation of the **<img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> leveraging replicated selection** and the **allele frequency differences** are contained in the `selection_signature_mapping.R` script. 
### 4.1 <img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> leveraging replicated selection
The function below will calculate the <img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> between all four subpopulations as: <br /> <br />
<img src="https://render.githubusercontent.com/render/math?math=F_{ST}=\frac{s^2}{\mu(p)*(1-\mu(p))%2B(\frac{s^2}{4})}"> 
<br /> <br /> according to [Weir and Cockerham, 1984](https://doi.org/10.1111/j.1558-5646.1984.tb05657.x). <br /> <br />

The function also calculates the <img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> only between the subpopulations selected in same direction and only between the subpopulations selected in opposite directions. This values are required for the calculation of the FDR for selection. <br /> <br />            
 
Additionally this function is calculating the absolute sum of <img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> between the subpopulations selected in opposite directions. This value can be used to exclude markers with a significant <img src="https://render.githubusercontent.com/render/math?math=F_{ST}">, but which diverge between subpopulation selected in the same direction. This could have happend due to selection for specific environmental conditions. Therefore, we use this adjustment to identify truly selected sites. <br /> <br />

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
```

### 4.2 Allele frequency differences
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
```                    
## 5 Significance thresholds
The calculation of significance thresholds and the plotting functions are contained in the `Significance_thresholds_and_plotting.R`. The `Filtering_for_coverage_average_RD_missingness.R` and `selection_signature_mapping.R script` can or should be run on a inactive Linux session on a high-throughput computing device. These scripts are usually run on extremly large data sets (raw sequence data or large VCF files). The `Significance_thresholds_and_plotting.R` is usually run on a much smaller data set, since many markers were removed in the filtering procedure. Furthermore, when windows font types want to be used,  the script needs to be run on a windows device. <br />   
### 5.1 Based on the empirical distribution
The significance thresholds based on the empirical distribution, were calculated by taking the 99.9th and 99.99th percentile of the empirical distribution of the <img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> and absolute allele frequency differences:
```{r}
FST_sig_thres_1 <- quantile(FST_values_od_cor$Fst, probs = 0.9999, na.rm = TRUE)
FST_sig_thres_2 <- quantile(FST_values_od_cor$Fst, probs = 0.999, na.rm = TRUE)    
AFD_sig_thres_1 <- quantile(allele_freq_diff$DiffStat, probs = 0.9999, na.rm = TRUE)
AFD_sig_thres_2 <- quantile(allele_freq_diff$DiffStat, probs = 0.999, na.rm = TRUE)
```          
The significance threshold is stored, so it can be used later directly for plotting. <br />
### 5.2 Based on drift simulations 
The significance thresholds based on drift simulations were calculated in the `Simulation_of_drift.R` and then only retrieved from this script. The simulation of drift is described below and the script is also available in the repository.
```{r}
AFD_drift_sim_sig_thres <- 0.4686
FST_drift_sim_sig_thres <-  0.3540692
```    
### 5.3 Simulation of Drift
The simulation of drift was conducted by using the `DriftSimulator.R` from [Beissinger (2021)](http://beissingerlab.github.io/Software/). The `DriftSimulator.R` script was run with a drift simulation script from [Kumar et al., 2021](https://academic.oup.com/pcp/article/62/7/1199/6279219), which enables the implementation of the drift simulator over a large set of markers. The script, which enables the simulation of drift over a large set of markers is available as `Run_drift_simulator.R`. We ran the `DriftSimulator.R` over a data set consisting out of 4,226,822 simulated SNP markers <br /> <br />
```{r}

```
The drift simulator samples the initial allele frequencies for 10 cycles, with a population size of 20,000 individuals with equal number of males and females to simulate random mating under "normal" conditions. Afterwards drift is simulated with a lower number of males and females to simulate the field conditions subdivision into subpopulations which led to limited spatial pollen distribution for three generations. Finally, sampling of 96 individuals out of 5000 for genotyping is simulated as well and included into the drift simulator [Turner et al., 2011](http://www.genetics.org/content/suppl/2012/03/30/genetics.112.139337.DC1). Additionally, the marker coverage was also sampled, the minimal marker coverage was set to at least 40 out of 96 observations at a marker, so that the marker coverage was always sampled between 40 to 96 observations per marker. The script can be repeatably run. We ran the script for 10 simulations. <br /> <br /> <img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> and allele frequency differences were calculated for all simulated markers, which corresponded in our case to 42,000,000 markers. We choosed then among the 42,000,000 markers the most extreme value as significance threshold.

### 5.4 Based on the FDR for selection
**FDR for selection for all possible values of the statistics**
The function `calculate_FDR_for_selection()` generates a table to show all posible values of a statistic and the number of observed markers diverged between subpopulation selected in the same and opposite directions at a certain value. The FDR for selection is received by dividing the number of observed markers diverged between subpopulation selected in the same direction by the number of observed markers diverged between subpopulation selected in opposite directions. <br /> <br />
```{r}

```
Even though, the significance threshold based on the FDR for selection is already printed by the previous function, the following function returns the threshold so it can be used directly in the Manhatten plot or to check the overlap between all the different statistics. 
```{r}


```
## 6 Sauron plot
The Sauron plot from [Turner and Miller (2012)](http://www.genetics.org/content/suppl/2012/03/30/genetics.112.139337.DC1) can be used to visualize the scope of drift and the scope of selection. Sauron plots of genetic differentiation are created by plotting the FST statistics (A) and the allele frequency differences (B) observed in the subpopulations selected in the same direction (blue) and in opposite directions (red) against each other at each SNP marker. The transparent red colored edges correspond to a false discovery rate (FDR) for selection <10%. The diverged markers observed in the subpopulations selected in the same direction (blue) are provoked by drift and other factors, but not by selection. The diverged markers observed in the subpopulations selected in opposite directions (red) are provoked by drift, other factors and selection. Therefore the observations which exceed the cloud of blue points are expected to be provoked only by selection. The y- and x-axis correspond to the range of <img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> and allele frequency differences.  
<img src="https://user-images.githubusercontent.com/63467079/149146525-ce94e222-dff8-4ad4-8dcb-ad14f7530032.png" width="700" height="350">
The Sauron plots for the allele frequency differences and <img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> are created by two different functions, since they plotting areas are different from each other. The "Sauron plot" for the <img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> does not even look like the eye of [Sauron](https://twitter.com/strnr/status/457201981007073280) anymore. The Sauron plot for the allele frequency difference is created by the `create_sauron_plot_allele_freq_diff()` function, which is shown above:
```{r}
create_sauron_plot_allele_freq_diff <- function(data_table,
                                                sig_threshold,
                                                stat_opposite_dir1,
                                                stat_opposite_dir2,
                                                stat_same_dir1,
                                                stat_same_dir2){
  windowsFonts(my = windowsFont('Calibri'))
  create_FDR_below_10_per_neg_area <- function(x){
    sig_threshold[1]+-1*x
  }
  y_coord_value <- create_FDR_below_10_per_neg_area(x = seq(-1,0,0.01))
  dt_coords_neg_area <- cbind(y_coord_value,seq(-1,0,0.01))
  dt_coords_neg_area <- as.data.frame(dt_coords_neg_area)
  dt_coords_neg_area$z <- rep(sig_threshold[1],nrow(dt_coords_neg_area))
  colnames(dt_coords_neg_area) <- c("x_coord_value","y_coord_value","threshold_value")
  dt_coords_neg_area <- dt_coords_neg_area[which(dt_coords_neg_area$y_coord_value > sig_threshold[1]),]
  #### Negative area of the Sauron plot
  create_FDR_below_10_per_pos_area <- function(x){
    sig_threshold[2]+-1*x
  }
  y_coord_value <- create_FDR_below_10_per_pos_area(x = 1-seq(0,1,0.01))
  dt_coords_pos_area <- cbind(y_coord_value,1-seq(0,1,0.01))
  dt_coords_pos_area <- as.data.frame(dt_coords_pos_area)
  dt_coords_pos_area$z <- rep(sig_threshold[2],nrow(dt_coords_pos_area))
  colnames(dt_coords_pos_area) <- c("x_coord_value","y_coord_value","threshold_value")
  dt_coords_pos_area <- dt_coords_pos_area[which(dt_coords_pos_area$y_coord_value < sig_threshold[2]),]
  sauron_1 <- ggplot()+
    geom_area(data = dt_coords_pos_area, aes(x=x_coord_value, y=threshold_value), alpha = 0.1, fill = "tomato")+
    geom_area(data = dt_coords_pos_area, aes(x=x_coord_value, y=y_coord_value), fill = "white")+
    geom_area(data = dt_coords_neg_area, aes(x=x_coord_value, y=threshold_value), alpha = 0.1, fill = "tomato")+
    geom_area(data = dt_coords_neg_area, aes(x=x_coord_value, y=y_coord_value),  fill = "white")+
    geom_point(data = data_table, aes(x = stat_same_dir1,y = stat_same_dir1),shape=1, colour = "royalblue3")+
    geom_point(data = data_table, aes(x = stat_opposite_dir1,y = stat_opposite_dir2),shape=1, colour = "firebrick")+
    theme(text = element_text(size =14, family = "my"),
          panel.background = element_rect(fill = "white"),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_text(size =14, family = "my", colour = "royalblue3", face = "bold"),
          axis.title.x = element_text(size =14, family = "my", colour = "royalblue3", face = "bold"),
          axis.line.x = element_line(color = "black"),
          axis.text = element_text(size =14, family = "my", colour = "black"))+
    labs(x=paste("\t","\t","\t","\t","\t","\t","\t","\t","Low pheno 1 vs Low pheno 2"),
         y=paste("\t","\t","\t","High pheno 1 vs High pheno 2"))+
    ylim(-1,1)+
    xlim(-1,1)
  d11 <-  ggplot_build(sauron_1)$data[[1]]
  empty_plot <- ggplot()+
    theme(panel.background = element_rect(fill = "white"))
  get_xaxis<-function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    legend <- tmp$grobs[[12]]
    return(legend)
  }
  xaxis_sauron <- get_xaxis(sauron_1)
  get_yaxis<-function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    legend <- tmp$grobs[[13]]
    return(legend)
  }
  yaxis_sauron <- get_yaxis(sauron_1)
  sauron_2 <- ggplot()+
    geom_area(data = dt_coords_pos_area, aes(x=x_coord_value, y=threshold_value), alpha = 0.1, fill = "tomato")+
    geom_area(data = dt_coords_pos_area, aes(x=x_coord_value, y=y_coord_value), fill = "white")+
    geom_area(data = dt_coords_neg_area, aes(x=x_coord_value, y=threshold_value), alpha = 0.1, fill = "tomato")+
    geom_area(data = dt_coords_neg_area, aes(x=x_coord_value, y=y_coord_value),  fill = "white")+
    geom_point(data = data_table, aes(x = stat_opposite_dir1,y = stat_opposite_dir2),shape=1, colour = "firebrick")+
    geom_point(data = data_table, aes(x = stat_same_dir1,y = stat_same_dir2),shape=1, colour = "royalblue3")+
    theme(text = element_text(size =14, family = "my"),
          panel.background = element_rect(fill = "white"),
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                          colour = "grey81"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "grey100"),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_text(size =14, family = "my", colour = "firebrick", face = "bold"),
          axis.title.x = element_text(size =14, family = "my", colour = "firebrick", face = "bold"),
          axis.line.x = element_line(color = "black"),
          axis.text = element_text(size =14, family = "my", colour = "black"))+
    labs(x=  "Low pheno 1 vs High pheno 1",
         y= "Low pheno 2 vs High pheno 2",
         tag = "B")+
    geom_abline(intercept = sig_threshold[1], slope = -1, size = 1, alpha = 0.2, colour = "firebrick")+
    geom_abline(intercept = sig_threshold[2], slope = -1, size = 1, alpha = 0.2, colour = "firebrick")+
    xlim(sig_threshold[1],sig_threshold[2])+
    ylim(sig_threshold[1],sig_threshold[2])
  sauron_2
  sauron_plot_all_afd <- grid.arrange(yaxis_sauron,sauron_2, xaxis_sauron,empty_plot,
                                      ncol=2, nrow = 2, 
                                      layout_matrix = rbind(c(1,2),c(4,3)),
                                      widths = c(0.2, 4), heights = c(4,0.2))
  return(sauron_plot_all_afd)
}
create_sauron_plot_allele_freq_diff(data_table = allele_freq_diff,
                                    sig_threshold = AFD_FDR_for_sel_sig_thres,
                                    stat_opposite_dir1 = allele_freq_diff$Short1_vs_Tall1,
                                    stat_opposite_dir2 = allele_freq_diff$Short2_vs_Tall2,
                                    stat_same_dir1 = allele_freq_diff$Short1_vs_Short2,
                                    stat_same_dir2 = allele_freq_diff$Tall1_vs_Tall2)
```
## 7 Manhatten plots
In these [Manhatten plot](https://en.wikipedia.org/wiki/Manhattan_plot#:~:text=A%20Manhattan%20plot%20is%20a%20type%20of%20scatter,genome-wide%20association%20studies%20%28GWAS%29%20to%20display%20significant%20SNPs.) the positions of the markers are plotted against the statistical value observed at this marker. Here the Manhatten plots are only shown for Chromosome 3, which showed in our study the most interesting results. The different significance thresholds are also plotted. <br /> 
**Differentiation along Chromosome 3 expressed by the <img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> (A) and absolute allele frequency differences (B) 
![2021_new_combined_manhatten_plots](https://user-images.githubusercontent.com/63467079/149943665-fdedbd0f-cb06-44f4-a2d2-c3ffd53b258a.png)
<br /> 
Because we plot the different signifcance thresholds and we observed quite much recombination the Manhatten plots were created for the different chromosomes separately. Therefore the data set first needs to divided into the markers observed at the different chromosomes. This can be easly done by a function from the `data.table` package.  
```{r}
Fst_tab <- as.data.table(FST_values_od_cor)
Fst_Chr_1 <- Fst_tab[Chromosome=="chr1",]
Fst_Chr_2 <- Fst_tab[Chromosome=="chr2",]
Fst_Chr_3 <- Fst_tab[Chromosome=="chr3",]
Fst_Chr_4 <- Fst_tab[Chromosome=="chr4",]
Fst_Chr_5 <- Fst_tab[Chromosome=="chr5",]
Fst_Chr_6 <- Fst_tab[Chromosome=="chr6",]
Fst_Chr_7 <- Fst_tab[Chromosome=="chr7",]
Fst_Chr_8 <- Fst_tab[Chromosome=="chr8",]
Fst_Chr_9 <- Fst_tab[Chromosome=="chr9",]
Fst_Chr_10 <- Fst_tab[Chromosome=="chr10",]          
```
The combined Manhatten plots for both statistics are created per chromosome separately by the following function:
```{r}
          
```
In case one wants to plot only one statistic in the Manhatten plot, this function creates the only one Manhatten plot:
```{r}
create_manhatten_plot_per_chr_one_stat <- function(data,
                                          statistic, 
                                          SNP_position,
                                          sign_thres_emp_dis_999,
                                          sign_thres_emp_dis_9999,
                                          sign_thres_drift_sim,
                                          sign_thres_FDR){
  windowsFonts(my = windowsFont('Calibri'))
  chr_plot <- ggplot(data = data)+
    geom_point(aes(x = SNP_position,y = statistic), colour = "black")+
    theme(text = element_text(size =12, family = "my"),
          panel.background = element_rect(fill = "white"),
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                          colour = "grey81"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "grey100"),
          panel.border = element_blank(),
          axis.line.y = element_line(colour = "black"),
          axis.line.x = element_line(colour = "black"),
          axis.ticks = element_blank(),
          axis.title.y = element_text(size =12, family = "my", colour = "black", face = "bold"),
          axis.text.y = element_text(size =12, family = "my", colour = "black"),
          axis.title.x = element_text(size =12, family = "my", colour = "black", face = "bold"),
          axis.text.x = element_blank(),
          legend.text = element_text(size =10, family = "my", colour = "black"),
          legend.title = element_text(size =10, family = "my", colour = "black", face = "bold"))+
    labs(y = "Fst value",x = "Position", tag = "A")+
  geom_hline(aes(yintercept = sign_thres_emp_dis_999, linetype=str_wrap("99.9th percentile",10)), colour="darkred", size = 1)+
    geom_hline(aes(yintercept = sign_thres_emp_dis_9999, linetype=str_wrap("99.99th",10)), colour="tomato", size = 1)+
    geom_hline(aes(yintercept = drift_simulation_based_sig_thresold, linetype=str_wrap("drift simulations",10)), colour="purple", size = 1)+
    geom_hline(aes(yintercept = FDR_of_selection_significance_threshold, linetype=str_wrap("FDR",10)), colour="gold", size = 1)+
    scale_linetype_manual(name = str_wrap("Significance threshold",10), values = c(1,1,1,1), 
                          guide = guide_legend(override.aes = list(color = c("darkred", "tomato",
                                                                             "purple","gold"))))
  return(chr_plot)
}
manhatten_plot <- create_manhatten_plot_per_chr_one_stat(data = Fst_Chr_1,
                              statistic = Fst_Chr_1$Fst,
                              SNP_position = Fst_Chr_1$Position,
                              sign_thres_emp_dis_999 = sign_thres_emp_dis_999,
                              sign_thres_emp_dis_9999 = sign_thres_emp_dis_9999,
                              sign_thres_drift_sim = sign_thres_drift_sim,
                              sign_thres_FDR = sign_thres_FDR)       
```
