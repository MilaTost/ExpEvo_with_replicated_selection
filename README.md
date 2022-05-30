# Experimental evolution with replicated selection
## Table of contents
[0 Introduction](https://github.com/milaleonie/ExpEvo_with_replicated_selection/blob/main/README.md#introduction) <br />
[1 Phenotypic data analysis](https://github.com/milaleonie/ExpEvo_with_replicated_selection/blob/main/README.md#1-phenotypic-data-analysis) <br />
2 Pipeline for the analysis of GBS data adapted from [Wickland et al. 2013](https://github.com/dpwickland/GB-eaSy) <br />
[3 Filtering](https://github.com/milaleonie/ExpEvo_with_replicated_selection/blob/main/README.md#filtering) <br />
&emsp;[3.1 Filtering for individual samples with low coverage](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#21-filtering-for-individual-samples-with-low-coverage) <br />
&emsp;[3.2 Filtering for read depth per sample](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#22-filtering-for-read-depth-per-sample) <br />
&emsp;[3.3 Filtering for missingness](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#23-filtering-for-missingness) <br />
&emsp;[3.4 Removal of non-diallelic and non-polymorphic markers](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#24-removal-of-non-diallelic-and-non-polymorphic-markers) <br />
[4 Selection signature mapping with FST leveraging replicated selection](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#3-selection-signature-mapping) <br />

[5 Significance thresholds](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#4-significance-thresholds) <br />
&emsp; [5.1 Based on the empirical distribution](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#41-based-on-the-empiric-distribution) <br />
&emsp; [5.2 Based on drift simulations](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#42-based-on-drift-simulations) <br />
&emsp; [5.3 Simulation of drift](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#43-simulation-of-drift) <br />
&emsp; [5.4 Based on the FDRfS](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#44-based-on-the-fdr-for-selection)<br />
[6 Sauron plot](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#5-sauron-plot) <br />
[7 Manhatten plots](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#6-other-plotting-scripts) <br />
[8 LD decay](https://github.com/milaleonie/Selection_signature_mapping_with_replicated_selection/blob/main/README.md#8-ld-decay)

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
Short_plants_2016 <- data[PlantHeight_group == "Selected for short plant height" & year == "2016",]
Tall_plants_2016 <- data[PlantHeight_group == "Selected for tall plant height" & year == "2016",]

Short_plants_2020 <- data[PlantHeight_group == "Selected for short plant height" & year == "2020",]
Tall_plants_2020 <- data[PlantHeight_group == "Selected for tall plant height" & year == "2020",]

t.test(Short_plants_2020$PlantHeight,Short_plants_2016$PlantHeight,
       alternative = "less",
       paired = TRUE)
t.test(Tall_plants_2020$PlantHeight,Tall_plants_2016$PlantHeight,
       alternative = "greater",
       paired = TRUE)
t.test(Short_plants_2020$PlantHeight,Tall_plants_2020$PlantHeight,
       alternative = "less",
       paired = TRUE)
```
Additionally we also used the trait measurments from the base population and all generations of selection to show the decrease and increase in the selected trait. In our case the selected trait was plant height. <br /> <br /> 
**Measured plant height in all years in the subpopulations selected for short plant (green) and tall plant height (purple).**
<img src="https://user-images.githubusercontent.com/63467079/150105713-a27b5365-4822-483e-abe4-6fff10332bc7.png" width="600" height="360"> <br />
The computation of the t-test statistic and the script for plotting the measured phenotypes across all years are available in the `Phenotypic_data_analysis.R` script. 
       
## 2 Pipeline for the analysis of GBS data adapted from [Wickland et al. 2013](https://github.com/dpwickland/GB-eaSy)
For the analysis of our raw reads from paired-end genotyping-by-sequencing (GBS) with ApeKI according to [Elshire et al. 2011](https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0019379&type=printable), we used the GB-eaSy pipeline from [Wickland et al. 2013](https://github.com/dpwickland/GB-eaSy). <br />
The pipeline consists out of several steps, which comprises:
- [Demultiplexing and trimming of the adapter sequence](https://github.com/dpwickland/GB-eaSy#step-2-demultiplex-raw-reads)
- [Alignement to the reference genome](https://github.com/dpwickland/GB-eaSy#step-3-align-to-reference)
- [Create a list of sorted bam files](https://github.com/dpwickland/GB-eaSy#step-4-create-list-of-sorted-bam-files)
- [Generate pileup and SNP calling](https://github.com/dpwickland/GB-eaSy#step-5-and-6-generate-pileup-and-call-snps)
- Filtering for quality parameters <br /> <br />

This bash-script was run on every sequenced plate separately, since every sequenced plate had it's own adapters and the same set of barcodes was used for all sequenced plates in our case. The *Demultiplexing* and *Alignement to the reference genome* steps were conducted by the `GB_eaSy_demultiplexing_and_alignemnt_S1.bash`,`GB_eaSy_demultiplexing_and_alignemnt_S2.bash`,`GB_eaSy_demultiplexing_and_alignemnt_S3.bash` and `GB_eaSy_demultiplexing_and_alignemnt_S4.bash` scripts. <br /> <br />
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
## 4 Scan for selection signatures
Our scan for selection was based on the <img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> leveraging replicated selection. 
The function for the calculation of the **<img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> leveraging replicated selection** is contained in the `selection_signature_mapping.R` script. 
<br /> <br /> 
The function below will calculate the <img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> between all possible and non-redundant comparisons between the subpopulations selected in the same and opposite directions at each marker as: <br /> <br />
<img src="https://render.githubusercontent.com/render/math?math=F_{ST} = \frac{s^2}{\mu(p)*(1-\mu(p))%2B(\frac{s^2}{2})}"> 
<br /> <br /> according to [Weir and Cockerham, 1984](https://doi.org/10.1111/j.1558-5646.1984.tb05657.x). <br /> <br />

The function also calculates the <img src="https://render.githubusercontent.com/render/math?math=F_{ST} Sum"> between the subpopulations selected in same and opposite directions. This values are required for the calculation of the false discovery rate for selection (FDRfS). In this statistic, observations which only occured in one comparison are excluded. Those observations might have been caused by drift. Selection is a repeatable force, so that we should be able to observe the same pattern in both comparisons. <br /> <br />            
 
```{r}
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
cat("The data set is loaded.","\n")
time_1 <- Sys.time()
print(time_1 - time_0)
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
```
          
## 5 Significance thresholds
Significance thresholds for selection were calculated three ways: 1) based on the empirical distribution; 2) based on drift simulations; and 3) based on the false discovery rate for selection (FDRfS) [Turner and Miller (2012)](http://www.genetics.org/content/suppl/2012/03/30/genetics.112.139337.DC1). <br /> 
The calculation of significance thresholds and the plotting functions are contained in the `Significance_thresholds_and_plotting.R`. The `Filtering_for_coverage_average_RD_missingness.R` and `Selection_signature_mapping.R` script can or should be run on a inactive Linux session on a high-throughput computing device. These scripts are usually run on extremly large data sets (raw sequence data or large VCF files). The `Significance_thresholds_and_plotting.R` is usually run on a much smaller data set, since many markers were removed in the filtering procedure. Furthermore, when windows font types want to be used,  the script needs to be run on a windows device. <br />   
### 5.1 Based on the empirical distribution
The significance thresholds based on the empirical distribution, were calculated by taking the 99.9th and 99.99th percentile of the empirical distribution of the <img src="https://render.githubusercontent.com/render/math?math=F_{ST}">:
```{r}
FST_sig_thres_1 <- quantile(FST_values_od_cor$Fst, probs = 0.9999, na.rm = TRUE)
FST_sig_thres_2 <- quantile(FST_values_od_cor$Fst, probs = 0.999, na.rm = TRUE)    
```          
The significance threshold is stored, so it can be used later directly for plotting. <br />
### 5.2 Based on drift simulations 
The significance thresholds based on drift simulations were calculated in the `Simulation_of_drift.R` and then only retrieved from this script. The simulation of drift is described below and the script is also available in the repository.
  
### 5.3 Simulation of Drift
The simulation of drift was conducted by using the `DriftSimulator.R` from [Beissinger (2021)](http://beissingerlab.github.io/Software/). The `DriftSimulator.R` script was run with a drift simulation script similar to the one from [Kumar et al., 2021](https://academic.oup.com/pcp/article/62/7/1199/6279219), which enables the implementation of the drift simulator over a large set of markers. The script, which enables the simulation of drift over a large set of markers is available as `Run_drift_simulator_over_many_markers.R`. The `Run_drift_simulator_over_many_markers.R` script contains a function which simulates drift acting a single locus. We ran 5,000,000 simulations. For each simulation, initial allele frequencies were set based on allele frequency spectrum observed in generation 0 [Gyawali et al., 2019](https://pubmed.ncbi.nlm.nih.gov/31590656/) . Drift was simulated with a population size of 5000 individuals with 250 female and 5000 male individuals for three generations. We assumed that every ear contributed ~500 kernels. Every kernel could have been pollinated by one of the male parents, which resulted in much higher harvest than 5000 kernels. Therefore, 5000 kernels were randomly drawn from the entire harvest to represent the seeds planted for the next generation. In the third generation, 96 individuals out of 5000 were sampled to represent the individuals that were actually genotyped (Turner et al., 2011; Kumar et al., 2021). Variable marker coverage was also simulated; marker coverage was sampled from a uniform distribution between 40 and 96 observations per marker to match our filtering process [Turner et al., 2011](http://www.genetics.org/content/suppl/2012/03/30/genetics.112.139337.DC1). Additionally, the marker coverage was also sampled, the minimal marker coverage was set to at least 40 out of 96 observations at a marker, so that the marker coverage was always sampled between 40 to 96 observations per marker [Turner et al., 2011](http://www.genetics.org/content/suppl/2012/03/30/genetics.112.139337.DC1). <br /> <br />. The drift simulator was run as: <br />
          
```{r}
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
real_data <- read.table("/usr/users/mtost/wd_test_GB_easy_rerun/R_analysis_20162020/GB1307_Allele_and_GT_freq_per_pop_after_MAF005_filt.txt",
                        header = TRUE)
initial_Short1 <- real_data$freq_P_Short1Shoe0
initial_Short2 <- real_data$freq_P_Short2Shoe0
initial_Tall1 <- real_data$freq_P_Tall1Shoe0
initial_Tall2 <- real_data$freq_P_Tall2Shoe0

total_pop <- 5000
males_pop <- 5000
females_pop <- 250
sim <- 5000000
cycles <- 3
kernels_per_plant <- 500

initial_Short1_new <- sample(initial_Short1, sim, replace = TRUE)
initial_Short2_new <- sample(initial_Short2, sim, replace = TRUE)
initial_Tall1_new <- sample(initial_Short1, sim, replace = TRUE)
initial_Tall2_new <- sample(initial_Short1, sim, replace = TRUE)

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

### RUN ----------------------------------------------------------------------
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
```
<img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> values were calculated for all simulated markers, which corresponded in our case to 5,000,000 simulations. We summed those up and choosed the the 99.9999th percentile of the emprirical distribution of all observations as significance threshold, similar to [Kumar et al., 2021](https://academic.oup.com/pcp/article/62/7/1199/6279219).  

### 5.4 Based on the FDRfS
**FDRfS for all possible values of the statistics**
The function `calculate_FDR_for_selection()` generates a table to show all posible values of a statistic and the number of observed markers diverged between subpopulation selected in the same and opposite directions at a certain value. The FDR for selection is received by dividing the number of observed markers diverged between subpopulation selected in the same direction by the number of observed markers diverged between subpopulation selected in opposite directions. <br />
The calculation of the FDR for selection is also demonstrated with the following table which shows for different statistics the number of markers diverged between subpopulations selected in the same and opposite directions and the corresponding FDR for selections:          
|Statistic|Markers diverged same direction|Markers diverged opposite directions|FDR for selection|
|---------|-------------------------------|------------------------------------|-----------------|
|0	|1417796	                      |1421826	                   |0.9972           |
|0.01	|1104313	                      |1115771	                   |0.9897           |
|0.02	|400204	                      |415091	                             |0.9641           |
|0.03	|243600	                      |259984	                             |0.9370           |
|0.04	|182635	                      |198752	                             |0.9189           |
|0.05	|142278	                      |159055	                             |0.8945           |
|0.06	|112868	                      |128888	                             |0.8757           |
|0.07	|90337	                      |105548	                             |0.8559           |
|0.08	|72791	                      |87109	                             |0.8356           |
|...	|...	                      |...	                             |...              |
|0.67	|1	                      |41	                             |0.0244           |
|0.671	|1	                      |41	                             |0.0244           |
|0.672	|1	                      |41	                             |0.0244           |
|0.673	|1	                      |41	                             |0.0244           |
|0.674	|1	                      |41	                             |0.0244           |
|0.675	|1	                      |40	                             |0.0250           |
|0.676	|1	                      |40	                             |0.0250           |
|0.677	|1	                      |39	                             |0.0256           |
|...	|...	                      |...	                             |...              |
|0.752	|1	                      |19	                             |0.0526           |
|0.753	|0		            |19	                             |0                |

**A table like this is automatically generated by the** `calculate_FDR_for_selection()` **function.**
          
```{r}
calculate_FDR_for_selection_based_on_sumFST <- function(stat_opposite_dir1,
                                                         stat_opposite_dir2,
                                                         stat_same_dir1,
                                                         stat_same_dir2,
                                                         statistic, 
                                                         FDR){
  stat_opposite_dir1 <- as.numeric(stat_opposite_dir1)
  stat_opposite_dir2 <- as.numeric(stat_opposite_dir2)
  stat_same_dir1 <- as.numeric(stat_same_dir1)
  stat_same_dir2 <- as.numeric(stat_same_dir2)
    dt <- matrix(data = seq(0,1,0.001), nrow = length(seq(0,1,0.001)), ncol = 3)
    for (i in 1:nrow(dt)) {
      dt[i,2] <- length(which(stat_same_dir1 + stat_same_dir2 >=dt[i,1]))
      dt[i,3] <- length(which(stat_opposite_dir1 + stat_opposite_dir2 >=dt[i,1]))
    }
  dt <- as.data.frame(dt)
  colnames(dt) <- c("Statistic","Markers diverged same dir","Markers diverged opposite dir")
  dt$FDR_for_selection <- dt[,2]/dt[,3]
  if(statistic == "FST"){
    FDR_in_per <- FDR*100
    FDR_below5per <- dt[which(dt$FDR_for_selection < FDR),]
    sig_threshold <- FDR_below5per$Statistic[which.min(FDR_below5per$Statistic)]
    cat("The significance threshold based on the FDR for selection <",FDR_in_per,"%","\n",
        "corresponds to an statistic of",sig_threshold, "\n",
        "This threshold was exceeded by:",max(dt[which(dt$FDR_for_selection < 0.05),2]),
        "markers observed between subpopulations selected in the same direction and",
        max(dt[which(dt$FDR_for_selection < 0.05),3]),
        "markers observed between subpopulations selected in opposite directions")
  }
  return(dt)
}
dt_FDR_for_selection_sum_FST <- calculate_FDR_for_selection_based_on_sumFST(stat_opposite_dir1 = FST_values_od_cor$FST_value_opposite_dir1,
                                                                              stat_opposite_dir2 = FST_values_od_cor$FST_value_opposite_dir2,
                                                                              stat_same_dir1 = FST_values_od_cor$FST_value_same_dir1,
                                                                              stat_same_dir2 = FST_values_od_cor$FST_value_same_dir2,
                                                                              statistic = "FST",
                                                                              FDR = 0.05)
```
Even though, the significance threshold based on the FDR for selection is already printed by the previous function, the following function returns the threshold so it can be used directly in the Manhatten plot or to check the overlap between all the different statistics. The FDRfS can be choosen in the function. In our case we choosed a FDRfS < 5%.
```{r}
get_sig_thresh_FDR_for_selection_based_on_sumFST <- function(stat_opposite_dir1,
                                                        stat_opposite_dir2,
                                                        stat_same_dir1,
                                                        stat_same_dir2,
                                                        statistic, 
                                                        FDR){
  stat_opposite_dir1 <- as.numeric(stat_opposite_dir1)
  stat_opposite_dir2 <- as.numeric(stat_opposite_dir2)
  stat_same_dir1 <- as.numeric(stat_same_dir1)
  stat_same_dir2 <- as.numeric(stat_same_dir2)
  dt <- matrix(data = seq(0,1,0.001), nrow = length(seq(0,1,0.001)), ncol = 3)
  for (i in 1:nrow(dt)) {
    dt[i,2] <- length(which(stat_same_dir1 + stat_same_dir2 >=dt[i,1]))
    dt[i,3] <- length(which(stat_opposite_dir1 + stat_opposite_dir2 >=dt[i,1]))
  }
  dt <- as.data.frame(dt)
  colnames(dt) <- c("Statistic","Markers diverged same dir","Markers diverged opposite dir")
  dt$FDR_for_selection <- dt[,2]/dt[,3]
  if(statistic == "FST"){
    FDR_in_per <- FDR*100
    FDR_below5per <- dt[which(dt$FDR_for_selection < FDR),]
    sig_threshold <- FDR_below5per$Statistic[which.min(FDR_below5per$Statistic)]
    cat("The significance threshold based on the FDR for selection <",FDR_in_per,"%","\n",
        "corresponds to an statistic of",sig_threshold, "\n",
        "This threshold was exceeded by:",max(dt[which(dt$FDR_for_selection < 0.05),2]),
        "markers observed between subpopulations selected in the same direction and",
        max(dt[which(dt$FDR_for_selection < 0.05),3]),
        "markers observed between subpopulations selected in opposite directions")
  }
  return(sig_threshold)
}
sig_threshold_sum_FST <- get_sig_thresh_FDR_for_selection_based_on_sumFST(stat_opposite_dir1 = FST_values_od_cor$FST_value_opposite_dir1,
                                                                            stat_opposite_dir2 = FST_values_od_cor$FST_value_opposite_dir2,
                                                                            stat_same_dir1 = FST_values_od_cor$FST_value_same_dir1,
                                                                            stat_same_dir2 = FST_values_od_cor$FST_value_same_dir2,
                                                                            statistic = "FST",
                                                                            FDR = 0.05)
```
## 6 Sauron plot       
<img src="https://user-images.githubusercontent.com/63467079/170998919-de260c76-c32f-4267-a85a-186fd2e857d3.png", width="700" height="350">
This plot depicts how the FDRfS was computed and provides a visualization of the scope of drift and selection. The Sauron plot comes from [Turner and Miller (2012)](http://www.genetics.org/content/suppl/2012/03/30/genetics.112.139337.DC1), but it can be also created for the <img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> statisitc. Sauron plot of genetic differentiation for FSTSum observed between the subpopulations selected in the same direction (blue) and in opposite directions (red). Each dot represents one SNP. The transparent red colored edges correspond to a false discovery rate (FDR) for selection < 5%. The y- and x-axis correspond to the range of <img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> values. 
The "Sauron plot" for the <img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> does not even look like the eye of [Sauron](https://twitter.com/strnr/status/457201981007073280) anymore, as the "Sauron plot" from [Turner and Miller (2012)](http://www.genetics.org/content/suppl/2012/03/30/genetics.112.139337.DC1) did. The Sauron plot for the allele frequency difference is created by the `create_sauron_plot_FST()` function, which is shown above:
```{r}
create_sauron_plot_FST <- function(data_table,
                                   sig_threshold,
                                   stat_opposite_dir1,
                                   stat_opposite_dir2,
                                   stat_same_dir1,
                                   stat_same_dir2){
  windowsFonts(my = windowsFont('Calibri'))
  create_FDR_area <- function(x){
    sig_threshold+-1*x
  }
  y_coord_value <- create_FDR_area(x = seq(0,1,0.01))
  dt_coords_area <- cbind(y_coord_value,seq(0,1,0.01))
  dt_coords_area <- as.data.frame(dt_coords_area)
  dt_coords_area$z <- rep(sig_threshold,nrow(dt_coords_area))
  colnames(dt_coords_area) <- c("x_coord_value","y_coord_value","threshold_value")
  dt_coords_threshold_area <- dt_coords_area[which(dt_coords_area$y_coord_value < sig_threshold),]
  dt_coords_area_rest1 <- cbind(seq(sig_threshold,1,0.01),seq(sig_threshold,1,0.01))
  dt_coords_area_rest1 <- as.data.frame(dt_coords_area_rest1)
  colnames(dt_coords_area_rest1) <- c("x_coord_value","y_coord_value")
  sauron_1 <- ggplot()+
    geom_area(data = dt_coords_area, aes(x=x_coord_value, y=threshold_value), alpha = 0.4, fill = "tomato")+
    geom_area(data = dt_coords_area, aes(x=x_coord_value, y=y_coord_value), fill = "white")+
    geom_area(data = dt_coords_area_rest1, aes(x=x_coord_value, y=y_coord_value),  orientation = "x", alpha = 0.4, fill = "tomato") +
    geom_area(data = dt_coords_area_rest1, aes(x=x_coord_value, y=y_coord_value), orientation = "y", alpha = 0.4, fill = "tomato") +
    geom_point(data = data_table, aes(x = stat_same_dir1,y = stat_same_dir2),shape=1, colour = "royalblue3")+
    geom_point(data = data_table, aes(x = stat_opposite_dir1,y = stat_opposite_dir2),shape=1, colour = "firebrick")+
    theme(text = element_text(size =font_size, family = "my"),
          panel.background = element_rect(fill = "white"),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_text(size = font_size, family = "my", colour = "royalblue3", face = "bold"),
          axis.title.x = element_text(size = font_size, family = "my", colour = "royalblue3", face = "bold"),
          axis.line.x = element_line(color = "black"),
          axis.text = element_text(size = font_size, family = "my", colour = "black"))+
    labs(x=paste("\t","\t","\t","\t","\t","\t","\t","\t","Short 1 vs Short 2"),
         y=paste("\t","\t","\t","Tall 1 vs Tall 2"))+
    ylim(0,0.75)+
    xlim(0,0.75)
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
    geom_area(data = dt_coords_area, aes(x=x_coord_value, y=threshold_value), alpha = 0.4, fill = "tomato")+
    geom_area(data = dt_coords_area, aes(x=x_coord_value, y=y_coord_value), fill = "white")+
    geom_area(data = dt_coords_area_rest1, aes(x=x_coord_value, y=y_coord_value),  orientation = "x", alpha = 0.4, fill = "tomato") +
    geom_area(data = dt_coords_area_rest1, aes(x=x_coord_value, y=y_coord_value), orientation = "y", alpha = 0.4, fill = "tomato") +
    geom_point(data = data_table, aes(x = stat_opposite_dir1,y = stat_opposite_dir2), shape = 1,  colour = "firebrick")+
    geom_point(data = data_table, aes(x = stat_same_dir1,y = stat_same_dir2),shape = 1, colour = "royalblue3")+
    theme(text = element_text(size = font_size, family = "my"),
          panel.background = element_rect(fill = "white"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.position = "bottom",
          axis.title.y = element_text(size = font_size, family = "my", colour = "firebrick", face = "bold"),
          axis.title.x = element_text(size = font_size, family = "my", colour = "firebrick", face = "bold"),
          axis.line.x = element_line(color = "black"),
          axis.text = element_text(size = font_size, family = "my", colour = "black"))+
    labs(x=  "Short 1 vs Short 2",
         y= "Tall 2 vs Tall 2")+
    geom_abline(intercept = sig_threshold, slope = -1, size = 1, alpha = 0.2, colour = "firebrick4")+
    xlim(0,0.75)+
    ylim(0,0.75)
  sauron_2
  sauron_plot_all_FST<- grid.arrange(yaxis_sauron,sauron_2, xaxis_sauron,empty_plot,
                                     ncol=2, nrow = 2, 
                                     layout_matrix = rbind(c(1,2),c(4,3)),
                                     widths = c(0.2, 4), heights = c(4,0.2))
  return(sauron_plot_all_FST)
}
Sauron_plot_FST <- create_sauron_plot_FST(data_table = FST_values_od_cor,
                                          sig_threshold = sig_threshold_sum_FST ,
                                          stat_opposite_dir1 = FST_values_od_cor$FST_value_opposite_dir1,
                                          stat_opposite_dir2 = FST_values_od_cor$FST_value_opposite_dir2,
                                          stat_same_dir1 = FST_values_od_cor$FST_value_same_dir1,
                                          stat_same_dir2 = FST_values_od_cor$FST_value_same_dir2)
Sauron_plot_FST
```

The two Sauron plots can be combined:
```{r}
combined_sauron_plots <- grid.arrange(Sauron_plot_AFD, Sauron_plot_FST, 
                                    ncol=2, nrow = 1,
                                    labels = c("A","B"))
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
create_manhatten_plot_per_chr_both_stat <- function(data_1,
                                                    data_2,
                                                    statistic_1, 
                                                    statistic_2,
                                                    name_of_statistic_1,
                                                    name_of_statistic_2,
                                                    SNP_position,
                                                    sign_thres_emp_dis_999,
                                                    sign_thres_emp_dis_9999,
                                                    sign_thres_drift_sim,
                                                    sign_thres_FDR){
  windowsFonts(my = windowsFont('Calibri'))
  Stat1_chr_plot <- ggplot(data = data_1)+
    geom_point(aes(x = SNP_position_1,y = statistic_1), colour = "black")+
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
    labs(y = name_of_statistic_1,x = "Position", tag = "A")+
    geom_hline(yintercept = sign_thres_emp_dis_999, colour="darkred", size = 1)+
    geom_hline(yintercept = sign_thres_emp_dis_9999, colour="tomato", size = 1)+
    geom_hline(yintercept = drift_simulation_based_sig_thresold, colour="purple", size = 1)+
    geom_hline(yintercept = FDR_of_selection_significance_threshold, colour="gold", size = 1)
  Stat2_chr_plot <- ggplot(data = data_2)+
    geom_point(aes(x = SNP_position_2,y = statistic_2), colour = "black")+
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
    labs(y = name_of_statistic_2,x = "Position", tag = "B")+
    geom_hline(aes(yintercept = sign_thres_emp_dis_999, linetype=str_wrap("99.9th percentile",10)), colour="darkred", size = 1)+
    geom_hline(aes(yintercept = sign_thres_emp_dis_9999, linetype=str_wrap("99.99th",10)), colour="tomato", size = 1)+
    geom_hline(aes(yintercept = drift_simulation_based_sig_thresold, linetype=str_wrap("drift simulations",10)), colour="purple", size = 1)+
    geom_hline(aes(yintercept = FDR_of_selection_significance_threshold, linetype=str_wrap("FDR",10)), colour="gold", size = 1)+
    scale_linetype_manual(name = str_wrap("Significance threshold",10), values = c(1,1,1,1), 
                          guide = guide_legend(override.aes = list(color = c("darkred", "tomato",
                                                                             "purple","gold"))))
  combined_Manhattenplots <- grid.arrange(Stat1_chr_plot,Stat2_chr_plot,
                                      ncol=2, nrow = 1,
                                      widths = c(8, 8.5), heights = c(4))
  return(combined_Manhattenplots)
}
manhatten_plot <- create_manhatten_plot_per_chr_both_stat(data_1 = Fst_Chr_3,
                                                          data_2 = AFD_Chr_3,
                                                          statistic_1 = Fst_Chr_3$Fst_value, 
                                                          statistic_2 = AFD_Chr_3$Abs_allele_freq_diff,
                                                          name_of_statistic_1 = "Fst",
                                                          name_of_statistic_2 = "|allele frequency difference|",
                                                          SNP_position_1 = Fst_Chr_3$Fst_value,
                                                          SNP_position_2 = AFD_Chr_3$Abs_allele_freq_diff,
                                                          Stat1_sign_thres_emp_dis_999 = FST_emp_dis_999_perc,
                                                          Stat1_sign_thres_emp_dis_9999 = FST_emp_dis_9999_perc,
                                                          Stat1_sign_thres_drift_sim = FST_drift_sim_sig_thres,
                                                          Stat1_sign_thres_FDR = FDR_for_sel_sig_thres,
                                                          Stat2_sign_thres_emp_dis_999 = AFD_emp_dis_999_perc,
                                                          Stat2_sign_thres_emp_dis_9999 = AFD_emp_dis_9999_perc,
                                                          Stat2_sign_thres_drift_sim = AFD_drift_sim_sig_thres,
                                                          Stat2_sign_thres_FDR = abs_AFD_FDR_for_sel_sig_thres)          
```
In case one wants to plot only one statistic in the Manhatten plot, this function creates the only one Manhatten plot:
```{r}
create_manhatten_plot_per_chr_one_stat <- function(data,
                                                   name_of_statistic,
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
    labs(y = name_of_statistic,x = "Position")+
    geom_hline(aes(yintercept = sign_thres_emp_dis_999, linetype=str_wrap("99.9th percentile",10)), colour="darkred", size = 1)+
    geom_hline(aes(yintercept = sign_thres_emp_dis_9999, linetype=str_wrap("99.99th",10)), colour="tomato", size = 1)+
    geom_hline(aes(yintercept = sign_thres_drift_sim, linetype=str_wrap("drift simulations",10)), colour="purple", size = 1)+
    geom_hline(aes(yintercept = sign_thres_FDR, linetype=str_wrap("FDR",10)), colour="gold", size = 1)+
    scale_linetype_manual(name = str_wrap("Significance threshold",10), values = c(1,1,1,1), 
                          guide = guide_legend(override.aes = list(color = c("darkred", "tomato",
                                                                             "purple","gold"))))
  return(chr_plot)
}
FST_manhatten_plot <- create_manhatten_plot_per_chr_one_stat(data = Fst_Chr_3,
                                                             name_of_statistic = "Fst",
                                                             statistic = Fst_Chr_3$Fst_value,
                                                             SNP_position = Fst_Chr_3$Position,
                                                             sign_thres_emp_dis_999 = FST_emp_dis_999_perc,
                                                             sign_thres_emp_dis_9999 = FST_emp_dis_9999_perc,
                                                             sign_thres_drift_sim = FST_drift_sim_sig_thres,
                                                             sign_thres_FDR = FST_FDR_for_sel_sig_thres)   
AFD_manhatten_plot <- create_manhatten_plot_per_chr_one_stat(data = AFD_Chr_3,
                                                             name_of_statistic = "|allele frequency difference|",
                                                             statistic = AFD_Chr_3$Abs_allele_freq_diff,
                                                             SNP_position = AFD_Chr_3$Position,
                                                             sign_thres_emp_dis_999 = AFD_emp_dis_999_perc,
                                                             sign_thres_emp_dis_9999 = AFD_emp_dis_9999_perc,
                                                             sign_thres_drift_sim = AFD_drift_sim_sig_thres,
                                                             sign_thres_FDR = abs_AFD_FDR_for_sel_sig_thres) 
```
         
## 8 LD decay
The extent of linkage disequilibrium (LD) was estimated based on all 4,029,092 SNP markers with the PLINK software  v1.90 as squared correlation between markers as R2 [(Purcell et al., 2007)](https://www.cell.com/ajhg/fulltext/S0002-9297(07)61352-4). The script for the computation of LD decay between markers is available as `LD_decay_calculation_with_plink.bash`.
![2021_V3_LD_decay_ggplot](https://user-images.githubusercontent.com/63467079/150421503-dcebed82-9f2c-45f0-9457-f8db16320ec7.png)
