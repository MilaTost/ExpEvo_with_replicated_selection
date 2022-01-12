# Selection_signature_mapping_scripts_with_replicated_selection
## Table of contents
#### 0 Introduction
#### 1 Pipeline for the analysis of GBS data adapted from [Wickland et al. 2013](https://github.com/dpwickland/GB-eaSy)
#### 2 Filtering
#### 3 <img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> with replicated selection
#### 4 Allele frequency differences
#### 5 Significance threshold based on the FDR for selection
#### 6 Significance thresholds based on drift simulations 
#### 7 Significance thresholds based on the empiric distribution
#### 8 Simulation of Drift
#### 9 Sauron plot
#### 10 Other plotting scripts


## Introduction
This repository contains scripts for selection signature mapping with replicated selection.
Usually in experimental evolution studies a random-mating population is selected in divergent directions for one
trait of interest. The highly divergent subpopulations are then scanned for signatures of selection.
A previous problem in these studies was the determination of significance thresholds. 
One solution to this problem is the usage of **replicated selection**. 
This was introduced by Turner and Miller (2012), but only tested on *Drospohila melanogaster*. <br /> <br />
In this approach two subpopulations are selected in the same direction. 
By that we differentiate between changes caused by selection and changes caused by other factors like drift.
The **replicated selection** is used to calculate a *FDR for selection*. 
Turner and Miller (2012) applied this significance threshold only on the
statistic based on the allele frequency differences. 
Whereas we applied the this new significance threshold to the commonly used 
          <img src="https://render.githubusercontent.com/render/math?math=F_{ST}">
statistic.

## Pipeline for the analysis of GBS data adapted from [Wickland et al. 2013](https://github.com/dpwickland/GB-eaSy)


## Filtering
After the VCF file was created with the GB-eaSy pipeline from [Wickland et al. 2013](https://github.com/dpwickland/GB-eaSy) and filtered for some quality parameters we did some additional filtering in R. 
#### The filtering steps included:
- Filtering for individual samples with low coverage
- Filtering for read depth per sample
- Filtering for missingness
- Removal of non-diallelic and non-polymorphic markers 

#### Loading of required packages and the VCF file:
Most of the function come from the vcfR package from [Knaus and Grünwald 2018](https://github.com/knausb/vcfR#:~:text=VcfR%20is%20an%20R%20package%20intended%20to%20allow,rapidly%20read%20from%20and%20write%20to%20VCF%20files.)
Most calculations and filtering steps over the entire set of markers were run by using the `mclapply()` function from the `parallel()` package.
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
The last command also tells you, how many markers the VCF file comprises.

#### Generate a plot to check the quality parameters and if the previous filtering worked:

```{r}
chrom <- create.chromR(name="Supercontig", vcf=vcf_S4, verbose=FALSE)
png(file='101_Quality_plot_distribution_before_filtering_vcf.png',width = 480, height = 480, units = "px", pointsize=12)
par(
  family='sans',
  cex.axis=1,
  bg='white')
chromoqc(chrom, dp.alpha = 66)
dev.off()
```

#### Filtering for individual samples with low coverage
In these steps the missing observation per individual sample are calculated.
```{r}
gt <- extract.gt(vcf_S4, element = "GT", as.numeric=TRUE)
coverage_ind <- function(i){
  length(which(!is.na(gt[,i])))
}
coverage_ind_tot <- mapply(coverage_ind,1:ncol(gt))
coverage_ind_rel <- coverage_ind_tot/NROW(gt)
NAs_per_ind <- cbind(coverage_ind_tot,coverage_ind_rel)
rownames(NAs_per_ind) <- colnames(gt)
NAs_per_ind <- as.data.frame(NAs_per_ind)
```
The actually filtering is done by these commands:
```{r}
cat("Before filtering the vcf_S4 file contained:",ncol(vcf_S4@gt)-1,"individuals.","\n")
index_ind <- rownames(NAs_per_ind[which(NAs_per_ind$coverage_ind_rel < 0.1),])
look_up_index <- function(i){
  which(colnames(gt)==index_ind[i])
}
index_vcf_S4 <- unlist(mclapply(1:length(index_ind),look_up_index))
vcf_S4@gt <- vcf_S4@gt[,-index_vcf_S4]
cat("The vcf_S4 file still contains:",ncol(vcf_S4@gt)-1,"individuals.","\n")
```

#### Filtering for read depth per sample
In these steps the average read depth per sample is calculated:
```{r}
cat("VCF file contains:",nrow(vcf_S4@fix),"markers before filtering for read depth.","\n")
dp <- extract.info(vcf_S4, element = "DP", as.numeric=TRUE)
an <- extract.info(vcf_S4, element = "AN", as.numeric=TRUE)
obs <- (an/2)
average_depth <- dp/obs
average_depth <- as.data.frame(average_depth)
rownames(average_depth) <- paste(vcf_S4@fix[,1],vcf_S4@fix[,2], sep = "_")
average_depth$SNP_ID <- paste(vcf_S4@fix[,1],vcf_S4@fix[,2], sep = "_")
average_depth$tot_DP <- dp
average_depth$observations <- obs
average_depth$observations <- as.numeric(average_depth$observations)
average_depth$tot_DP <- as.numeric(average_depth$tot_DP)
average_depth$average_depth <- as.numeric(average_depth$average_depth)

write.table(average_depth,"/usr/users/mtost/GB_easy_results_analysis_wd/GB1002_average_read_depth.txt", sep = "  ", row.names = TRUE,
            quote = FALSE)
```
The actually filtering is done by these commands:
```{r}
index_depth <- which(average_depth$average_depth > 1 & average_depth$average_depth < 10)
vcf_S4@gt <- vcf_S4@gt[index_depth,]
vcf_S4@fix <- vcf_S4@fix[index_depth,]
```


#### Filtering for missingness
```{r}
dp <- extract.gt(vcf_S4, element = "DP", as.numeric=TRUE)
gt <- extract.gt(vcf_S4, element = "GT", as.numeric=TRUE)
gt_pop <- gt
colnames(gt_pop) <- str_sub(colnames(gt),1,9)
Shoepag_1 <- gt_pop[,which(colnames(gt_pop)=="Shoepag_1")]
Shoepag_2 <- gt_pop[,which(colnames(gt_pop)=="Shoepag_2")]
Shoepag_3 <- gt_pop[,which(colnames(gt_pop)=="Shoepag_3")]
Shoepag_4 <- gt_pop[,which(colnames(gt_pop)=="Shoepag_4")]

get_missing_obs_pop1 <- function(i){
  length(which(is.na(Shoepag_1[i,])))
}
mis_obs_pop1 <- unlist(mclapply(1:nrow(Shoepag_1), get_missing_obs_pop1))

get_missing_obs_pop2 <- function(i){
  length(which(is.na(Shoepag_2[i,])))
}
mis_obs_pop2 <- unlist(mclapply(1:nrow(Shoepag_2), get_missing_obs_pop2))

get_missing_obs_pop3 <- function(i){
  length(which(is.na(Shoepag_3[i,])))
}
mis_obs_pop3 <- unlist(mclapply(1:nrow(Shoepag_3), get_missing_obs_pop3))

get_missing_obs_pop4 <- function(i){
  length(which(is.na(Shoepag_4[i,])))
}
mis_obs_pop4 <- unlist(mclapply(1:nrow(Shoepag_4), get_missing_obs_pop4))

NAs_per_marker_per_pop <- cbind(mis_obs_pop1,mis_obs_pop2,mis_obs_pop3,mis_obs_pop4)
NAs_per_marker_per_pop <- as.data.frame(NAs_per_marker_per_pop)
NAs_per_marker_per_pop <- cbind(rownames(gt_pop),NAs_per_marker_per_pop)
colnames(NAs_per_marker_per_pop) <- c("SNP_ID","mis_obs_pop1","mis_obs_pop2","mis_obs_pop3","mis_obs_pop4")
NAs_per_marker_per_pop <- as.data.frame(NAs_per_marker_per_pop)

rm(gt_pop,Shoepag_1,Shoepag_2,Shoepag_3,Shoepag_4)

### Calculate the range of missing markers -------------------------------------
################################################################################
get_range_of_missingness <- function(i){
  abs(range(NAs_per_marker_per_pop[i,2:5])[1]-range(NAs_per_marker_per_pop[i,2:5])[2])
}

markers_diff <- unlist(mclapply(1:nrow(NAs_per_marker_per_pop),get_range_of_missingness))
NAs_per_marker_per_pop$range_of_missingness <- markers_diff

cat("Calculation of missing observations per markers is done.","\n")

```

#### Removal of non-diallelic and non-polymorphic markers
```{r}

```

## <img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> with replicated selection


## Allele frequency differences


## Significance threshold based on the FDR for selection


## Significance thresholds based on drift simulations 


## Significance thresholds based on the empiric distribution


## Simulation of Drift


## Sauron plot

![GB1006_Sauron_plots_combined_values_2021_12_17](https://user-images.githubusercontent.com/63467079/149146525-ce94e222-dff8-4ad4-8dcb-ad14f7530032.png)

## Other plotting scripts

