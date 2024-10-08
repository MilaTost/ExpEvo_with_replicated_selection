# Experimental evolution with replicated selection
## Table of contents
[Introduction](https://github.com/MilaTost/ExpEvo_with_replicated_selection/blob/main/README.md#0-introduction) <br />
[Phenotypic data analysis](https://github.com/MilaTost/ExpEvo_with_replicated_selection/blob/main/README.md#1-phenotypic-data-analysis) <br />
&emsp;[Truncation selection thresholds](https://github.com/MilaTost/ExpEvo_with_replicated_selection/blob/main/README.md#11-truncation-selection-thresholds) <br />
&emsp;[Plant height measurements](https://github.com/MilaTost/ExpEvo_with_replicated_selection/blob/main/README.md#12-plant-height-measurements) <br />
&emsp;[Realized heritability from the Breeder's equation](https://github.com/MilaTost/ExpEvo_with_replicated_selection/blob/main/README.md#13-realized-heritability-from-the-breeders-equation) <br />
2 Pipeline for the analysis of GBS data adapted from [Wickland et al. 2013](https://github.com/dpwickland/GB-eaSy) <br />
[Filtering](https://github.com/MilaTost/ExpEvo_with_replicated_selection/blob/main/README.md#3-filtering) <br />
[Allele frequency estimation](https://github.com/MilaTost/ExpEvo_with_replicated_selection/blob/main/README.md#4-allele-frequency-estimation) <br />
[Estimation of the effective population size](https://github.com/MilaTost/ExpEvo_with_replicated_selection/blob/main/README.md#5-estimation-of-the-effective-population-size) <br />
[LD decay](https://github.com/MilaTost/ExpEvo_with_replicated_selection/blob/main/README.md#6-ld-decay) <br />
[Scan for selection signatures mapping with FST leveraging replicated selection](https://github.com/MilaTost/ExpEvo_with_replicated_selection/blob/main/README.md#7-scan-for-selection-signatures) <br />
[Test for selection: "Sauron plot"](https://github.com/MilaTost/ExpEvo_with_replicated_selection/blob/main/README.md#8-test-for-selection-sauron-plot) <br />
[Window based analysis](https://github.com/MilaTost/ExpEvo_with_replicated_selection/tree/main#9-window-based-analysis) <br />
[Comparison of significance thresholds and candidate gene identification](https://github.com/MilaTost/ExpEvo_with_replicated_selection/blob/main/README.md#10-comparison-of-significance-thresholds-and-candidate-gene-identification) <br />
[Haplotype estimation and haplotype block calculation](https://github.com/MilaTost/ExpEvo_with_replicated_selection/tree/main#11-haplotype-estimation-and-haplotype-block-calculation) <br />
[Generation of the LDheatmap](https://github.com/MilaTost/ExpEvo_with_replicated_selection/tree/main#11-haplotype-estimation-and-haplotype-block-calculation) <br /> <br />
## Introduction
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
statistic. <br />
**Experimental design of an experimental evolution study with replicated selection** <br />
<img width="500" alt="Figure_1" src="https://github.com/MilaTost/ExpEvo_with_replicated_selection/assets/63467079/558bfbe8-5844-451a-ab2c-229f72488a92">

## Phenotypic data analysis
### Truncation selection thresholds
The truncation thresholds of the 5% tallest or 5% shortest plants within the populations can be calculated like this:
```{r}
quantile(Short_Population_1$PlantHeight, 0.05, na.rm = TRUE)
quantile(Short_Population_2$PlantHeight, 0.05, na.rm = TRUE)
quantile(Tall_Population_1$PlantHeight, 0.05, na.rm = TRUE)
quantile(Tall_Population_2$PlantHeight, 0.05, na.rm = TRUE)
```
### Plant height measurements
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
**Measured plant height in all years in the subpopulations selected for short plant (green) and tall plant height (purple).** <br />
<img width="600" alt="Figure_2" src="https://github.com/MilaTost/ExpEvo_with_replicated_selection/assets/63467079/596ca155-6c47-4b77-9e65-607d1825f9c2"> <br />
The computation of the t-test statistic and the script for plotting the measured phenotypes across all years are available in the `Phenotypic_data_analysis.R` script.
### Realized heritability from the Breeder's equation
The realized heritability was calcaluated according to (Lush, 1937). This is done by the `Realized_heritability_calculations.R` script. To calculate the realized heritability the selection differential and the response to selection are calculated between the different generations of selection. In the following example the realized heritability was calculated between generation 0 and generation 1 of the subpopulations short plants 1. <br />

```{r}
calculated_realized_herit <- function(Data_Gen1,
                                      Data_Gen2){
  mean_Gen_1 <- mean(Data_Gen1$PlantHeight, na.rm = TRUE)
  mean_Gen_2 <- mean(Data_Gen2$PlantHeight, na.rm = TRUE)
  Selected_prop <- quantile(Data_Gen1$PlantHeight, 0.05, na.rm = TRUE)
  Selection_Differential <- Selected_prop - mean_Gen_1
  Response_to_Selection <- mean_Gen_2 - mean_Gen_1
  herit_BS <- as.numeric(Response_to_Selection)/as.numeric(Selection_Differential)
  return(herit_BS)
}
calculated_realized_herit(Data_Gen1 = Short_1_2016,
                          Data_Gen2 = Short_1_2017)
```

## Pipeline for the analysis of GBS data adapted from [Wickland et al. 2013](https://github.com/dpwickland/GB-eaSy)
For the analysis of our raw reads from paired-end genotyping-by-sequencing (GBS) with ApeKI according to [Elshire et al. 2011](https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0019379&type=printable), we used the GB-eaSy pipeline from [Wickland et al. 2013](https://github.com/dpwickland/GB-eaSy). <br />
The pipeline consists out of several steps, which comprises:
- [Demultiplexing and trimming of the adapter sequence](https://github.com/dpwickland/GB-eaSy#step-2-demultiplex-raw-reads)
- [Alignement to the reference genome](https://github.com/dpwickland/GB-eaSy#step-3-align-to-reference)
- [Create a list of sorted bam files](https://github.com/dpwickland/GB-eaSy#step-4-create-list-of-sorted-bam-files)
- [Generate pileup and SNP calling](https://github.com/dpwickland/GB-eaSy#step-5-and-6-generate-pileup-and-call-snps)
- Filtering for quality parameters <br /> <br />

This bash-script was run on every sequenced plate separately, since every sequenced plate had it's own adapters and the same set of barcodes was used for all sequenced plates in our case. The *Demultiplexing* and *Alignement to the reference genome* steps were conducted by the `GB_eaSy_demultiplexing_and_alignemnt_S1.bash`,`GB_eaSy_demultiplexing_and_alignemnt_S2.bash`,`GB_eaSy_demultiplexing_and_alignemnt_S3.bash` and `GB_eaSy_demultiplexing_and_alignemnt_S4.bash` scripts. <br /> <br />
The *Generate Pileup*,*SNP calling* and *Filtering for quality parameters* were run for all sequencing runs together. The *Filtering for quality parameters* was changed a lot and differed from the GB-eaSy pipline introduced by [Wickland et al. 2013](https://github.com/dpwickland/GB-eaSy). The  *Generate Pileup*,*SNP calling* and *Filtering for quality parameters* steps were conducted by the `GB_eaSy_pileup_SNP_calling_filtering.bash` script.

## Filtering
After the VCF file was created with the GB-eaSy pipeline from [Wickland et al. 2013](https://github.com/dpwickland/GB-eaSy) and filtered for some quality parameters we did some additional filtering in R. All  listed *Filterung* steps were conducted by the `Filtering_for_coverage_average_RD_missingness.R` script. The Rscript contains the *Filtering* functions for the following filtering steps.  
### The filtering steps:
- Filtering for individual samples with low coverage
- Filtering for read depth per sample
- Filtering for missingness
- Removal of non-diallelic and non-polymorphic markers
<br /> <br />
The Rscript contains many functions from the `vcfR` package from [Knaus and Grünwald 2018](https://github.com/knausb/vcfR#:~:text=VcfR%20is%20an%20R%20package%20intended%20to%20allow,rapidly%20read%20from%20and%20write%20to%20VCF%20files.).
<br />

## Allele frequency estimation
The allele frequency estimation was done with the `Allele_frequency_estimation.R` rscript.<br/> <br />

## Estimation of the effective population size  
The effective population size was calculated using the known demographic parameters of the populations and based off of the formula:
          <img width="120" alt="effective_pop_size" src="https://user-images.githubusercontent.com/63467079/171614297-2936b14e-3c6c-4914-b0b5-6e6cca2c5fc7.png"> <br /> where ![formula](https://render.githubusercontent.com/render/math?math=N_{m}) is the number of male and ![formula](https://render.githubusercontent.com/render/math?math=N_{f}) the number of female individuals (Crow and Kimura, 1970). In every generation the ~5% shortest or tallest plants were harvested, amounting to ~250 plants. Under the assumption of complete random mating, all 5000 plants could contribute as male parents for the next generation, the effective population size was calculated with the `Estimation_of_effective_population_size.R` rscript.<br /> <br />

In field trials, the assumption of random mating is violated because of the effects of assortative mating due to varying flowering dates and limited spatial pollen dispersion (Allard, 1999). Therefore, we evaluated the flowering dates of the 96 randomly chosen plants from each subpopulation to approximate the number of simultaneously flowering tassels and silks. It was assumed that silks remain receptive up to 5 days after silk emergence (Nieh et al., 2014), so that we calculated the number of flowering tassels during this time interval and projected it onto the entire subpopulation. This calculation of this is also available in the  `Estimation_of_effective_population_size.R` rscript. <br /> <br />

## LD decay
The extent of linkage disequilibrium (LD) was estimated based on 100,000 SNP markers with TASSEL v5 using the squared correlation between markers, ![formula](https://render.githubusercontent.com/render/math?math=R^{2}), over a window size of 2,000 bp (Bradbury et al., 2007). <br />
#### Filtering and thinning out for LD decay calculation
We filtered the individuals out with a marker coverage below `0.6`. For this analysis, we required that every SNP marker was observed at least 80 times in each subpopulation. This filtering resulted in 1,243,604 SNP markers. Additionally, we thinned the data set by randomly sampling 10,000 markers per chromosome. <br />  
#### LD decay calculation
The calculation can be found in the script `Calc_LD_decay_calc_with_TASSEL_based_on_thinned_out_VCF.bash`. <br />
#### LD decay plotting
LD decay was modeled with a nonlinear regression model as expected value <img width="120" alt="Formula_Github" src="https://github.com/MilaTost/ExpEvo_with_replicated_selection/assets/63467079/dedff33d-bf1f-40e4-9191-9dbfc33b846b"> with N as the number of individuals at each site and C as the recombination coefficient between sites (Remington et al., 2001). The calculation and plotting script is available as `Plot_LD_decay_based_on_TASSEL_output.R`. We also created a pairwise LD heatmap with the R package LDheatmap (Shin et al., 2006). This calculation is contained in the `Create_Haplotype_figure.R` script. <br />

## 7 Scan for selection signatures
Our scan for selection was based on the <img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> leveraging replicated selection. <img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> is calculated as: <br /> <br />
<img width="150" alt="FST_formula" src="https://user-images.githubusercontent.com/63467079/171010281-cf3648d2-1baf-4fc3-8e23-83d3dbd792ca.png">
<br /> <br /> according to [Weir and Cockerham, 1984](https://doi.org/10.1111/j.1558-5646.1984.tb05657.x). <br />

#### Calculation of FSTSum
The sum of FST was calculated between the two non-redundant comparisons of divergently selected subpopulations (Short 1 vs Tall 1; Short 2 vs Tall 2) and subpopulation selected in the same direction (Short 1 vs Short; Tall 1 vs Tall 2). These values are required for the calculation of the false discovery rate for selection (FDRfS). For plotting we only use the FSTSum value calculated between divergently selected subpopulations (Short 1 vs Tall 1; Short 2 vs Tall 2). The results are presented [here](https://github.com/MilaTost/ExpEvo_with_replicated_selection/edit/main/README.md#manhattan-plot). <br />       

## Window-based analysis
We implemented a window-based analysis to assess linked selection. A cubic smoothing spline was applied to the single-marker based FSTSum values. From the fitted spline, inflection points were calculated and used as window boundaries for newly defined regions. This analysis was performed with the R package GenWin from [Beissinger et al., 2015](https://gsejournal.biomedcentral.com/articles/10.1186/s12711-015-0105-9). Briefly, the GenWin package returns window boundaries that we used to define regions with linked markers that may be analyzed together ([Beissinger et al., 2015). The calculation of window boundaries with the `splineAnalyze()` function from the [GenWin package](https://cran.r-project.org/web/packages/GenWin/index.html) are available in the `Window_based_analysis.R`. The results are presented [here](https://github.com/MilaTost/ExpEvo_with_replicated_selection/edit/main/README.md#manhattan-plot). <br />

## Calculation of significance thresholds
Significance thresholds for selection were calculated three ways: 1) based on the empirical distribution; 2) based on drift simulations; and 3) based on the false discovery rate for selection (FDRfS) [Turner and Miller (2012)](http://www.genetics.org/content/suppl/2012/03/30/genetics.112.139337.DC1). <br />
The calculation of significance thresholds and the plotting functions are contained in the `Significance_thresholds_and_plotting.R`. <br /> <br />
### Based on the empirical distribution
The significance thresholds based on the empirical distribution, were calculated by taking the 99.9th and 99.99th percentile of the empirical distribution of the <img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> :
```{r}
FST_sig_thres_1 <- quantile(FST_values_od_cor$Fst, probs = 0.9999, na.rm = TRUE)
FST_sig_thres_2 <- quantile(FST_values_od_cor$Fst, probs = 0.999, na.rm = TRUE)    
```          
The significance threshold is stored, so it can be used later directly for plotting. <br /> <br />

### Based on the FDRfS
**FDRfS for all possible values of the statistics**
The function `calculate_FDR_for_selection()` generates a table to show all posible values of a statistic and the number of observed markers diverged between subpopulation selected in the same and opposite directions at a certain value. The FDR for selection is received by dividing the number of observed markers diverged between subpopulation selected in the same direction by the number of observed markers diverged between subpopulation selected in opposite directions. <br />

The different thresholds are compared here: <br />
<img src="C:\Users\mtost\Documents\2024_Shoepeg_resubmission_V3\Plots\final plots\S42.png" width="600">
**Observations between the subpopulations selected in opposite directions (A) and the same direction (B) with the significance thresholds based on the 99.99th percentile of the empirical distribution (light purple), based on the 99.9th percentile of the empirical distribution (dark purple), and by the false discovery rate for selection (FDRfS) (yellow).** <br /> <br />
The plotting functions are available in the `Plot_Manhattan_plot_with_WStat.R` script. <br /> <br />

## Haplotype estimation and haplotype block calculation
The data were phased using fastPHASE version 1.4.8 [(Sheet and Stephens, 2006)](https://stephenslab.uchicago.edu/software.html#fastphase) for further dissection. We phased the genomic data of the region from 9.437 to 10.457 Mb on chromosome 3 with 10 iterations of the expectation-maximization (EM) algorithm (Sheet and Stephens, 2006). The model was supplied with labels indicating the different subpopulations. The scripts are available as `Run_fastPHASE_random9MB_region_on_chr3.bash`. Before we ran fastPHASE, we prepared the data with the Rscript `Prepare_data_for_fastPHASE.R`. Haplotype blocks were calculated with the R package HaploBlocker [Pook et al., 2019](https://github.com/tpook92/HaploBlocker). The results look like this: <br />
<img src="https://github.com/MilaTost/ExpEvo_with_replicated_selection/assets/63467079/76e648e6-fd51-41ab-bc79-6a9a201b0c80.png" width="600"> <br />
The haplotype block calculation and plotting of the blocks are contained in the `Create_Haplotype_figure.R` script. <br /> <br />

## Generation of the LDheatmap
We also created a pairwise LD heatmap with the R package LDheatmap [Shin et al., 2006] (https://sfustatgen.github.io/LDheatmap/index.html) to supplement our haplotype investigation. We looked at LD across putatively selected regions in the different subpopulations. The pairwise LD heatmap shows strong LD in the subpopulations selected for short plant height across the entire region from 9.437 to 10.457 Mb on chromosome 3. <br />
<img src="https://github.com/MilaTost/ExpEvo_with_replicated_selection/assets/63467079/d61e86d6-e710-48cf-b892-92f940ff505e.png" width="600">
