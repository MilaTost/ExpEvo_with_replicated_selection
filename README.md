# Experimental evolution with replicated selection
## Table of contents
[0 Introduction](https://github.com/MilaTost/ExpEvo_with_replicated_selection/blob/main/README.md#0-introduction) <br />
[1 Phenotypic data analysis](https://github.com/MilaTost/ExpEvo_with_replicated_selection/blob/main/README.md#1-phenotypic-data-analysis) <br />
&emsp;[1.1 Truncation selection thresholds](https://github.com/MilaTost/ExpEvo_with_replicated_selection/blob/main/README.md#11-truncation-selection-thresholds) <br />
&emsp;[1.2 Plant height measurements](https://github.com/MilaTost/ExpEvo_with_replicated_selection/blob/main/README.md#12-plant-height-measurements) <br />
&emsp;[1.3 Realized heritability from the Breeder's equation](https://github.com/MilaTost/ExpEvo_with_replicated_selection/blob/main/README.md#13-realized-heritability-from-the-breeders-equation) <br />
2 Pipeline for the analysis of GBS data adapted from [Wickland et al. 2013](https://github.com/dpwickland/GB-eaSy) <br />
[3 Filtering](https://github.com/MilaTost/ExpEvo_with_replicated_selection/blob/main/README.md#3-filtering) <br />
&emsp;[3.1 Filtering for individual samples with low coverage](https://github.com/MilaTost/ExpEvo_with_replicated_selection/blob/main/README.md#31-filtering-for-individual-samples-with-low-coverage) <br />
&emsp;[3.2 Filtering for read depth per sample](https://github.com/MilaTost/ExpEvo_with_replicated_selection/blob/main/README.md#32-filtering-for-read-depth-per-sample) <br />
&emsp;[3.3 Filtering for missingness](https://github.com/MilaTost/ExpEvo_with_replicated_selection/blob/main/README.md#33-filtering-for-missingness) <br />
&emsp;[3.4 Removal of non-diallelic and non-polymorphic markers](https://github.com/milaleonie/ExpEvo_with_replicated_selection/blob/main/README.md#34-removal-of-non-diallelic-and-non-polymorphic-markers) <br />
[4 Allele frequency estimation](https://github.com/MilaTost/ExpEvo_with_replicated_selection/blob/main/README.md#4-allele-frequency-estimation) <br />
[5 Estimation of the effective population size](https://github.com/MilaTost/ExpEvo_with_replicated_selection/blob/main/README.md#5-estimation-of-the-effective-population-size) <br />
[6 LD decay](https://github.com/MilaTost/ExpEvo_with_replicated_selection/blob/main/README.md#6-ld-decay) <br />
[7 Scan for selection signatures mapping with FST leveraging replicated selection](https://github.com/MilaTost/ExpEvo_with_replicated_selection/blob/main/README.md#7-scan-for-selection-signatures) <br />
[8 Test for selection: "Sauron plot"](https://github.com/MilaTost/ExpEvo_with_replicated_selection/blob/main/README.md#8-test-for-selection-sauron-plot) <br />
[9 Window based analysis](https://github.com/MilaTost/ExpEvo_with_replicated_selection/tree/main#9-window-based-analysis) <br />
[10 Comparison of significance thresholds and candidate gene identification](https://github.com/MilaTost/ExpEvo_with_replicated_selection/blob/main/README.md#10-comparison-of-significance-thresholds-and-candidate-gene-identification) <br />
&emsp; [10.1 Based on the empirical distribution](https://github.com/MilaTost/ExpEvo_with_replicated_selection/blob/main/README.md#101-based-on-the-empirical-distribution) <br />
&emsp; [10.2 Based on drift simulations](https://github.com/MilaTost/ExpEvo_with_replicated_selection/blob/main/README.md#102-based-on-drift-simulations) <br />
&emsp; [Simulation of drift](https://github.com/MilaTost/ExpEvo_with_replicated_selection/blob/main/README.md#simulation-of-drift) <br />
&emsp; [10.3 Based on the FDRfS](https://github.com/MilaTost/ExpEvo_with_replicated_selection/tree/main#103-based-on-the-fdrfs) <br />
[11 Haplotype estimation and haplotype block calculation](https://github.com/MilaTost/ExpEvo_with_replicated_selection/tree/main#11-haplotype-estimation-and-haplotype-block-calculation) <br />
[12 Generation of the LDheatmap](https://github.com/MilaTost/ExpEvo_with_replicated_selection/tree/main#11-haplotype-estimation-and-haplotype-block-calculation) <br />
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
statistic. <br /> 
**Experimental design of an experimental evolution study with replicated selection** <br />
<img width="500" alt="Figure_1" src="https://github.com/MilaTost/ExpEvo_with_replicated_selection/assets/63467079/558bfbe8-5844-451a-ab2c-229f72488a92">

## 1 Phenotypic data analysis
### 1.1 Truncation selection thresholds
The truncation thresholds of the 5% tallest or 5% shortest plants within the populations can be calculated like this:
```{r}
quantile(Short_Population_1$PlantHeight, 0.05, na.rm = TRUE)
quantile(Short_Population_2$PlantHeight, 0.05, na.rm = TRUE)
quantile(Tall_Population_1$PlantHeight, 0.05, na.rm = TRUE)
quantile(Tall_Population_2$PlantHeight, 0.05, na.rm = TRUE)
```
### 1.2 Plant height measurements
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
### 1.3 Realized heritability from the Breeder's equation
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
The Rscript contains many functions from the `vcfR` package from [Knaus and Grünwald 2018](https://github.com/knausb/vcfR#:~:text=VcfR%20is%20an%20R%20package%20intended%20to%20allow,rapidly%20read%20from%20and%20write%20to%20VCF%20files.).
<br />

## 4 Allele frequency estimation 
The allele frequency estimation was done with the `Allele_frequency_estimation.R` rscript.<br/> <br />
          
## 5 Estimation of the effective population size  
The effective population size was calculated using the known demographic parameters of the populations and based off of the formula: 
          <img width="120" alt="effective_pop_size" src="https://user-images.githubusercontent.com/63467079/171614297-2936b14e-3c6c-4914-b0b5-6e6cca2c5fc7.png"> <br /> where ![formula](https://render.githubusercontent.com/render/math?math=N_{m}) is the number of male and ![formula](https://render.githubusercontent.com/render/math?math=N_{f}) the number of female individuals (Crow and Kimura, 1970). In every generation the ~5% shortest or tallest plants were harvested, amounting to ~250 plants. Under the assumption of complete random mating, all 5000 plants could contribute as male parents for the next generation, the effective population size was calculated with the `Estimation_of_effective_population_size.R` rscript.<br /> <br />
         
In field trials, the assumption of random mating is violated because of the effects of assortative mating due to varying flowering dates and limited spatial pollen dispersion (Allard, 1999). Therefore, we evaluated the flowering dates of the 96 randomly chosen plants from each subpopulation to approximate the number of simultaneously flowering tassels and silks. It was assumed that silks remain receptive up to 5 days after silk emergence (Nieh et al., 2014), so that we calculated the number of flowering tassels during this time interval and projected it onto the entire subpopulation. This calculation of this is also available in the  `Estimation_of_effective_population_size.R` rscript. <br /> <br />

## 6 LD decay
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
The sum of FST was calculated between the two non-redundant comparisons of divergently selected subpopulations (Short 1 vs Tall 1; Short 2 vs Tall 2) and subpopulation selected in the same direction (Short 1 vs Short; Tall 1 vs Tall 2). These values are required for the calculation of the false discovery rate for selection (FDRfS). For plotting we only use the FSTSum value calculated between divergently selected subpopulations (Short 1 vs Tall 1; Short 2 vs Tall 2). <br />        

## 8 Visulization of selection and drift: "Sauron plot"      
<img src="https://user-images.githubusercontent.com/63467079/171002013-a510cf2f-5150-4808-bff8-b3a5d92d40b4.png" width="400" height="400"> <br /> <br /> 
This so-called "Sauron plot" depicts how the FDRfS was computed and provides a visualization of the scope of drift and selection. The Sauron plot comes from [Turner and Miller (2012)](http://www.genetics.org/content/suppl/2012/03/30/genetics.112.139337.DC1), but it can be also created for the <img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> statisitc. Sauron plot of genetic differentiation for FSTSum observed between the subpopulations selected in the same direction (blue) and in opposite directions (red). Each dot represents one SNP. The transparent red colored edges correspond to a false discovery rate (FDR) for selection < 5%. The y- and x-axis correspond to the range of <img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> values. 
The "Sauron plot" for the <img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> does not even look like the eye of [Sauron](https://twitter.com/strnr/status/457201981007073280) anymore, as the "Sauron plot" from [Turner and Miller (2012)](http://www.genetics.org/content/suppl/2012/03/30/genetics.112.139337.DC1) did. The Sauron plot for the allele frequency difference is created by the `create_sauron_plot_FST()` function, which is available in the `FstSum_value_calculation.R` script. <br /> <br />

## 9 Window-based analysis
We implemented a window-based analysis to assess linked selection. A cubic smoothing spline was applied to the single-marker based FSTSum values. From the fitted spline, inflection points were calculated and used as window boundaries for newly defined regions. This analysis was performed with the R package GenWin from [Beissinger et al., 2015](https://gsejournal.biomedcentral.com/articles/10.1186/s12711-015-0105-9). Briefly, the GenWin package returns window boundaries that we used to define regions with linked markers that may be analyzed together ([Beissinger et al., 2015). The calculation of window boundaries with the `splineAnalyze()` function from the [GenWin package](https://cran.r-project.org/web/packages/GenWin/index.html) are available in the `Window_based_analysis.R`. <br /> 

## 10 Calculation of significance thresholds
Significance thresholds for selection were calculated three ways: 1) based on the empirical distribution; 2) based on drift simulations; and 3) based on the false discovery rate for selection (FDRfS) [Turner and Miller (2012)](http://www.genetics.org/content/suppl/2012/03/30/genetics.112.139337.DC1). <br />
The calculation of significance thresholds and the plotting functions are contained in the `Significance_thresholds_and_plotting.R`. <br /> <br />
### 10.1 Based on the empirical distribution
The significance thresholds based on the empirical distribution, were calculated by taking the 99.9th and 99.99th percentile of the empirical distribution of the <img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> :
```{r}
FST_sig_thres_1 <- quantile(FST_values_od_cor$Fst, probs = 0.9999, na.rm = TRUE)
FST_sig_thres_2 <- quantile(FST_values_od_cor$Fst, probs = 0.999, na.rm = TRUE)    
```          
The significance threshold is stored, so it can be used later directly for plotting. <br /> <br />
### 10.2 Based on drift simulations 
The significance thresholds based on drift simulations were calculated in the `Simulation_of_drift.R` and then only retrieved from this script. The simulation of drift is described below and the script is also available in the repository.<br /> <br />
  
### Simulation of Drift
The simulation of drift was conducted by using the `DriftSimulator.R` from [Beissinger (2021)](http://beissingerlab.github.io/Software/). The `DriftSimulator.R` script was run with a drift simulation script similar to the one from [Kumar et al., 2021](https://academic.oup.com/pcp/article/62/7/1199/6279219), which enables the implementation of the drift simulator over a large set of markers. The script, which enables the simulation of drift over a large set of markers is available as `Run_drift_simulator_over_many_markers.R`. The `Run_drift_simulator_over_many_markers.R` script contains a function which simulates drift acting a single locus. We ran 5,000,000 simulations. For each simulation, initial allele frequencies were set based on allele frequency spectrum observed in generation 0 [Gyawali et al., 2019](https://pubmed.ncbi.nlm.nih.gov/31590656/) . Drift was simulated with a population size of 5000 individuals with 250 female and 5000 male individuals for three generations. We assumed that every ear contributed ~500 kernels. Every kernel could have been pollinated by one of the male parents, which resulted in much higher harvest than 5000 kernels. Therefore, 5000 kernels were randomly drawn from the entire harvest to represent the seeds planted for the next generation. In the third generation, 96 individuals out of 5000 were sampled to represent the individuals that were actually genotyped (Turner et al., 2011; Kumar et al., 2021). Variable marker coverage was also simulated; marker coverage was sampled from a uniform distribution between 40 and 96 observations per marker to match our filtering process [Turner et al., 2011](http://www.genetics.org/content/suppl/2012/03/30/genetics.112.139337.DC1). Additionally, the marker coverage was also sampled, the minimal marker coverage was set to at least 40 out of 96 observations at a marker, so that the marker coverage was always sampled between 40 to 96 observations per marker [Turner et al., 2011](http://www.genetics.org/content/suppl/2012/03/30/genetics.112.139337.DC1). <br /> <br /> The drift simulator is available as `DriftSimulator.R` script. <br /> <br />
<img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> values were calculated for all simulated markers, which corresponded in our case to 500,000 simulations. We summed those up and choosed the the 99.9999th percentile of the emprirical distribution of all observations as significance threshold, similar to [Kumar et al., 2021](https://academic.oup.com/pcp/article/62/7/1199/6279219). <br /> <br />

### 10.3 Based on the FDRfS
**FDRfS for all possible values of the statistics**
The function `calculate_FDR_for_selection()` generates a table to show all posible values of a statistic and the number of observed markers diverged between subpopulation selected in the same and opposite directions at a certain value. The FDR for selection is received by dividing the number of observed markers diverged between subpopulation selected in the same direction by the number of observed markers diverged between subpopulation selected in opposite directions. <br />
The calculation of the FDR for selection is also demonstrated with the following table which shows for different statistics the number of markers diverged between subpopulations selected in the same and opposite directions and the corresponding FDR for selections:          
|Statistic|Markers diverged same direction|Markers diverged opposite directions|FDR for selection|
|---------|-------------------------------|------------------------------------|-----------------|
|0	|1417796	                      |1421826	                   |0.9972           |
|0.01	|1104313	                      |1115771	                   |0.9897           |
|0.02	|400204	                      |415091	                             |0.9641           |
|0.03	|243600	                      |259984	                             |0.9370           |
|...	|...	                      |...	                             |...              |
|0.674	|1	                      |41	                             |0.0244           |
|0.675	|1	                      |40	                             |0.0250           |
|0.676	|1	                      |40	                             |0.0250           |
|0.677	|1	                      |39	                             |0.0256           |
|...	|...	                      |...	                             |...              |
|0.752	|1	                      |19	                             |0.0526           |
|0.753	|0		            |19	                             |0                |

**A table like this is automatically generated by the** `calculate_FDR_for_selection()` **function.** <br /> <br /> 
Even though, the significance threshold based on the FDR for selection is already printed by the previous function, the following function returns the threshold so it can be used directly in the Manhatten plot or to check the overlap between all the different statistics. The FDRfS can be choosen in the function. In our case we choosed a FDRfS < 5%. <br /> <br /> 

## Manhattan plot
In the [Manhatten plot](https://en.wikipedia.org/wiki/Manhattan_plot#:~:text=A%20Manhattan%20plot%20is%20a%20type%20of%20scatter,genome-wide%20association%20studies%20%28GWAS%29%20to%20display%20significant%20SNPs.) the positions of the markers are plotted against the <img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> value observed at this marker. <br /> <br /> 
<img src="https://github.com/MilaTost/ExpEvo_with_replicated_selection/assets/63467079/663ee6e3-7891-4fd9-9aef-7f8a470a7284.png" width="600">
**Observations between the subpopulations selected in opposite directions (A) and the same direction (B) and expressed by Wstat values across regions (C, D) with the significance thresholds based on the 99.9th percentile of the empirical distribution (light purple), based on the 99.99th percentile of the empirical distribution (dark purple), drift simulations (pink), and by the false discovery rate for selection (FDRfS) (yellow).** <br /> <br /> 
The plotting functions are available in the `Plot_Manhattan_plot_with_WStat_and_FSTSum.R` script. <br /> <br /> 

## 11 Haplotype estimation and haplotype block calculation
The data were phased using fastPHASE version 1.4.8 [(Sheet and Stephens, 2006)](https://stephenslab.uchicago.edu/software.html#fastphase) for further dissection. We phased the genomic data of the region from 9.437 to 10.457 Mb on chromosome 3 with 10 iterations of the expectation-maximization (EM) algorithm (Sheet and Stephens, 2006). The model was supplied with labels indicating the different subpopulations. The scripts are available as `Run_fastPHASE_random9MB_region_on_chr3.bash`. Before we ran fastPHASE, we prepared the data with the Rscript `Prepare_data_for_fastPHASE.R`. Haplotype blocks were calculated with the R package HaploBlocker [Pook et al., 2019](https://github.com/tpook92/HaploBlocker). The results look like this: <br />
<img src="https://github.com/MilaTost/ExpEvo_with_replicated_selection/assets/63467079/76e648e6-fd51-41ab-bc79-6a9a201b0c80.png" width="600"> <br /> 
The haplotype block calculation and plotting of the blocks are contained in the `Create_Haplotype_figure.R` script. <br /> <br />

## Generation of the LDheatmap














