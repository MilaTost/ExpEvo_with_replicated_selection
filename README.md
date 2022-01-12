# Selection_signature_mapping_scripts_with_replicated_selection
## Table of contents
#### 1 [Pipeline for the analysis of GBS data adapted from][Wickland et al. 2013](https://github.com/dpwickland/GB-eaSy)
#### 2 Filtering
#### 3 <img src="https://render.githubusercontent.com/render/math?math=F_{ST}"> with replicated selection
#### 4 Allele frequency differences
#### 5 Significance threshold based on the FDR for selection
#### 6 Significance thresholds based on drift simulations and the empiric distribution
#### 7 Simulation of Drift
#### 8 Sauron plot
#### 9 Other plotting scripts


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

