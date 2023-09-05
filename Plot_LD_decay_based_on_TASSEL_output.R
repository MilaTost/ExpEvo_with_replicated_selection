# Calculations for expected values of r2 under drift equilibrium [Hill and Weir (1988)
# implemented in Remington, et al. (2001)]
rm(list=ls())
# Calculate and save LD statistics from TASSEL using sliding window of 50 SNPs and maf = 0.05
# read in table with NaN and distance calculated from TASSEL or wherever
cat("Installation of the packages starts!", "\n")
install.packages("ggplot2",
                 lib = "/home/uni08/mtost/R/x86_64-pc-linux-gnu-library/4.1",
                 repos='http://cran.rstudio.com/')
library(ggplot2)
### Test plotting on local machine ---------------------------------------------
#setwd("C:/Users/mtost/Documents/Shoepeg_paper_final_versions/Scripts/Test_dir/")
#ld <- read.delim("small_ld_output_TASSEL.txt",stringsAsFactors = FALSE,header=TRUE, sep = "\t") 
#windowsFonts(my = windowsFont('Calibri'))
###############################################################################
### 		For more detail, please see following reference
###############################################################################
## Remington, D. L., Thornsberry, J. M., Matsuoka, Y.,
## Wilson, L. M., Whitt, S. R., Doebley, J., ... & Buckler, E. S. (2001). Structure of linkage disequilibrium
## and phenotypic associations in the maize genome. 
## Proceedings of the national academy of sciences, 98(20), 11479-11484. 
## https://doi.org/10.1073/pnas.201394398
###############################################################################
setwd("/usr/users/mtost/Shoepeg_resubmission_new_analysis/Plot_LD_decay_with_TASSEL_output/")
result_dir <- "/usr/users/mtost/Shoepeg_resubmission_new_analysis/Results/"
# import TASSEL LD output file
ld <- read.delim("2022_GB10_Shoepeg_after_filtering.ld.win2000.txt",stringsAsFactors = FALSE,header=TRUE, sep = "\t") 
cat("The data was read in!","\n")
##remove sites that have NaN for distance or r2
ld_sub <- ld[ld$R.2 != "NaN",]
ld_sub$dist <- as.numeric(ld_sub$Dist_bp)
ld_sub2 <- ld_sub[ld_sub$dist != "NaN",]
ld_sub2$rsq <- ld_sub2$R.2

file <- ld_sub2[,c(1,2,7,8,15:19)]

# C values range from about 0.5 to 2, start with 0.1
Cstart <- c(C=0.1)

# fit a non linear model using the arbitrary C value, 
# N is the number of the genotypes that have the SNP site
modelC <- nls(rsq ~ ( (10+C*dist)/( (2+C*dist) * (11+C*dist) ) ) * 
                ( 1+( (3+C*dist) * (12+12*C*dist+(C*dist)^2) ) / ( 2*N*(2+C*dist) * (11+C*dist) ) ), 
              data=file, start=Cstart, control=nls.control(maxiter=100))

# extract rho, the recombination parameter in 4Nr
rho <- summary(modelC)$parameters[1]

# feed in the new value of rho to obtain LD values adjusted for their distances along the chromosome/genome
newrsq <- ( (10+rho*file$dist) / ( (2+rho*file$dist) * (11+rho*file$dist) ) ) *
  ( 1 + ( (3+rho * file$dist) * (12+12*rho*file$dist + (rho*file$dist)^2) ) / 
      (2*file$N*(2+rho*file$dist) * (11+rho*file$dist) ) )

newfile <- data.frame(file$dist, newrsq)

maxld <- max(newfile$newrsq,na.rm=TRUE) #using max LD value from adjusted data
halfdecay <- maxld*0.5
halfdecaydist <- newfile$file.dist[which.min(abs(newfile$newrsq-halfdecay))]
cat("The maxld is ", maxld, "\n")
cat("The halfdecay is at a distance", halfdecaydist, "\n")
ld_01 <- 0.1
ld_01_decaydist <- newfile$file.dist[which.min(abs(newfile$newrsq-ld_01))]
cat("The decay distance of R2 = 0.1 is", ld_01_decaydist, "\n")
newfile <- newfile[order(newfile$file.dist),]
after_2000b <- file[file$dist <= 2000, ]
mean_decay_after_2000b <- mean(as.numeric(after_2000b$rsq), na.rm = TRUE)
cat("The mean decay after 2000 bp is", mean_decay_after_2000b, "\n")
after_1000b <- file[file$dist <= 1000, ]
mean_decay_after_1000b <- mean(as.numeric(after_1000b$rsq), na.rm = TRUE)
# plotting the values
xlim_file <- c(min(file$dist, na.rm = TRUE), max(file$dist, na.rm = TRUE))
ylim_file <- c(min(file$rsq, na.rm = TRUE), max(file$rsq, na.rm = TRUE))
xlim_newfile <- c(min(newfile$file.dist, na.rm = TRUE), max(newfile$file.dist, na.rm = TRUE))
ylim_newfile <- c(min(newfile$newrsq, na.rm = TRUE), max(newfile$newrsq, na.rm = TRUE))
cat("Ylim of the file is", ylim_file, "\n")
cat("Xlim of the file is", xlim_file, "\n")
cat("Ylim of the newfile is", ylim_newfile, "\n")
cat("Xlim of the newfile is", xlim_newfile, "\n")
cat("Plotting starts","\n")
newfile$newrsq <- round(newfile$newrsq, 2)
newfile$file.dist
plot_r2 <- ggplot()+
  geom_point(data = file, 
             aes(x = dist, y = rsq))+
  #geom_line(data = newfile, 
            #aes(x = file.dist, y = newrsq), colour = "red", linewidth = 1)+
  geom_smooth(data = newfile, 
              aes(x = file.dist, y = newrsq))+
  theme(text = element_text(size = 14, #family = "my"
                            ),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_line(color = "black"),
        strip.text.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black", face = "bold"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black", face = "bold"),
        axis.text.x = element_text(size = 14, colour = "black"),
        legend.text = element_text(size = 14, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"))+
  geom_hline(aes(yintercept = 0.1), colour = "palegreen3")+
  geom_vline(aes(xintercept = ld_01_decaydist), colour = "palegreen3")+
  labs(x = "Distance between markers in bp",
       y = "LD decay")+
  ylim(0,1)+
  xlim(1,104759819)
plot_r2
cat("Plotting is finsished", "\n")
cat("Saving as .png starts", "\n")
#result_dir <- "C:/Users/mtost/Documents/Shoepeg_paper_final_versions/Figures/New_figures/"
ggsave(paste0(result_dir, Sys.Date(),"_V9_LD_decay_ggplot_TASSEL_Calibri_font.png"),
       plot_r2,
       height = 10,
       width = 10,
       dpi = 1800)
q()