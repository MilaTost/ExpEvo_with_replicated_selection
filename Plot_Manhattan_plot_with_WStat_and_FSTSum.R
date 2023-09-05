# FDRfS significance thresholds based on Wstat based on the FST OD and SD ------
rm(list = ls())
#### Loading of the packages ---------------------------------------------------
#install.packages("ggplot2", dependencies = TRUE, repos = "http://cran.us.r-project.org")
#remove.packages("vctrs")
# this worked to fix ggplot
library(data.table)
library(stringr)
library(foreach)
library(parallel)
library(tidyverse)
library(doMC)
library(gridExtra)
library(ggplot2)
cat("Packages are loaded.","\n")
cat("Core registration starts.","\n")
cores <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK'))
cores <- detectCores(all.tests = FALSE, logical = TRUE)
registerDoMC(cores)
#### Load the data -------------------------------------------------------------
###### Load the WStat values ---------------------------------------------------
setwd("C:/Users/mtost/Documents/Shoepeg_paper_final_versions/Scripts/Results/")
WStat_OD_M3 <- fread("2022-12-19_Method_3_window_boundaries_based_on_FSTSum_OD_NA_excluded.txt")
WStat_SD_M3 <- fread("2022-12-19_Method_3_window_boundaries_based_on_FSTSum_SD_NA_excluded.txt")
dt_FDRfS_M3 <- fread("2022-12-19_Method_3_dt_FDRfS_WStat_calc_on_FSTSum.txt")
###### Load the FST values -----------------------------------------------------
setwd("C:/Users/mtost/Documents/Shoepeg_paper_final_versions/Scripts/Results/")
FST_values_od_cor <- fread("2022_02_16_GB1005_corrected_dt_Fst_values.txt",
                           skip = 1)
#FST_values_od_cor <- fread("2023-02-28_FST_values.txt",
                           #skip = 1)
dim(FST_values_od_cor)
#rm(FST_values_od_cor)
colnames_for_dt <- read.table("2022_02_16_GB1005_corrected_dt_Fst_values.txt",
                              nrows = 1, header = TRUE)
#colnames_for_dt <- read.table("2023-02-28_FST_values.txt",
                              #nrows = 1, header = TRUE)
FST_values_od_cor <- FST_values_od_cor[,2:ncol(FST_values_od_cor)]
colnames(FST_values_od_cor) <- colnames(colnames_for_dt)
FST_values_od_cor <- as.data.frame(FST_values_od_cor)
#### Calculate the FSTSum ------------------------------------------------------
calc_FSTSum_NAs_excluded_entire_dt <- function(comparison_1,
                                               comparison_2){
  calc_FSTSum_NAs_excluded <- function(i){
    if(!is.na(comparison_1[i]) & !is.na(comparison_2[i])){
      FSTSum_value <- comparison_1[i] + comparison_2[i]
    }
    if(is.na(comparison_1[i]) & !is.na(comparison_2[i])){
      FSTSum_value <- NA
    }
    if(!is.na(comparison_1[i]) & is.na(comparison_2[i])){
      FSTSum_value <- NA
    }
    return(FSTSum_value) 
  }
  dt_FSTSum_values <- mclapply(1:length(comparison_1), calc_FSTSum_NAs_excluded)
  dt_FSTSum_values <- unlist(dt_FSTSum_values)
  return(dt_FSTSum_values)
}
FSTSum_OD <- calc_FSTSum_NAs_excluded_entire_dt(comparison_1 = FST_values_od_cor$FST_value_opposite_dir1,
                                                comparison_2 = FST_values_od_cor$FST_value_opposite_dir2)
FSTSum_SD <- calc_FSTSum_NAs_excluded_entire_dt(comparison_1 = FST_values_od_cor$FST_value_same_dir1,
                                                comparison_2 = FST_values_od_cor$FST_value_same_dir2)
FST_values_od_cor$FSTSum_OD <- FSTSum_OD
FST_values_od_cor$FSTSum_SD <- FSTSum_SD

rm(FSTSum_SD, FSTSum_OD)
###### Load the parameters -----------------------------------------------------
sig_threshold_WStat <- 35
dt_FDRfS_M3
#### Get the regions with higher WStat values ----------------------------------
region_high_WStat_OD <- WStat_OD_M3[which(WStat_OD_M3$Wstat > sig_threshold_WStat),]
region_high_WStat_SD <- WStat_SD_M3[which(WStat_SD_M3$Wstat > sig_threshold_WStat),]
#### Get the region under selection --------------------------------------------
(9953933-9911000)/1000
-9923000

#### Plot the entire genome ----------------------------------------------------
method_GenWin = 3
### How many markers would have been detected? ---------------------------------
OD1_sig_thres_dis_999_perc <- quantile(FST_values_od_cor$FST_value_opposite_dir1, probs = 0.999, na.rm = TRUE)  
OD1_sig_thres_dis_9999_perc <- quantile(FST_values_od_cor$FST_value_opposite_dir1, probs = 0.9999, na.rm = TRUE)
OD2_sig_thres_dis_999_perc <- quantile(FST_values_od_cor$FST_value_opposite_dir2, probs = 0.999, na.rm = TRUE)  
OD2_sig_thres_dis_9999_perc <- quantile(FST_values_od_cor$FST_value_opposite_dir2, probs = 0.9999, na.rm = TRUE)
OD_sig_thres_dis_999_perc <- quantile(FST_values_od_cor$FSTSum_OD, probs = 0.999, na.rm = TRUE)  
OD_sig_thres_dis_9999_perc <- quantile(FST_values_od_cor$FSTSum_OD, probs = 0.9999, na.rm = TRUE)
mean(length(which(FST_values_od_cor$FST_value_opposite_dir1 > OD1_sig_thres_dis_9999_perc)),
length(which(FST_values_od_cor$FST_value_opposite_dir2 > OD2_sig_thres_dis_9999_perc)))
mean(length(which(FST_values_od_cor$FST_value_opposite_dir1 > OD1_sig_thres_dis_999_perc)),
length(which(FST_values_od_cor$FST_value_opposite_dir2 > OD2_sig_thres_dis_999_perc)))
length(which(FST_values_od_cor$FSTSum_OD > OD_sig_thres_dis_999_perc))
length(which(FST_values_od_cor$FSTSum_OD > OD_sig_thres_dis_9999_perc))
sig_thres_FDRfS <- 0.662
length(which(FST_values_od_cor$FST_sum_OD > sig_thres_FDRfS))
length(which(FST_values_od_cor$FST_sum_SD > sig_thres_FDRfS))
env_sel_dt <- FST_values_od_cor[which(FST_values_od_cor$FST_value_same_dir1 > 0.3 & FST_values_od_cor$FST_value_same_dir2 < 0.3),]
env_sel_dt <- FST_values_od_cor[which(FST_values_od_cor$FSTSum_SD > sig_thres_FDRfS),]

dim(env_sel_dt)
env_sel_dt <- FST_values_od_cor[which(FST_values_od_cor$FST_value_same_dir1 < 0.331 & FST_values_od_cor$FST_value_same_dir2 > 0.331),]
dim(env_sel_dt)
#chr3:185529884..185529961
dim(env_sel_dt)
env_sel_dt
write.table(env_sel_dt,"C:/Users/mtost/Documents/Shoepeg_paper_final_versions/Scripts/Results/2023_02_27_Regions_potentially_under_sel_for_spec_env_cond.txt",sep = "  ", row.names = TRUE,
            quote = FALSE,
            col.names = TRUE)
#chr1:216085673..218921975
#chr8:172105113..172210878
#chr9:142919977..142982657
#chr9:143121351..143151744
#chr7:28654680..28654707
### Plot FSTSum ----------------------------------------------------------------
sig_thres_dis_999_perc <- quantile(FST_values_od_cor$FSTSum_OD, probs = 0.999, na.rm = TRUE)  
sig_thres_dis_9999_perc <- quantile(FST_values_od_cor$FSTSum_OD, probs = 0.9999, na.rm = TRUE)
sig_thres_drift_sim <-  0.2043
sig_thres_FDRfS <- 0.662
FST_values_od_cor$Chromosome <- as.factor(fct_inorder(FST_values_od_cor$Chromosome))
FST_values_od_cor$SNP_ID <- as.factor(fct_inorder(FST_values_od_cor$SNP_ID))
str(FST_values_od_cor)
min_y <- 0
max_y <- 1
length(which(FST_values_od_cor$FSTSum_OD > sig_thres_drift_sim))
plot_FSTSum_for_genome_OD <- function(data,
                                      tag_of_plot,
                                      min_y,
                                      max_y){
  font_size <- 14
  windowsFonts(my = windowsFont('Calibri'))
  my_pal_col <- (c("chr1" = "firebrick4", "chr2" = "firebrick1",
                   "chr3" = "firebrick4", "chr4" = "firebrick1",
                   "chr5" = "firebrick4", "chr6" = "firebrick1", 
                   "chr7" = "firebrick4", "chr8" = "firebrick1",
                   "chr9" = "firebrick4", "chr10" = "firebrick1"))
  plot_FSTSum <- ggplot()+
    geom_point(data = FST_values_od_cor,
               aes(x = fct_inorder(SNP_ID),
                   y = FSTSum_OD,
                   colour = Chromosome))+
    theme(text = element_text(size = font_size, colour = "black", family = "my"),
          panel.background = element_rect(fill = "white"),
          panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                          colour = "grey81"),
          panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                          colour = "grey100"),
          panel.border = element_blank(),
          axis.line.y = element_line(colour = "black"),
          axis.line.x = element_line(colour = "black"),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = font_size, colour = "black"),
          axis.title.x = element_text(size = font_size, colour = "black", face = "bold"),
          axis.title.y = element_text(size = font_size, colour = "black", face = "bold"),
          legend.position = "none")+
    labs(tag = tag_of_plot, x = "Position (bp)", y = "FSTSum")+
    geom_hline(aes(yintercept = sig_thres_FDRfS), linetype="solid", colour="gold")+
    geom_hline(aes(yintercept = sig_thres_drift_sim), colour="deeppink", linetype="dotted",)+
    geom_hline(aes(yintercept = as.numeric(sig_thres_dis_9999_perc)), linetype="dotted", colour="purple1")+
    geom_hline(aes(yintercept = as.numeric(sig_thres_dis_999_perc)), linetype="dotted", colour="purple4")+
    scale_color_manual(values = my_pal_col)+
    ylim(min_y, max_y)
  return(plot_FSTSum)
}
plot_FSTSum_OD <- plot_FSTSum_for_genome_OD(data = FST_values_od_cor,
                                              tag_of_plot = "A",
                                              min_y = min_y, 
                                              max_y = max_y)
plot_FSTSum_OD
plot_FSTSum_for_genome_SD <- function(data,
                                       tag_of_plot,
                                       min_y,
                                       max_y){
  font_size <- 14
  windowsFonts(my = windowsFont('Calibri'))
  my_pal_col <- (c("chr1" = "royalblue4", "chr2" = "royalblue2",
                   "chr3" = "royalblue4", "chr4" = "royalblue2",
                   "chr5" = "royalblue4", "chr6" = "royalblue2", 
                   "chr7" = "royalblue4", "chr8" = "royalblue2",
                   "chr9" = "royalblue4", "chr10" = "royalblue2"))
  plot_FSTSum <- ggplot()+
    geom_point(data = FST_values_od_cor,
               aes(x = fct_inorder(SNP_ID),
                   y = FSTSum_SD,
                   colour = Chromosome))+
    theme(text = element_text(size = font_size, colour = "black", family = "my"),
          panel.background = element_rect(fill = "white"),
          panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                          colour = "grey81"),
          panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                          colour = "grey100"),
          panel.border = element_blank(),
          axis.line.y = element_line(colour = "black"),
          axis.line.x = element_line(colour = "black"),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = font_size, colour = "black"),
          axis.title.x = element_text(size = font_size, colour = "black", face = "bold"),
          axis.title.y = element_text(size = font_size, colour = "black", face = "bold"),
          legend.position = "none")+
    geom_hline(aes(yintercept = sig_thres_FDRfS), linetype="solid", colour="gold")+
    geom_hline(aes(yintercept = sig_thres_drift_sim), colour="deeppink", linetype="dotted",)+
    geom_hline(aes(yintercept = as.numeric(sig_thres_dis_9999_perc)), linetype="dotted", colour="purple1")+
    geom_hline(aes(yintercept = as.numeric(sig_thres_dis_999_perc)), linetype="dotted", colour="purple4")+
    labs(tag = tag_of_plot, x = "Position (bp)", y = "FSTSum")+
    scale_color_manual(values = my_pal_col)+
    ylim(min_y, max_y)
  return(plot_FSTSum)
}
plot_FSTSum_SD <- plot_FSTSum_for_genome_SD(data = FST_values_od_cor,
                                              tag_of_plot = "C",
                                              min_y = min_y, 
                                              max_y = max_y)
result_dir <- "C:/Users/mtost/Documents/Shoepeg_paper_final_versions/Figures/New_figures/"
png(paste0(result_dir, Sys.Date(),"_FSTSum_plot_combined_V5.png"), width = 6, height = 8, units="in", res = 1800)
grid.arrange(plot_FSTSum_OD,
             plot_FSTSum_SD,
             nrow = 2,
             ncol = 1, 
             heights = c(4,4),
             widths = 6)
dev.off()
rm(FST_values_od_cor)
##### Plot the WStat calculated based on FSTSum --------------------------------
sig_threshold_WStat <- 35
WStat_OD <- as.data.frame(WStat_OD_M3)
WStat_SD <- as.data.frame(WStat_SD_M3)
rm(WStat_SD_M3, WStat_OD_M3)
sig_thres_dis_999_perc <- quantile(WStat_OD$Wstat, probs = 0.999, na.rm = TRUE)  
sig_thres_dis_9999_perc <- quantile(WStat_OD$Wstat, probs = 0.9999, na.rm = TRUE)
length(which(WStat_OD$Wstat > sig_thres_dis_999_perc))
length(which(WStat_OD$Wstat > sig_thres_dis_9999_perc))
windowsFonts(my = windowsFont('Calibri'))

all_WStat_values <- c(WStat_OD$Wstat,
                      WStat_SD$Wstat)
min_y <- min(all_WStat_values, na.rm = TRUE)
max_y <- max(all_WStat_values, na.rm = TRUE)
rm(all_WStat_values)
WStat_OD$Chromosome <- as.factor(WStat_OD$Chromosome)
WStat_OD$Window_ID_Start <- paste0(WStat_OD$Chromosome,"_", WStat_OD$WindowStart)
WStat_OD$Window_ID_Start <- as.factor(fct_inorder(WStat_OD$Window_ID_Start))
WStat_SD$Chromosome <- as.factor(WStat_SD$Chromosome)
WStat_SD$Window_ID_Start <- paste0(WStat_SD$Chromosome,"_", WStat_SD$WindowStart)
WStat_SD$Window_ID_Start <- as.factor(fct_inorder(WStat_SD$Window_ID_Start))
WStat_OD$Chromosome <- as.factor(WStat_OD$Chromosome)
WStat_SD$Chromosome <- as.factor(WStat_SD$Chromosome)
plot_WStatSum_for_genome_OD <- function(data,
                                        tag_of_plot,
                                        put_selected_regions_OD_Reg_Start,
                                        put_selected_regions_OD_Reg_Stop,
                                        FDRfS_threshold,
                                        min_y,
                                        max_y){
  font_size <- 14
  my_pal_col <- (c("1" = "firebrick4","2" = "firebrick1",
                   "3" = "firebrick4","4" = "firebrick1",
                   "5" = "firebrick4","6" = "firebrick1", 
                   "7" = "firebrick4","8" = "firebrick1",
                   "9" = "firebrick4","10" = "firebrick1"))
  plot_WStatSum <- ggplot()+
    geom_line(data = data,
              aes(x = Window_ID_Start, y = Wstat, 
                  colour = Chromosome, group = Chromosome))+
    theme(text = element_text(size = font_size, colour = "black", family = "my"),
          panel.background = element_rect(fill = "white"),
          panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                          colour = "grey81"),
          panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                          colour = "grey100"),
          panel.border = element_blank(),
          axis.line.y = element_line(colour = "black"),
          axis.line.x = element_line(colour = "black"),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(colour = "black"),
          axis.title.x = element_text(colour = "black", face = "bold"),
          axis.title.y = element_text(colour = "black", face = "bold"),
          legend.position = "none")+
    geom_hline(aes(yintercept = FDRfS_threshold), colour = "gold2")+
    geom_hline(aes(yintercept = as.numeric(sig_thres_dis_9999_perc)), linetype="dotted", colour="purple1")+
    geom_hline(aes(yintercept = as.numeric(sig_thres_dis_999_perc)), linetype="dotted", colour="purple4")+
    labs(tag = tag_of_plot, x = "Position (bp)", y = "WStat")+
    scale_colour_manual(values = my_pal_col)+
    ylim(min_y, max_y)
  return(plot_WStatSum)
}
plot_OD <- plot_WStatSum_for_genome_OD(data = WStat_OD,
                                       FDRfS_threshold = sig_threshold_WStat,
                                       tag_of_plot = "B",
                                       min_y = min_y, 
                                       max_y = max_y)
plot_OD
plot_WStatSum_for_genome_SD <- function(data,
                                        tag_of_plot,
                                        FDRfS_threshold,
                                        min_y,
                                        max_y){
  font_size <- 14
  my_pal_col <- (c("1" = "royalblue4", "2" = "royalblue2",
                   "3" = "royalblue4", "4" = "royalblue2",
                   "5" = "royalblue4", "6" = "royalblue2", 
                   "7" = "royalblue4", "8" = "royalblue2",
                   "9" = "royalblue4", "10" = "royalblue2"))
  plot_WStatSum <- ggplot()+
    geom_line(data = data,
              aes(x = Window_ID_Start, y = Wstat, 
                  colour = Chromosome, group = Chromosome))+
    theme(text = element_text(size = font_size, family = "my", colour = "black"),
          panel.background = element_rect(fill = "white"),
          panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                          colour = "grey81"),
          panel.grid.major = element_line(linewidth = 0.5, linetype = 'solid',
                                          colour = "grey100"),
          panel.border = element_blank(),
          axis.line.y = element_line(colour = "black"),
          axis.line.x = element_line(colour = "black"),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(colour = "black"),
          axis.title.y = element_text(colour = "black", face = "bold"),
          axis.title.x = element_text(colour = "black", face = "bold"),
          legend.position = "none")+
    geom_hline(aes(yintercept = FDRfS_threshold), colour = "gold2")+
    geom_hline(aes(yintercept = as.numeric(sig_thres_dis_9999_perc)), linetype="dotted", colour="purple1")+
    geom_hline(aes(yintercept = as.numeric(sig_thres_dis_999_perc)), linetype="dotted", colour="purple4")+
    labs(tag = tag_of_plot, x = "Position (bp)", y = "WStat")+
    scale_color_manual(values = my_pal_col)+
    ylim(min_y, max_y)
  return(plot_WStatSum)
}
plot_SD <- plot_WStatSum_for_genome_SD(data = WStat_SD,
                                       FDRfS_threshold = sig_threshold_WStat,
                                       tag_of_plot = "D",
                                       min_y = min_y, 
                                       max_y = max_y)
plot_SD
png(paste0(result_dir, Sys.Date(),"_Wstat_applied_on_FSTSum_plot_combined_v6.png"), width = 4, height = 8, units="in", res = 1800)
grid.arrange(plot_OD,
             plot_SD,
             nrow = 2,
             ncol = 1, 
             heights = c(4,4),
             widths = 4)
dev.off()
### Create the legend ----------------------------------------------------------
font_size <- 12
sig_thres_dis_999_perc <- quantile(FST_values_od_cor$FSTSum_OD, probs = 0.999, na.rm = TRUE)  
sig_thres_dis_9999_perc <- quantile(FST_values_od_cor$FSTSum_OD, probs = 0.9999, na.rm = TRUE)
sig_thres_drift_sim <-  0.2043
sig_thres_FDRfS <- 0.662
legend_chr_plot <- ggplot(data = FST_values_od_cor)+
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
        axis.title.y = element_text(size =font_size-2, family = "my", colour = "black", face = "bold"),
        axis.text.y = element_text(size =font_size-2, family = "my", colour = "black"),
        axis.title.x = element_text(size =font_size-2, family = "my", colour = "black", face = "bold"),
        axis.text.x = element_blank(),
        legend.text = element_text(size =font_size-2, family = "my", colour = "black"),
        legend.title = element_text(size =font_size-2, family = "my", colour = "black", face = "bold"),
        legend.key.size = unit(1.75, "cm"),
        legend.key.width = unit(0.5,"cm"),
        legend.key = element_rect(fill = "white"))+
  labs(y = "FST",x = "Position")+
  geom_hline(aes(yintercept = sig_thres_FDRfS, linetype=str_wrap("FDRfS",12)), colour="gold", linewidth = 1)+
  geom_hline(aes(yintercept = sig_thres_dis_9999_perc, linetype=str_wrap("Outlier threshold: 99.99th percentile",12)), colour="purple1", linewidth = 1)+
  geom_hline(aes(yintercept = sig_thres_dis_999_perc, linetype =str_wrap("Outlier threshold: 99.9th percentile",12)), colour="purple4", linewidth = 1)+
  geom_hline(aes(yintercept = sig_thres_drift_sim, linetype=str_wrap("Simulations of drift",12)), colour="deeppink", linewidth = 1)+
  scale_linetype_manual(name = str_wrap("Significance threshold",10), values = c(1,1,1,1),
                        guide = guide_legend(override.aes = list(colour = c("gold","purple1","purple4","deeppink"),
                                                                 linetype = c("solid","dotted","dotted","dotted"))))
legend_chr_plot
get_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 
legend <- get_legend(legend_chr_plot)
result_dir <- "C:/Users/mtost/Documents/Shoepeg_paper_final_versions/Figures/New_figures/"
ggsave(paste0(result_dir, Sys.Date(),"_Wstat_applied_on_FSTSum_plot_legend_V5.png"),
       legend,
       height = 4, width = 2,
       dpi =1800)
#### Combine the manhattan plots -----------------------------------------------
combined_WStat_Manhattan_plots <- grid.arrange(plot_OD,
                                               plot_SD,
                                               nrow = 1,
                                               ncol = 2, 
                                               heights = c(4),
                                               widths = c(4.1,4))
result_dir <- "C:/Users/mtost/Documents/Shoepeg_paper_final_versions/Figures/New_figures/"
# V1: Plot candidate gene region -----------------------------------------------
#### Load the data -------------------------------------------------------------
setwd("C:/Users/mtost/Documents/Shoepeg_paper_final_versions/Scripts/Results/")
WStat_OD_M3 <- fread("2022-12-19_Method_3_window_boundaries_based_on_FSTSum_OD_NA_excluded.txt")
WStat_SD_M3 <- fread("2022-12-19_Method_3_window_boundaries_based_on_FSTSum_SD_NA_excluded.txt")
dt_FDRfS_M3 <- fread("2022-12-19_Method_3_dt_FDRfS_WStat_calc_on_FSTSum.txt")
result_dir <- "C:/Users/mtost/Documents/Shoepeg_paper_final_versions/Figures/New_figures/"
setwd("C:/Users/mtost/Documents/Masterarbeit/Data_analysis/Final_data_sets_for_paper/")
candidate_genes <- read.table("candidate_gene_region_and_markers.txt",
                              header = TRUE)
#### Get the regions with higher WStat values ----------------------------------
sig_threshold_WStat <- 35
region_high_WStat_OD <- WStat_OD_M3[which(WStat_OD_M3$Wstat > sig_threshold_WStat),]
region_high_WStat_SD <- WStat_SD_M3[which(WStat_SD_M3$Wstat > sig_threshold_WStat),]
round(region_high_WStat_OD$Wstat,2)
region_high_WStat_OD$SNPcount
#### Load the plotting parameters ----------------------------------------------
windowsFonts(my = windowsFont('Calibri'))
font_size <- 12
#### Prepare the data for plotting ---------------------------------------------
FST_values_od_cor <- as.data.table(FST_values_od_cor)
FST_values_Chr3 <- FST_values_od_cor[Chromosome == "chr3"]
WStat_SD_M3_Chr3 <- WStat_SD_M3[Chromosome == "3"]
WStat_OD_M3_Chr3 <- WStat_OD_M3[Chromosome == "3"]
which(WStat_OD_M3$Wstat > sig_threshold_WStat)
WStat_SD_M3[which(WStat_SD_M3$Wstat > sig_threshold_WStat),]
WStat_OD_M3[which(WStat_OD_M3$Wstat > sig_threshold_WStat),]
rm(WStat_OD_M3, WStat_SD_M3, FST_values_od_cor)
rm(WStat_OD, WStat_SD)
max_y_WStat <- max(WStat_OD_M3_Chr3$Wstat)
min_y_WStat <- min(WStat_OD_M3_Chr3$Wstat)
#### Start plotting ------------------------------------------------------------
##### Method 3; FSTSum SD ------------------------------------------------------
candidate_gene_region_FSTSum_M3_SD <- ggplot()+
  geom_point(data = FST_values_Chr3, aes(x = Position, y = mean_FST_SD_NA_excluded), colour = "black")+
  theme(text = element_text(size =font_size, family = "my"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size =font_size, colour = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size =font_size-6, colour = "black"))+
  labs(y = "FSTSum",x = "Position", tag = "B")+
  annotate("text", x = 9953933, y = 0.75, label = "Marker found by Gyawali et al. 2019", colour = "darkorange", family = "my")+
  annotate("segment", x = 9953933, y = 1, xend = 9953933, yend = 0.5,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")), 
           colour = "darkorange")+
  annotate("text", x = 10062219, y = 0.8, label = "IAA8", colour = "tomato2", family = "my")+
  annotate("segment", x = 10062219, y = 1, xend = 10065595, yend = 0.5,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")),
           colour = "tomato2")+
  annotate("text", x = 10440993, y = 0.85, label = "Dwarf1", colour = "tomato2", family = "my")+
  annotate("segment", x = 10440993, y = 1, xend = 10443340, yend = 0.5,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")),
           colour = "tomato2")+
  xlim(min(candidate_genes$Start)-500000,max(candidate_genes$End)+500000)+
  ylim(0,1)
candidate_gene_region_FSTSum_M3_SD
##### Method 3; FSTSum OD ------------------------------------------------------
candidate_gene_region_FSTSum_M3_OD <- ggplot()+
  geom_point(data = FST_values_Chr3, aes(x = Position, y = mean_FST_OD_NA_excluded), colour = "black")+
  theme(text = element_text(size =font_size, family = "my"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size =font_size, colour = "black"),
        axis.title.y = element_text(size =font_size, colour = "black", face = "bold"),
        axis.text.x = element_text(size =font_size-6, colour = "black"),
        legend.text = element_text(size =font_size, colour = "black"),
        legend.title = element_text(size =font_size, colour = "black", face = "bold"))+
  labs(y = "FstSum", x = "Position", tag = "A")+
  annotate("text", x = 9953933, y = 0.75, label = "Marker found by Gyawali et al. 2019", colour = "darkorange", family = "my")+
  annotate("segment", x = 9953933, y = 1, xend = 9953933, yend = 0.5,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")), 
           colour = "darkorange")+
  annotate("text", x = 10062219, y = 0.8, label = "IAA8", colour = "tomato2", family = "my")+
  annotate("segment", x = 10062219, y = 1, xend = 10065595, yend = 0.5,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")),
           colour = "tomato2")+
  annotate("text", x = 10440993, y = 0.85, label = "Dwarf1", colour = "tomato2", family = "my")+
  annotate("segment", x = 10440993, y = 1, xend = 10443340, yend = 0.5,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")),
           colour = "tomato2")+
  annotate("text", x = region_high_WStat_OD$WindowStart[1], y = 0.9, label = "Significant region from Wstat approach", colour = "darkblue", family = "my")+
  annotate("rect", xmin = region_high_WStat_OD$WindowStart[1], xmax = region_high_WStat_OD$WindowStop[1], ymin = 0, ymax = 1, fill = "darkblue", alpha = 0.6)+
  xlim(min(candidate_genes$Start)-500000,max(candidate_genes$End)+500000)+
  ylim(0,1)
candidate_gene_region_FSTSum_M3_OD
##### Method 3; WStat OD ------------------------------------------------------
candidate_gene_region_WStat_M3_SD <- ggplot()+
  annotate("rect", xmin = WStat_SD_M3_Chr3$WindowStart, xmax = WStat_SD_M3_Chr3$WindowStop, ymin = 0, ymax = WStat_SD_M3_Chr3$Wstat, fill = "darkmagenta", alpha = 0.6)+
  theme(text = element_text(size =font_size, family = "my"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size =font_size, colour = "black"),
        axis.title.x = element_text(size =font_size, colour = "black", face = "bold"),
        axis.text.x = element_text(size =font_size-6, colour = "black"))+
  labs(y = "WStat",x = "Position", tag = "D")+
  annotate("text", x = 9953933, y = 0.75*52, label = "Marker found by Gyawali et al. 2019", colour = "darkorange", family = "my")+
  annotate("segment", x = 9953933, y = 1*52, xend = 9953933, yend = 0.5*52,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")), 
           colour = "darkorange")+
  annotate("text", x = 10062219, y = 0.8*52, label = "IAA8", colour = "tomato2", family = "my")+
  annotate("segment", x = 10062219, y = 1*52, xend = 10065595, yend = 0.5*52,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")),
           colour = "tomato2")+
  annotate("text", x = 10440993, y = 0.85*52, label = "Dwarf1", colour = "tomato2", family = "my")+
  annotate("segment", x = 10440993, y = 1*52, xend = 10443340, yend = 0.5*52,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")),
           colour = "tomato2")+
  xlim(min(candidate_genes$Start)-500000,max(candidate_genes$End)+500000)
candidate_gene_region_WStat_M3_SD
##### Method 3; Wstat OD ------------------------------------------------------
candidate_gene_region_WStat_M3_OD <- ggplot()+
  annotate("rect", xmin = WStat_OD_M3_Chr3$WindowStart, xmax = WStat_OD_M3_Chr3$WindowStop, ymin = 0, ymax = WStat_OD_M3_Chr3$Wstat, fill = "darkmagenta", alpha = 0.6)+
  theme(text = element_text(size = font_size, family = "my"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.title.x = element_text(size =font_size, colour = "black", face = "bold"),
        axis.title.y = element_text(size =font_size, colour = "black", face = "bold"),
        axis.text.y = element_text(size =font_size, colour = "black"),
        axis.text.x = element_text(size =font_size-2, colour = "black", angle = 45))+
  labs(y = "WStat",x = "Position", tag = "C")+
  annotate("text", x = 9953933, y = 0.75*max_y_WStat, label = "Marker found by Gyawali et al. 2019", colour = "darkorange", family = "my")+
  annotate("segment", x = 9953933, y = 1*max_y_WStat, xend = 9953933, yend = 0.5*max_y_WStat,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")), 
           colour = "darkorange")+
  annotate("text", x = 10062219, y = 0.8*max_y_WStat, label = "IAA8", colour = "tomato2", family = "my")+
  annotate("segment", x = 10062219, y = 1*max_y_WStat, xend = 10065595, yend = 0.5*max_y_WStat,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")),
           colour = "tomato2")+
  annotate("text", x = 10440993, y = 0.85*max_y_WStat, label = "Dwarf1", colour = "tomato2", family = "my")+
  annotate("segment", x = 10440993, y = 1*max_y_WStat, xend = 10443340, yend = 0.5*max_y_WStat,
           arrow = arrow(type = "closed", length = unit(0.02, "npc")),
           colour = "tomato2")+
  xlim(min(candidate_genes$Start)-500000,max(candidate_genes$End)+500000)
candidate_gene_region_WStat_M3_OD
#### Combine the plots ---------------------------------------------------------
M3_combined_candidate_gene_regions <- grid.arrange(candidate_gene_region_FSTSum_M3_OD,
                                                   candidate_gene_region_FSTSum_M3_SD,
                                                   candidate_gene_region_WStat_M3_OD,
                                                   candidate_gene_region_WStat_M3_SD,
                                                   nrow = 2,
                                                   ncol = 2, 
                                                   heights = c(6,6),
                                                   widths = c(6,6.5),
                                                   layout_matrix = rbind(c(1,2),
                                                                         c(3,4)))
setwd(result_dir)
ggsave("2022_12_20_FSTSum_and_WStat_based_on_FSTSum_M3_M4.png",
       M3_combined_candidate_gene_regions,
       height = 5, width = 10,
       dpi =1800)
# V2: Visualize haplotypes in the candidate gene region ------------------------
# Visualize haplotypes in the WStat region -------------------------------------
### Load haplotypes ------------------------------------------------------------
setwd("C:/Users/mtost/Documents/Shoepeg_paper_final_versions/Scripts/Results/")
Haplotypes <- fread("2022-12-30_haplotypes_of_CHR3_candidate_gene_regions.txt",
                    skip = 1)
Haplotypes_colnames <- read.table("2022-12-30_haplotypes_of_CHR3_candidate_gene_regions.txt",
                                  nrows = 1,
                                  header = TRUE)
colnames(Haplotypes) <- c("Row", colnames(Haplotypes_colnames))
dim(Haplotypes)
Haplotypes <- as.data.frame(Haplotypes)
Haplotypes[1:2, 1:5]
SNP_IDs <- Haplotypes[,2]
#### Get the regions -----------------------------------------------------------
start_SNPs <- 9911000
end_SNPs <- 9923000
dist <- 1500
marker_positions <- str_sub(SNP_IDs, 6)
which(marker_positions <= start_SNPs+dist & marker_positions >= start_SNPs-dist)
which(marker_positions <= end_SNPs+dist & marker_positions >= end_SNPs-dist)
SNP_IDs <- Haplotypes$Marker
pos_index_start <- 1562
pos_index_end <- 1616
Haplotypes_WStat <- Haplotypes[pos_index_start:pos_index_end, ]
SNP_IDs_WStat <- Haplotypes_WStat$Marker
#### Get different subpopulations ----------------------------------------------
index_Shoepeg_1 <- which(str_sub(colnames(Haplotypes_WStat), 1, 9)=="Shoepag_1")
index_Shoepeg_2 <- which(str_sub(colnames(Haplotypes_WStat), 1, 9)=="Shoepag_2")
index_Shoepeg_3 <- which(str_sub(colnames(Haplotypes_WStat), 1, 9)=="Shoepag_3")
index_Shoepeg_4 <- which(str_sub(colnames(Haplotypes_WStat), 1, 9)=="Shoepag_4")
Haplotypes_Shoepeg_1_WStat <- Haplotypes_WStat[,index_Shoepeg_1]
Haplotypes_Shoepeg_2_WStat <- Haplotypes_WStat[,index_Shoepeg_2]
Haplotypes_Shoepeg_3_WStat <- Haplotypes_WStat[,index_Shoepeg_3]
Haplotypes_Shoepeg_4_WStat <- Haplotypes_WStat[,index_Shoepeg_4]
rownames(Haplotypes_Shoepeg_1_WStat) <- SNP_IDs_WStat
rownames(Haplotypes_Shoepeg_2_WStat) <- SNP_IDs_WStat
rownames(Haplotypes_Shoepeg_3_WStat) <- SNP_IDs_WStat
rownames(Haplotypes_Shoepeg_4_WStat) <- SNP_IDs_WStat
#### Prepare the data for the geom_tile() plot ---------------------------------
prepare_data_for_tile_plotting <- function(genotype_data){
  y <- rownames(genotype_data)
  x <- colnames(genotype_data)
  rownames(genotype_data) <- NULL
  colnames(genotype_data) <- NULL
  transpose_marker_by_marker <- function(i){
    z <- t(genotype_data[i,])
    new_data <- cbind(x, rep(y[i], length(x)), z)
    colnames(new_data) <- c("x", "y", "z")
    new_data <- as.data.frame(new_data)
    return(new_data)
  }
  entire_dt <- foreach(i = 1:nrow(genotype_data), .combine = rbind) %do% transpose_marker_by_marker(i)
  return(entire_dt)
}
data_new_Shoepeg_1_WStat <- prepare_data_for_tile_plotting(genotype_data = Haplotypes_Shoepeg_1_WStat)
data_new_Shoepeg_2_WStat <- prepare_data_for_tile_plotting(genotype_data = Haplotypes_Shoepeg_2_WStat)
data_new_Shoepeg_3_WStat <- prepare_data_for_tile_plotting(genotype_data = Haplotypes_Shoepeg_3_WStat)
data_new_Shoepeg_4_WStat <- prepare_data_for_tile_plotting(genotype_data = Haplotypes_Shoepeg_4_WStat)
data_new_Shoepeg_1[1:100,]
# Heatmap ----------------------------------------------------------------------
my_pal_col_1 <- (c("A/G" = "darkorange1",
                   "G/A" = "darkorange1",
                   "T/A" = "sienna",
                   "A/T" = "sienna",
                   "A/C" = "mediumpurple",
                   "C/A" = "mediumpurple",
                   "C/T" = "cyan 2",
                   "T/C" = "cyan 2",
                   "G/T" = "olivedrab2",
                   "T/G" = "olivedrab2",
                   "C/G" = "springgreen",
                   "G/C" = "springgreen",
                   "."   = "grey21",
                   "A/A" = "tomato1",
                   "G/G" = "gold1",
                   "C/C" = "skyblue3", 
                   "T/T" = "forestgreen"))
#str(data_new_Shoepeg_1)
#data_new_Shoepeg_1 <- as.data.frame(data_new_Shoepeg_1)
#data_new_Shoepeg_1$y <- str_sub(data_new_Shoepeg_1$y, 6)
#data_new_Shoepeg_1$x <- str_sub(data_new_Shoepeg_1$x, 9)
#data_new_Shoepeg_1$x <- as.numeric(data_new_Shoepeg_1$x)
#data_new_Shoepeg_1$y <- as.numeric(data_new_Shoepeg_1$y)
#data_new_Shoepeg_1$x <- as.factor(data_new_Shoepeg_1$x)
#data_new_Shoepeg_1$z <- as.factor(data_new_Shoepeg_1$z)
#### Reduce the Haplotypes -----------------------------------------------------
windowsFonts(my = windowsFont('Calibri'))
font_size <- 10
#my_pal_col_1 <- (c("." = "grey21",
                  #"A/A" = "tomato1",
                   #"G/G" = "gold1",
                   #"C/C" = "skyblue3", 
                   #"T/T" = "springgreen"))
reduce_haplotypes <- function(data){
  data$z[which(data$z == "A/G")] <- "."
  data$z[which(data$z == "G/A")] <- "."
  data$z[which(data$z == "T/A")] <- "."
  data$z[which(data$z == "A/T")] <- "."
  data$z[which(data$z == "A/C")] <- "."
  data$z[which(data$z == "C/A")] <- "."
  data$z[which(data$z == "C/T")] <- "."
  data$z[which(data$z == "T/C")] <- "."
  data$z[which(data$z == "G/T")] <- "."
  data$z[which(data$z == "T/G")] <- "."
  data$z[which(data$z == "C/G")] <- "."
  data$z[which(data$z == "G/C")] <- "."
  return(data)
}
#data_new_Shoepeg_1 <- reduce_haplotypes(data = data_new_Shoepeg_1)
#data_new_Shoepeg_2 <- reduce_haplotypes(data = data_new_Shoepeg_2)
#data_new_Shoepeg_3 <- reduce_haplotypes(data = data_new_Shoepeg_3)
#data_new_Shoepeg_4 <- reduce_haplotypes(data = data_new_Shoepeg_4)
#### Prepare the data ----------------------------------------------------------
df_shoe1_WStat <- data.frame(
  y = factor(data_new_Shoepeg_1_WStat$x),
  x = data_new_Shoepeg_1_WStat$y,
  z = factor(data_new_Shoepeg_1_WStat$z)
)
df_shoe2_WStat <- data.frame(
  y = factor(data_new_Shoepeg_2_WStat$x),
  x = data_new_Shoepeg_2_WStat$y,
  z = factor(data_new_Shoepeg_2_WStat$z)
)
df_shoe3_WStat <- data.frame(
  y = factor(data_new_Shoepeg_3_WStat$x),
  x = data_new_Shoepeg_3_WStat$y,
  z = factor(data_new_Shoepeg_3_WStat$z)
)
df_shoe4_WStat <- data.frame(
  y = factor(data_new_Shoepeg_4_WStat$x),
  x = data_new_Shoepeg_4_WStat$y,
  z = factor(data_new_Shoepeg_4_WStat$z)
)
#### Plot Shoepeg 1 
regions_of_chr_Shoepeg_1 <- ggplot(df_shoe1_WStat, aes(y = y, x = x, fill = z, colour = z)) + 
  geom_tile()+
  theme(text = element_text(size = font_size, family = "my"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = font_size, family = "my", colour = "black"),
        axis.title.x = element_blank())+
  labs(y = str_wrap("Individual samples in short 1",12), x = "SNP Positions",
       tag = "E")+
  scale_fill_manual(values = my_pal_col_1)+
  scale_color_manual(values = my_pal_col_1)
regions_of_chr_Shoepeg_1
#### Plot Shoepeg 2 
regions_of_chr_Shoepeg_2 <- ggplot(df_shoe2_WStat, aes(y = y, x = x, fill = z, colour = z)) + 
  geom_tile()+
  theme(text = element_text(size = font_size, family = "my"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = font_size, family = "my", colour = "black"),
        axis.title.x = element_text(size = font_size, family = "my", colour = "black"))+
  labs(y = str_wrap("Individual samples in tall 2",12), x = "SNP Positions")+
  scale_fill_manual(values = my_pal_col_1)+
  scale_color_manual(values = my_pal_col_1)
regions_of_chr_Shoepeg_2
#### Plot Shoepeg 3
regions_of_chr_Shoepeg_3 <- ggplot(df_shoe3_WStat, aes(y = y, x = x, fill = z, colour = z)) + 
  geom_tile()+
  theme(text = element_text(size = font_size, family = "my"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = font_size, family = "my", colour = "black"),
        axis.title.x = element_blank())+
  labs(y = str_wrap("Individual samples in tall 1",12), x = "SNP Positions")+
  scale_fill_manual(values = my_pal_col_1)+
  scale_color_manual(values = my_pal_col_1)
regions_of_chr_Shoepeg_3
#### Plot Shoepeg 4
regions_of_chr_Shoepeg_4 <- ggplot(df_shoe4_WStat, aes(y = y, x = x, fill = z, colour = z)) + 
  geom_tile()+
  theme(text = element_text(size = font_size, family = "my"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = font_size, family = "my", colour = "black"),
        axis.title.x = element_blank())+
  labs(y = str_wrap("Individual samples in short 2",12), x = "SNP Positions")+
  scale_fill_manual(values = my_pal_col_1)+
  scale_color_manual(values = my_pal_col_1)
regions_of_chr_Shoepeg_4
### Get the legend
regions_of_chr_legend <- ggplot(df_shoe4_WStat, aes(y = y, x = x, fill = z, colour = z)) + 
  geom_tile()+
  theme(text = element_text(size = font_size, family = "my"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        legend.title = element_blank(),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = font_size, family = "my", colour = "black", face = "bold"),
        axis.title.x = element_text(size = font_size, family = "my", colour = "black", face = "bold"))+
  labs(y = "Individual samples", x = "SNP Positions")+
  scale_fill_manual(values = my_pal_col_1)+
  scale_color_manual(values = my_pal_col_1)
get_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 
legend <- get_legend(regions_of_chr_legend)
regions_of_chr_all_Shoepeg_pop <- grid.arrange(regions_of_chr_Shoepeg_1,
                                               regions_of_chr_Shoepeg_4,
                                               regions_of_chr_Shoepeg_3,
                                               regions_of_chr_Shoepeg_2, 
                                               legend,
                                               nrow = 4,
                                               ncol = 2, 
                                               heights = c(2, 2, 2, 2.1),
                                               widths = c(6,0.5),
                                               layout_matrix = cbind(c(1,2,3,4),
                                                                     c(5)))
ggsave("C:/Users/mtost/Documents/Shoepeg_paper_final_versions/Figures/New_figures/Haplotypes_of_candidate_gene_region_all_Shoepeg_pop_all_GTs_WStat_reg.png",
       regions_of_chr_all_Shoepeg_pop,
       height = 6, width = 10,
       dpi =1800)
# Visualize haplotypes in the FSTSum region -------------------------------------
### Load haplotypes ------------------------------------------------------------
setwd("C:/Users/mtost/Documents/Shoepeg_paper_final_versions/Scripts/Results/")
Haplotypes <- fread("2022-12-30_haplotypes_of_CHR3_candidate_gene_regions.txt",
                    skip = 1)
Haplotypes_colnames <- read.table("2022-12-30_haplotypes_of_CHR3_candidate_gene_regions.txt",
                                  nrows = 1,
                                  header = TRUE)
colnames(Haplotypes) <- c("Row", colnames(Haplotypes_colnames))
dim(Haplotypes)
Haplotypes <- as.data.frame(Haplotypes)
Haplotypes[1:2, 1:5]
SNP_IDs <- Haplotypes[,2]
#### Get the regions -----------------------------------------------------------
start_SNPs <- 10062671-100
end_SNPs <- 10062671+100
dist <- 1500
marker_positions <- str_sub(SNP_IDs, 6)
which(marker_positions <= start_SNPs+dist & marker_positions >= start_SNPs-dist)
which(marker_positions <= end_SNPs+dist & marker_positions >= end_SNPs-dist)
SNP_IDs <- Haplotypes$Marker
pos_index_start <- 1900
pos_index_end <- 1928
Haplotypes_FSTSum <- Haplotypes[pos_index_start:pos_index_end, ]
SNP_IDs_FSTSum <- Haplotypes_FSTSum$Marker
dim(Haplotypes_FSTSum)
#### Get different subpopulations ----------------------------------------------
index_Shoepeg_1 <- which(str_sub(colnames(Haplotypes_FSTSum), 1, 9)=="Shoepag_1")
index_Shoepeg_2 <- which(str_sub(colnames(Haplotypes_FSTSum), 1, 9)=="Shoepag_2")
index_Shoepeg_3 <- which(str_sub(colnames(Haplotypes_FSTSum), 1, 9)=="Shoepag_3")
index_Shoepeg_4 <- which(str_sub(colnames(Haplotypes_FSTSum), 1, 9)=="Shoepag_4")
Haplotypes_Shoepeg_1_FSTSum <- Haplotypes_FSTSum[,index_Shoepeg_1]
Haplotypes_Shoepeg_2_FSTSum <- Haplotypes_FSTSum[,index_Shoepeg_2]
Haplotypes_Shoepeg_3_FSTSum <- Haplotypes_FSTSum[,index_Shoepeg_3]
Haplotypes_Shoepeg_4_FSTSum <- Haplotypes_FSTSum[,index_Shoepeg_4]
rownames(Haplotypes_Shoepeg_1_FSTSum) <- SNP_IDs_FSTSum
rownames(Haplotypes_Shoepeg_2_FSTSum) <- SNP_IDs_FSTSum
rownames(Haplotypes_Shoepeg_3_FSTSum) <- SNP_IDs_FSTSum
rownames(Haplotypes_Shoepeg_4_FSTSum) <- SNP_IDs_FSTSum
#### Prepare the data for the geom_tile() plot ---------------------------------
prepare_data_for_tile_plotting <- function(genotype_data){
  y <- rownames(genotype_data)
  x <- colnames(genotype_data)
  rownames(genotype_data) <- NULL
  colnames(genotype_data) <- NULL
  transpose_marker_by_marker <- function(i){
    z <- t(genotype_data[i,])
    new_data <- cbind(x, rep(y[i], length(x)), z)
    colnames(new_data) <- c("x", "y", "z")
    new_data <- as.data.frame(new_data)
    return(new_data)
  }
  entire_dt <- foreach(i = 1:nrow(genotype_data), .combine = rbind) %do% transpose_marker_by_marker(i)
  return(entire_dt)
}
data_new_Shoepeg_1_FSTSum <- prepare_data_for_tile_plotting(genotype_data = Haplotypes_Shoepeg_1_FSTSum)
data_new_Shoepeg_2_FSTSum <- prepare_data_for_tile_plotting(genotype_data = Haplotypes_Shoepeg_2_FSTSum)
data_new_Shoepeg_3_FSTSum <- prepare_data_for_tile_plotting(genotype_data = Haplotypes_Shoepeg_3_FSTSum)
data_new_Shoepeg_4_FSTSum <- prepare_data_for_tile_plotting(genotype_data = Haplotypes_Shoepeg_4_FSTSum)
data_new_Shoepeg_1_FSTSum[1:100,]
# Heatmap ----------------------------------------------------------------------
my_pal_col_1 <- (c("A/G" = "darkorange1",
                   "G/A" = "darkorange1",
                   "T/A" = "sienna",
                   "A/T" = "sienna",
                   "A/C" = "mediumpurple",
                   "C/A" = "mediumpurple",
                   "C/T" = "cyan 2",
                   "T/C" = "cyan 2",
                   "G/T" = "olivedrab2",
                   "T/G" = "olivedrab2",
                   "C/G" = "springgreen",
                   "G/C" = "springgreen",
                   "."   = "grey21",
                   "A/A" = "tomato1",
                   "G/G" = "gold1",
                   "C/C" = "skyblue3", 
                   "T/T" = "forestgreen"))
#str(data_new_Shoepeg_1)
#data_new_Shoepeg_1 <- as.data.frame(data_new_Shoepeg_1)
#data_new_Shoepeg_1$y <- str_sub(data_new_Shoepeg_1$y, 6)
#data_new_Shoepeg_1$x <- str_sub(data_new_Shoepeg_1$x, 9)
#data_new_Shoepeg_1$x <- as.numeric(data_new_Shoepeg_1$x)
#data_new_Shoepeg_1$y <- as.numeric(data_new_Shoepeg_1$y)
#data_new_Shoepeg_1$x <- as.factor(data_new_Shoepeg_1$x)
#data_new_Shoepeg_1$z <- as.factor(data_new_Shoepeg_1$z)
#### Reduce the Haplotypes -----------------------------------------------------
windowsFonts(my = windowsFont('Calibri'))
font_size <- 10
#my_pal_col_1 <- (c("." = "grey21",
#"A/A" = "tomato1",
#"G/G" = "gold1",
#"C/C" = "skyblue3", 
#"T/T" = "springgreen"))
reduce_haplotypes <- function(data){
  data$z[which(data$z == "A/G")] <- "."
  data$z[which(data$z == "G/A")] <- "."
  data$z[which(data$z == "T/A")] <- "."
  data$z[which(data$z == "A/T")] <- "."
  data$z[which(data$z == "A/C")] <- "."
  data$z[which(data$z == "C/A")] <- "."
  data$z[which(data$z == "C/T")] <- "."
  data$z[which(data$z == "T/C")] <- "."
  data$z[which(data$z == "G/T")] <- "."
  data$z[which(data$z == "T/G")] <- "."
  data$z[which(data$z == "C/G")] <- "."
  data$z[which(data$z == "G/C")] <- "."
  return(data)
}
#data_new_Shoepeg_1 <- reduce_haplotypes(data = data_new_Shoepeg_1)
#data_new_Shoepeg_2 <- reduce_haplotypes(data = data_new_Shoepeg_2)
#data_new_Shoepeg_3 <- reduce_haplotypes(data = data_new_Shoepeg_3)
#data_new_Shoepeg_4 <- reduce_haplotypes(data = data_new_Shoepeg_4)
#### Prepare the data ----------------------------------------------------------
df_shoe1_FSTSum <- data.frame(
  y = factor(data_new_Shoepeg_1_FSTSum$x),
  x = data_new_Shoepeg_1_FSTSum$y,
  z = factor(data_new_Shoepeg_1_FSTSum$z)
)
df_shoe2_FSTSum <- data.frame(
  y = factor(data_new_Shoepeg_2_FSTSum$x),
  x = data_new_Shoepeg_2_FSTSum$y,
  z = factor(data_new_Shoepeg_2_FSTSum$z)
)
df_shoe3_FSTSum <- data.frame(
  y = factor(data_new_Shoepeg_3_FSTSum$x),
  x = data_new_Shoepeg_3_FSTSum$y,
  z = factor(data_new_Shoepeg_3_FSTSum$z)
)
df_shoe4_FSTSum <- data.frame(
  y = factor(data_new_Shoepeg_4_FSTSum$x),
  x = data_new_Shoepeg_4_FSTSum$y,
  z = factor(data_new_Shoepeg_4_FSTSum$z)
)
#### Plot Shoepeg 1 
regions_of_chr_Shoepeg_1 <- ggplot(df_shoe1_FSTSum, aes(y = y, x = x, fill = z, colour = z)) + 
  geom_tile()+
  theme(text = element_text(size = font_size, family = "my"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = font_size, family = "my", colour = "black"),
        axis.title.x = element_blank())+
  labs(y = str_wrap("Individual samples in short 1",12), x = "SNP Positions",
       tag = "F")+
  scale_fill_manual(values = my_pal_col_1)+
  scale_color_manual(values = my_pal_col_1)
regions_of_chr_Shoepeg_1
#### Plot Shoepeg 2 
regions_of_chr_Shoepeg_2 <- ggplot(df_shoe2_FSTSum, aes(y = y, x = x, fill = z, colour = z)) + 
  geom_tile()+
  theme(text = element_text(size = font_size, family = "my"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = font_size, family = "my", colour = "black"),
        axis.title.x = element_text(size = font_size, family = "my", colour = "black"))+
  labs(y = str_wrap("Individual samples in tall 2",12), x = "SNP Positions")+
  scale_fill_manual(values = my_pal_col_1)+
  scale_color_manual(values = my_pal_col_1)
regions_of_chr_Shoepeg_2
#### Plot Shoepeg 3
regions_of_chr_Shoepeg_3 <- ggplot(df_shoe3_FSTSum, aes(y = y, x = x, fill = z, colour = z)) + 
  geom_tile()+
  theme(text = element_text(size = font_size, family = "my"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = font_size, family = "my", colour = "black"),
        axis.title.x = element_blank())+
  labs(y = str_wrap("Individual samples in tall 1",12), x = "SNP Positions")+
  scale_fill_manual(values = my_pal_col_1)+
  scale_color_manual(values = my_pal_col_1)
regions_of_chr_Shoepeg_3
#### Plot Shoepeg 4
regions_of_chr_Shoepeg_4 <- ggplot(df_shoe4_FSTSum, aes(y = y, x = x, fill = z, colour = z)) + 
  geom_tile()+
  theme(text = element_text(size = font_size, family = "my"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = font_size, family = "my", colour = "black"),
        axis.title.x = element_blank())+
  labs(y = str_wrap("Individual samples in short 2",12), x = "SNP Positions")+
  scale_fill_manual(values = my_pal_col_1)+
  scale_color_manual(values = my_pal_col_1)
regions_of_chr_Shoepeg_4
### Get the legend
regions_of_chr_legend <- ggplot(df_shoe4_FSTSum, aes(y = y, x = x, fill = z, colour = z)) + 
  geom_tile()+
  theme(text = element_text(size = font_size, family = "my"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        legend.title = element_blank(),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = font_size, family = "my", colour = "black", face = "bold"),
        axis.title.x = element_text(size = font_size, family = "my", colour = "black", face = "bold"))+
  labs(y = "Individual samples", x = "SNP Positions")+
  scale_fill_manual(values = my_pal_col_1)+
  scale_color_manual(values = my_pal_col_1)
get_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 
legend <- get_legend(regions_of_chr_legend)
regions_of_chr_all_Shoepeg_pop <- grid.arrange(regions_of_chr_Shoepeg_1,
                                               regions_of_chr_Shoepeg_4,
                                               regions_of_chr_Shoepeg_3,
                                               regions_of_chr_Shoepeg_2, 
                                               legend,
                                               nrow = 4,
                                               ncol = 2, 
                                               heights = c(2, 2, 2, 2.1),
                                               widths = c(3,0.5),
                                               layout_matrix = cbind(c(1,2,3,4),
                                                                     c(5)))
ggsave("C:/Users/mtost/Documents/Shoepeg_paper_final_versions/Figures/New_figures/Haplotypes_of_candidate_gene_region_all_Shoepeg_pop_all_GTs_iAA8_region_V2.png",
       regions_of_chr_all_Shoepeg_pop,
       height = 6, width = 4,
       dpi =1800)
#### Get the legend for the plot -----------------------------------------------
ChickWeight <- as.data.table(ChickWeight)
ChickWeight <- ChickWeight[Chick == 1,]
ChickWeight$Chick <- "WStat"
ChickWeight$Diet <- "FSTSum"
legend_plot <- ggplot()+
  geom_point(data = ChickWeight, aes(x = Time, y = weight, colour = Diet))+
  geom_bar(data = ChickWeight, stat="identity", position="dodge", aes(fill = Chick, y = weight,
                                                                      x = Time))+
  theme(text = element_text(size =font_size, family = "my"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.ticks = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size =font_size, colour = "black"),
        axis.title.x = element_text(size =font_size, colour = "black", face = "bold"),
        axis.text.x = element_text(size =font_size, colour = "black"),
        legend.key = element_blank(),
        legend.text = element_text(size =font_size, colour = "black"),
        legend.title = element_blank(),
        legend.position = "bottom")+
  scale_fill_manual(values = c("WStat" = "tomato2"))+
  scale_colour_manual(values = c("FSTSum" = "black"))
get_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
}
legend_candidate_gene_region <- get_legend(legend_plot)
