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
setwd("/YOUR/OWN/PATH/")
WStat_OD_M3 <- fread("2022-12-19_Method_3_window_boundaries_based_on_FSTSum_OD_NA_excluded.txt")
WStat_SD_M3 <- fread("2022-12-19_Method_3_window_boundaries_based_on_FSTSum_SD_NA_excluded.txt")
dt_FDRfS_M3 <- fread("2022-12-19_Method_3_dt_FDRfS_WStat_calc_on_FSTSum.txt")

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
