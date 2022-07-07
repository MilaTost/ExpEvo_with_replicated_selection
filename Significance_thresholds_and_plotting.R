### Load the packages and clear the global environment -------------------------
rm(list = ls())
library(data.table)
library(stringr)
library(ggplot2)
library(gridExtra)
library(foreach)
library(parallel)
library(tidyverse)
### Load the data --------------------------------------------------------------
FST_values_od_cor <- read.table("C:/Users/mtost/Documents/Masterarbeit/Data_analysis/Final_data_sets_for_paper/2022_02_16_GB1005_corrected_dt_Fst_values.txt")
# Generate the significance thresholds based on the empirical distribution ---
sig_thres_dis_999_perc <- quantile(FST_values_od_cor$Sum_FST, probs = 0.999, na.rm = TRUE)  
sig_thres_emp_dis_9999_perc <- quantile(FST_values_od_cor$Sum_FST, probs = 0.9999, na.rm = TRUE)
# Load the significance threshold based on the drift simulations ---------------
sig_thres_drift_sim <-  0.7620932
### The significance thresholds based on drift simulations were 
### calculated in the drift simulation scripts 
### See Drift_simulations.R

# Generate significance threshold based on the FDRfS ----------------------------------
### FDR for selection for all possible values of the statistics
# How many markers were observed at a certain FST value?
# We generate a table to show the different possible FST values
# difference at the number of observed markers at these.
## Calculate the FDR for selection based on the sum FST ------------------------
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
write.table(dt_FDR_for_selection_sum_FST,"C:/Users/mtost/Documents/Masterarbeit/GB10_output/GB1006_FDR_for_the_different_sum_FST.txt",sep = "  ", row.names = TRUE,
            quote = FALSE)

### Select a significance threshold and save it --------------------------------
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
sig_thres_FDRfS <- get_sig_thresh_FDR_for_selection_based_on_sumFST(stat_opposite_dir1 = FST_values_od_cor$FST_value_opposite_dir1,
                                                                          stat_opposite_dir2 = FST_values_od_cor$FST_value_opposite_dir2,
                                                                          stat_same_dir1 = FST_values_od_cor$FST_value_same_dir1,
                                                                          stat_same_dir2 = FST_values_od_cor$FST_value_same_dir2,
                                                                          statistic = "FST",
                                                                          FDR = 0.05)
# Test for selection: Sauron plot --------------------------------------------
create_sauron_plot_FST <- function(data_table,
                                   font_size,
                                   font_type,
                                   sig_threshold,
                                   stat_opposite_dir1,
                                   stat_opposite_dir2,
                                   stat_same_dir1,
                                   stat_same_dir2){
  windowsFonts(my = windowsFont(font_type))
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
    theme(text = element_text(size = font_size, family = "my"),
          panel.background = element_rect(fill = "white"),
          axis.line = element_line(colour = "black"),
          axis.title.y = element_text(size = font_size, family = "my", colour = "royalblue3", face = "bold"),
          axis.title.x = element_text(size = font_size, family = "my", colour = "royalblue3", face = "bold"),
          axis.line.x = element_line(color = "black"),
          axis.text = element_text(size = font_size-2, family = "my", colour = "black"))+
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
          axis.text = element_text(size = font_size-2, family = "my", colour = "black"))+
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
                                          font_size = 16,
                                          font_type = 'Calibri',
                                          sig_threshold = sig_thres_FDRfS,
                                          stat_opposite_dir1 = FST_values_od_cor$FST_value_opposite_dir1,
                                          stat_opposite_dir2 = FST_values_od_cor$FST_value_opposite_dir2,
                                          stat_same_dir1 = FST_values_od_cor$FST_value_same_dir1,
                                          stat_same_dir2 = FST_values_od_cor$FST_value_same_dir2)
Sauron_plot_FST
ggsave("C:/Users/mtost/Documents/Shoepeg_final_versions/Figures/Figure_test.png",
       Sauron_plot_FST,
       height = 6, width = 6,
       dpi =1200)
# Manhatten plot -------------------------------------------------------------
## Load the significance thresholds --------------------------------------------
sig_thres_dis_999_perc <- quantile(FST_values_od_cor$Sum_FST, probs = 0.999, na.rm = TRUE)  
sig_thres_emp_dis_9999_perc <- quantile(FST_values_od_cor$Sum_FST, probs = 0.9999, na.rm = TRUE)
sig_thres_drift_sim <-  0.7620932
sig_thres_FDRfS <- 0.662

## Prepare the data set, so that the chromosomes are correctly represented -----
library(stringr)
FST_values_od_cor$Chromosome <- str_sub(FST_values_od_cor$Chromosome,4,5)
FST_values_od_cor$Chromosome <- as.factor(FST_values_od_cor$Chromosome)
str(FST_values_od_cor)
FST_values_od_cor$Chromosome <- as.character(FST_values_od_cor$Chromosome)
### Plotting function Manhatten Fst plot ---------------------------------------
## Create Manhatten plots for presentations ------------------------------------
create_manhatten_plot_per_chr_one_stat <- function(data,
                                                   name_of_statistic,
                                                   font_size,
                                                   sig_thres_emp_dis_999,
                                                   sig_thres_emp_dis_9999,
                                                   sig_thres_drift_sim,
                                                   sig_thres_FDR){
  windowsFonts(my = windowsFont('Calibri'))
  my_pal_col <- (c("1" = "pink","2" = "pink3","3" = "pink",
                   "4" = "pink3","5" = "pink","6" = "pink3", 
                   "7" = "pink","8" = "pink3","9" = "pink",
                   "10" = "pink3"))
  Fst_tab_small <- data[1:10,]
  legend_chr_plot <- ggplot(data = Fst_tab_small)+
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
          legend.title = element_text(size =font_size-2, family = "my", colour = "black", face = "bold"))+
    labs(y = "FST",x = "Position")+
    geom_hline(aes(yintercept = sig_thres_drift_sim, linetype=str_wrap("1) Drift Sim",12)), colour="purple", size = 1)+
    geom_hline(aes(yintercept = sig_thres_FDR, linetype=str_wrap("2) FDRfS",12)), colour="gold", size = 1)+
    geom_hline(aes(yintercept = sig_thres_emp_dis_9999, linetype=str_wrap("3) 99.99th",12)), colour="darkred", size = 1)+
    geom_hline(aes(yintercept = sig_thres_emp_dis_999, linetype=str_wrap("4) 99.9th",12)), colour="tomato", size = 1)+
    scale_linetype_manual(name = str_wrap("Significance threshold",10), values = c(1,1,1,1),
                          guide = guide_legend(override.aes = list(color = c("purple","gold","darkred","tomato"),
                                                                   linetype = c("solid","solid","dotted","dotted"))))
  get_legend <- function(a.gplot){ 
    tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
    legend <- tmp$grobs[[leg]] 
    legend
  } 
  legend <- get_legend(legend_chr_plot)
  chr_plot <- ggplot(data =  data, aes(x =  fct_inorder(SNP_ID), y = Sum_FST))+
    geom_point(aes(x = fct_inorder(SNP_ID), y = Sum_FST, colour = Chromosome, group = Chromosome))+
    theme(text = element_text(size =font_size, family = "my"),
          panel.background = element_rect(fill = "white"),
          panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                          colour = "grey81"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                          colour = "grey100"),
          panel.border = element_blank(),
          axis.line.y = element_line(colour = "black"),
          axis.line.x = element_line(colour = "black"),
          axis.ticks = element_blank(),
          axis.title.y = element_text(size =font_size, family = "my", colour = "black", face = "bold"),
          axis.text.y = element_text(size =font_size-2, family = "my", colour = "black"),
          axis.title.x = element_text(size =font_size, family = "my", colour = "black", face = "bold"),
          axis.text.x = element_blank(),
          legend.text = element_blank(),
          legend.title = element_blank())+
    labs(y = "FstSum",x = "Position")+
    geom_hline(aes(yintercept = sig_thres_FDR), linetype="solid", colour="gold", size = 1)+
    geom_hline(aes(yintercept = sig_thres_drift_sim), linetype="solid", colour="purple", size = 1)+
    geom_hline(aes(yintercept = as.numeric(sig_thres_emp_dis_9999)), linetype="dotted", colour="darkred", size = 1)+
    geom_hline(aes(yintercept = as.numeric(sig_thres_emp_dis_999)), linetype="dotted", colour="tomato", size = 1)+
    scale_color_manual(values=my_pal_col)+
    guides(color = "none", linetype = "none")
  Manhatten_plot <- grid.arrange(chr_plot, legend,
                                 ncol=2, nrow = 1, 
                                 layout_matrix = rbind(c(1,2)),
                                 widths = c(17, 3), heights = c(4))
  return(Manhatten_plot)
}
FST_manhatten_plot <- create_manhatten_plot_per_chr_one_stat(data = FST_values_od_cor,
                                                             name_of_statistic = "Fst",
                                                             font_size = 16,
                                                             sig_thres_emp_dis_999 = sig_thres_emp_dis_9999_perc,
                                                             sig_thres_emp_dis_9999 = sig_thres_emp_dis_9999_perc,
                                                             sig_thres_drift_sim = sig_thres_drift_sim,
                                                             sig_thres_FDR = sig_thres_FDRfS)
ggsave("C:/Users/mtost/Documents/Shoepeg_final_versions/Figures/Figure_5.png",
       FST_manhatten_plot,
       height = 6, width = 10,
       dpi =2000)
# Plot candidate gene region -------------------------------------------------
FST_values_od_cor <- as.data.table(FST_values_od_cor)
Fst_Chr_3 <-  FST_values_od_cor[Chromosome=="3",]
ymax_Fst <- 1
sig_thres_FDRfS <- 0.662
## Plot candidate gene region for only my markers ------------------------------
windowsFonts(my = windowsFont("Calibri"))
font_size <- 16
## Load the data set with the candidate genes 
candidate_genes <- read.table("C:/Users/mtost/Documents/Masterarbeit/Data_analysis/Final_data_sets_for_paper/candidate_gene_region_and_markers.txt",
                              header = TRUE)
candidate_gene_region <- ggplot(data =  Fst_Chr_3, aes(x = Position, y = mean_FST_OD_NA_excluded))+
  geom_point(aes(x = Position, y = Sum_FST))+
  theme(text = element_text(size =font_size, family = "my"),
        panel.background = element_rect(fill = "white"),
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                        colour = "grey81"),
        panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                        colour = "grey100"),
        panel.border = element_blank(),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size =font_size, family = "my", colour = "black", face = "bold"),
        axis.text.y = element_text(size =font_size-2, family = "my", colour = "black"),
        axis.title.x = element_text(size =font_size, family = "my", colour = "black", face = "bold"),
        axis.text.x = element_text(size =font_size-2, family = "my", colour = "black"),
        legend.text = element_text(size =font_size-2, family = "my", colour = "black"),
        legend.title = element_text(size =font_size-2, family = "my", colour = "black", face = "bold"))+
  labs(y = "FstSum",x = "Position")+
  geom_hline(aes(yintercept = sig_thres_FDRfS), linetype="solid", colour="gold", size = 1)+
  geom_vline(aes(xintercept = 9953933), linetype="solid", colour="darkorange", size = 1)+
  annotate("text", x = 9953933, y = 0.75, label = "Marker found by Gyawali et al. 2019", colour = "darkorange", family = "my")+
  geom_vline(aes(xintercept = 10179485), linetype="solid", colour="darkmagenta", size = 1)+
  annotate("text", x = 10179485, y = 0.70, label = "Marker found by Peiffer et al. 2014", colour = "darkmagenta", family = "my")+
  annotate("rect", xmin=10062219, xmax=10065595, ymin = 0, ymax = ymax_Fst, fill = "red", alpha = 0.6)+
  annotate("text", x = 10062219, y = 0.8, label = "iAA8", colour = "tomato2", family = "my")+
  annotate("rect", xmin=10440993, xmax=10443340, ymin = 0, ymax = ymax_Fst, fill = "red", alpha = 0.6)+
  annotate("text", x = 10440993, y = 0.85, label = "Dwarf1", colour = "tomato2", family = "my")+
  xlim(min(candidate_genes$Start)-1000000,max(candidate_genes$End)+1000000)+
  ylim(0,ymax_Fst)
candidate_gene_region
ggsave("C:/Users/mtost/Documents/Shoepeg_final_versions/Figures/Figure_6_V2.png",
       candidate_gene_region,
       height = 4, width = 10,
       dpi =2400)
## How many markers exceeded the different significance thresholds ---------
sig_thres_dis_999_perc <- quantile(FST_values_od_cor$Sum_FST, probs = 0.999, na.rm = TRUE)  
sig_thres_emp_dis_9999_perc <- quantile(FST_values_od_cor$Sum_FST, probs = 0.9999, na.rm = TRUE)
sig_thres_drift_sim <-  0.7620932
## The significance thresholds based on drift simulations were 
## calculated in the drift simulation scripts
sig_thres_FDRfS <- 0.662
## The significance thresholds based on drift simulations were 
## calculated in the drift simulation scripts
##### Infos about the markers which passed this threshold ----------------------
markers_passed_treshold <- function(data, threshold, values_of_the_statistic, SNP_IDs){
  data <- as.data.frame(data)
  dt_markers <- data[which(values_of_the_statistic >= threshold),which(colnames(data)==SNP_IDs)]
  cat(length(dt_markers),"markers exceed the significance threshold based on this significance threshold.")
  return(dt_markers)
}
###### How many markers exceeded the different significance thresholds ---------
Nr_markers_passed_sig_thres_dis_999_perc <- markers_passed_treshold(data = FST_values_od_cor, 
                                                                    threshold = sig_thres_dis_999_perc,
                                                                    values_of_the_statistic = FST_values_od_cor$Sum_FST,
                                                                    SNP_IDs = "SNP_ID")
Nr_markers_passed_sig_thres_dis_9999_perc <- markers_passed_treshold(data = FST_values_od_cor, 
                                                                     threshold = sig_thres_dis_9999_perc,
                                                                     values_of_the_statistic = FST_values_od_cor$Sum_FST,
                                                                     SNP_IDs = "SNP_ID")
Nr_markers_passed_sig_thres_drift_sim <- markers_passed_treshold(data = FST_values_od_cor, 
                                                                 threshold = sig_thres_drift_sim,
                                                                 values_of_the_statistic = FST_values_od_cor$Sum_FST,
                                                                 SNP_IDs = "SNP_ID")
Nr_markers_passed_sig_thres_FDRfS <- markers_passed_treshold(data = FST_values_od_cor, 
                                                             threshold = sig_thres_FDRfS,
                                                             values_of_the_statistic = FST_values_od_cor$Sum_FST,
                                                             SNP_IDs = "SNP_ID")
### Create a data table with significant markers -------------------------------
Nr_markers_passed_sig_thres_FDRfS <- as.data.frame(Nr_markers_passed_sig_thres_FDRfS)
write.table(Nr_markers_passed_sig_thres_FDRfS,"C:/Users/mtost/Documents/Masterarbeit/Data_analysis/Final_data_sets_for_paper/2022_04_25_FST_FDRfS_significant_markers.txt",sep = "  ", row.names = TRUE,
            quote = FALSE)