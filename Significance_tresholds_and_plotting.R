### Load the packages and clear the global environment -------------------------
rm(list = ls())
library(data.table)
library(stringr)
library(ggplot2)
library(gridExtra)
library(foreach)
### Load the data --------------------------------------------------------------
allele_freq_diff <- read.table("C:/Users/mtost/Documents/Masterarbeit/Data_analysis/Final_data_sets_for_paper/2022_01_20_GB1006_Allele_freq_diff.txt")

FST_values_od_cor <- read.table("C:/Users/mtost/Documents/Masterarbeit/GB10_output/GB1007_Fst_values.txt")
head(FST_values_od_cor)
head(allele_freq_diff)
### Significance thresholds based on the empirical distribution ----------------
FST_emp_dis_999_perc <- quantile(FST_values_od_cor$FST_value, probs = 0.999, na.rm = TRUE)  
FST_emp_dis_9999_perc <- quantile(FST_values_od_cor$FST_value, probs = 0.9999, na.rm = TRUE)
AFD_emp_dis_999_perc <- quantile(allele_freq_diff$Abs_allele_freq_diff, probs = 0.999, na.rm = TRUE)
AFD_emp_dis_9999_perc <- quantile(allele_freq_diff$Abs_allele_freq_diff, probs = 0.9999, na.rm = TRUE)

##### Infos about the markers which passed this threshold ----------------------
markers_passed_treshold <- function(data, threshold, values_of_the_statistic, SNP_IDs){
  data <- as.data.frame(data)
  dt_markers <- data[which(values_of_the_statistic > threshold),which(colnames(data)==SNP_IDs)]
  cat(length(dt_markers),"markers exceed the significance threshold based on this significance threshold.")
  return(dt_markers)
}
FST_emp_dis_999_perc_sig_markers <- markers_passed_treshold(data = FST_values_od_cor, 
                                                            threshold = FST_emp_dis_999_perc,
                                                            values_of_the_statistic = FST_values_od_cor$Fst_value,
                                                            SNP_IDs = "SNP_ID")
FST_emp_dis_9999_perc_sig_markers <- markers_passed_treshold(data = FST_values_od_cor, 
                                                            threshold = FST_emp_dis_9999_perc,
                                                            values_of_the_statistic = FST_values_od_cor$Fst_value,
                                                            SNP_IDs = "SNP_ID")
AFD_emp_dis_999_perc_sig_markers <- markers_passed_treshold(data = allele_freq_diff, 
                                                            threshold = AFD_emp_dis_999_perc,
                                                            values_of_the_statistic = allele_freq_diff$Abs_allele_freq_diff,
                                                            SNP_IDs = "SNP_ID")
AFD_emp_dis_9999_perc_sig_markers <- markers_passed_treshold(data = allele_freq_diff, 
                                                            threshold = AFD_emp_dis_9999_perc,
                                                            values_of_the_statistic = allele_freq_diff$Abs_allele_freq_diff,
                                                            SNP_IDs = "SNP_ID")
### Significance thresholds based on drift simulations -------------------------
AFD_drift_sim_sig_thres <- 0.4686
FST_drift_sim_sig_thres <-  0.3540692
## The significance thresholds based on drift simulations were 
## calculated in the drift simulation scripts

### Significance thresholds based on the FDR for selection ---------------------
## FDR for selection for all possible values of the statistics -----------------
# How many markers were observed at a certain absolute 
# allele frequency difference.
# We generate a table to show the different possible absolute allele frequency 
# difference at the number of observed markers at these.
calculate_FDR_for_selection <- function(stat_opposite_dir1,
                                        stat_opposite_dir2,
                                        stat_same_dir1,
                                        stat_same_dir2,
                                        statistic){
  stat_opposite_dir1 <- as.numeric(stat_opposite_dir1)
  stat_opposite_dir2 <- as.numeric(stat_opposite_dir2)
  stat_same_dir1 <- as.numeric(stat_same_dir1)
  stat_same_dir2 <- as.numeric(stat_same_dir2)
  if(statistic == "abs_AFD"){
    dt <- matrix(data = seq(0,2,0.01), nrow = length(seq(0,2,0.01)), ncol = 3)
    for (i in 1:nrow(dt)) {
      dt[i,2] <- length(which(stat_same_dir1 + stat_same_dir2>=dt[i,1]))
      dt[i,3] <- length(which(stat_opposite_dir1 + stat_opposite_dir2>=dt[i,1]))
    }
  }
  if(statistic == "AFD"){
    dt_pos <- matrix(data = seq(0,2,0.01), nrow = length(seq(0,2,0.01)), ncol = 3)
    for (i in 1:nrow(dt_pos)) {
      dt_pos[i,2] <- length(which(stat_same_dir1 + stat_same_dir2>=dt_pos[i,1]))
      dt_pos[i,3] <- length(which(stat_opposite_dir1 + stat_opposite_dir2>=dt_pos[i,1]))
    }
    dt_neg <- matrix(data = seq(-2,0,0.01), length(seq(0,2,0.01)), ncol = 3)
    for (i in 1:nrow(dt_neg)) {
      dt_neg[i,2] <- length(which(stat_same_dir1 + stat_same_dir2<=dt_neg[i,1]))
      dt_neg[i,3] <- length(which(stat_opposite_dir1 + stat_opposite_dir2<=dt_neg[i,1]))
    }
    dt <- rbind(dt_neg,dt_pos)
    }
  if(statistic == "FST"){
    dt <- matrix(data = seq(0,1,0.01), nrow = length(seq(0,1,0.01)), ncol = 3)
    for (i in 1:nrow(dt)) {
      dt[i,2] <- length(which(abs(stat_same_dir1 + stat_same_dir2)>=dt[i,1]))
      dt[i,3] <- length(which(abs(stat_opposite_dir1 + stat_opposite_dir2)>=dt[i,1]))
    }
  }
  dt <- as.data.frame(dt)
  colnames(dt) <- c("Statistic","Markers diverged same dir","Markers diverged opposite dir")
  dt$FDR_for_selection <- dt[,2]/dt[,3]
  if(statistic == "FST" | statistic == "abs_AFD"){
    dt <- dt[!is.infinite(dt[,3]),]
    dt <- dt[!is.na(dt[,3]),]
    sig_threshold <- min(dt[which(dt$FDR_for_selection < 0.1),1])
    cat("The significance threshold based on the FDR for selection < 10%","\n",
        "corresponds to an statistic of",sig_threshold, "\n")
  }
  if(statistic == "AFD"){
    dt <- dt[!is.infinite(dt[,3]),]
    dt <- dt[!is.na(dt[,3]),]
    neg_area <- dt[which(dt$Statistic < 0),]
    pos_area <- dt[which(dt$Statistic > 0),]
    neg_sig_threshold <- max(neg_area[which(neg_area$FDR_for_selection < 0.1),1])
    pos_sig_threshold <- min(pos_area[which(pos_area$FDR_for_selection < 0.1),1])
    cat("The significance threshold based on the FDR for selection < 10%","\n",
        "corresponds to an allele frequency difference of",
        neg_sig_threshold,"and",pos_sig_threshold,"\n")
  }
  return(dt)
}
dt_FDR_for_selection_AFD <- calculate_FDR_for_selection(stat_opposite_dir1 = allele_freq_diff$Low1_vs_High1,
                                                        stat_opposite_dir2 = allele_freq_diff$Low2_vs_High2,
                                                        stat_same_dir1 = allele_freq_diff$Low1_vs_Low2,
                                                        stat_same_dir2 = allele_freq_diff$High1_vs_High2,
                                                        statistic = "AFD")
write.table(dt_FDR_for_selection_AFD,"C:/Users/mtost/Documents/Masterarbeit/GB10_output/GB1006.1_FDR_for_the_different_AFD.txt",sep = "  ", row.names = TRUE,
            quote = FALSE)
dt_FDR_for_selection_abs_AFD <- calculate_FDR_for_selection(stat_opposite_dir1 = allele_freq_diff$Low1_vs_High1,
                                                        stat_opposite_dir2 = allele_freq_diff$Low2_vs_High2,
                                                        stat_same_dir1 = allele_freq_diff$Low1_vs_Low2,
                                                        stat_same_dir2 = allele_freq_diff$High1_vs_High2,
                                                        statistic = "abs_AFD")
write.table(dt_FDR_for_selection_abs_AFD,"C:/Users/mtost/Documents/Masterarbeit/GB10_output/GB1006.2_FDR_for_the_different_abs_AFD.txt",sep = "  ", row.names = TRUE,
            quote = FALSE)
dt_FDR_for_selection_FST <- calculate_FDR_for_selection(stat_opposite_dir1 = FST_values_od_cor$FST_value_opposite_dir1,
                                                        stat_opposite_dir2 = FST_values_od_cor$FST_value_opposite_dir2,
                                                        stat_same_dir1 = FST_values_od_cor$FST_value_same_dir1,
                                                        stat_same_dir2 = FST_values_od_cor$FST_value_same_dir2,
                                                        statistic = "FST")
write.table(dt_FDR_for_selection_FST,"C:/Users/mtost/Documents/Masterarbeit/GB10_output/GB1006.3_FDR_for_the_different_abs_FST.txt",sep = "  ", row.names = TRUE,
            quote = FALSE)
### Select a significance threshold
get_FDR_for_selection_sign_thres <- function(stat_opposite_dir1,
                                             stat_opposite_dir2,
                                             stat_same_dir1,
                                             stat_same_dir2,
                                             statistic){
  stat_opposite_dir1 <- as.numeric(stat_opposite_dir1)
  stat_opposite_dir2 <- as.numeric(stat_opposite_dir2)
  stat_same_dir1 <- as.numeric(stat_same_dir1)
  stat_same_dir2 <- as.numeric(stat_same_dir2)
  if(statistic == "abs_AFD"){
    dt <- matrix(data = seq(0,2,0.01), nrow = length(seq(0,2,0.01)), ncol = 3)
    for (i in 1:nrow(dt)) {
      dt[i,2] <- length(which(stat_same_dir1 + stat_same_dir2>=dt[i,1]))
      dt[i,3] <- length(which(stat_opposite_dir1 + stat_opposite_dir2>=dt[i,1]))
    }
  }
  if(statistic == "AFD"){
    dt_pos <- matrix(data = seq(0,2,0.01), nrow = length(seq(0,2,0.01)), ncol = 3)
    for (i in 1:nrow(dt_pos)) {
      dt_pos[i,2] <- length(which(stat_same_dir1 + stat_same_dir2>=dt_pos[i,1]))
      dt_pos[i,3] <- length(which(stat_opposite_dir1 + stat_opposite_dir2>=dt_pos[i,1]))
    }
    dt_neg <- matrix(data = seq(-2,0,0.01), length(seq(0,2,0.01)), ncol = 3)
    for (i in 1:nrow(dt_neg)) {
      dt_neg[i,2] <- length(which(stat_same_dir1 + stat_same_dir2<=dt_neg[i,1]))
      dt_neg[i,3] <- length(which(stat_opposite_dir1 + stat_opposite_dir2<=dt_neg[i,1]))
    }
    dt <- rbind(dt_neg,dt_pos)
  }
  if(statistic == "FST"){
    dt <- matrix(data = seq(0,1,0.01), nrow = length(seq(0,1,0.01)), ncol = 3)
    for (i in 1:nrow(dt)) {
      dt[i,2] <- length(which(abs(stat_same_dir1 + stat_same_dir2)>=dt[i,1]))
      dt[i,3] <- length(which(abs(stat_opposite_dir1 + stat_opposite_dir2)>=dt[i,1]))
    }
  }
  dt <- as.data.frame(dt)
  colnames(dt) <- c("Statistic","Markers diverged same dir","Markers diverged opposite dir")
  dt$FDR_for_selection <- dt[,2]/dt[,3]
  if(statistic == "FST" | statistic == "abs_AFD"){
    dt <- dt[!is.infinite(dt[,3]),]
    dt <- dt[!is.na(dt[,3]),]
    sig_threshold <- max(dt[which(dt$FDR_for_selection < 0.1),1])
    cat("The significance threshold based on the FDR for selection < 10%","\n",
        "corresponds to an statistic of",sig_threshold, "\n")
  }
  if(statistic == "AFD"){
    dt <- dt[!is.infinite(dt[,3]),]
    dt <- dt[!is.na(dt[,3]),]
    neg_area <- dt[which(dt$Statistic < 0),]
    pos_area <- dt[which(dt$Statistic > 0),]
    neg_sig_threshold <- max(neg_area[which(neg_area$FDR_for_selection < 0.1),1])
    pos_sig_threshold <- min(pos_area[which(pos_area$FDR_for_selection < 0.1),1])
    cat("The significance threshold based on the FDR for selection < 10%","\n",
        "corresponds to an allele frequency difference of",
        neg_sig_threshold,"and",pos_sig_threshold,"\n")
    sig_threshold <- c(neg_sig_threshold,pos_sig_threshold)
  }
  return(sig_threshold)
}
FST_FDR_for_sel_sig_thres <- get_FDR_for_selection_sign_thres(stat_opposite_dir1 = FST_values_od_cor$FST_value_opposite_dir1,
                                                              stat_opposite_dir2 = FST_values_od_cor$FST_value_opposite_dir2,
                                                              stat_same_dir1 = FST_values_od_cor$FST_value_same_dir1,
                                                              stat_same_dir2 = FST_values_od_cor$FST_value_same_dir2,
                                                              statistic = "FST")
AFD_FDR_for_sel_sig_thres <- get_FDR_for_selection_sign_thres(stat_opposite_dir1 = allele_freq_diff$Low1_vs_High1,
                                                              stat_opposite_dir2 = allele_freq_diff$Low2_vs_High2,
                                                              stat_same_dir1 = allele_freq_diff$Low1_vs_Low2,
                                                              stat_same_dir2 = allele_freq_diff$High1_vs_High2,
                                                              statistic = "AFD")
abs_AFD_FDR_for_sel_sig_thres <- get_FDR_for_selection_sign_thres(stat_opposite_dir1 = allele_freq_diff$Low1_vs_High1,
                                                                  stat_opposite_dir2 = allele_freq_diff$Low2_vs_High2,
                                                                  stat_same_dir1 = allele_freq_diff$Low1_vs_Low2,
                                                                  stat_same_dir2 = allele_freq_diff$High1_vs_High2,
                                                              statistic = "abs_AFD")
### Sauron plot ----------------------------------------------------------------
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
Sauron_plot_AFD <- create_sauron_plot_allele_freq_diff(data_table = allele_freq_diff,
                                                       sig_threshold = AFD_FDR_for_sel_sig_thres,
                                                       stat_opposite_dir1 = allele_freq_diff$Low1_vs_High1,
                                                       stat_opposite_dir2 = allele_freq_diff$Low2_vs_High2,
                                                       stat_same_dir1 = allele_freq_diff$Low1_vs_High1,
                                                       stat_same_dir2 = allele_freq_diff$Low1_vs_High2)
### Create the Sauron plot for the FST statistic -------------------------------
create_sauron_plot_FST <- function(data_table,
                                   sig_threshold,
                                   stat_opposite_dir1,
                                   stat_opposite_dir2,
                                   stat_same_dir1,
                                   stat_same_dir2){
  if(.Platform$OS.type == "windows") {
    windowsFonts(my = windowsFont('Calibri'))}
  else {my <- "Helvetica-Narrow"}
  create_FDR_area <- function(x){
    sig_threshold+-1*x
  }
  y_coord_value <- create_FDR_area(x = seq(-1,0,0.01))
  dt_coords_area <- cbind(y_coord_value,seq(-1,0,0.01))
  dt_coords_area <- as.data.frame(dt_coords_area)
  dt_coords_area$z <- rep(sig_threshold,nrow(dt_coords_area))
  colnames(dt_coords_area) <- c("x_coord_value","y_coord_value","threshold_value")
  dt_coords_area <- dt_coords_area[which(dt_coords_area$y_coord_value > sig_threshold),]
  sauron_1 <- ggplot()+
    geom_area(data = dt_coords_area, aes(x=x_coord_value, y=threshold_value), alpha = 0.1, fill = "tomato")+
    geom_area(data = dt_coords_area, aes(x=x_coord_value, y=y_coord_value), fill = "white")+
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
    ylim(0,1)+
    xlim(0,1)
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
    geom_area(data = dt_coords_area, aes(x=x_coord_value, y=threshold_value), alpha = 0.1, fill = "tomato")+
    geom_area(data = dt_coords_area, aes(x=x_coord_value, y=y_coord_value), fill = "white")+
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
    xlim(0,sig_threshold)+
    ylim(0,sig_threshold)
  sauron_2
  sauron_plot_all_FST<- grid.arrange(yaxis_sauron,sauron_2, xaxis_sauron,empty_plot,
                                      ncol=2, nrow = 2, 
                                      layout_matrix = rbind(c(1,2),c(4,3)),
                                      widths = c(0.2, 4), heights = c(4,0.2))
  return(sauron_plot_all_FST)
}

Sauron_plot_FST <- create_sauron_plot_FST(data_table = allele_freq_diff,
                                          sig_threshold = AFD_FDR_for_sel_sig_thres,
                                          stat_opposite_dir1 = allele_freq_diff$Short1_vs_Tall1,
                                          stat_opposite_dir2 = allele_freq_diff$Short2_vs_Tall2,
                                          stat_same_dir1 = allele_freq_diff$Short1_vs_Short2,
                                          stat_same_dir2 = allele_freq_diff$Tall1_vs_Tall2)

Sauron_plot_AFD
# Combine Sauron plots ---------------------------------------------------------
combined_sauron_plots <- grid.arrange(Sauron_plot_AFD, Sauron_plot_FST, 
                                    ncol=2, nrow = 1,
                                    labels = c("A","B"))
## this needs to be tested
### Manhatten plot -------------------------------------------------------------
# create Manhatten plots per Chromosome
# Use significance thresholds
# create a legend for the different significance thresholds
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

allele_freq_diff <- as.data.table(allele_freq_diff)
AFD_Chr_1 <- allele_freq_diff[Chromosome=="chr1",]
AFD_Chr_2 <- allele_freq_diff[Chromosome=="chr2",]
AFD_Chr_3 <- allele_freq_diff[Chromosome=="chr3",]
AFD_Chr_4 <- allele_freq_diff[Chromosome=="chr4",]
AFD_Chr_5 <- allele_freq_diff[Chromosome=="chr5",]
AFD_Chr_6 <- allele_freq_diff[Chromosome=="chr6",]
AFD_Chr_7 <- allele_freq_diff[Chromosome=="chr7",]
AFD_Chr_8 <- allele_freq_diff[Chromosome=="chr8",]
AFD_Chr_9 <- allele_freq_diff[Chromosome=="chr9",]
AFD_Chr_10 <- allele_freq_diff[Chromosome=="chr10",]

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
