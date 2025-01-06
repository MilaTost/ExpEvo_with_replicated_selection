rm(list=ls())
# Plot haplotypes --------------------------------------------------------------
library(compiler)
library(dplyr)
library(devtools)
library(data.table)
library(foreach)
library(stringr)
library(ggplot2)
library(forcats)
library(gridExtra)
#### Load new data of the 3 significant regions --------------------------------
# General data sets
setwd(input_dir)
real_subpop_labels <- readLines("new_labels.txt", warn = FALSE)
real_subpop_labels <- str_split(real_subpop_labels, " ")
real_subpop_labels <- unlist(real_subpop_labels)
subpop_names <- vector(length = length(real_subpop_labels))
subpop_names[str_which(real_subpop_labels, "1")] <- "Shoepeg_1"
subpop_names[str_which(real_subpop_labels, "2")] <- "Shoepeg_2"
subpop_names[str_which(real_subpop_labels, "3")] <- "Shoepeg_3"
subpop_names[str_which(real_subpop_labels, "4")] <- "Shoepeg_4"
subpop_names
# Region on chr3 ---------------------------------------------------------------
Haplotypes_chr3 <- readLines("Clusters_chr3_pos_9411000_to_10423000_hapguess_switch.out", warn = FALSE)
start_SNPs_WStat <- 9494350
end_SNPs_WStat <- 9529950
start_SNPs_WStat <- 10062550
end_SNPs_WStat <- 10063550
Haplotypes_colnames <- read.table("2024-04-24_marker_positions_random_region_9411000_to_10423000.txt",
                                  header = TRUE)
marker_positions <- Haplotypes_colnames$x
#### Prepare the data ----------------------------------------------------------
prepare_haplotype_sequences <- function(Haplotypes,
                                        marker_positions,
                                        subpop_labels){
  Haplotypes <- Haplotypes[22:length(Haplotypes)]
  labels_of_subpop <- seq(1, length(Haplotypes)-2, 3)
  Labels_of_subpop <- Haplotypes[labels_of_subpop]
  haplotype_1 <- seq(2, length(Haplotypes)-1, 3)
  Haplotype_Allele_1 <- Haplotypes[haplotype_1]
  haplotype_2 <- seq(3, length(Haplotypes), 3)
  Haplotype_Allele_2 <- Haplotypes[haplotype_2]
  prepare_haplotype_sequences <- function(Haplotypes){
    prepare_haplotype_sequences_split <- function(i){
      splitted_haplo <- str_split(Haplotypes[i], " ")
      unlisted_split_haplo <- unlist(splitted_haplo)
      return(unlisted_split_haplo)
    }
    dt_Haplotypes <- foreach(i = 1:length(Haplotypes), .combine = cbind) %do% prepare_haplotype_sequences_split(i)
    return(dt_Haplotypes)
  }
  dt_Haplotypes_A1 <- prepare_haplotype_sequences(Haplotypes = Haplotype_Allele_1)
  dt_Haplotypes_A2 <- prepare_haplotype_sequences(Haplotypes = Haplotype_Allele_2)
  diploid_seq <- paste0(dt_Haplotypes_A1, dt_Haplotypes_A2)
  dt_diploid_seq <- matrix(data = diploid_seq,
                           nrow = nrow(dt_Haplotypes_A1),
                           ncol = ncol(dt_Haplotypes_A1),
                           byrow = FALSE)
  dt_Haplotypes_A1 <- cbind(marker_positions, dt_Haplotypes_A1)
  dt_Haplotypes_A2 <- cbind(marker_positions, dt_Haplotypes_A2)
  colnames(dt_Haplotypes_A1) <- c("Marker", subpop_labels)
  colnames(dt_Haplotypes_A2) <- c("Marker", subpop_labels)
  dt_diploid_seq <- cbind(marker_positions, dt_diploid_seq)
  Haplotypes <- list("Allele_1" = dt_Haplotypes_A1,
                     "Allele_2" = dt_Haplotypes_A2,
                     "Subpopulation_labels" = Labels_of_subpop,
                     "DT_Diploid_sequences" = dt_diploid_seq)
  rm(dt_Haplotypes_A1, dt_Haplotypes_A2)
  return(Haplotypes)
}
new_Haplotypes <- prepare_haplotype_sequences(Haplotypes = Haplotypes_chr3,
                                              marker_positions = marker_positions,
                                              subpop_labels = subpop_names)
new_Haplotypes <- prepare_haplotype_sequences(Haplotypes = Haplotypes_chr7,
                                              marker_positions = marker_positions,
                                              subpop_labels = subpop_names)
dt_Haplotypes_A1 <- as.data.frame(new_Haplotypes$Allele_1)
dt_Haplotypes_A2 <- as.data.frame(new_Haplotypes$Allele_2)
### Prepare the data -----------------------------------------------------------
prepare_data <- function(data_1,
                         data_2,
                         allele_1,
                         allele_2){
  data_1 <- as.matrix(data_1)
  data_2 <- as.matrix(data_2)
  prep_data <- function(data){
    prepare_data_per_element <- function(i){
      t_col <- t(data[i,])
      new_col<- t_col[2:length(t_col)]
      new_new_col <- cbind(colnames(data)[2:ncol(data)], rep(data[i,1], length(new_col)), new_col)
      new_new_col <- as.data.frame(new_new_col)
      return(new_new_col)
    }
    new_data <- foreach(i = 1:nrow(data), .combine = rbind) %do% prepare_data_per_element(i)
    return(new_data)
  }
  new_data_1 <- prep_data(data = data_1)
  new_data_2 <- prep_data(data = data_2)
  # everything is okay here
  final_data_1 <- cbind(new_data_1, rep(allele_1, nrow(new_data_1)))
  final_data_2 <- cbind(new_data_2, rep(allele_2, nrow(new_data_2)))
  final_data_1 <- as.data.frame(final_data_1)
  final_data_2 <- as.data.frame(final_data_2)
  colnames(final_data_1) <- c("Individual", "Marker", "Ref_encoding", "Allele")
  colnames(final_data_2) <- c("Individual", "Marker", "Ref_encoding", "Allele")
  final_data <- rbind(final_data_1, final_data_2)
  colnames(final_data) <- c("Individual", "Marker", "Ref_encoding", "Allele")
  return(final_data)
}
dt_Haplotypes <- prepare_data(data_1 = dt_Haplotypes_A1,
                              data_2 = dt_Haplotypes_A2,
                              allele_1 = "A1",
                              allele_2 = "A2")
dt_Haplotypes$Marker <- as.factor(dt_Haplotypes$Marker)
dt_Haplotypes$Ind_Allele <- paste0(dt_Haplotypes$Individual, "_",dt_Haplotypes$Allele)
dt_Haplotypes$Population <- 0
index_col_pop <- which(colnames(dt_Haplotypes) == "Population")
dt_Haplotypes[which(dt_Haplotypes$Individual == "Shoepeg_1"), index_col_pop] <- "Short_1"
dt_Haplotypes[which(dt_Haplotypes$Individual == "Shoepeg_4"), index_col_pop] <- "Short_2"
dt_Haplotypes[which(dt_Haplotypes$Individual == "Shoepeg_3"), index_col_pop] <- "Tall_1"
dt_Haplotypes[which(dt_Haplotypes$Individual == "Shoepeg_2"), index_col_pop] <- "Tall_2"
### Load plotting parameters ---------------------------------------------------
##### Load colours
fill_A <- "midnightblue"
fill_T <- "tomato3"
fill_G <- "lightgoldenrod"
fill_T <- "cyan3"
fill_A <- "purple"
fill_T <- "yellow"
fill_G <- "orange"
fill_C <- "blue"
# New colours
#fill_A <- "red4"
#fill_T <- "springgreen"
#fill_G <- "skyblue4"
#fill_C <- "goldenrod1"
my_pal_fill <- c("A" = fill_A,
                 "T" = fill_T,
                 "G" = fill_G,
                 "C" = fill_C)
my_pal_colour <- c("A" = fill_A,
                   "T" = fill_T,
                   "G" = fill_G,
                   "C" = fill_C)
#### Prepare the data
dt_Haplotypes$Population <- as.factor(dt_Haplotypes$Population)
dt_Haplotypes <- dt_Haplotypes[order(dt_Haplotypes$Population),]
dt_Haplotypes$X <- rep(0, nrow(dt_Haplotypes))
dt_Haplotypes$Marker
levels(dt_Haplotypes$Population)
dt_Haplotypes_Short_1 <- dt_Haplotypes[dt_Haplotypes$Population == "Short_1",]
dt_Haplotypes_Tall_1 <- dt_Haplotypes[dt_Haplotypes$Population == "Tall_1",]
dt_Haplotypes_Short_2 <- dt_Haplotypes[dt_Haplotypes$Population == "Short_2",]
dt_Haplotypes_Tall_2 <- dt_Haplotypes[dt_Haplotypes$Population == "Tall_2",]
dt_Haplotypes_Short_1$Ref_encoding
### Reduce figure in size ------------------------------------------------------
vector_marker <- str_split(dt_Haplotypes$Marker, "_")
dt_marker <- matrix(data = unlist(vector_marker),
                    nrow = length(dt_Haplotypes$Marker),
                    byrow = TRUE)
bp_marker <- as.numeric(dt_marker[, 2])
#CHR3: 9911000	9923000 --------------------------------------------------------
#region_start <- 9911000
#region_end <- 9923000
#CHR3: 9494350	9529950 --------------------------------------------------------
region_start <- 9494350
region_end <- 9529950
#CHR3: 10062550	10063550 -------------------------------------------------------
region_start <- 10062550
region_end <- 10063550
unique(bp_marker[which(round((bp_marker)/1000000, 1) == round(region_start/1000000, 1))])
unique(bp_marker[which(round((bp_marker)/100000, 1) == round(region_end/100000, 1))])
10056542
10064651
#CHR3: 9494350	9529950
start_pos <- 9475947
end_pos <- 9529994
#CHR3: 10062550	10063550
start_pos <- 10062407
end_pos <- 10064547
#CHR7: 42687000	42725000 -------------------------------------------------------
region_start <- 42687000
region_end <- 42725000
region_start <- 42684450
region_end <- 42706850
unique(bp_marker[which(round((bp_marker)/1000000, 1) == round(region_start/1000000, 1))])
unique(bp_marker[which(round((bp_marker)/1000000, 1) == round(region_end/1000000, 1))])
start_pos <- 42684426
end_pos <-   42735280
start_pos <- 42684450
end_pos <- 42706850
### Prepare the data -----------------------------------------------------------
new_dt <- cbind(dt_Haplotypes[, 1:2], bp_marker, dt_Haplotypes[, 3:ncol(dt_Haplotypes)])
new_dt$bp_marker <- as.numeric(new_dt$bp_marker)
new_dt <- new_dt[order(new_dt$bp_marker),]
index_start_vec <- which(new_dt$bp_marker < start_pos)
index_start <- index_start_vec[length(index_start_vec)]
index_end_vec <- which(new_dt$bp_marker > end_pos)
index_end <- index_end_vec[1]
#new_dt <- new_dt[which(new_dt$bp_marker == end_pos)[length(which(new_dt$bp_marker == end_pos ))],]
dt_Haplotypes_Short_1 <- new_dt[which(new_dt$Population == "Short_1"),]
dt_Haplotypes_Tall_1 <- new_dt[which(new_dt$Population == "Tall_1"),]
dt_Haplotypes_Short_2 <- new_dt[which(new_dt$Population == "Short_2"),]
dt_Haplotypes_Tall_2 <- new_dt[which(new_dt$Population == "Tall_2"),]
nrow(dt_Haplotypes_Short_1)
nrow(dt_Haplotypes_Short_2)
nrow(dt_Haplotypes_Tall_1)
nrow(dt_Haplotypes_Tall_2)
### Plotting -------------------------------------------------------------------
##### Load other parameters
windowsFonts(my = windowsFont('Calibri'))
font_size <- 14
panel_background_colour <- "white"
plot_font_colour <- "black"
plot_background_colour <- "white"
theme_Mila <- theme(panel.background = element_rect(fill = panel_background_colour, colour = plot_background_colour),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.border =  element_rect(fill = NA, colour = panel_background_colour),
                    plot.background = element_rect(fill = plot_background_colour, colour = plot_background_colour),
                    axis.line = element_line(colour = plot_font_colour),
                    axis.line.x = element_line(color = plot_font_colour),
                    axis.ticks.y = element_blank(),
                    axis.ticks.x = element_blank(),
                    legend.position = "none",
                    legend.title = element_text(colour = plot_font_colour, size = font_size, family = "my", face = "bold"),
                    legend.text = element_text(colour = plot_font_colour, size = font_size, family = "my"),
                    plot.title = element_text(colour = plot_font_colour, size = font_size, family = "my"),
                    #axis.title.x = element_text(colour = plot_font_colour, size = font_size, family = "my", face = "bold"),
                    #axis.text.x = element_text(vjust = 0.5, size = font_size-2, family = "my", colour = plot_font_colour, angle = 45),
                    axis.text.x = element_blank(),
                    axis.text.y = element_blank(),
                    axis.title.x = element_blank(),
                    #axis.text.y = element_text(size = font_size-10, family = "my", colour = plot_font_colour),
                    axis.title.y = element_text(colour = plot_font_colour, size = font_size, family = "my", face = "bold"),
                    plot.tag = element_blank())
HAP_plot_Short1 <- ggplot(dt_Haplotypes_Short_1) +
  geom_tile(aes(x = fct_inorder(Marker), y = fct_inorder(Ind_Allele),
                colour = Ref_encoding, fill = Ref_encoding), lwd = 2)+
  theme_Mila +
  scale_fill_manual(values = my_pal_fill)+
  scale_colour_manual(values = my_pal_colour)+
  labs(y = "Short 1", x = "Marker position")
HAP_plot_Short1
HAP_plot_Short2 <- ggplot(dt_Haplotypes_Short_2) +
  geom_tile(aes(x = fct_inorder(Marker), y = Ind_Allele,
                colour = Ref_encoding, fill = Ref_encoding), lwd = 2)+
  theme_Mila +
  scale_fill_manual(values = my_pal_fill)+
  scale_colour_manual(values = my_pal_colour)+
  labs(y = "Short 2", x = "Marker position")
HAP_plot_Tall1 <- ggplot(dt_Haplotypes_Tall_1) +
  geom_tile(aes(x = fct_inorder(Marker), y = fct_inorder(Ind_Allele),
                colour = Ref_encoding, fill = Ref_encoding), lwd = 2)+
  theme_Mila +
  scale_fill_manual(values = my_pal_fill)+
  scale_colour_manual(values = my_pal_colour)+
  labs(y = "Tall 1", x = "Marker position")
theme_Mila <- theme(panel.background = element_rect(fill = panel_background_colour, colour = plot_background_colour),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.border =  element_rect(fill = NA, colour = panel_background_colour),
                    plot.background = element_rect(fill = plot_background_colour, colour = plot_background_colour),
                    axis.line = element_line(colour = plot_font_colour),
                    axis.line.x = element_line(color = plot_font_colour),
                    axis.ticks.y = element_blank(),
                    axis.ticks.x = element_blank(),
                    legend.position = "none",
                    legend.title = element_text(colour = plot_font_colour, size = font_size, family = "my", face = "bold"),
                    legend.text = element_text(colour = plot_font_colour, size = font_size, family = "my"),
                    plot.title = element_text(colour = plot_font_colour, size = font_size, family = "my"),
                    axis.title.x = element_text(colour = plot_font_colour, size = font_size, family = "my", face = "bold"),
                    #axis.text.x = element_text(vjust = 0.5, size = font_size-2, family = "my", colour = plot_font_colour, angle = 45),
                    axis.text.x = element_blank(),
                    axis.text.y = element_blank(),
                    #axis.text.y = element_text(size = font_size-10, family = "my", colour = plot_font_colour),
                    axis.title.y = element_text(colour = plot_font_colour, size = font_size, family = "my", face = "bold"),
                    plot.tag = element_blank())
HAP_plot_Tall2 <- ggplot(dt_Haplotypes_Tall_2) +
  geom_tile(aes(x = fct_inorder(Marker), y = fct_inorder(Ind_Allele),
                colour = Ref_encoding, fill = Ref_encoding), lwd = 2)+
  theme_Mila +
  scale_fill_manual(values = my_pal_fill)+
  scale_colour_manual(values = my_pal_colour)+
  labs(y = "Tall 2", x = "Marker position")
# Extract the legend
theme_Mila <- theme(panel.background = element_rect(fill = panel_background_colour, colour = plot_background_colour),
                    panel.grid.minor = element_blank(),
                    panel.grid.major = element_blank(),
                    panel.border =  element_rect(fill = NA, colour = panel_background_colour),
                    plot.background = element_rect(fill = plot_background_colour, colour = plot_background_colour),
                    axis.line = element_line(colour = plot_font_colour),
                    axis.line.x = element_line(color = plot_font_colour),
                    axis.ticks.y = element_blank(),
                    axis.ticks.x = element_blank(),
                    legend.position = "bottom",
                    legend.title = element_blank(),
                    legend.text = element_text(colour = plot_font_colour, size = font_size, family = "my"),
                    plot.title = element_text(colour = plot_font_colour, size = font_size, family = "my"),
                    axis.title.x = element_text(colour = plot_font_colour, size = font_size, family = "my", face = "bold"),
                    #axis.text.x = element_text(vjust = 0.5, size = font_size-2, family = "my", colour = plot_font_colour, angle = 45),
                    axis.text.x = element_blank(),
                    axis.text.y = element_blank(),
                    #axis.text.y = element_text(size = font_size-10, family = "my", colour = plot_font_colour),
                    axis.title.y = element_text(colour = plot_font_colour, size = font_size, family = "my", face = "bold"),
                    plot.tag = element_blank())
HAP_plot_legend <- ggplot(dt_Haplotypes_Short_1) +
  geom_tile(aes(x = fct_inorder(Marker), y = Ind_Allele,
                colour = Ref_encoding, fill = Ref_encoding), lwd = 2)+
  theme_Mila +
  scale_fill_manual(values = my_pal_fill)+
  scale_colour_manual(values = my_pal_colour)+
  labs(y = "Short 1", x = "Position")
get_legend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
legend_HAP <- get_legend(HAP_plot_legend)
# Combine the plot
library(gridExtra)
HAP_plots <- grid.arrange(HAP_plot_Short1,
                           HAP_plot_Short2,
                           HAP_plot_Tall1,
                           HAP_plot_Tall2,
                           legend_HAP,
                           ncol = 1, nrow = 5,
                           layout_matrix = rbind(1,2,3,4,5),
                           widths = 10, heights = c(4,4,4,5,0.5))
# Save the plot
result_dir <- "C:/Users/mtost/Documents/2024_Shoepeg_paper_V4/New figures/"
ggsave(paste0(result_dir, Sys.Date(),"_Haplotypes_chr3_9.4_to_to_10.0635Mb.png"),
       HAP_plots,
       device = png,
       height = 6,
       width = 12,
       dpi = 900)
ggsave(paste0(result_dir, Sys.Date(),"_Haplotypes_chr3_10.0625_to_10.0635Mb.png"),
       HAP_plots,
       device = png,
       height = 6,
       width = 12,
       dpi = 900)
ggsave(paste0(result_dir, Sys.Date(),"_Haplotypes_chr7_42187000_to_43225000.png"),
       HAP_plots,
       device = png,
       height = 6,
       width = 12,
       dpi = 900)
ggsave(paste0(result_dir, Sys.Date(),"_Haplotype_plot_chr5_only_WStat_region.png"),
       HAP_plots,
       device = png,
       height = 6,
       width = 12,
       dpi = 900)
ggsave(paste0(result_dir, Sys.Date(),"_NEW_Haplotype_plot_chr5.png"),
       HAP_plots,
       device = png,
       height = 6,
       width = 12,
       dpi = 900)
ggsave(paste0(result_dir, Sys.Date(),"NEW_Haplotype_plot_chr7.png"),
       HAP_plots,
       device = png,
       height = 6,
       width = 12,
       dpi = 900)
ggsave(paste0(result_dir, Sys.Date(),"_Haplotype_plot_chr7_only_WStat_region.png"),
       HAP_plots,
       device = png,
       height = 6,
       width = 12,
       dpi = 900)
