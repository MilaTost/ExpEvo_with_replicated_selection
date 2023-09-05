# Create Haplotype-LD-decay figure for 9.2 Mb region ---------------------------
#### Loading of the packages ---------------------------------------------------
rm(list=ls())
cat("The packages are getting loaded.","\n")
library(compiler)
library(devtools)
library(data.table)
library(foreach)
library(stringr)
library(RandomFieldsUtils)
library(ggplot2)
library(HaploBlocker)
library(forcats)
library(gridExtra)
#### Load the plotting parameters ----------------------------------------------
windowsFonts(my = windowsFont('Calibri'))
font_size <- 18
# Reconstruction of Haplotype blocks -------------------------------------------
#### Load haplotypes -----------------------------------------------------------
##### Region 1 -----------------------------------------------------------------
setwd("C:/Users/mtost/Documents/Shoepeg_paper_final_versions/Scripts/Results/")
Haplotypes_colnames_reg1 <- read.table("2023-04-04_marker_positions_random_region_chr1 299047470 300064986.txt",
                                       header = TRUE)
chromosome_nr <- "chr1"
marker_positions_reg1 <- Haplotypes_colnames_reg1$x
SNP_IDs_reg1 <- marker_positions_reg1
tail(marker_positions_reg1)
Haplotypes_reg1 <- readLines("2023_04_03_Results_of_random_region_chr1_pos_299047470_300064986_hapguess_switch.out", warn = FALSE)
start_SNPs <- marker_positions_reg1[1]
end_SNPs <- marker_positions_reg1[length(marker_positions_reg1)]
(end_SNPs - start_SNPs)/1000000
length(marker_positions_reg1)
##### Region 2 -----------------------------------------------------------------
setwd("C:/Users/mtost/Documents/Shoepeg_paper_final_versions/Scripts/Results/")
Haplotypes_colnames_reg2 <- read.table("2023-04-04_marker_positions_random_region_chr7 177034868 178052384.txt",
                                       header = TRUE)
chromosome_nr <- "chr7"
marker_positions_reg2 <- Haplotypes_colnames_reg2$x
SNP_IDs_reg2 <- marker_positions_reg2
Haplotypes_reg2 <- readLines("2023_04_03_Results_of_random_region_chr7_pos_177034868_178052384_hapguess_switch.out", warn = FALSE)
start_SNPs <- marker_positions_reg2[1]
end_SNPs <- marker_positions_reg2[length(marker_positions_reg2)]
(end_SNPs - start_SNPs)/1000000
length(marker_positions_reg2)
dt_Haplotypes_colnames_reg <- read.table("2022-12-30_haplotypes_of_CHR3_candidate_gene_regions.txt",
                                         nrows = 1,
                                         header = TRUE)
Haplotypes_colnames_reg <- colnames(dt_Haplotypes_colnames_reg)
##### Candidate gene region Gen0 -----------------------------------------------
setwd("C:/Users/mtost/Documents/Shoepeg_paper_final_versions/Scripts/Results/")
Haplotypes_colnames_Gen0 <- read.table("2023-04-04_marker_positions_cand_gene_reg_Gen0.txt",
                                       header = TRUE)
marker_positions_Gen0 <- Haplotypes_colnames_Gen0$x
SNP_IDs_Gen0 <- marker_positions_Gen0
Haplotypes_Gen0 <- readLines("2023-04-03_Results_Shoepeg_cand_gene_region_Gen0_pos_9435824_10453340_hapguess_switch.out", warn = FALSE)
dt_Haplotypes_colnames_Gen0 <- fread("colnames_Shoepeg_Gen0.txt",
                                     header = TRUE)
individual_names <- colnames(dt_Haplotypes_colnames_Gen0)[11:ncol(dt_Haplotypes_colnames_Gen0)]
Haplotypes_colnames_Gen0 <- c("Position", individual_names[str_which(individual_names, "Shoe_0")])
length(Haplotypes_colnames_Gen0)
#### Prepare the data for the block calculation --------------------------------
prepare_haplotype_sequences <- function(Haplotypes,
                                        marker_positions,
                                        Haplotype_names){
  number <- length(22:length(Haplotypes))
  Haplotypes <- Haplotypes[22:length(Haplotypes)]
  # ncol of dt_Haplotypes_A1 is wrong, is 368, should be 366
  labels_of_subpop <- seq(1, number-2, 3)
  Labels_of_subpop <- Haplotypes[labels_of_subpop]
  haplotype_1 <- seq(2, number-1, 3)
  Haplotype_Allele_1 <- Haplotypes[haplotype_1]
  haplotype_2 <- seq(3, number, 3)
  Haplotype_Allele_2 <- Haplotypes[haplotype_2]
  prepare_haplotype_sequences <- function(Haplotypes_dt){
    
    prepare_haplotype_sequences_split <- function(i){
      splitted_haplo <- str_split(Haplotypes_dt[i], " ")
      unlisted_split_haplo <- unlist(splitted_haplo)
      return(unlisted_split_haplo)
    }
    dt_Haplotypes <- foreach(i = 1:length(Haplotypes_dt), .combine = cbind) %do% prepare_haplotype_sequences_split(i)
    return(dt_Haplotypes)
  }
  dt_Haplotypes_A1 <- prepare_haplotype_sequences(Haplotypes_dt = Haplotype_Allele_1)
  dt_Haplotypes_A2 <- prepare_haplotype_sequences(Haplotypes_dt = Haplotype_Allele_2)
  diploid_seq <- paste0(dt_Haplotypes_A1, dt_Haplotypes_A2)
  dt_diploid_seq <- matrix(data = diploid_seq, 
                           nrow = nrow(dt_Haplotypes_A1),
                           ncol = ncol(dt_Haplotypes_A1),
                           byrow = FALSE)
  dt_Haplotypes_A1 <- cbind(marker_positions, dt_Haplotypes_A1)
  dt_Haplotypes_A2 <- cbind(marker_positions, dt_Haplotypes_A2)
  dt_diploid_seq <- cbind(marker_positions, dt_diploid_seq)
  colnames(dt_Haplotypes_A1) <- Haplotype_names
  colnames(dt_Haplotypes_A2) <- Haplotype_names
  colnames(dt_diploid_seq) <- Haplotype_names
  Haplotypes <- list("Allele_1" = dt_Haplotypes_A1,
                     "Allele_2" = dt_Haplotypes_A2,
                     "Subpopulation_labels" = Labels_of_subpop,
                     "DT_Diploid_sequences" = dt_diploid_seq)
  rm(dt_Haplotypes_A1, dt_Haplotypes_A2)
  return(Haplotypes)
}
##### Region 1
new_Haplotypes_reg1 <- prepare_haplotype_sequences(Haplotypes = Haplotypes_reg1,
                                                   Haplotype_names = Haplotypes_colnames_reg,
                                                   marker_positions = marker_positions_reg1)
dt_Haplotypes_A1_reg1 <- as.data.frame(new_Haplotypes_reg1$Allele_1)
dt_Haplotypes_A2_reg1 <- as.data.frame(new_Haplotypes_reg1$Allele_2)
##### Region 2
new_Haplotypes_reg2 <- prepare_haplotype_sequences(Haplotypes = Haplotypes_reg2,
                                                   Haplotype_names = Haplotypes_colnames_reg,
                                                   marker_positions = marker_positions_reg2)
dt_Haplotypes_A1_reg2 <- as.data.frame(new_Haplotypes_reg2$Allele_1)
dt_Haplotypes_A2_reg2 <- as.data.frame(new_Haplotypes_reg2$Allele_2)
##### Gen0
new_Haplotypes_Gen0 <- prepare_haplotype_sequences(Haplotypes = Haplotypes_Gen0,
                                                   Haplotype_names = Haplotypes_colnames_Gen0,
                                                   marker_positions = marker_positions_Gen0)
dt_Haplotypes_A1_Gen0 <- as.data.frame(new_Haplotypes_Gen0$Allele_1)
dt_Haplotypes_A2_Gen0 <- as.data.frame(new_Haplotypes_Gen0$Allele_2)
#Haplotypes <- as.data.frame(new_Haplotypes$DT_Diploid_sequences)
### this is not how you do it!!!
rm(new_Haplotypes_reg1,
   new_Haplotypes_reg2,
   new_Haplotypes_Gen0)
##### Prepare the data set for block calculation for diploids ------------------
Haplotypes_reg1 <- cbind(dt_Haplotypes_A1_reg1, dt_Haplotypes_A2_reg1[,-1])
Haplotypes_reg2 <- cbind(dt_Haplotypes_A1_reg2, dt_Haplotypes_A2_reg2[,-1])
Haplotypes_Gen0 <- cbind(dt_Haplotypes_A1_Gen0, dt_Haplotypes_A2_Gen0[,-1])
##### something weird is happening
#### Divide data set into subgroups --------------------------------------------
##### Region 1
subpopulations <- str_sub(colnames(Haplotypes_reg1),1,9)
Haplotypes_subpop_short_1_reg1 <- Haplotypes_reg1[,which(subpopulations == "Shoepag_1")]
Haplotypes_subpop_short_2_reg1 <- Haplotypes_reg1[,which(subpopulations == "Shoepag_4")]
Haplotypes_subpop_tall_1_reg1 <- Haplotypes_reg1[,which(subpopulations == "Shoepag_3")]
Haplotypes_subpop_tall_2_reg1 <- Haplotypes_reg1[,which(subpopulations == "Shoepag_2")]
##### Region 2
Haplotypes_subpop_short_1_reg2 <- Haplotypes_reg2[,which(subpopulations == "Shoepag_1")]
Haplotypes_subpop_short_2_reg2 <- Haplotypes_reg2[,which(subpopulations == "Shoepag_4")]
Haplotypes_subpop_tall_1_reg2 <- Haplotypes_reg2[,which(subpopulations == "Shoepag_3")]
Haplotypes_subpop_tall_2_reg2 <- Haplotypes_reg2[,which(subpopulations == "Shoepag_2")]
##### Gen0
Haplotypes_subpop_short_1_Gen0 <- Haplotypes_Gen0[ ,which(str_detect(colnames(Haplotypes_Gen0),"Short_1"))]
Haplotypes_subpop_short_2_Gen0 <- Haplotypes_Gen0[ ,which(str_detect(colnames(Haplotypes_Gen0),"Short_2"))]
Haplotypes_subpop_tall_1_Gen0 <- Haplotypes_Gen0[ ,which(str_detect(colnames(Haplotypes_Gen0),"Tall_1"))]
Haplotypes_subpop_tall_2_Gen0 <- Haplotypes_Gen0[ ,which(str_detect(colnames(Haplotypes_Gen0),"Tall_2"))]
#### Include also the likelihoods of haplotypes ---------------------------------
#### Do the block calculation ---------------------------------------------------
##### Region 1
blocklist_short_1_reg1 <- block_calculation(dhm = Haplotypes_subpop_short_1_reg1,
                                            adaptive_mode = TRUE)
blocklist_short_2_reg1 <- block_calculation(dhm = Haplotypes_subpop_short_2_reg1,
                                            adaptive_mode = TRUE)
blocklist_tall_1_reg1 <- block_calculation(dhm = Haplotypes_subpop_tall_1_reg1,
                                           adaptive_mode = TRUE)
blocklist_tall_2_reg1 <- block_calculation(dhm = Haplotypes_subpop_tall_2_reg1,
                                           adaptive_mode = TRUE)
calc_blocks_short1_reg1 <- blocklist_startend(blocklist_short_1_reg1)
block_length_short1_reg1 <- calc_blocks_short1_reg1[ ,2] - calc_blocks_short1_reg1[ ,1]
mean(block_length_short1_reg1)

calc_blocks_short2_reg1 <- blocklist_startend(blocklist_short_2_reg1)
block_length_short2_reg1 <- calc_blocks_short2_reg1[ ,2] - calc_blocks_short2_reg1[ ,1]
mean(block_length_short2_reg1)

calc_blocks_tall1_reg1 <- blocklist_startend(blocklist_tall_1_reg1)
block_length_tall1_reg1 <- calc_blocks_tall1_reg1[ ,2] - calc_blocks_tall1_reg1[ ,1]
mean(block_length_tall1_reg1)

calc_blocks_tall2_reg1 <- blocklist_startend(blocklist_tall_2_reg1)
block_length_tall2_reg1 <- calc_blocks_tall2_reg1[ ,2] - calc_blocks_tall2_reg1[ ,1]
mean(block_length_tall2_reg1)
##### Region 2
blocklist_short_1_reg2 <- block_calculation(dhm = Haplotypes_subpop_short_1_reg2,
                                            adaptive_mode = TRUE)
blocklist_short_2_reg2 <- block_calculation(dhm = Haplotypes_subpop_short_2_reg2,
                                            adaptive_mode = TRUE)
blocklist_tall_1_reg2 <- block_calculation(dhm = Haplotypes_subpop_tall_1_reg2,
                                           adaptive_mode = TRUE)
blocklist_tall_2_reg2 <- block_calculation(dhm = Haplotypes_subpop_tall_2_reg2,
                                           adaptive_mode = TRUE)
calc_blocks_short1_reg2 <- blocklist_startend(blocklist_short_1_reg2)
block_length_short1_reg2 <- calc_blocks_short1_reg2[ ,2] - calc_blocks_short1_reg2[ ,1]
mean(block_length_short1_reg2)

calc_blocks_short2_reg2 <- blocklist_startend(blocklist_short_2_reg2)
block_length_short2_reg2 <- calc_blocks_short2_reg2[ ,2] - calc_blocks_short2_reg2[ ,1]
mean(block_length_short2_reg2)

calc_blocks_tall1_reg2 <- blocklist_startend(blocklist_tall_1_reg2)
block_length_tall1_reg2 <- calc_blocks_tall1_reg2[ ,2] - calc_blocks_tall1_reg2[ ,1]
mean(block_length_tall1_reg2)

calc_blocks_tall2_reg2 <- blocklist_startend(blocklist_tall_2_reg2)
block_length_tall2_reg2 <- calc_blocks_tall2_reg2[ ,2] - calc_blocks_tall2_reg2[ ,1]
mean(block_length_tall2_reg2)
##### Gen0
blocklist_short_1_Gen0 <- block_calculation(dhm = Haplotypes_subpop_short_1_Gen0,
                                            adaptive_mode = TRUE)
blocklist_short_2_Gen0 <- block_calculation(dhm = Haplotypes_subpop_short_2_Gen0,
                                            adaptive_mode = TRUE)
blocklist_tall_1_Gen0 <- block_calculation(dhm = Haplotypes_subpop_tall_1_Gen0,
                                           adaptive_mode = TRUE)
blocklist_tall_2_Gen0 <- block_calculation(dhm = Haplotypes_subpop_tall_2_Gen0,
                                           adaptive_mode = TRUE)
calc_blocks_short1_Gen0 <- blocklist_startend(blocklist_short_1_Gen0)
block_length_short1_Gen0 <- calc_blocks_short1_Gen0[ ,2] - calc_blocks_short1_Gen0[ ,1]
mean(block_length_short1_Gen0)

calc_blocks_short2_Gen0 <- blocklist_startend(blocklist_short_2_Gen0)
block_length_short2_Gen0 <- calc_blocks_short2_Gen0[ ,2] - calc_blocks_short2_Gen0[ ,1]
mean(block_length_short2_Gen0)

calc_blocks_tall1_Gen0 <- blocklist_startend(blocklist_tall_1_Gen0)
block_length_tall1_Gen0 <- calc_blocks_tall1_Gen0[ ,2] - calc_blocks_tall1_Gen0[ ,1]
mean(block_length_tall1_Gen0)

calc_blocks_tall2_Gen0 <- blocklist_startend(blocklist_tall_2_Gen0)
block_length_tall2_Gen0 <- calc_blocks_tall2_Gen0[ ,2] - calc_blocks_tall2_Gen0[ ,1]
mean(block_length_tall2_Gen0)
#### Do block calculation with subgroups for entire data set -------------------
#### Start Haplotype block plotting ---------------------------------------------
#blocklist_all_subpops <- block_calculation(dhm = Haplotypes[1:240,],
#adaptive_mode = TRUE,
#subgroups = list(index_short1,
#index_short2,
#index_tall1,
#index_tall2))
# does not work out on all individuals due to computational limitations
#blocklist_all_subpops
#### Prepare data for plotting -------------------------------------------------
get_haplotype_data_set <- function(blocklist,
                                   SNP_IDs){
  extract_haplotypes <- function(i){
    names_avail_samples <- paste0("Haplotype", "_", blocklist[[i]][[6]])
    data <- cbind(names_avail_samples, rep(paste("block", i)))
    return(data)
  }
  dt_haplotypes <- foreach(i = 1:length(lengths(blocklist)), .combine = rbind) %do% extract_haplotypes(i)
  dt_haplotypes <- as.data.frame(dt_haplotypes)
  colnames(dt_haplotypes) <- c("Haplotype","Block")
  calc_blocks <- blocklist_startend(blocklist)
  block_length <- calc_blocks[ ,2] - calc_blocks[ ,1]
  calc_blocks_new <- cbind(calc_blocks, block_length)
  calc_blocks <- calc_blocks_new
  get_block_length <- function(i){
    block_length_new <- as.numeric(block_length[which(unique(dt_haplotypes$Block)[i] == rownames(calc_blocks))])
    block_index <- which(unique(dt_haplotypes$Block)[i] == dt_haplotypes$Block)
    block_length_vec <- rep(block_length_new, length(block_index))
    #cat("For Block", i, block_length_vec, "\n")
    return(block_length_vec)
  }
  block_sizes <- foreach(i = 1:length(rownames(calc_blocks)), .combine = c) %do% get_block_length(i)
  dt_haplotypes$Block_Length <- block_sizes
  translate_block_length <- function(i){
    SNPs_in_Block <- SNP_IDs[calc_blocks[i,1]:calc_blocks[i,2]]
    window_start <- SNPs_in_Block[1]
    window_end <- SNPs_in_Block[length(SNPs_in_Block)]
    window_boundaries <- cbind(window_start, window_end)
    block_index <- which(unique(dt_haplotypes$Block)[i] == dt_haplotypes$Block)
    rep(window_boundaries, length(block_index))
    window_boundaries_vec <- rep(window_boundaries, length(block_index))
    return(window_boundaries_vec)
  }
  dt_window_boundaries <- foreach(i = 1:length(rownames(calc_blocks)), .combine = c) %do% translate_block_length(i)
  dt_window_boundaries <- matrix(data = dt_window_boundaries, 
                                 nrow = nrow(dt_haplotypes),
                                 ncol = 2, byrow = TRUE)
  dt_haplotypes$block_start <- dt_window_boundaries[,1]
  dt_haplotypes$block_end <- dt_window_boundaries[,2]
  block_start_in_bp <- as.numeric(str_sub(dt_haplotypes$block_end, 6))
  block_end_in_bp <- as.numeric(str_sub(dt_haplotypes$block_start, 6))
  dt_haplotypes$block_start_in_bp <- block_start_in_bp
  dt_haplotypes$block_end_in_bp <- block_end_in_bp
  return(dt_haplotypes)
}
##### Region 1
haplotypes_dt_short1_reg1 <- get_haplotype_data_set(blocklist = blocklist_short_1_reg1,
                                                    SNP_IDs = SNP_IDs_reg1)
haplotypes_dt_short2_reg1 <- get_haplotype_data_set(blocklist = blocklist_short_2_reg1,
                                                    SNP_IDs = SNP_IDs_reg1)
haplotypes_dt_tall1_reg1 <- get_haplotype_data_set(blocklist = blocklist_tall_1_reg1,
                                                   SNP_IDs = SNP_IDs_reg1)
haplotypes_dt_tall2_reg1 <- get_haplotype_data_set(blocklist = blocklist_tall_2_reg1,
                                                   SNP_IDs = SNP_IDs_reg1)
##### Region 2
haplotypes_dt_short1_reg2 <- get_haplotype_data_set(blocklist = blocklist_short_1_reg2,
                                                    SNP_IDs = SNP_IDs_reg2)
haplotypes_dt_short2_reg2 <- get_haplotype_data_set(blocklist = blocklist_short_2_reg2,
                                                    SNP_IDs = SNP_IDs_reg2)
haplotypes_dt_tall1_reg2 <- get_haplotype_data_set(blocklist = blocklist_tall_1_reg2,
                                                   SNP_IDs = SNP_IDs_reg2)
haplotypes_dt_tall2_reg2 <- get_haplotype_data_set(blocklist = blocklist_tall_2_reg2,
                                                   SNP_IDs = SNP_IDs_reg2)
#haplotypes_dt_all_subpop <- get_haplotype_data_set(blocklist = blocklist_all_subpops,
                                                   #SNP_IDs = SNP_IDs)
##### Gen0
haplotypes_dt_short1_Gen0 <- get_haplotype_data_set(blocklist = blocklist_short_1_Gen0,
                                                    SNP_IDs = SNP_IDs_Gen0)
haplotypes_dt_short2_Gen0 <- get_haplotype_data_set(blocklist = blocklist_short_2_Gen0,
                                                    SNP_IDs = SNP_IDs_Gen0)
haplotypes_dt_tall1_Gen0 <- get_haplotype_data_set(blocklist = blocklist_tall_1_Gen0,
                                                   SNP_IDs = SNP_IDs_Gen0)
haplotypes_dt_tall2_Gen0 <- get_haplotype_data_set(blocklist = blocklist_tall_2_Gen0,
                                                   SNP_IDs = SNP_IDs_Gen0)
##### Coding region of the gene iAA8 
start_SNPs_iAA8 <- 10060719
end_SNPs_iAA8 <- 10067095
##### Coding region of the gene Dwarf1 
start_SNPs_Dwarf1 <- 10440993
end_SNPs_Dwarf1 <- 10443340
###### How many blocks were observed? ------------------------------------------
####### Generation 0
haplotypes_dt_short1_Gen0$Block <- as.factor(haplotypes_dt_short1_Gen0$Block)
cat("Number of Blocks in Short 1:", length(levels(haplotypes_dt_short1_Gen0$Block)), "blocks","\n")
cat("On average, we observed:", round(mean(table(haplotypes_dt_short1_Gen0$Block)),0), "haplotypes", "\n")

haplotypes_dt_short2_Gen0$Block <- as.factor(haplotypes_dt_short2_Gen0$Block)
cat("Number of Blocks in Short 2:", length(levels(haplotypes_dt_short2_Gen0$Block)), "blocks","\n")
cat("On average, we observed:", round(mean(table(haplotypes_dt_short2_Gen0$Block)),0), "haplotypes", "\n")

haplotypes_dt_tall1_Gen0$Block <- as.factor(haplotypes_dt_tall1_Gen0$Block)
cat("Number of Blocks in Tall 1:", length(levels(haplotypes_dt_tall1_Gen0$Block)), "blocks","\n")
cat("On average, we observed:", round(mean(table(haplotypes_dt_tall1_Gen0$Block)),0), "haplotypes", "\n")

haplotypes_dt_tall2_Gen0$Block <- as.factor(haplotypes_dt_tall2_Gen0$Block)
cat("Number of Blocks in Tall 2:", length(levels(haplotypes_dt_tall2_Gen0$Block)), "blocks","\n")
cat("On average, we observed:", round(mean(table(haplotypes_dt_tall2_Gen0$Block)),0), "haplotypes", "\n")
####### Random Region 1
haplotypes_dt_short1_reg1$Block <- as.factor(haplotypes_dt_short1_reg1$Block)
cat("Number of Blocks in Short 1:", length(levels(haplotypes_dt_short1_reg1$Block)), "blocks","\n")
cat("On average, we observed:", round(mean(table(haplotypes_dt_short1_reg1$Block)),0), "haplotypes", "\n")

haplotypes_dt_short2_reg1$Block <- as.factor(haplotypes_dt_short2_reg1$Block)
cat("Number of Blocks in Short 2:", length(levels(haplotypes_dt_short2_reg1$Block)), "blocks","\n")
cat("On average, we observed:", round(mean(table(haplotypes_dt_short2_reg1$Block)),0), "haplotypes", "\n")

haplotypes_dt_tall1_reg1$Block <- as.factor(haplotypes_dt_tall1_reg1$Block)
cat("Number of Blocks in Tall 1:", length(levels(haplotypes_dt_tall1_reg1$Block)), "blocks","\n")
cat("On average, we observed:", round(mean(table(haplotypes_dt_tall1_reg1$Block)),0), "haplotypes", "\n")

haplotypes_dt_tall2_reg1$Block <- as.factor(haplotypes_dt_tall2_reg1$Block)
cat("Number of Blocks in Tall 2:", length(levels(haplotypes_dt_tall2_reg1$Block)), "blocks","\n")
cat("On average, we observed:", round(mean(table(haplotypes_dt_tall2_reg1$Block)),0), "haplotypes", "\n")
####### Random Region 2
haplotypes_dt_short1_reg2$Block <- as.factor(haplotypes_dt_short1_reg2$Block)
cat("Number of Blocks in Short 1:", length(levels(haplotypes_dt_short1_reg2$Block)), "blocks","\n")
cat("On average, we observed:", round(mean(table(haplotypes_dt_short1_reg2$Block)),0), "haplotypes", "\n")

haplotypes_dt_short2_reg2$Block <- as.factor(haplotypes_dt_short2_reg2$Block)
cat("Number of Blocks in Short 2:", length(levels(haplotypes_dt_short2_reg2$Block)), "blocks","\n")
cat("On average, we observed:", round(mean(table(haplotypes_dt_short2_reg2$Block)),0), "haplotypes", "\n")

haplotypes_dt_tall1_reg2$Block <- as.factor(haplotypes_dt_tall1_reg2$Block)
cat("Number of Blocks in Tall 1:", length(levels(haplotypes_dt_tall1_reg2$Block)), "blocks","\n")
cat("On average, we observed:", round(mean(table(haplotypes_dt_tall1_reg2$Block)),0), "haplotypes", "\n")

haplotypes_dt_tall2_reg2$Block <- as.factor(haplotypes_dt_tall2_reg2$Block)
cat("Number of Blocks in Tall 2:", length(levels(haplotypes_dt_tall2_reg2$Block)), "blocks","\n")
cat("On average, we observed:", round(mean(table(haplotypes_dt_tall2_reg2$Block)),0), "haplotypes", "\n")
###### Dwarf1
####### Generation 0
######## Short 1
haplotypes_dt_short1_Gen0$block_start <- as.numeric(haplotypes_dt_short1_Gen0$block_start)
index_1 <- which(haplotypes_dt_short1_Gen0$block_start < start_SNPs_Dwarf1)
index_2 <- which(haplotypes_dt_short1_Gen0$block_end < end_SNPs_Dwarf1)
dt <- haplotypes_dt_short1_Gen0[index_1[length(index_1)]:index_2[length(index_2)],]
dt$Block <- as.factor(dt$Block)
cat("In total we have", length(levels(dt$Block)), "blocks.", "\n", levels(dt$Block), "\n")
cat("In Block", levels(dt$Block)[1], "we observed", length(which(dt$Block == levels(dt$Block)[1])), "haplotypes.", "\n",
    "In Block", levels(dt$Block)[2], "we observed", length(which(dt$Block == levels(dt$Block)[2])), "haplotypes.", "\n")
dt[which(dt$Block == levels(dt$Block)[1]),]
dt[which(dt$Block == levels(dt$Block)[2]),]
######## Short 2
haplotypes_dt_short2_Gen0$block_start <- as.numeric(haplotypes_dt_short2_Gen0$block_start)
index_1 <- which(haplotypes_dt_short2_Gen0$block_start < start_SNPs_Dwarf1)
index_2 <- which(haplotypes_dt_short2_Gen0$block_end < end_SNPs_Dwarf1)
dt <- haplotypes_dt_short2_Gen0[index_1[length(index_1)]:index_2[length(index_2)],]
dt$Block <- as.factor(dt$Block)
cat("In total we have", length(levels(dt$Block)), "blocks.", "\n", levels(dt$Block), "\n")
cat("In Block", levels(dt$Block)[1], "we observed", length(which(dt$Block == levels(dt$Block)[1])), "haplotypes.", "\n",
    "In Block", levels(dt$Block)[2], "we observed", length(which(dt$Block == levels(dt$Block)[2])), "haplotypes.", "\n",
    "In Block", levels(dt$Block)[3], "we observed", length(which(dt$Block == levels(dt$Block)[3])), "haplotypes.", "\n")
dt[which(dt$Block == levels(dt$Block)[1]),]
dt[which(dt$Block == levels(dt$Block)[2]),]
dt[which(dt$Block == levels(dt$Block)[3]),]
######## Tall 1
haplotypes_dt_tall1_Gen0$block_start <- as.numeric(haplotypes_dt_tall1_Gen0$block_start)
index_1 <- which(haplotypes_dt_tall1_Gen0$block_start < start_SNPs_Dwarf1)
index_2 <- which(haplotypes_dt_tall1_Gen0$block_end < end_SNPs_Dwarf1)
dt <- haplotypes_dt_tall1_Gen0[index_1[length(index_1)]:index_2[length(index_2)],]
dt$Block <- as.factor(dt$Block)
cat("In total we have", length(levels(dt$Block)), "blocks.", "\n", levels(dt$Block), "\n")
cat("In Block", levels(dt$Block)[1], "we observed", length(which(dt$Block == levels(dt$Block)[1])), "haplotypes.", "\n",
    "In Block", levels(dt$Block)[2], "we observed", length(which(dt$Block == levels(dt$Block)[2])), "haplotypes.", "\n",
    "In Block", levels(dt$Block)[3], "we observed", length(which(dt$Block == levels(dt$Block)[3])), "haplotypes.", "\n")
dt[which(dt$Block == levels(dt$Block)[1]),]
dt[which(dt$Block == levels(dt$Block)[2]),]
dt[which(dt$Block == levels(dt$Block)[3]),]
######## Tall 2
haplotypes_dt_tall2_Gen0$block_start <- as.numeric(haplotypes_dt_tall2_Gen0$block_start)
index_1 <- which(haplotypes_dt_tall2_Gen0$block_start < start_SNPs_Dwarf1)
index_2 <- which(haplotypes_dt_tall2_Gen0$block_end < end_SNPs_Dwarf1)
dt <- haplotypes_dt_tall2_Gen0[index_1[length(index_1)]:index_2[length(index_2)],]
dt$Block <- as.factor(dt$Block)
cat("In total we have", length(levels(dt$Block)), "blocks.", "\n", levels(dt$Block), "\n")
cat("In Block", levels(dt$Block)[1], "we observed", length(which(dt$Block == levels(dt$Block)[1])), "haplotypes.", "\n",
    "In Block", levels(dt$Block)[2], "we observed", length(which(dt$Block == levels(dt$Block)[2])), "haplotypes.", "\n",
    "In Block", levels(dt$Block)[3], "we observed", length(which(dt$Block == levels(dt$Block)[3])), "haplotypes.", "\n")
dt[which(dt$Block == levels(dt$Block)[1]),]
dt[which(dt$Block == levels(dt$Block)[2]),]
dt[which(dt$Block == levels(dt$Block)[3]),]
###### iAA8 gene
distance <- 10000

#### Prepare the data for plotting ---------------------------------------------
#### Plotting ------------------------------------------------------------------
#### V1: Plot haplotypes 9.2 Mb region -----------------------------------------
#### Load the parameters 
##### Coding region of the gene iAA8 
start_SNPs_iAA8 <- 10060719
end_SNPs_iAA8 <- 10067095
##### Coding region of the gene Dwarf1 
start_SNPs_Dwarf1 <- 10440993
end_SNPs_Dwarf1 <- 10443340
##### Start plotting 
font_size <- 18
plot_Haplotype_Blocks <- function(blocklist_dt,
                                  marker_positions,
                                  Subpopulation_name){
  blocklist_dt$Haplotype <- as.factor(fct_inorder(blocklist_dt$Haplotype))
  Haplo_in_the_middle <- paste0("Haplotype_",length(levels(blocklist_dt$Haplotype))/2)
  Haplo_above_that <- paste0("Haplotype_",(length(levels(blocklist_dt$Haplotype))/2)+40)
  Haplo_slightly_above_that <- paste0("Haplotype_",(length(levels(blocklist_dt$Haplotype))/2)+30)
  font_size <- font_size
  windowsFonts(my = windowsFont('Calibri'))
  blocklist_dt$Block <- as.factor(blocklist_dt$Block)
  n <- length(levels(blocklist_dt$Block))
  nr_colours <- round(n/6,0)+1
  my_pal_col_1 <- (c(rep(c("tomato2","cyan4","darkseagreen3","dodgerblue1","pink1","mediumpurple"), nr_colours)))
  Haplotype_blocks_plot <- ggplot(blocklist_dt, aes(xmin = block_start, xmax = block_end,
                                                    ymin = Haplotype, ymax = Haplotype)) +
    geom_rect(data = blocklist_dt,
              aes(xmin = block_start, xmax = block_end,
                  ymin = Haplotype, ymax = Haplotype,
                  fill = Block, colour = Block))+
    theme(text = element_text(size = font_size, family = "my"),
          panel.background = element_rect(fill = "grey42"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth=1),
          legend.title = element_blank(),
          legend.position = "none",
          axis.line.y = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          #axis.title.x = element_text(size = font_size, family = "my", colour = "black"),
          #axis.text.x = element_text(size = font_size, family = "my", colour = "black"),
          axis.title.y = element_text(size = font_size, family = "my", colour = "black"))+
    scale_fill_manual(values = my_pal_col_1)+
    scale_color_manual(values = my_pal_col_1)+
    geom_vline(aes(xintercept = 10060719), colour = "firebrick", linewidth = 0.5)+
    geom_vline(aes(xintercept = 10067095), colour = "firebrick", linewidth = 0.5)+
    geom_vline(aes(xintercept = 10440993), colour = "firebrick", linewidth = 0.5)+
    geom_vline(aes(xintercept = 10443340), colour = "firebrick", linewidth = 0.5)+
    xlim(marker_positions[1], marker_positions[length(marker_positions)])+
    labs(title = Subpopulation_name, x = "Position",
         y = "Haplotypes")
  Haplotype_blocks_plot
  return(Haplotype_blocks_plot)
}
plot_Haplotype_Blocks_tall2 <- function(blocklist_dt,
                                        Subpopulation_name,
                                        marker_positions){
  blocklist_dt$Haplotype <- as.factor(fct_inorder(blocklist_dt$Haplotype))
  blocklist_dt$block_start <- as.numeric(blocklist_dt$block_start)
  blocklist_dt$block_end <- as.numeric(blocklist_dt$block_end)
  marker_positions <- as.numeric(marker_positions)
  font_size <- font_size
  windowsFonts(my = windowsFont('Calibri'))
  blocklist_dt$Block <- as.factor(blocklist_dt$Block)
  n <- length(levels(blocklist_dt$Block))
  nr_colours <- round(n/6,0)+1
  my_pal_col_1 <- (c(rep(c("tomato2","cyan4","darkseagreen3","dodgerblue1","pink1","mediumpurple"), nr_colours)))
  Haplotype_blocks_plot <- ggplot(blocklist_dt, aes(xmin = block_start, xmax = block_end,
                                                    ymin = Haplotype, ymax = Haplotype)) +
    geom_rect(data = blocklist_dt,
              aes(xmin = block_start, xmax = block_end,
                  ymin = Haplotype, ymax = Haplotype,
                  fill = Block, colour = Block))+
    theme(text = element_text(size = font_size, family = "my"),
          panel.background = element_rect(fill = "grey42"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA, linewidth=1),
          legend.title = element_blank(),
          legend.position = "none",
          axis.line.y = element_blank(),
          axis.line.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_text(size = font_size, family = "my", colour = "black"),
          axis.text.x = element_text(size = font_size-2, family = "my", colour = "black", angle = 75, hjust = 1),
          axis.title.y = element_text(size = font_size, family = "my", colour = "black"))+
    scale_fill_manual(values = my_pal_col_1)+
    scale_color_manual(values = my_pal_col_1)+
    labs(title = Subpopulation_name, x = "Position",
         y = "Haplotypes")+
    geom_vline(aes(xintercept = 10060719), colour = "firebrick", linewidth = 0.5)+
    geom_vline(aes(xintercept = 10067095), colour = "firebrick", linewidth = 0.5)+
    geom_vline(aes(xintercept = 10440993), colour = "firebrick", linewidth = 0.5)+
    geom_vline(aes(xintercept = 10443340), colour = "firebrick", linewidth = 0.5)+
    xlim(marker_positions[1], marker_positions[length(marker_positions)])
  Haplotype_blocks_plot
  return(Haplotype_blocks_plot)
}
###### Region 1
Short_1_plot_reg1 <- plot_Haplotype_Blocks(blocklist_dt = haplotypes_dt_short1_reg1,
                                           marker_positions = marker_positions_reg1,
                                           Subpopulation_name = "Short 1")
Short_1_plot_reg1
Short_2_plot_reg1 <- plot_Haplotype_Blocks(blocklist_dt = haplotypes_dt_short2_reg1,
                                           marker_positions = marker_positions_reg1,
                                           Subpopulation_name = "Short 2")
Short_2_plot_reg1
Tall_1_plot_reg1 <- plot_Haplotype_Blocks(blocklist_dt = haplotypes_dt_tall1_reg1,
                                          marker_positions = marker_positions_reg1,
                                          Subpopulation_name = "Tall 1")
Tall_1_plot_reg1
Tall_2_plot_reg1 <- plot_Haplotype_Blocks_tall2(blocklist_dt = haplotypes_dt_tall2_reg1,
                                           Subpopulation_name = "Tall 2",
                                           marker_positions = marker_positions_reg1)
Tall_2_plot_reg1
setwd("C:/Users/mtost/Documents/Shoepeg_paper_final_versions/Scripts/Results/")
png("2023-04-04-Haplotypes_random_region_1.png", width = 16, height = 8.25, units="in", res = 1800)
grid.arrange(Short_1_plot_reg1,
             Short_2_plot_reg1,
             Tall_1_plot_reg1,
             Tall_2_plot_reg1,
             nrow = 4,
             ncol = 1, 
             heights = c(4,4,4,6),
             widths = 8)
dev.off()
combined_haplotype_block_plot <- grid.arrange(Short_1_plot,
                                              Short_2_plot,
                                              Tall_1_plot,
                                              Tall_2_plot,
                                              nrow = 4,
                                              ncol = 1, 
                                              heights = c(4,4,4,6),
                                              widths = 8)
###### Region 2
Short_1_plot_reg2 <- plot_Haplotype_Blocks(blocklist_dt = haplotypes_dt_short1_reg2,
                                           marker_positions = marker_positions_reg2,
                                           Subpopulation_name = "Short 1")
Short_1_plot_reg2
Short_2_plot_reg2 <- plot_Haplotype_Blocks(blocklist_dt = haplotypes_dt_short2_reg2,
                                           marker_positions = marker_positions_reg2,
                                           Subpopulation_name = "Short 2")
Short_2_plot_reg2
Tall_1_plot_reg2 <- plot_Haplotype_Blocks(blocklist_dt = haplotypes_dt_tall1_reg2,
                                          marker_positions = marker_positions_reg2,
                                          Subpopulation_name = "Tall 1")
Tall_1_plot_reg2
Tall_2_plot_reg2 <- plot_Haplotype_Blocks_tall2(blocklist_dt = haplotypes_dt_tall2_reg2,
                                                Subpopulation_name = "Tall 2",
                                                marker_positions = marker_positions_reg2)
Tall_2_plot_reg2
setwd("C:/Users/mtost/Documents/Shoepeg_paper_final_versions/Scripts/Results/")
png("2023-04-04-Haplotypes_random_region_2.png", width = 16, height = 8.25, units="in", res = 1800)
grid.arrange(Short_1_plot_reg2,
             Short_2_plot_reg2,
             Tall_1_plot_reg2,
             Tall_2_plot_reg2,
             nrow = 4,
             ncol = 1, 
             heights = c(4,4,4,6),
             widths = 8)
dev.off()
combined_haplotype_block_plot_reg2 <- grid.arrange(Short_1_plot_reg2,
                                                   Short_2_plot_reg2,
                                                   Tall_1_plot_reg2,
                                                   Tall_2_plot_reg2,
                                                   nrow = 4,
                                                   ncol = 1, 
                                                   heights = c(4,4,4,6),
                                                   widths = 8)
###### Gen0
Short_1_plot_Gen0 <- plot_Haplotype_Blocks(blocklist_dt = haplotypes_dt_short1_Gen0,
                                           marker_positions = marker_positions_Gen0,
                                           Subpopulation_name = "Short 1")
Short_1_plot_Gen0
Short_2_plot_Gen0 <- plot_Haplotype_Blocks(blocklist_dt = haplotypes_dt_short2_Gen0,
                                           marker_positions = marker_positions_Gen0,
                                           Subpopulation_name = "Short 2")
Short_2_plot_Gen0
Tall_1_plot_Gen0 <- plot_Haplotype_Blocks(blocklist_dt = haplotypes_dt_tall1_Gen0,
                                          marker_positions = marker_positions_Gen0,
                                          Subpopulation_name = "Tall 1")
Tall_1_plot_Gen0
Tall_2_plot_Gen0 <- plot_Haplotype_Blocks_tall2(blocklist_dt = haplotypes_dt_tall2_Gen0,
                                                Subpopulation_name = "Tall 2",
                                                marker_positions = marker_positions_Gen0)
Tall_2_plot_Gen0
### Calculate block length -----------------------------------------------------
mean(Short_1_plot_Gen0$data$block_end-Short_1_plot_Gen0$data$block_start)
mean(Short_2_plot_Gen0$data$block_end-Short_2_plot_Gen0$data$block_start)
mean(Tall_1_plot_Gen0$data$block_end-Tall_1_plot_Gen0$data$block_start)
mean(Tall_2_plot_Gen0$data$block_end-Tall_2_plot_Gen0$data$block_start)
### Plot only coding region ----------------------------------------------------
marker_positions_Gen0_iAA8 <- marker_positions_Gen0[which(marker_positions_Gen0 > start_SNPs_iAA8-1000000 & 
                                                          marker_positions_Gen0 < end_SNPs_iAA8+1000000)]
haplotypes_dt_tall2_Gen0_iAA8 <- haplotypes_dt_tall2_Gen0[which(haplotypes_dt_tall2_Gen0$block_start > start_SNPs_iAA8-1000000 & 
                                                                haplotypes_dt_tall2_Gen0$block_end < end_SNPs_iAA8+1000000), ]
plot_Haplotype_Blocks_tall2(blocklist_dt = haplotypes_dt_tall2_Gen0_iAA8,
                            Subpopulation_name = "Tall 2",
                            marker_positions = marker_positions_Gen0_iAA8)
setwd("C:/Users/mtost/Documents/Shoepeg_paper_final_versions/Scripts/Results/")
png("2023-04-04-Haplotypes_candidate_gene_region_Generation0.png", width = 16, height = 8.25, units="in", res = 1800)
grid.arrange(Short_1_plot_Gen0,
             Short_2_plot_Gen0,
             Tall_1_plot_Gen0,
             Tall_2_plot_Gen0,
             nrow = 4,
             ncol = 1, 
             heights = c(4,4,4,6),
             widths = 8)
dev.off()
combined_haplotype_block_plot_reg2 <- grid.arrange(Short_1_plot_reg2,
                                                   Short_2_plot_reg2,
                                                   Tall_1_plot_reg2,
                                                   Tall_2_plot_reg2,
                                                   nrow = 4,
                                                   ncol = 1, 
                                                   heights = c(4,4,4,6),
                                                   widths = 8)


## Add LD-Decay plot -----------------------------------------------------------
setwd("C:/Users/mtost/Documents/Shoepeg_paper_final_versions/Scripts/Results/")
SNP_data <- fread("2023-03-10_haplotypes_of_CHR3_candidate_gene_regions.txt",
                       skip = 1)
names_SNP_data <- read.table("2023-03-10_haplotypes_of_CHR3_candidate_gene_regions.txt",
                       header = TRUE)
SNP_data <- as.data.frame(SNP_data)
colnames(SNP_data) <- colnames(names_SNP_data)
rownames(SNP_data) <- rownames(names_SNP_data)
### NEW APPROACH ---------------------------------------------------------------
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("snpStats")
# Install the latest development version from GitHub with
library(LDheatmap)
library(vcfR)
require(snpStats)
# Load the data ----------------------------------------------------------------
# Matrix of genotypes, coded as 0, 1 2 copies of an index allele
SNP_data <- SNP_data[, -1]
dim(SNP_data)
str(SNP_data)
SNP_data <- as.matrix(SNP_data)
gdat <- matrix(SNP_data, ncol = 378)
gdat
dim(gdat)
gdat <- as(gdat,"SnpMatrix")
dim(gdat)
SNPs <- unique(round(as.numeric(str_sub(rownames(SNP_data), 6))/100)*100)
length(SNPs)
#> object has no names - using numeric order for row/column names
rgb.palette <- colorRampPalette(rev(c("pink", "tomato","firebrick")), space = "rgb")
LD_plot <- LDheatmap(gdat,genetic.distances=SNPs, color=rgb.palette(5),
                     add.key = TRUE,
                     add.map = FALSE,
                     SNP.name = TRUE, flip = TRUE)
LDmat <- as.data.frame(LD_plot$LDmatrix)
colnames(LDmat) <- SNPs

setwd("C:/Users/mtost/Documents/Shoepeg_paper_final_versions/Scripts/Results/")
png("2023-03-10-LD_decay_plot.png", width = 16, height = 8.25, units="in", res = 1800)
LDheatmap(gdat,genetic.distances=SNPs, color=rgb.palette(5),
          add.key = TRUE,
          add.map = FALSE,
          SNP.name = TRUE, flip = TRUE)
dev.off()
### Create own legend for data -------------------------------------------------
windowsFonts(my = windowsFont('Calibri'))
rgb.palette <- colorRampPalette(rev(c("firebrick", "tomato","pink")), space = "rgb")
font_size <- 12
R2_values <- sample(seq(0.1,1,0.2), 1200, replace = TRUE)
x <- paste0("SNP", seq(1,1200))
y <- paste0("SNP", seq(1201,2400))
data <- cbind(x, y, R2_values)
data <- as.data.frame(data)
colnames(data) <- c("x", "y", "R2_values")
legend_chr_plot <- ggplot(data = data,
                          aes(y = y, x = x))+
  geom_tile(aes(fill = R2_values), colour = "white")+
  theme(text = element_text(size = 12, family = "my"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size =font_size, family = "my", colour = "black", face = "bold"),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size =font_size, family = "my", colour = "black", face = "bold"),
        axis.text.x = element_blank(),
        legend.text = element_text(size =font_size, family = "my", colour = "black"),
        legend.title = element_text(size =font_size, family = "my", colour = "black", face = "bold"),
        legend.position = "right")+
  labs(y = "LD", x = "Position")+
  scale_fill_manual(name = expression(paste(R^{2})),
                    values = rgb.palette(5),
                    labels = c("0 - 0.2", 
                               "0.2 - 0.4",
                               "0.4 - 0.6",
                               "0.6 - 0.8",
                               "0.8 - 1.0"))
legend_chr_plot
get_legend <- function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  legend
} 
legend <- get_legend(legend_chr_plot)
setwd("C:/Users/mtost/Documents/Shoepeg_paper_final_versions/Scripts/Results/")
png(paste0(Sys.Date(),"_legend_for_LD_plot.png"), width = 2, height = 8.25, units="in", res = 1800)
ggplot(data = data,
       aes(y = y, x = x))+
  geom_tile(aes(fill = R2_values), colour = "white")+
  theme(text = element_text(size = 12, family = "my"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        axis.line.y = element_line(colour = "black"),
        axis.line.x = element_line(colour = "black"),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size =font_size, family = "my", colour = "black", face = "bold"),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size =font_size, family = "my", colour = "black", face = "bold"),
        axis.text.x = element_blank(),
        legend.text = element_text(size =font_size, family = "my", colour = "black"),
        legend.title = element_text(size =font_size, family = "my", colour = "black", face = "bold"),
        legend.position = "right")+
  labs(y = "LD", x = "Position")+
  scale_fill_manual(name = expression(paste(R^{2})),
                    values = rgb.palette(5),
                    labels = c("0 - 0.2", 
                               "0.2 - 0.4",
                               "0.4 - 0.6",
                               "0.6 - 0.8",
                               "0.8 - 1.0"))
dev.off()
