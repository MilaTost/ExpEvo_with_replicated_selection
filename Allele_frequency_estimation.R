########### Allele frequency calculation ---------------------------------------
### Get the allele frequencies per population ----------------------------------
calculate_the_allele_freq <- function(data,
                                      population_name_1,
                                      population_name_2,
                                      population_name_3,
                                      population_name_4){
  time_0 <- Sys.time()
  cat("Calculation of the allele frequency spectrum started.","\n")
  gt <- extract.gt(data, element = "GT", as.numeric = FALSE)
  gt_file <- gt
  pop_1 <- gt_file[,which(str_detect(colnames(gt_file),population_name_1))]
  pop_2 <- gt_file[,which(str_detect(colnames(gt_file),population_name_2))]
  pop_3 <- gt_file[,which(str_detect(colnames(gt_file),population_name_3))]
  pop_4 <- gt_file[,which(str_detect(colnames(gt_file),population_name_4))]
  dp <- extract.gt(data, element = "DP")
  SNP_IDs <- rownames(dp)
  get_allele_freq_per_marker_per_pop <- function(data, population_name){
    get_allele_freq_per_marker <- function(i){
      N_alleles_total <- 2*length(which(!is.na(data[i,])))
      N_ALT_alleles <- 2*length(which(data[i,] == "1|1")) + 2*length(which(data[i,] == "1/1"))
      N_REF_alleles <- 2*length(which(data[i,] == "0|0")) + 2*length(which(data[i,] == "0/0"))
      N_HET_alleles <- length(which(data[i,] == "1|0")) + length(which(data[i,] == "1/0")) +
        length(which(data[i,] == "0|1")) + length(which(data[i,] == "0/1"))
      N_GT_tot <- length(which(!is.na(data[i,])))
      N_P_GT <- length(which(data[i,] == "0|0")) + length(which(data[i,] == "0/0"))
      N_Q_GT <- length(which(data[i,] == "1|1")) + length(which(data[i,] == "1/1"))
      freq_REF <- (N_REF_alleles + N_HET_alleles)/N_alleles_total
      freq_ALT <- (N_ALT_alleles + N_HET_alleles)/N_alleles_total
      freq_HET <- N_HET_alleles/N_GT_tot
      freq_P <- N_P_GT/N_GT_tot
      freq_Q <- N_Q_GT/N_GT_tot
      Exp_Het <- 2*N_P_GT/N_GT_tot*N_Q_GT/N_GT_tot
      marker_info <- cbind(freq_REF,freq_ALT,freq_HET,freq_P,freq_Q,N_GT_tot,Exp_Het)
      return(marker_info)
    }
    af_per_Marker_pop <- unlist(mclapply(1:nrow(data),get_allele_freq_per_marker))
    marker_tab_pop <- matrix(af_per_Marker_pop, ncol = 7, byrow = TRUE)
    marker_tab_pop <- as.data.frame(marker_tab_pop)
    colnames(marker_tab_pop) <- c(paste("freq_REF",population_name,sep="_"),paste("freq_ALT",population_name,sep="_"),
                                  paste("freq_HET",population_name,sep="_"),paste("freq_P",population_name,sep="_"),
                                  paste("freq_Q",population_name,sep="_"),paste("N_GT_tot",population_name,sep="_"),
                                  paste("Exp_HET",population_name,sep="_"))
    rownames(marker_tab_pop) <- SNP_IDs
    return(marker_tab_pop)
  }
  pop_1_dt <- get_allele_freq_per_marker_per_pop(data=pop_1, population_name = population_name_1)
  pop_2_dt <- get_allele_freq_per_marker_per_pop(data=pop_2, population_name = population_name_2)
  pop_3_dt <- get_allele_freq_per_marker_per_pop(data=pop_3, population_name = population_name_3)
  pop_4_dt <- get_allele_freq_per_marker_per_pop(data=pop_4, population_name = population_name_4)
  marker_table_per_pop <- cbind(pop_1_dt,pop_2_dt,
                                pop_3_dt,pop_4_dt)
  marker_table_per_pop <- as.data.frame(marker_table_per_pop)
  cat("Pop 1: For",length(which(marker_table_per_pop[,3] < marker_table_per_pop[,7])),
      "markers the observed heterozygosity was higher than the expected heterozygosity!","\n")
  cat("Pop 2: For",length(which(marker_table_per_pop[,10] < marker_table_per_pop[,14])),
      "markers the observed heterozygosity was higher than the expected heterozygosity!","\n")
  cat("Pop 3: For",length(which(marker_table_per_pop[,17] < marker_table_per_pop[,21])),
      "markers the observed heterozygosity was higher than the expected heterozygosity!","\n")
  cat("Pop 4: For",length(which(marker_table_per_pop[,24] < marker_table_per_pop[,28])),
      "markers the observed heterozygosity was higher than the expected heterozygosity!","\n")
  cat("Pop 1:",length(which(marker_table_per_pop[,1] < marker_table_per_pop[,2])),"times the REF allele was the minor allele","\n")
  cat("Pop 1:",length(which(marker_table_per_pop[,1] > marker_table_per_pop[,2])),"times the REF allele was the major allele","\n")
  cat("Pop 2:",length(which(marker_table_per_pop[,8] < marker_table_per_pop[,9])),"times the REF allele was the minor allele","\n")
  cat("Pop 2:",length(which(marker_table_per_pop[,8] > marker_table_per_pop[,9])),"times the REF allele was the major allele","\n")
  cat("Pop 3:",length(which(marker_table_per_pop[,15] < marker_table_per_pop[,16])),"times the REF allele was the minor allele","\n")
  cat("Pop 3:",length(which(marker_table_per_pop[,15] > marker_table_per_pop[,16])),"times the REF allele was the major allele","\n")
  cat("Pop 4:",length(which(marker_table_per_pop[,21] < marker_table_per_pop[,22])),"times the REF allele was the minor allele","\n")
  cat("Pop 4:",length(which(marker_table_per_pop[,21] > marker_table_per_pop[,22])),"times the REF allele was the major allele","\n")
  cat("Calculation of the allele frequency spectrum is done.","\n")
  time_1 <- Sys.time()
  cat("This was done within",print(time_1 - time_0),"\n")
  return(marker_table_per_pop)
}
marker_table_per_pop <- calculate_the_allele_freq_diff(data = vcf_S4, 
                                                       population_name_1 = "Shoepag_1",
                                                       population_name_2 = "Shoepag_2",
                                                       population_name_3 = "Shoepag_3",
                                                       population_name_4 = "Shoepag_4")
write.table(marker_table_per_pop,"2022_01_17_GB1004_Allele_and_GT_freq_per_pop.txt", sep = "  ", row.names = TRUE,
            quote = FALSE)