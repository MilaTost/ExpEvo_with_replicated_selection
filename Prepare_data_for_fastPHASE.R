# Prepare data to visualize haplotypes -----------------------------------------
#### Installing the packages ---------------------------------------------------
cat("Package installation starts.","\n")
install.packages("data.table",
                 lib = "/home/uni08/mtost/R/x86_64-pc-linux-gnu-library/4.1",
                 repos='http://cran.rstudio.com/')
install.packages("parallel",
                 lib = "/home/uni08/mtost/R/x86_64-pc-linux-gnu-library/4.1",
                 repos='http://cran.rstudio.com/')
install.packages("stringr",
                 lib = "/home/uni08/mtost/R/x86_64-pc-linux-gnu-library/4.1",
                 repos='http://cran.rstudio.com/')
install.packages("doMC",
                 lib = "/home/uni08/mtost/R/x86_64-pc-linux-gnu-library/4.1",
                 repos='http://cran.rstudio.com/')
install.packages("vcfR",
                 lib = "/home/uni08/mtost/R/x86_64-pc-linux-gnu-library/4.1",
                 repos='http://cran.rstudio.com/',
                 INSTALL_opts = '--no-lock')
#### Load required packages ----------------------------------------------------
library(stringr)
library(data.table)
library(doMC)
library(vcfR)
library(parallel)
#### Load data sets ------------------------------------------------------------
cores <- as.integer(Sys.getenv('SLURM_CPUS_PER_TASK'))
cores <- detectCores(all.tests = FALSE, logical = TRUE)
registerDoMC(cores)
time_0 <- Sys.time()
cat("The data set loading starts")
# vcf_S4 <- read.vcfR("YOUR_DATA_after_filtering.vcf.gz", verbose = FALSE)
#setwd("/path/to/your/own/working/directory/")
vcf_S4 <- read.vcfR("/usr/users/mtost/wd_GBeasy_rerun/Results/2022_GB10_Shoepeg_after_filtering.vcf.gz", verbose = FALSE)
cat("The data set is loaded.","\n")
time_1 <- Sys.time()
print(time_1 - time_0)
### Test script on local device ------------------------------------------------
#rm(list = ls())
#setwd("C:/Users/mtost/Documents/Shoepeg_paper_final_versions/Scripts/Test_dir/")
#vcf_S4 <- read.vcfR("Shoepeg_small.vcf.gz", verbose = FALSE)
#### Subset the data set for CHR3 ----------------------------------------------
gt <- extract.gt(vcf_S4, element = "GT", return.alleles = TRUE)
index_chr3 <- which(str_detect(rownames(gt),"chr3"))
gt <- gt[index_chr3,]
#### Significant region from FSTSum approach -----------------------------------
start_SNPs <- 9435824
end_SNPs <- 10355846
#### Start the subsetting ------------------------------------------------------
SNP_IDs <- rownames(gt)
marker_positions <- str_sub(SNP_IDs, 6)
pos_index_start <- which(marker_positions == start_SNPs)
pos_index_end <- which(marker_positions == end_SNPs)
gt_chr3_cand_reg <- gt[pos_index_start:pos_index_end, ]
gt <- gt_chr3_cand_reg
rm(gt_chr3_cand_reg)
#### How many markers do we have? ----------------------------------------------
cat("The region starts at position", start_SNPs, "and ends at", end_SNPs, "\n",
    "and contains", nrow(gt), "markers.", "\n")
### Save the new col- and rownames ---------------------------------------------
rownames_of_gt <- rownames(gt)
colnames_of_gt <- colnames(gt)
### Load the parameters --------------------------------------------------------
nr_markers <- nrow(gt)
nr_individuals <- ncol(gt)
### Get infos about the input file ---------------------------------------------
cat("The input file starts at marker position", rownames_of_gt[1], "\n",
    "and ends at marker position", rownames_of_gt[length(rownames_of_gt)], "\n")
cat("The number of markers is:", nr_markers, "\n")
cat("The number of individuals is:", nr_individuals, "\n")
### Prepare the data set -------------------------------------------------------
return_proper_VCF <- function(data){
  replace_single_dots <- function(i){
    gt_per_ind <- data[,i]
    index_dots <- which(gt_per_ind == ".")
    gt_per_ind[index_dots] <- "?/?"
    return(gt_per_ind)
  }
  new_gt <- mclapply(1:ncol(data), replace_single_dots)
  return(new_gt)
}
new_gt <- return_proper_VCF(data = gt)
### Construct entire input file and write it directly --------------------------
result_directory <- "/usr/users/mtost/Shoepeg_resubmission_new_analysis/Results/"
setwd(result_directory)
sink("2023_01_12_Shoepeg_9MB_region_of_candidate_genes.inp")
for (i in 1:nr_individuals) {
  return_seq_per_ind <- function(data,
                                 ind,
                                 nr_markers,
                                 nr_individuals){
    ind <- data[[ind]]
    chr3_sequence <- str_split(ind, "/")
    chr3_sequence <- unlist(chr3_sequence)
    haplotype_set_1 <- seq(1, (nr_markers)*2, 2)
    haplotype_set_2 <- seq(2, (nr_markers)*2, 2)
    haplotype_dt_1 <- chr3_sequence[haplotype_set_1]
    haplotype_dt_2 <- chr3_sequence[haplotype_set_2]
    cat("Shoepeg", i, "\n",
        haplotype_dt_1, "\n",
        haplotype_dt_2, "\n")
  }
  ind <- return_seq_per_ind(data = new_gt,
                            ind = i,
                            nr_markers = nr_markers,
                            nr_individuals = nr_individuals)
}
sink()
2733/44
