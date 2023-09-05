#!/bin/bash
#SBATCH -p medium
#SBATCH -t 48:00:00 
#SBATCH -J fastPHASE_9MB_candidate_on_chr3
#SBATCH -c 8
#SBATCH -N 1
#SBATCH --mem=200G
#SBATCH -o 2023_02_15_fastPHASE_9MB_candidate_on_chr3.out
#SBATCH -e 2023_02_15_fastPHASE_9MB_candidate_on_chr3.err
#SBATCH --mail-type=ALL

./fastPHASE -T10 -unew_labels.txt -o2023_02_15_Results_of_9MB_candidate_on_chr3 /usr/users/mtost/Shoepeg_resubmission_new_analysis/Results/2023-02-14_Shoepeg_9.2Mb_region_pos_9435824_10443340.inp

