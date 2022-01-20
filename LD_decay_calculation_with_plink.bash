#!/bin/bash
#SBATCH -p fat
#SBATCH -t 48:00:00
#SBATCH -J LD_decay_calc_with_plink
#SBATCH -c 8
#SBATCH -N 1
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mila.tost@agr.uni-goettingen.de
#SBATCH -o LD_decay_calc_with_plink.out
#SBATCH -e LD_decay_calc_with_plink.err

module load plink/1.90
plink --r2 --vcf /usr/users/mtost/wd_GBeasy_rerun/Results/2020_Shoepag_results_GB10.vcf --ld-window-r2 0.0



