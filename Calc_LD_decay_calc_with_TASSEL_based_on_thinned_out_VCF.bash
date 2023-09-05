#!/bin/bash
#SBATCH -p medium
#SBATCH -t 48:00:00 
#SBATCH -J LD_decay_TASSEL
#SBATCH -c 8
#SBATCH -N 1
#SBATCH --mem=400G
#SBATCH -o 2023_01_04_LD_decay_calc_TASSEL_thinned_out_VCF.out
#SBATCH -e 2023_01_04_LD_decay_calc_TASSEL_thinned_out_VCF.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=milaleoniexia@gmail.com

module load openjdk/11.0.2
/usr/users/mtost/software_TASSEL/tassel-5-standalone/run_pipeline.pl -Xms512m -Xmx400g -fork1 -vcf /usr/users/mtost/wd_GBeasy_rerun/Results/2022-12-23_GB10_Shoepeg_strong_filt_and_thinning_for_LD_decay_calc.vcf.gz -ld -ldWinSize 10000 -export 2023_GB10_Shoepeg_after_filtering_thinned_out.ld.win10000.txt
