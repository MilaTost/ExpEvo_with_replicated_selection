#!/bin/bash
#SBATCH -p fat
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mila.tost@agr.uni-goettingen.de
#SBATCH --mem=200G
#SBATCH -t 48:00:00

module load openjdk/11.0.2
module load bcftools/1.8
module load samtools/1.9
module load bwa/0.7.12

set -e
set -x

. /usr/users/mtost/wd_test_GB_easy_rerun/GB_parameters_mpileup_SNPcall_filt.txt

#start GNU parallel to run each region (e.g. chromosome, scaffold) on separate CPU core
parallel  --gnu --max-procs 4 --keep-order "\

## Create list of sorted BAM files
ls -1 Intermediate_files/2.bam_alignments/*.sorted_bam > Intermediate_files/2.bam_alignments/samples_list.txt
##GENERATE PILEUP 
bcftools mpileup --regions {} --output-type z --skip-indels --annotate AD,DP --fasta-ref $REF_GENOME --min-MQ 30 --min-BQ 30  --no-version -b Intermediate_files/2.bam_alignments/samples_list.txt -o Intermediate_files/3.mpileup/mpileup_{}.vcf.gz;\
##CALL SNPS
bcftools call --multiallelic-caller --variants-only --no-version Intermediate_files/3.mpileup/mpileup_{}.vcf.gz | sed -e 's|$(pwd)\/||g' -e 's/Intermediate_files\/2\.bam_alignments\///g' -e  's/\.R.\.fastq.sorted_bam//g'  > Intermediate_files/4.Raw_SNPs/raw_SNPs_{}.vcf;\
" ::: `grep ">" $REF_GENOME | cut -d ' ' -f1 | sed 's/>//g'`

## Merge only chromosome vcf files 
ls -v Intermediate_files/4.Raw_SNPs/*_chr*.vcf > Intermediate_files/4.Raw_SNPs/file_names.txt

##COMBINE SNPS FROM ALL REGIONS
bcftools concat --regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10 -f Intermediate_files/4.Raw_SNPs/file_names.txt > Results/2016_2020_Shoepeg_results_NO_FILT.vcf

##Filter for DP
bcftools filter -i 'DP>160 & DP<5000' Results/2020_Shoepeg_results_NO_FILT.vcf > Results/2020_Shoepeg_results_DP_FILT.vcf
