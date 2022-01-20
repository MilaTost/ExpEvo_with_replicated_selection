#!/bin/bash
#SBATCH -p fat
#SBATCH -N 1
#SBATCH -c 16
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
. /usr/users/mtost/wd_GB_easy/GB_eaSy_Shoepag3_parameters.txt


##DEMULTIPLEX READS
if [ -z $RAW_READS_R2 ] #if RAW_READS_R2 variable not set in parameters file, then demultiplex single-end reads
	then
		java -jar $GBSX --Demultiplexer -t $NUM_CORES -f1 $RAW_READS_R1 -i $BARCODES_FILE -gzip true -ca $ADAPTER_SEQ -minsl $MIN_LENGTH -o Intermediate_files/1.Demultiplexed_reads; rm Intermediate_files/1.Demultiplexed_reads/*undetermined*
		
elif [ -n $RAW_READS_R2 ] #if RAW_READS_R2 variable set in parameters file, then demultiplex paired-end reads
	then
		java -jar $GBSX --Demultiplexer -t $NUM_CORES -f1 $RAW_READS_R1 -f2 $RAW_READS_R2 -i $BARCODES_FILE -gzip true -ca $ADAPTER_SEQ -minsl $MIN_LENGTH -o Intermediate_files/1.Demultiplexed_reads; rm Intermediate_files/1.Demultiplexed_reads/*undetermined*
fi

##ALIGN TO REFERENCE
if [ -z $RAW_READS_R2 ] #if RAW_READS_R2 variable not set in parameters file, then align single-end reads to reference genome
	then
		parallel --max-procs $NUM_CORES --keep-order --xapply "bwa mem  $REF_GENOME {} | samtools sort -o Intermediate_files/2.bam_alignments/{/.}.sorted_bam; samtools index Intermediate_files/2.bam_alignments/{/.}.sorted_bam" ::: `ls Intermediate_files/1.Demultiplexed_reads/*.R1.fastq.gz` 
		
elif [ -n $RAW_READS_R2 ] #if RAW_READS_R2 variable set in parameters file, then align paired-end reads to reference genome
	then
		parallel --max-procs $NUM_CORES --keep-order --xapply "bwa mem  $REF_GENOME {1} {2} | samtools sort -o Intermediate_files/2.bam_alignments/{1/.}.sorted_bam; samtools index Intermediate_files/2.bam_alignments/{1/.}.sorted_bam" ::: `ls Intermediate_files/1.Demultiplexed_reads/*.R1.fastq.gz` ::: `ls Intermediate_files/1.Demultiplexed_reads/*.R2.fastq.gz`
fi

