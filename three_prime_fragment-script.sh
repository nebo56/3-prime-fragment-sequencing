#!/bin/bash -l
#$ -l h_vmem=16G
#$ -l tmem=16G
#$ -l h_rt=32:0:0
#$ -j y
#$ -S /bin/bash

# set paths
PYTHONPATH=/home/skgthab/programs/Python-2.7.5/bin
PERL5LIB=/home/skgthab/programs/ActivePerl-5.16.3.1603-x86_64-linux-glibc-2.3.5-296746/perl/bin
export PATH=$PYTHONPATH:$PATH
export PATH=$PERL5LIB:$PATH
export PATH=/home/skgthab/programs/tophat-2.0.9.Linux_x86_64:$PATH
export PATH=/home/skgthab/programs/bowtie2-2.1.0:$PATH
export PATH=/home/skgthab/programs/samtools-0.1.19:$PATH
export PATH=/home/skgthab/programs/fastx_toolkit0.0.13:$PATH
export PATH=/home/skgthab/programs/bedtools-2.17.0/bin:$PATH

#set the data, path, bowtie2 index file path, chromosome and genomic start position of the transcript
path=
data=$1
index=$2
scripts=
chromosome=
genomic_start_position=

########################
### 1. preprocessing ###
########################


# remove adapters
fastx_clipper -Q 33 -a AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG -i ${path}${data}.fq > ${path}${data}-clipped.fq

# convert fastq to fasta
fastq_to_fasta -Q 33 -n -i ${path}${data}-clipped.fq -o ${path}${data}.fasta

# delete everything up to the end of GCTGATGGCGATGAATGA
python ${scripts}remove_up_stream.py ${path}${data}.fasta GCTGATGGCGATGAATGA ${path}${data}-filtered.fasta

# remove and swap the random barcode to a new file
python ${scripts}swap_barcodes.py ${path}${data}-filtered.fasta ${path}${data}-noBarcodes.fasta ${path}${data}-Barcodes.fasta

##################
### 2. mapping ###
##################

# map the remaining fasta sequence
bowtie2-align -x /home/skgthab/bowtie-indexes/${index}/${index} -f ${path}${data}-noBarcodes.fasta -S ${path}${data}.sam

####################
### 3. filtering ###
####################

# filter out all reads with mismatches
samtools view -Sh ${path}${data}.sam | grep -e "^@" -e "XM:i:[0][^0-9]" > ${path}${data}-filtered.sam

# SAM to BED with collapsed read count by random barcodes
python ${scripts}SAMtoCollapsedSAMandBED.py ${path}${data}-filtered.sam ${path}${data}-Barcodes.fasta ${path}${data}-collapsed.sam ${path}${data}.bed

# add chromosome and starting genome position to BED
python ${scripts}setBEDpositions.py ${path}${data}.bed ${path}${data}-genome_wide.bed ${chromosome} ${genomic_start_position}

# create table with xnts per nt from the genome
python ${scripts}xnts_per_nt.py ${path}${data}.bed ${path}${data} /home/skgthab/bowtie-indexes/${index}/${index}.fa

# convert mapped reads SAM to BAM and BAI
samtools view -Sb ${path}${data}-collapsed.sam > ${path}${data}-collapsed.bam
samtools sort ${path}${data}-collapsed.bam ${path}${data}-collapsed-sorted
samtools index ${path}${data}-collapsed-sorted.bam ${path}${data}-collapsed-sorted.bam.bai

# clean
rm ${path}${data}-collapsed.sam
rm ${path}${data}-collapsed.bam
rm ${path}${data}-filtered.sam
rm ${path}${data}.sam
rm ${path}${data}-trimmed.fasta
rm ${path}${data}-filtered.fasta
rm ${path}${data}-noBarcodes.fasta
rm ${path}${data}.fasta

####################################
### 4. normalising and filtering ###
####################################

# use "3prime-fragment-plots-binning-filtering.R" script for additinal filtering and normalisation








