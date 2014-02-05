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

#set the data, path, bowtie2 index file path, chromosome of the transcript and a genomic start position of the transcript
data=
path=
index=
chromosome=
genomic_start_position=

# convert fastq to fasta
fastq_to_fasta -Q 33 -n -i ${path}${data}.fq -o ${path}${data}.fasta

# remove and swap the barcode to a new file
python ${path}swap_barcodes.py ${path}${data}.fasta ${path}${data}-noBarcodes.fasta ${path}${data}-Barcodes.fasta

#keep sequnces that contains GAGGCCATGCGTCGACTA allowing 4 missmatches and remove the adapter after
python ${path}filter_fasta.py ${path}${data}-noBarcodes.fasta GAGGCCATGCGTCGACTA ${path}${data}-filtered.fasta 4

# delete everything up to the end of GCACGA
python ${path}remove_up_stream.py ${path}${data}-filtered.fasta GCACGA ${path}${data}-trimmed.fasta

# map the remaining fasta sequence
bowtie2-align -x ${index} -f ${path}${data}-trimmed.fasta -S ${path}${data}.sam

#filter out all reads with mismatches
samtools view -Sh ${path}${data}.sam | grep -e "^@" -e "XM:i:[0][^0-9]" > ${path}${data}-filtered.sam

# SAM to BED with collapsed read count by random barcodes
python ${path}SAMtoCollapsedSAMandBED.py ${path}${data}-filtered.sam ${path}${data}-Barcodes.fasta ${path}${data}-collapsed.sam ${path}${data}.bed

# add chromosome and starting genome position to BED
python ${path}setBEDpositions.py ${path}${data}.bed ${path}${data}-genome_wide.bed ${chromosome} ${genomic_start_position}

# create table with xnts per nt from the genome
python ${path}xnts_per_nt.py ${path}${data}.bed ${path}${data} ${index}.fa

# convert mapped reads SAM to BAM and BAI
samtools view -Sb ${path}${data}-collapsed.sam > ${path}${data}-collapsed.bam
rm ${path}${data}-collapsed.sam
samtools sort ${path}${data}-collapsed.bam ${path}${data}-collapsed-sorted
rm ${path}${data}-collapsed.bam
samtools index ${path}${data}-collapsed-sorted.bam ${path}${data}-collapsed-sorted.bam.bai


