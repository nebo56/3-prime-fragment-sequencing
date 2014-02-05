3-prime-fragment-sequencing
===========================

Precisely determine pattern of 3‘ fragments

All the 3' fragment data was mapped with bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) alignment software to customize TPI sequence. For the data type  conversions we used fastx-tools (http://hannonlab.cshl.edu/fastx_toolkit/) and samtools (http://samtools.sourceforge.net/). For all the filtering and preprocessing we used custom scripts written in python. Final figures were normalised in R version 3.0.2 (http://www.r-project.org/ ) using ggplot2 (http://ggplot2.org/) package.

Pipeline:
1. preprocessing
- select reads that contains  GAGGCCATGCGTCGACTA sequence allowing 4 missmatches
- filter read sequences by deleting everything up to the end of  GCACGA

2. mapping
- map selected reads with bowtie2 mapping aligner on custom TPI genome allowing single hits only

3. filtering
- select reads with no miss matches
- count number of reads uniq random barcode for each position (these reads are parts of the PCR artifacts)
- select reads which are mapped on same strand only (anti sense reads are part of the background noise)

4. normalisation?
- 10 nucleotide binning 
- normalise counts with the number of all reads per sample
- remove the background noise by filtering each sample with wild type (TPI-WT)
 


