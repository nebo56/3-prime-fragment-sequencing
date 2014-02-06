All the 3' fragment data was mapped with bowtie2 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) alignment software to customize TPI sequence. For the data type  conversions we used fastx-tools (http://hannonlab.cshl.edu/fastx_toolkit/) and samtools (http://samtools.sourceforge.net/). For all the filtering and preprocessing we used custom scripts written in python. Final figures were normalised in R version 3.0.2 (http://www.r-project.org/ ) using ggplot2 (http://ggplot2.org/) package.

Pipeline (starts with 3_prime_fragment-script.sh):
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
 

Scripts:
- 3_prime_fragment-script.sh (Main script)
- filter_fasta.py (The script will remove everything up to the end of adapter sequence allowing N number of miss matches.)
- SAMtoCollapsedSAMandBED.py (The script will read .SAM file and write it to collapsed .SAM file ignoring "4" flag for strand and remove all reads with duplicated barcode for each position. Collapsed data will be also written to  a .BED file format.)
- remove_up_stream.py (The scrpt will remove everything up to the end of adapter sequence )
- setBEDpositions.py (The script will set the chromosome and extend BED positions.)
- swap_barcodes.py (The script will read fasta file and remove random barcode and experimental barcode from fasta. Random barcode will be saved to a new fasta file.)
- xnts_per_nt.py (Script will add a crosslink number from BED to every nt in the genome. Results will be written into 2 files seperated by strand of the binding.)
