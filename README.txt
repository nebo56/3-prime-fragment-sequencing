All the 3' fragment data was mapped with Bowtie2.1 (http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) alignment 
software to customize TPI FASTA sequence. For the data type  conversions we used fastx-tools 
(http://hannonlab.cshl.edu/fastx_toolkit/) and samtools (http://samtools.sourceforge.net/). For the filtering and 
preprocessing we used custom scripts written in Python. Final figures were normalised in R version 3.0.2 
(http://www.r-project.org/).

Pipeline description (3_prime_fragment-script.sh):
1. preprocessing
 - select reads that contains  GAGGCCATGCGTCGACTA sequence allowing 4 missmatches
 - filter read sequences by deleting everything up to the end of  GCACGA
2. mapping
 - map selected reads with bowtie2 mapping aligner on custom TPI genome allowing single hits only
3. filtering
 - select reads with no miss matches
 - count number of reads with uniq random barcode for each position (un-uniq reads are part of the PCR artifacts)
 - select reads which are mapped on the same strand only (anti sense reads are part of background noise)
4. normalising and filtering (R scripts)
 - use 10 and 30 nucleotide binning 
 - normalise counts with the number of all reads per sample
 - remove the background noise by filtering each TPI-mutant sample with TPI wild type (TPI-WT)


Scripts description:
 - 3_prime_fragment-script.sh
	Main script
 - filter_fasta.py
	The script will remove everything up to the end of adapter sequence allowing N number of miss matches.
 - SAMtoCollapsedSAMandBED.py 
	The script will read .SAM file and write it to collapsed .SAM file ignoring "4" flag for strand and remove all 
	reads with duplicated barcode for each position. Collapsed data will be also written to  a .BED file format.
 - remove_up_stream.py
	The scrpt will remove everything up to the end of adapter sequence.
 - setBEDpositions.py
	The script will set the chromosome and extend BED positions.
 - swap_barcodes.py
	The script will read fasta file and remove random barcode and experimental barcode from fasta. Random barcode 
	will be saved to a new fasta file.
 - xnts_per_nt.py
	Script will add a crosslink number from BED to every nt in the genome. Results will be written into 2 files 
	seperated by strand of the binding.
 - 3prime-fragment-plots-binning-filtering.R
	The script will import 3 prime fragments tables and added aditional columns with binned and normalised values. 
	Each one of them will be ploted in 1 nucleotide, 10 nt and 30 nt resolution. Set all paths to "*_same.tab" tables 
	from 3 prime fragments results.

