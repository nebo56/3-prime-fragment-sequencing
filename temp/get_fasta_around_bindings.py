'''
Created on Feb 17, 2014

@author: Nejc Haberman

Script will return a genome sequence for each binding in the region around crosslink site. If there are 10 cross links at one position the script will retrun 10 sequences around that region.
You also need to specify region in BED file in which you are interested in.
input:    BED (with crosslink positions), genome_sequnce, distance (NUM of nt around cross link), start_position (NUM region), end_position (NUM region)
output:   FASTA (fasta sequences aroud cross links)

special: we will ignore minus strand reads and chromosome because we are using it for custom genome TPI 3'framgents
'''

import sys

# Method will receive a sequence in string, position and distance around position. Method will return a sequnce from position-distance:position+distance
#def get_sequnce(seq, position, distance):


#Method will load a genome fasta sequnce and return it as one string
def load_genome(fname_in_genome):
    fin = open(fname_in_genome, "rt")
    header = fin.readline()
    genome = ""
    line = fin.readline()
    while line:
        genome += line.rstrip('\n').rstrip('\r')
        line = fin.readline()
    genome += line.rstrip('\n').rstrip('\r')
    fin.close()
    return genome
    

def get_fasta_around_bindings (fname_in_bed, fname_in_genome, distance, fasta_out, start_BED_position, end_BED_position):
    genome = load_genome(fname_in_genome)
    fin = open(fname_in_bed, "rt")
    fout = open(fasta_out, "w")
    header = fin.readline()
    line = fin.readline()
    
    while line:
        col = line.rstrip('\n').split('\t')
        position = int(col[1])
        cDNA = int(col[3])  #number of crosslinks
        strand = col[5]
        
        if (strand == '+') and (position >= start_BED_position) and (position <= end_BED_position): #we ignore minus strand and xnts that are out of BED region
            start = position - distance 
            if start < 0: start = 0
            end = position + distance + 1
            if end > (genome.__len__() - 1): end = genome.__len__() - 1
            fasta = genome[start:end]
            for i in range(0,cDNA):
                #fout.write('>' + line + '\n')
                fout.write(fasta + '\n')    
        line = fin.readline()
    fin.close()
    fout.close()
    

if sys.argv.__len__() == 7:
    fname_in_bed = sys.argv[1]
    fname_in_genome = sys.argv[2]
    distance = int(sys.argv[3])
    fasta_out = sys.argv[4]
    start_BED_position = int(sys.argv[5])
    end_BED_position = int(sys.argv[6])
    get_fasta_around_bindings (fname_in_bed, fname_in_genome, distance, fasta_out, start_BED_position, end_BED_position)
else:
    print("python get_fasta_around_bindings.py <input_BED_file> <input_genome_sequence> <NUM_distance> <output_FASTA_file> <start BED position> <end BED position>")



'''
path = "/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/test-fasta-motifs/Ule_Nov3-Library1-4/"
fname_in_bed = path + "Ule_Nov3_NoIndex_L007_R1_001._CCGG.bed"
fname_in_genome = path + "TPI-PTC48.fa"
distance = 20
fasta_out = path + "Ule_Nov3_NoIndex_L007_R1_001._CCGG-0_100-motifs.fasta"
start_BED_position = 164
end_BED_position = 264
get_fasta_around_bindings (fname_in_bed, fname_in_genome, distance, fasta_out, start_BED_position, end_BED_position)
'''
    
