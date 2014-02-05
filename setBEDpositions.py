'''
Created on Aug 20, 2013

@author: Nejc

The script will set the chromosome and extend bed position by the input.
'''

import sys

def flank_positions(fin_fname, fout_fname, chromosome, right_shift):
    fin = open(fin_fname, "rt")
    fout = open(fout_fname, "w")
    header = fin.readline()
    line = fin.readline()
    while line:
        col = line.rstrip('\n').rsplit('\t')
        pos1 = col[1]
        pos2 = col[2]
        cDNA = col[3]
        strand = col[5]
        pos1 = int(pos1) + int(right_shift)
        pos2 = pos1 + 1
        
        if strand == '-':
            cDNA = '-' + cDNA
            
        fout.write(chromosome + '\t' + str(pos1) + '\t' + str(pos2) + '\t' + cDNA + '\n')
        line = fin.readline()

if sys.argv.__len__() == 5:
    fin_fname = sys.argv[1]
    fout_fname = sys.argv[2]
    chromosome = sys.argv[3]
    right_shift = int(sys.argv[4])
    flank_positions(fin_fname, fout_fname, chromosome, right_shift)
else:
    #print str(sys.argv.__len__())
    print "error:\t4 arguments are needed\n" + '\n' +"example:\t $ python bed_setBEDpositions.py input_fname.bed output_fname.bed chr1 right_shiftNUM"
