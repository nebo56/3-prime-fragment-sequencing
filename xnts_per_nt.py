'''
Created on Jan 16, 2014

@author: Nejc Haberman

Script will add a crosslink number from BED to every nt in the genome. Results will be written into 2 files seperated by strand of the binding
'''
import sys

def load_xnts (fin_name):
    fin = open(fin_name, "rt")
    header = fin.readline()
    line = fin.readline()
    xnts = {}
    while line:
        col = line.rstrip('\n').rsplit('\t')
        pos = int(col[1])
        cDNA = col[3]
        strand = col[5]
        if strand == '-':
            cDNA = strand + cDNA
        xnts.setdefault(pos,{}).setdefault(strand,cDNA)
        line = fin.readline()
    fin.close()
    return xnts


def set_xnts_per_nt (fname_in, fname_out, fname_genome):
    xnts = load_xnts(fname_in)
    fin_genome = open(fname_genome, "rt")

    fout_same = open(fname_out + "_same.tab", "w")
    fout_anti = open(fname_out + "_anti.tab", "w")

    fout_same.write("position" + '\t' + "genome" + '\t' + "cDNA" + '\n')
    fout_anti.write("position" + '\t' + "genome" + '\t' + "cDNA" + '\n')

    header = fin_genome.readline()
    line = fin_genome.readline()
    pos = 0 #position
    while line:
        line = line.rstrip('\n')
        for i in range(0,line.__len__()-1):
            nt = line[i]
            same_xnt_pos = xnts.get(pos)
            if same_xnt_pos != None:
                same_xnt = same_xnt_pos.get('+')
                if same_xnt != None:
                        fout_same.write(str(pos) + '\t' + nt + '\t' + str(same_xnt) + '\n')
                else:
                        fout_same.write(str(pos) + '\t' + nt + '\t' + str(0) + '\n')
            else:
                fout_same.write(str(pos) + '\t' + nt + '\t' + str(0) + '\n')

            anti_xnt_pos = xnts.get(pos)
            if anti_xnt_pos != None:
                anti_xnt = anti_xnt_pos.get('-')
                if anti_xnt != None:
                        fout_anti.write(str(pos) + '\t' + nt + '\t' + str(anti_xnt) + '\n')
                else:
                        fout_anti.write(str(pos) + '\t' + nt + '\t' + str(0) + '\n')
            else:
                fout_anti.write(str(pos) + '\t' + nt + '\t' + str(0) + '\n')
            pos += 1
        line = fin_genome.readline()
        line = line.rstrip('\n')
    fout_same.close()
    fout_anti.close()


            
if sys.argv.__len__() == 4:
    fname_in = sys.argv[1]
    fname_out = sys.argv[2]
    fname_genome = sys.argv[3]
    set_xnts_per_nt (fname_in, fname_out, fname_genome)
else:
    print("python xnts_per_nt.py <input_file> <output_file> <genome_file>")

''' 
set_xnts_per_nt("/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/results/Ule_Nov3-Library1-4/TPI-PTC48/Ule_Nov3_NoIndex_L007_R1_001._CCGG.bed", "/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/results/Ule_Nov3-Library1-4/TPI-PTC48/Ule_Nov3_NoIndex_L007_R1_001._CCGG", "/media/skgthab/storage/UCL/2014.01.14@Niels-3primeFragment/genome-data/TPI-PTC48.fa")
'''
