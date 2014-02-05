'''
Created on Jan 14, 2014

@author: Nejc Haberman

The scrpt will remove everything up to the end of adapter sequence allowing N number of miss matches
'''
import sys

# method received a genome sequence, position on the genome and the comparing sequence. it returns the number off miss matches between the genome and import sequence 
def Hamming_check_0_or_1(genome, posn, sequence):
    errors = 0
    for i in range(0, sequence.__len__()):
        if genome[posn+i] != sequence[i]:
            errors += 1
    return errors 

# method will keep sequnces that contains adapter allowing 4 missmatches and remove the adapter after
def filter_fasta (fname_in_fa, adapter_seq, fname_out_fa, number_of_missmatches):

    fin = open(fname_in_fa, "rt")
    fout = open(fname_out_fa, "w")
   
    line = fin.readline()
    while line:
        if line[0] == '>':
            id = line
            seq = fin.readline()
            
            errors = Hamming_check_0_or_1(seq, 0, adapter_seq)  #number of missmatches
            if errors <= number_of_missmatches:
                new_seq = seq[adapter_seq.__len__()-1:] #remove adapter
                fout.write(id)
                fout.write(new_seq)
        line = fin.readline()
    fin.close()
    fout.close()


if sys.argv.__len__() == 5:
    fname_in_fa = sys.argv[1]
    fname_in_adapter = sys.argv[2]
    fname_out_fa = sys.argv[3]
    number_of_missmatches = int(sys.argv[4])
    filter_fasta(fname_in_fa, fname_in_adapter, fname_out_fa, number_of_missmatches)
else:
    print("python filter_fasta.py <input_fasta_file> <adapter_sequence> <output_fasta_file> <number_of_missmatches>")

