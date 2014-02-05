'''
Created on Jan 14, 2014

@author: Nejc Haberman

The scrpt will remove everything up to the end of adapter sequence 
'''
import sys

# method will filter sequence up to the end of the adaptor and save it as fasta
def filter_fasta (fname_in_fa, adapter_seq, fname_out_fa):

    fin = open(fname_in_fa, "rt")
    fout = open(fname_out_fa, "w")
   
    line = fin.readline()
    while line:
        if line[0] == '>':
            id = line
            seq = fin.readline()

            pos = seq.find(adapter_seq)
            if pos != -1:
                new_seq = seq[pos + adapter_seq.__len__():]
                fout.write(id)
                fout.write(new_seq)
        line = fin.readline()
    fin.close()
    fout.close()


if sys.argv.__len__() == 4:
    fname_in_fa = sys.argv[1]
    fname_in_adapter = sys.argv[2]
    fname_out_fa = sys.argv[3]
    filter_fasta(fname_in_fa, fname_in_adapter, fname_out_fa)
else:
    print("python remove_up_stream.py <input_fasta_file> <adapter_sequence> <output_fasta_file>")

