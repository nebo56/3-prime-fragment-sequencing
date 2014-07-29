'''
Created on Nov 28, 2013

@author: Nejc Haberman


The script will read fasta file and remove random barcode and AAA linker from fasta. Random barcode will be saved to a new fasta file and reads without AAA linker
will be saved in a separate fasta file.
'''


import sys

def swap_barcodes(fin_fasta, fout_fasta, fout_randomBarcodes):
    finFasta = open(fin_fasta, "rt")
    foutFasta = open(fout_fasta, "w")
    foutFastaNOaaa = open(fout_fasta + "-noAAA.fasta", "w")
    foutRandomBarcode = open(fout_randomBarcodes, "w")
    line = finFasta.readline()
    while line:
        if line[0] != '>':
            randomBarcode = line[0:6]
            aaa = line[6:10]
            seqRead = line[9:]
            if aaa == "AAA" or aaa == "aaa":    #save reads with AAA linker to a different FASTA format
                foutFastaNOaaa.write(id)
                foutFastaNOaaa.write(seqRead)
            else:
                foutFasta.write(id)
                foutRandomBarcode.write(id)
                foutFasta.write(seqRead)
                foutRandomBarcode.write(randomBarcode + '\n')
        else:   #write ID
            col = line.rstrip('\n').rsplit(' ')    #id is in the first column
            id = col[0] + '\n'
        line = finFasta.readline()
    finFasta.close()
    foutFasta.close()
    foutRandomBarcode.close()
    
if sys.argv.__len__() == 4:
    fin_fasta = sys.argv[1]
    fout_fasta = sys.argv[2]
    fout_randomBarcodes = sys.argv[3]
    swap_barcodes(fin_fasta, fout_fasta, fout_randomBarcodes)
else:
    print "you need 3 arguments to run the script"
    quit()
    
                                                                                                                                                                                                                                                               
