__author__='pier-luc'

from Bio import SeqIO
import sys

def formatFasta(inFile, outFile):
    output_handle = open(outFile, "w")
    for lines in open(inFile):
        if '>' in lines and lines[0] != '>':
            fasta_sym = lines.index(">")
            nucl = lines.find("nucleotides")
            output_handle.write(lines[0:fasta_sym-1]+"\n")
            output_handle.write(lines[fasta_sym:nucl+len("nucleotides")-1]+'\n')
            output_handle.write(lines[nucl+len("nucleotides"):])
        else:
            output_handle.write(lines)

    long_sequences = [] # Setup an empty list
    for record in SeqIO.parse(open(outFile, "rU"), "fasta"):
        if len(record.seq) > 500 :
            long_sequences.append(record)
 
    output_handle = open(outFile, "w")
    SeqIO.write(long_sequences, output_handle, "fasta")
    output_handle.close()


if __name__ == "__main__":
    formatFasta(sys.argv[1], sys.argv[2])


