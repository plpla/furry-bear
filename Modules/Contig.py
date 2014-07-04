__author__='pier-luc'
__date__ = "2014-06-09"

import Error

class Contig(object):
    """
    Class representing a contig: a reference sequence that may have multiple genes. Its main caracteristics are:
     sequence, sequence length and id.
    """
    def __init__(self, id, sequence):
        self.id = id
        self.sequence = sequence.upper()
        self.sequence_length = len(sequence)
        self.genes = [] #dict or list? Will be used to store genes references if needed.
        self.colored_kmers = 0
        self.length_in_kmers = 0

    def add_gene(self, gene):
        """
        Add a gene to the Contig.
        :param gene: A gene to add to the contig. More things should be checked... (len)
        """
        if type(gene).__name__ != "Gene":
            raise Error.LogicError("Contig:add_gene: Object is not good")
        if gene.contig_id != self.id:
            raise Error.LogicError("Contig:add_gene: Logic error: Not the same contig")
        else:
            self.genes.append(gene)
        #return True

    def get_tsv_for_identification(self):
        line = self.id+"\t"+str(self.length_in_kmers)+"\t"+str(self.colored_kmers)
        return line


if __name__ == "__main__":
    print("Running some test")
    c = Contig("Contig-23", "ACTGATCGACTGACTATCGTACGATCGA")
    print(type(c).__name__)
