__author__='pier-luc'
__date__= '2014-06-04'


from Bio.Seq import Seq
from Bio.Alphabet import IUPAC


class Gene(object):
    """
    This class represent a gene entry (in fact it's an orf) as given by prodigal.
    """
    def __init__(self):
        self.name = ""
        self.contig_id = ""
        self.start = 0
        self.end = 0
        self.strand = ''
        self.start_type = ""
        self.metadata = ""
        self.sequence = ""  #TODO:  removed but need to check sequence conformity with postion + getSeq function

    def set_info(self, name, contig_id, start, end, strand, start_type, metadata):
        self.name = name
        self.contig_id = contig_id
        self.start = start
        self.end = end
        self.strand = strand
        self.start_type = start_type
        self.metadata = metadata
        self.sequence = ""

    def get_string(self):
        gene = (self.name+" # " + str(self.start) + " # " + str(self.end) + " # " +self.strand + " # " + self.metadata +
                '\n' + self.sequence)
        #TODO: Output sequence in a fancy format
        return gene

    def translate_sequence(self):
        seq = Seq(self.sequence, IUPAC.unambiguous_dna)
        protein_seq = str(seq.translate(table=11))
        return protein_seq