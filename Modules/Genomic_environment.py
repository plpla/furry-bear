__author__='pier-luc'
__date__ = "2014-06-10"

import sys
import copy
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import math
import pysam

import Error
import Gene


class Genomic_environment(object):
    """
    Contains the genomic context of a position. Used to evaluate SNP proportions at the read level.
    Uses pysam but eventually, we could provide an alternative for samtools user that prefer to
    compute pileup outside python.
    """

    def __init__(self):
        self.contig = ""    #redundant info but usefull for parsing bam-file
        self.position = -1    #Biological positon on contig, not Pythonesque. IOW: starting at 1 not 0
        self.reads = {}   #reads dict where key is read id and value a dict of position with base call


    def populate(self, sam_file_name, minimum_alignment_score):
        if self.contig == "":
            RuntimeError("contig must be set before reading a bam file")
        if self.contig[0]==">":
            current_contig_to_analyse = self.contig.lstrip('>') #Necessary because there is no ">" in the bam file...
        else:
            current_contig_to_analyse = self.contig
        sys.stderr.write("Loading file %s\n" %sam_file_name)
        samfile = pysam.Samfile(sam_file_name, 'rb')
        if not samfile._hasIndex(): #if no index, we must build it
            samfile.close()
            sys.stderr.write("Building index for %s\n" % sam_file_name)
            pysam.index(sam_file_name)
            samfile = pysam.Samfile(sam_file_name, 'rb')
        if self.position-3 < 0:
            sys.stderr.write("%s position %s. I have problem computing this position\n" % (self.contig, self.position))
        for pileup_data in samfile.pileup(current_contig_to_analyse, max([0,self.position-3]), self.position+1):
            #print(str(self.position-3)+" "+str(pileup_data.pos)+" "+str(self.position+1))
            if self.position-3 <= pileup_data.pos <= self.position+1:
                #print('in')
                for pileup_read in pileup_data.pileups:
                    if not pileup_read.alignment.qname in self.reads:
                        self.reads[pileup_read.alignment.qname] = {}
                    if ord(pileup_read.alignment.qual[pileup_read.qpos])-33 > minimum_alignment_score:
                            self.reads[pileup_read.alignment.qname][int(pileup_data.pos+1)] = \
                                pileup_read.alignment.seq[pileup_read.qpos] #using biological position, not python.
        samfile.close()

    def group_environment(self):
        """
        Group environment by their sequences
        :return:A dictionary with sequences as key and number of occurences as value
        """
        sequences = {}
        for read in self.reads:
            seq = ""
            for pos in self.reads[read]:
                #print("Position:%s" %pos)
                seq += self.reads[read][pos]
            if not seq in sequences:
                sequences[seq] = 1
            else:
                sequences[seq] += 1
        return sequences

    def translate_all_frames(self):
        """
        Will translate every possible amino acid that are in the environment.
        :return: A tsv table to stdout
        """
        nucleotides_sequences = self.group_environment()
        total = 0
        for entry in nucleotides_sequences:
            total += nucleotides_sequences[entry]
        if total != len(self.reads):
            Error.LogicError("Genomic_environment: translate_all_frames: Total is different" +
                             " than depth")
        frames = []
        for i in range(0, 6):
            new_frame = {}
            frames.append(new_frame)
        for entry in nucleotides_sequences:
            #Kinda fuzzy... we must find the amino acid in all frames.
            #We will have a list of frame where each fram is a dict stored with aa.
            if len(entry) != 5:
                Error.error("Length of sequence: %s is not 5" % entry)
            #I dont see a way to avoid hardcoding :(
            #Frame 1
            print(entry)
            print(entry[0:3])
            #print(len(entry[0:3]))
            codon = Seq(entry[0:3], IUPAC.unambiguous_dna)
            aa = str(codon.translate(table=11))
            if not aa in frames[0]:
                frames[0][aa] = 0
            frames[0][aa] += nucleotides_sequences[entry]
            #frame 2
            print(entry[1:4])
            codon = Seq(entry[1:4], IUPAC.unambiguous_dna)
            aa = str(codon.translate(table=11))
            if not aa in frames[1]:
                frames[1][aa] = 0
            frames[1][aa] += nucleotides_sequences[entry]
            #Frame 3
            print(entry[2:5])
            codon = Seq(entry[2:5], IUPAC.unambiguous_dna)
            aa = str(codon.translate(table=11))
            if not aa in frames[2]:
                frames[2][aa] = 0
            frames[2][aa] += nucleotides_sequences[entry]

            reverse_entry = entry[::-1]

            #Frame 4 (reverse starting from right)
            print(reverse_entry[0:3])
            codon = Seq(reverse_entry[0:3], IUPAC.unambiguous_dna)
            aa = str(codon.translate(table=11))
            if not aa in frames[3]:
                frames[3][aa] = 0
            frames[3][aa] += nucleotides_sequences[entry]
            #Frame 5
            print(reverse_entry[1:4])
            codon = Seq(reverse_entry[1:4], IUPAC.unambiguous_dna)
            aa = str(codon.translate(table=11))
            if not aa in frames[4]:
                frames[4][aa] = 0
            frames[4][aa] += nucleotides_sequences[entry]
            #Frame 6
            print(reverse_entry[2:5])
            codon = Seq(reverse_entry[2:5], IUPAC.unambiguous_dna)
            aa = str(codon.translate(table=11))
            if not aa in frames[5]:
                frames[5][aa] = 0
            frames[5][aa] += nucleotides_sequences[entry]
        ##We output this in a cute way
        header = "Contig\tSNP position(contig)\tFrame\tAmino acid\tProportion\n"
        print(header),
        for frame in range(0, 6):
            frame_number = str(frame+1)
            for amino_acid in frames[frame]:
                line = ""
                line = self.contig+"\t"+str(self.position)+'\t'+frame_number+"\t"+amino_acid+"\t"+str(float(frames[frame][amino_acid]/total))
                print(line)


    def translate_using_sequence(self, contig):
        """
        Translate the SNPs at this position and give the proportion.
        :param contig: A Contig object with associated Gene(s)
        :return:A line to stdout
        """
        if type(contig).__name__ != "Contig":
            Error.error("Genomic_environment: translate_using_sequence: The object passed is not of type Contig")
        #For each gene in Contig... we search for the one that contain the position
        gene = Gene.__init__()
        for orf in contig.genes:
            if orf.start <= self.position <= orf.end:
                gene = orf
        if gene.sequence == "":
            Error.error("Genomic_environment: translate_using_sequence: There is no sequence in the selected gene")
        orf_seq = Seq(gene.sequence, IUPAC.unambiguous_dna)
        amino_acid_sequence = str(orf_seq.translate(table=11))
        orf_seq = orf_seq.tomutable()
        nucleotides_sequences = self.group_environment()
        total = len(self.reads)
        header = "Contig\tSNP position(contig)\tFrame\tAmino acid\tProportion\n"
        #We put - in the frame number column.
        print(header),
        for entry in nucleotides_sequences:
            mutation = entry[2]
            mutated_seq = copy.copy(orf_seq)
            mutated_seq[self.position-gene.start] = mutation
            mutated_prot_seq = str(mutated_seq.translate(table=11))
            changed_amino_acid_position = math.ceil(float(self.position-gene.start+1)/3)
            changed_amino_acid = mutated_prot_seq[changed_amino_acid_position-1]    #-1 cuz its a python string.
            line = (self.contig+"\t"+self.position+"\t-\t"+changed_amino_acid+"\t" +
                    str(float(nucleotides_sequences[entry])/total))
            print(line)

##Output problem: dict entries are not sorted... so seq is not good!
    def output_by_reads(self):
        for read in self.reads:
            seq = ""
            current_position = self.position-2  #At this point we are in biological positions
            while current_position <= self.position+2:
                try:
                    seq += self.reads[read][current_position]
                except KeyError:
                    seq+= "N"   #We had an N if there is no nucleotide at the required position... it a test.
                current_position += 1
            print("%s\t%s\t%s" % (self.position, read, seq))
##TODO: might give a bad sequence because dict are not sorted.
    def output_by_sequence(self):
        sequences = self.group_environment()
        for entry in sequences:
            print("%s\t%s\t%s" % (self.position, entry, sequences[entry]))




#TODO: move this to a test. With multiple file... checking if coverage and basecount is similar to the one
#provided by GATK
if __name__ == "__main__":
    print("hi! lets run some test!")
    position = Genomic_environment()
    position.position = 166
    position.contig="contig-23"
    position.populate("/rap/nne-790-ab/projects/pplante/Metagenome_snp_pipeline_2014-06-20/P1J0_SortedBamAlignments/" +
                      "P1J0_Lane1_R2_1.fastq.bz2.sam.gz.sorted.bam", 3)
    position.populate("/rap/nne-790-ab/projects/pplante/Metagenome_snp_pipeline_2014-06-20/P1J0_SortedBamAlignments/" +
                      "P1J0_Lane1_R1_1.fastq.bz2.sam.gz.sorted.bam", 3)
    position.populate("/rap/nne-790-ab/projects/pplante/Metagenome_snp_pipeline_2014-06-20/P1J0_SortedBamAlignments/" +
                      "P1J0_Lane1_R2_2.fastq.bz2.sam.gz.sorted.bam", 3)
    position.populate("/rap/nne-790-ab/projects/pplante/Metagenome_snp_pipeline_2014-06-20/P1J0_SortedBamAlignments/" +
                      "P1J0_Lane1_R1_2.fastq.bz2.sam.gz.sorted.bam", 3)
    position.populate("/rap/nne-790-ab/projects/pplante/Metagenome_snp_pipeline_2014-06-20/P1J0_SortedBamAlignments/" +
                      "P1J0_Lane1_R2_3.fastq.bz2.sam.gz.sorted.bam", 3)
    position.populate("/rap/nne-790-ab/projects/pplante/Metagenome_snp_pipeline_2014-06-20/P1J0_SortedBamAlignments/" +
                      "P1J0_Lane1_R1_3.fastq.bz2.sam.gz.sorted.bam", 3)

    position.populate("/rap/nne-790-ab/projects/pplante/Metagenome_snp_pipeline_2014-06-20/P1J0_SortedBamAlignments/" +
                      "P1J0_Lane1_R1_4.fastq.bz2.sam.gz.sorted.bam", 3)
    position.populate("/rap/nne-790-ab/projects/pplante/Metagenome_snp_pipeline_2014-06-20/P1J0_SortedBamAlignments/" +
                      "P1J0_Lane1_R2_4.fastq.bz2.sam.gz.sorted.bam", 3)

    position.populate("/rap/nne-790-ab/projects/pplante/Metagenome_snp_pipeline_2014-06-20/P1J0_SortedBamAlignments/" +
                      "P1J0_Lane1_R1_5.fastq.bz2.sam.gz.sorted.bam", 3)
    position.populate("/rap/nne-790-ab/projects/pplante/Metagenome_snp_pipeline_2014-06-20/P1J0_SortedBamAlignments/" +
                      "P1J0_Lane1_R2_5.fastq.bz2.sam.gz.sorted.bam", 3)

    position.populate("/rap/nne-790-ab/projects/pplante/Metagenome_snp_pipeline_2014-06-20/P1J0_SortedBamAlignments/" +
                      "P1J0_Lane1_R1_6.fastq.bz2.sam.gz.sorted.bam", 3)
    position.populate("/rap/nne-790-ab/projects/pplante/Metagenome_snp_pipeline_2014-06-20/P1J0_SortedBamAlignments/" +
                      "P1J0_Lane1_R2_6.fastq.bz2.sam.gz.sorted.bam", 3)

    position.populate("/rap/nne-790-ab/projects/pplante/Metagenome_snp_pipeline_2014-06-20/P1J0_SortedBamAlignments/" +
                      "P1J0_Lane1_R1_7.fastq.bz2.sam.gz.sorted.bam", 3)
    position.populate("/rap/nne-790-ab/projects/pplante/Metagenome_snp_pipeline_2014-06-20/P1J0_SortedBamAlignments/" +
                      "P1J0_Lane1_R2_7.fastq.bz2.sam.gz.sorted.bam", 3)

    position.populate("/rap/nne-790-ab/projects/pplante/Metagenome_snp_pipeline_2014-06-20/P1J0_SortedBamAlignments/" +
                      "P1J0_Lane1_R1_8.fastq.bz2.sam.gz.sorted.bam", 3)
    position.populate("/rap/nne-790-ab/projects/pplante/Metagenome_snp_pipeline_2014-06-20/P1J0_SortedBamAlignments/" +
                      "P1J0_Lane1_R2_8.fastq.bz2.sam.gz.sorted.bam", 3)

    position.populate("/rap/nne-790-ab/projects/pplante/Metagenome_snp_pipeline_2014-06-20/P1J0_SortedBamAlignments/" +
                      "P1J0_Lane1_R1_9.fastq.bz2.sam.gz.sorted.bam", 3)
    position.populate("/rap/nne-790-ab/projects/pplante/Metagenome_snp_pipeline_2014-06-20/P1J0_SortedBamAlignments/" +
                      "P1J0_Lane1_R2_9.fastq.bz2.sam.gz.sorted.bam", 3)

    position.populate("/rap/nne-790-ab/projects/pplante/Metagenome_snp_pipeline_2014-06-20/P1J0_SortedBamAlignments/" +
                      "P1J0_Lane1_R1_10.fastq.bz2.sam.gz.sorted.bam", 3)
    position.populate("/rap/nne-790-ab/projects/pplante/Metagenome_snp_pipeline_2014-06-20/P1J0_SortedBamAlignments/" +
                      "P1J0_Lane1_R2_10.fastq.bz2.sam.gz.sorted.bam", 3)

    position.populate("/rap/nne-790-ab/projects/pplante/Metagenome_snp_pipeline_2014-06-20/P1J0_SortedBamAlignments/" +
                      "P1J0_Lane1_R1_11.fastq.bz2.sam.gz.sorted.bam", 3)
    position.populate("/rap/nne-790-ab/projects/pplante/Metagenome_snp_pipeline_2014-06-20/P1J0_SortedBamAlignments/" +
                      "P1J0_Lane1_R2_11.fastq.bz2.sam.gz.sorted.bam", 3)

    print("Number of reads: %s" % len(position.reads))
    #for i in position.reads:
    #    print(i)
    #    for a in position.reads[i]:
    #        print("\t Pos: %s\t Nucl %s" % (a, position.reads[i][a]))
    position.output_by_reads()
    position.output_by_sequence()
    position.translate_all_frames()


















