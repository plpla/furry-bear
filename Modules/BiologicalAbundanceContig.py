#!/usr/bin/python

from Contig import *
from ContigIdentification import *

class BiologicalAbundanceContigObject(Contig):
    "Class representing a contig and the information about it from RayMeta"
    def __init__(self, name="", length_in_kmer = 0, colored_kmers = 0):
        Contig.__init__(self, name, length_in_kmer, colored_kmers)
        self.contigIdentifications = []
        self.LCA_id = ""
        self.LCA_name = ""

    def addNewContigIdentification(self, Name, Length, Matches):
        newContigIdentification = ContigIdentificationObject(Name, Length, Matches)
        self.contigIdentifications.append(newContigIdentification)
    #print("addNewContigIdentification:"+str(len(self.contigIdentifications)))

    def calculatePLvalues(self):
        for entry in self.contigIdentifications:
            entry.calculatePLvalue(Contig.length_in_kmers, Contig.colored_kmers)

    def selectBestIdentifications(self, number_of_best_id):
        #New step: removePL-value of 0;
        non_nul_contig_id = []
        #print("Select before non-nul");
        #print(len(self.contigIdentifications));
        for entry in self.contigIdentifications:
            entry.calculatePLvalue(Contig.length_in_kmers, Contig.colored_kmers)
            #print(entry.getPLvalue())
            #print(ContigObject.getLengthInKmer(self))
            #print(ContigObject.getColoredKmers(self))
            if entry.getPLvalue() > 0:
                non_nul_contig_id.append(entry)
        self.contigIdentifications = non_nul_contig_id
        #print("Select after select non-=null")
        #print(len(self.contigIdentifications))
        if len(self.contigIdentifications) <= 0:
            return
        else:
            self.contigIdentifications.sort(key=lambda x: x.PLvalue, reverse=True)
            #print("Select after sort");
            #print(len(self.contigIdentifications));
            if len(self.contigIdentifications) > number_of_best_id:
                self.contigIdentifications = self.contigIdentifications[0:number_of_best_id]

    def removeNulIdentification(self):
        contig_id_modified = []
        for entry in self.contigIdentifications:
            if entry.getPLvalue() > 0:
                contig_id_modified.append(entry)
        return contig_id_modified

    #	def removeHumanIdentifications(self):
    #		for entry in self.contigIdentifications:
    #			if("Homo sapiens chromosome" in entry.getSequenceName()):

    def showTSV(self):
        #self.calculatePLvalues();
        if len(self.contigIdentifications) != 0:
            #print("Longueur contigIdentifications:"+str(len(self.contigIdentifications)));
            for entry in self.contigIdentifications:
                line = Contig.showTSV(self)
                line = (line+"\t"+entry.getSequenceName()+"\t"+str(entry.getSequenceLengthInKmer())+"\t" +
                        str(entry.getMatchesInContig())+"\t"+str(entry.getPLvalue()))
                print(line)

