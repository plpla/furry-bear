#! /usr/bin/python

import sys

import Utility
import Error

#cPickle is pickle but using C instead of pure python. It is fastest ans since there is more than 2M entries in the dict,
#I recommand using cPickle.
try:
    import cPickle as pickle
except:
    import pickle

class GenomesToTaxon():
    def __init__(self):
        self.Converter = {}

    def prepareConverter(self, fileIn):
        Utility.isValid(fileIn)
        sys.stderr.write("counting lines to prepare converter\n")
        numOfLines = Utility.countLines(fileIn)
        sys.stderr.write(str(numOfLines)+" to read\n")
        readed = 0
        for line in open(fileIn):
            genome = int(line.split()[0])
            taxon = int(line.split()[1])
            self.Converter[genome] = taxon
            readed += 1
            if readed % 100000 == 0:
                sys.stderr.write(str(readed)+" lines readed on "+str(numOfLines)+"\n")

    def dumpConverter(self, fileOut):
        sys.stderr.write("Writing converter to file\n")
        f = open(fileOut,"wb")
        pickle.dump(self.Converter,f, protocol=2)
        f.close()
        sys.stderr.write("Converter dumped\n")

    def loadConverter(self, file_name):
        sys.stderr.write("Loading converter from file\n")
        if Utility.isValid(file_name):
            f = open(file_name, "rb")
            self.Converter = pickle.load(f)
            f.close()
            sys.stderr.write("Converted loaded\n")
        else:
            Error.error("Converter.bin can not be opened. You should produce or reproduce it using prepare.py");

    def convertToTaxon(self, genome):
        if self.genomeIsValid(genome):
            return self.Converter[genome]

    def getConverter(self):
        return self.Converter

    def genomeIsValid(self, genome):
        sys.stderr.write("Searching for genome id:"+str(genome)+"\n")
        detect = False
        """
        try:
            if(self.getConverter()[genome]):
                detect=1;
        except:
            detect=0;
        """
        if genome in self.getConverter():
            detect = True
        else:
            detect = False
        return detect


#test main      
if __name__ == "__main__":
    if len(sys.argv) == 2:
        file_name = sys.argv[1]
        converter = GenomesToTaxon()
        converter.prepareConverter(file_name, "Data/Converter.bin")
    newConverter = GenomesToTaxon()
    newConverter.loadConverter("Data/Converter.bin")

