__author__='pier-luc'

from Bio import SeqIO


def formatFasta(inFile, outFile):
    """
    Sometimes, there are format errror in fasta files. That should correct them
    :param inFile: A fasta files that might contains format error
    :param outFile: A fasta file with no format error
    :return: Nothing
    """
    SeqIO.convert(inFile, "fasta", outFile, "fasta")


def checkFiles(listOfFile):
    for i in listOfFile:
        if not(isValid(i)):
            Error.error("File"+str(i)+"can't be opened")

def isValid(file_name):
    valid = True
    try:
        a = open(file_name, 'r')
        a.close()
    except:
        valid = False
    return valid

def countLines(file_name):
    lines = 0
    try:
        for i in open(file_name):
            lines += 1
    except:
        lines = 0
    if lines == 0:
        Error.error("File "+str(file_name)+" is empty or can not be opened")
    return lines

############################################################################################
'''
Functions that will be relocated/deleted/replaced.
Temporary file
'''


import subprocess
import BiologicalAbundanceContig
import Error
import sys


def testValidity():
    if len(sys.argv) == 3:
        try:
            directory = sys.argv[1]
            numberOfBestMatch = int(sys.argv[2])
        except:
            print(__doc__)
            sys.exit(1)
    elif len(sys.argv) == 4:
        try:
            directory = sys.argv[1]
            pathToContigsIdFile = sys.argv[3]
            tempFile = open(pathToContigsIdFile, 'r')
            tempFile.close()
            numberOfBestMatch = int(sys.argv[2])
        except:
            print(__doc__)
            sys.exit(1)
    else:
        print(__doc__)
        sys.exit(1)
    try:
        tempFile = open(directory+"/_DeNovoAssembly/Contigs.tsv", 'r')
        tempFile.close()
    except:
        Error.error("Unable to open: "+directory+"/_DeNovoAssembly/Contigs.tsv")
        sys.exit(1)

def readPathsFile(pathToContigsIdFile):
    #TODO: verifier que si les fichiers ont ou non le ContigIdentifications.tsv a la fin.
    contigIdentificationsFiles = []
    try:
        tempfile = open(pathToContigsIdFile, 'r')
        tempfile.close()
    except:
        Error.error("Unable to open: "+pathToContigsIdFile)
        sys.exit(1)
    for lines in open(pathToContigsIdFile, 'r'):
        contigIdentificationsFiles.append(lines[0:len(lines)-1])
    return contigIdentificationsFiles

def getPathsFromDirectory(directory):
    command = "find " + directory + "| grep ContigIdentifications.tsv"
    files = subprocess.check_output(command, shell=True)
    contigIdentificationsFiles = files.split('\n')
    for entry in contigIdentificationsFiles:
        if entry == "" or '\n' in entry:
            contigIdentificationsFiles.remove(entry)
    return contigIdentificationsFiles

def readContigsTSVfile(pathToFile):
    directory = pathToFile + "/_DeNovoAssembly/Contigs.tsv"
    biologicalAbundanceContigs = {}
    for lines in open(directory, 'r'):
        if lines == "":
            break
        elif lines[0] == "#":
            continue
        else:
            name = lines.split()[0]
            length = int(lines.split()[2])
            coloredKmer = int(lines.split()[3])
            newBiologicalAbundanceContig = BiologicalAbundanceContig(name,length,coloredKmer)
            #print(name);
            biologicalAbundanceContigs[name] = newBiologicalAbundanceContig
    return biologicalAbundanceContigs

def readContigIdentificationFiles(files, biologicalAbundanceContigs):
    for line in open(files, 'r'):
        if line == "":
            break
        if line[0] == '#':
            continue
        else:
            try:

                contigName = line.split()[0]
                sequenceName = line.split('\t')[6]
                sequenceLength = int(line.split('\t')[7])
                matches = int(line.split('\t')[8])
                biologicalAbundanceContigs[contigName].addNewContigIdentification(sequenceName, sequenceLength, matches)

            except:
                message = ("File: " + files +
                           " can not be read or do not correspond to the usual patern. It has been discarded")
                Error.warning(message)

def showTSV(biologicalAbundanceContigs):
    header=("Contig name"+'\t'+"Contig length in k-mer"+'\t'+"Contig colored k-mer"+'\t'+"Sequence name"+'\t' +
            "Sequence length in k-mer"+"\t"+"Matches in sequence"+'\t'+"pl-value")
    print(header)
    for contigs in biologicalAbundanceContigs:
        biologicalAbundanceContigs[contigs].showTSV()
