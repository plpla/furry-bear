__author__ = 'Pier-Luc Plante'
__license__ = 'GPL'
__usage__ = ("This script will prepare everything in order to launch alignment with BWA and variant calling with GATK" +
             " on Colosse.\n In order to work, you must have 2 directory containing symbolic links to reads and " +
             "assembly directories.\n Require biopython.\n python PrepareAndAlign.py contigsDir readsDir \n" +
             "assembly_dir=True/False (optional, default=true)")

import os
import sys
from datetime import datetime

from Modules.ExperimentDesign import ExperimentDesign
from Modules.fastaFixer import formatFasta


def prepareFileTable(pathToContigs, pathToReads, have_assembly_dir=True):
    """
    This function shall be replaced as needed. It must create a file structured like this:
    READS_PATH=...
    CONTIGS_PATH=...
    ASSEMBLY_DIR= (is there a directory called "Assembly" after the sample directory un the contigs path.) True/False
    ContigDirectory  ReadsSample1Dir   ReadsSample2Dir ...

    Path shall be absolute?
    Current directory structure; Sample_PXJY in reads ans Sample_PXJY-Assembly in contigs
   TODO: put an exemple in the directory to make sure it is clear...
    """
    cwd=os.getcwd()
    os.chdir(cwd+"/"+pathToContigs)
    contigsDirectoryList = os.listdir(cwd+"/"+pathToContigs)
    os.chdir(cwd+"/"+pathToReads)
    readsDirectoryList = os.listdir(cwd+"/"+pathToReads)
    os.chdir(cwd)
    resultList = {}
    for assembly in contigsDirectoryList:
        sampleName = assembly.split()[0]
        j = sampleName.split('_')[1][3] #i.e. Sample_P3J0
        if j == '0':
            p = sampleName.split('_')[1][1]
            resultList[sampleName] = []
            for reads in readsDirectoryList:
                if reads.split('_')[1][1] == p:
                    resultList[sampleName].append(reads)
    outFile = open("fileTable", 'w')
    outFile.write("READS_PATH="+pathToReads+'\n')
    outFile.write("CONTIGS_PATH="+pathToContigs+'\n')
    if have_assembly_dir:
        outFile.write("ASSEMBLY_DIR=TRUE\n")
    else:
        outFile.write("ASSEMBLY_DIR=FALSE\n")
    for sample in resultList:
        line = sample
        resultList[sample].sort()
        for readsSet in resultList[sample]:
            line += '\t'+readsSet
        outFile.write(line+'\n')
    outFile.close()


def getDate():
    now = datetime.now().strftime("%Y%m%d")
    return now


def writeAlignmentScript(refDir, readsDir, refFileName, sampleName, assemblyName, groupID, assemblyDir):

    cwd = os.getcwd()
    extension = "-BlueTsunamie-"+getDate()
    outFile = open("JobsScripts/"+sampleName+extension+".sh", 'w')
    outFile.write("#!/bin/bash\n")
    outFile.write("#PBS -N "+sampleName+extension+".stdout\n")
    outFile.write("#PBS -o "+sampleName+extension+".stdout\n")
    outFile.write("#PBS -e "+sampleName+extension+".stderr\n")
    outFile.write("#PBS -A "+groupID+"\n")
    outFile.write("#PBS -l walltime=48:00:00\n")
    outFile.write("#PBS -l nodes=1:ppn=8\n")
    outFile.write("#PBS -q default\n")
    outFile.write('cd "${PBS_O_WORKDIR}"\n\n')
    outFile.write("source /rap/nne-790-ab/software/NGS-Pipelines/LoadModules.sh\n\n") #TODO use my own modules
    #Moved this section to the prepareAlignementReferences function.
    #outFile.write("module load nne-790-ab/mugqic\n")
    #outFile.write("module load mugqic/python/2.7.6\n\n")
    #outFile.write("mkdir Contigs_references\n\n")
    #if assemblyDir:
    #    outFile.write("python2.7 fastaFixer.py "+ refDir+assemblyName + '/Assembly/' + refFileName +
    #                  " Contigs_references/"+assemblyName+".fasta\n\n")
    #else:
    #    outFile.write("python2.7 fastaFixer.py "+ refDir+assemblyName + '/' + refFileName +
    #                  " Contigs_references/"+assemblyName+".fasta\n\n")

    outFile.write("./BashScripts/BlueTsunami Contigs_references/"+assemblyName+".fasta " + readsDir + sampleName + " 8 " +
                  sampleName + extension + '\n\n')
    outFile.write("./BashScripts/RunGATK_colosse.sh "+cwd+"/Contigs_references/"+assemblyName+".fasta " + cwd +
                      '/' + sampleName + extension+"/SortedBamAlignments " + sampleName + extension+"-GATK")
    outFile.close()


def prepareAlignmentScripts(design):
    try:
        os.mkdir("JobsScripts")
    except:
        sys.stderr.write("JobsScripts directory already exists. Will overwrite files if there is any colisions\n")
    for sampleGroup in design.design:
        for sample in design.design[sampleGroup]:
            writeAlignmentScript(design.contigs_path, design.reads_path, "Contigs-gt-500.fasta", sample, sampleGroup,
                                 "nne-790-ad", 1)


def createJobLaunchScript():
    fileList=os.listdir(os.getcwd()+"/JobsScripts/")
    outFile=open("JobsScripts/SubmitJobs.sh", 'w')
    outFile.write("#!/bin/bash\n")
    for f in fileList:
        if f!="SubmitJobs.sh":
            outFile.write("msub JobsScripts/"+f+'\n')
    outFile.close()


def prepareAlignmentReferences(design):
    try:
        os.mkdir("Contigs_references")
    except:
        sys.stderr.write("Contigs_references directory already exists. Will overwrite files if there is any colisions\n")
    for sample in design.design:
        if design.assembly_dir==True:
            formatFasta(design.contigs_path+'/'+sample+"/Assembly/Contigs.fasta", "Contigs_references/"+sample+".fasta")
        if design.assembly_dir==False:
            formatFasta(design.contigs_path+'/'+sample+"/Contigs.fasta", "Contigs_references/"+sample+".fasta")
    sys.stderr.write("Alignement references ready\n")


if __name__ == "__main__":
    if len(sys.argv!=3):
        print(__usage__)
        sys.exit(0)
    sys.stderr.write("Preparing the files needed for the alignment\n")
    have_assembly_dir=True
    if len(sys.argv) ==4 and "assembly_dir" in sys.argv[3]:
        if sys.argv[3].split('=')[1] == "True":
            have_assembly_dir = True
        if sys.argv[3].split('=')[1] == "False":
            have_assembly_dir = False
    prepareFileTable(sys.argv[1], sys.argv[2], have_assembly_dir)
    design = ExperimentDesign()
    design.loadDesignFromFile("fileTable")
    prepareAlignmentReferences(design)
    prepareAlignmentScripts(design)
    createJobLaunchScript()
    sys.stderr.write("You can now launch alignment using:\n")
    sys.stderr.write("sh JobsScripts/SubmitJobs.sh \n")


