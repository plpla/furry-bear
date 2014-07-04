__author__ = 'Pier-Luc Plante'
__license__ = 'GPL'


from datetime import datetime
import sys


class ExperimentDesign(object):
    """
    This class contains the experimental design. All needed informations to repeat the experiment is in there
    """
    def __init__(self):
        self.reads_path = ""
        self.contigs_path = ""
        self.design = {}
        self.assembly_dir = bool

    def loadDesignFromFile(self, fileTable = "fileTable"):
        date = datetime.now().strftime("%Y%m%d")
        f = open(fileTable, 'r')
        try:
            line=f.readline()
            if "READS_PATH" in line:
                sys.stderr.write("Reads path is ok\n")
                self.reads_path = line.split('=')[1].rstrip('\n')
            else:
                raise IOError("File is invalid")
            line = f.readline()
            if "CONTIGS_PATH" in line:
                sys.stderr.write("Contigs path is ok\n")
                self.contigs_path = line.split('=')[1].rstrip('\n')
            else:
                raise IOError("File is invalid")
            line = f.readline()
            if "ASSEMBLY_DIR=TRUE\n" == line:
                sys.stderr.write("Assembly directory found\n")
                self.assembly_dir = True
            elif "ASSEMBLY_DIR=FALSE\n" == line:
                self.assembly_dir = False
                sys.stderr.write("Assembly directory not found\n")
            else:
                raise IOError("File is invalid")
        except IOError as e:
            sys.stderr.write("I/O error({0}): {1}".format(e.errno, e.strerror))
            exit(0)
        design = {}
        while 1:
            line = f.readline()
            if line == "":
                break
            else:
                contig = line.split()[0]
                design[contig] = []
                for sample in line.split()[1:]:
                    design[contig].append(sample)
        self.design = design
        f.close()
        return design

