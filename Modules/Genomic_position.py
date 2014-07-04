__author__='pier-luc'
__date__ = "2014-06-05"


import Error
from Genomic_environment import Genomic_environment
import sys


class Genomic_position(object):
    #L,heritage would have been preferable... if one day im not too lazy
    def __init__(self):
        self.position = -1
        self.depth = 0
        self.nucl_a = 0
        self.nucl_c = 0
        self.nucl_g = 0
        self.nucl_t = 0
        self.environment = Genomic_environment()

    def set_environment_contig(self, contig_name):
        self.environment.contig = contig_name
        self.environment.position = self.position


    def update_depth(self):
        self.depth = self.nucl_t+self.nucl_g+self.nucl_c+self.nucl_a

    def set_nucl_proportions(self, line):
        self.nucl_a = int(line.split()[0].split(":")[1])
        self.nucl_c = int(line.split()[1].split(":")[1])
        self.nucl_g = int(line.split()[2].split(":")[1])
        self.nucl_t = int(line.split()[3].split(":")[1])
        self.update_depth()

    def fill_from_gatk_line(self, file_line):
        try:
            self.position = int(file_line.split('\t')[0].split(':')[1])
            self.nucl_a = int(file_line.split('\t')[4].split(' ')[0].split(':')[1])
            self.nucl_c = int(file_line.split('\t')[4].split(' ')[1].split(':')[1])
            self.nucl_g = int(file_line.split('\t')[4].split(' ')[2].split(':')[1])
            self.nucl_t = int(file_line.split('\t')[4].split(' ')[3].split(':')[1])
            self.environment.position = self.position
            self.update_depth()
        except:
            sys.stderr.write("fill_from_gatk_line: line is not standard:\n Problematic line: %s\n" % file_line)

    def get_most_abundant(self):
        """
        Will work in most cases except when the most abundant is shared by 2 nucl. In this case, only one is returned.
        TODO: Fix or workaround?
        :return:List ["Nucl", count]
        """
        most_abundant = max([self.nucl_a, self.nucl_c, self.nucl_g, self.nucl_t])
        if most_abundant == self.nucl_a:
            return ["A", self.nucl_a]
        elif most_abundant == self.nucl_t:
            return ["T", self.nucl_t]
        elif most_abundant == self.nucl_c:
            return ["C", self.nucl_c]
        elif most_abundant == self.nucl_g:
            return ["G", self.nucl_g]
        else:
            ArithmeticError("Genomic_position:get_most_abundant: Logic error: Can't find max.")

    def get_mutations_string(self):
        """
        No parsing here. If you want to select mutations, do it somewhere else!
        :return:String with this format: A:  C:  G:  T:"
        """
        proportion_a = str(self.nucl_a/self.depth)
        proportion_c = str(self.nucl_c/self.depth)
        proportion_g = str(self.nucl_g/self.depth)
        proportion_t = str(self.nucl_t/self.depth)
        line = "A:"+proportion_a+" C:"+proportion_c+" G:"+proportion_g+" T:"+proportion_t
        return line
