__author__='pier-luc'
__date__ = "2014-06-05"



import argparse
import sys


class Option_parser(object):
    def __init__(self, args):
        self.parser = argparse.ArgumentParser(description="Tool to analyse metagenomic assembly and"+
                                                          " profiling after RayMeta")
        self.subparsers = self.parser.add_subparsers(help='sub-command help')


        self.pipeline1 = self.subparsers.add_parser("pipe1", help="Find SNP located in ORF, select SNP base on " +
                                                                   "parameters, search context in BAM files and outp" +
                                                                   "ut a table")

        self.pipeline1.add_argument("-d", type=int, help="Minimum depth", default=50, required=False)
        self.pipeline1.add_argument("-n", type=int, help="Number of process for hard calculation", default=1,
                                     required=False)
        self.pipeline1.add_argument("-c", type=str,
                                     help="Contigs.fasta file",
                                     required=True)
        self.pipeline1.add_argument("-g", type=str,
                                     help="A tab file of 3 column containing the experiment design." +
                                          "First: Sample name, second: path to GATK file, "+
                                          "third: path to BAM file directory. The first line is the reference",
                                     required=True)
        self.pipeline1.add_argument("-r", type=float,
                                     help="Minimum mutation ratio",
                                     default=0.0, required=False)
        self.pipeline1.add_argument("-o", type=str, help="Output",
                                     default="stdout", required=False)
        self.pipeline1.add_argument("-p", type=str,
                                     help="File containing gene sequences for the contig file (produce by Prodigal",
                                     required=False)
        self.pipeline1.add_argument("-s",
                                     help="Remove synonymous mutations",
                                     type=bool, default=False, required=False)
        self.pipeline1.add_argument("-m",
                                     help="Calculate dn and ds. Will not output a table",
                                     type=bool, default=False, required=False)

        #self.parser.add_argument("")
        self.arguments = vars(self.parser.parse_args(args))

    def get_arguments(self):
        return self.arguments

    def get_parser(self):
        return self.parser


#Main test
if __name__ == "__main__":
    parser = Option_parser(sys.argv[1:])
    arg = parser.get_arguments()
    print(arg)
    print(arg["o"])




