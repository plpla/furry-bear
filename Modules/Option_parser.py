__author__='pier-luc'
__date__ = "2014-06-05"



import argparse
import sys


class Option_parser(object):
    def __init__(self, args):
        self.parser = argparse.ArgumentParser(description="Tool to analyse metagenomic assembly and"+
                                                          " profiling after RayMeta")
        self.subparsers = self.parser.add_subparsers(help='sub-command help')
        self.parser_prepare = self.subparsers.add_parser('prepare',
                                                         help=('Prepare data for the tree construction.' +
                                                         "Requiered to be run at least one time"))
        self.parser_find_contig_id = self.subparsers.add_parser('identify', help='Find the last common ancester')

        self.parser_snp = self.subparsers.add_parser("snp", help="Find and select SNPs in metagenome.")

        self.parser_prepare.add_argument('-t', type=str,
                                         help='The TreeOfLife-Edges.tsv', default=None, required=True)
        self.parser_prepare.add_argument('-n', type=str,
                                         help='The Taxon-Names.tsv file', default=None, required=True)
        self.parser_prepare.add_argument('-f', type=str,
                                         help='The GenomeToTaxon.tsv file', default=None, required=True)

        self.parser_find_contig_id.add_argument('-d', type=str,
                                                help='The biological abundance directory', required=True)
        self.parser_find_contig_id.add_argument('-c', type=open,
                                                help='The contig file (.fasta)', required=True)
        self.parser_find_contig_id.add_argument('-i', type=open,
                                                help='File containing a list of contigIdentification file',
                                                required=False)
        self.parser_find_contig_id.add_argument('-b', type=int,
                                                help='Maximum number of best match to consider. Impact compute time' +
                                                     "(default 10)",
                                                default=10,  required=False)

        self.parser_snp.add_argument("-d", type=int, help="Minimum depth", default=50, required=False)
        self.parser_snp.add_argument("-n", type=int, help="Number of process for hard calculation", default=2,
                                     required=False)
        self.parser_snp.add_argument("-c", type=str,
                                     help="Contigs.fasta file",
                                     required=True)
        self.parser_snp.add_argument("-g", type=str,
                                     help="A tab file of 3 column containing the experiment design." +
                                          "First: Sample name, second: path to GATK file, "+
                                          "third: path to BAM file directory. The first line is the reference",
                                     required=True)
        self.parser_snp.add_argument("-r", type=float,
                                     help="Minimum mutation ratio",
                                     default=0.0, required=False)
        self.parser_snp.add_argument("-o", type=str, help="Output",
                                     default="stdout", required=False)
        self.parser_snp.add_argument("-p", type=str,
                                     help="File containing gene sequences for the contig file (produce by Prodigal",
                                     required=False)
        self.parser_snp.add_argument("-s",
                                     help="Remove synonymous mutations",
                                     type=bool, default=False, required=False)
        self.parser_snp.add_argument("-m",
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




