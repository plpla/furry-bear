__author__='pier-luc'
__date__ = "2014-06-09"

import Error
import glob



class Aligned_contig(object):
    """
    Class that contains a Contigs, its alignement (multiple genomic positions).
    """
    def __init__(self, contig):
        if type(contig).__name__ != "Contig":
            Error.error("Aligned_contig: __init__: The object passed is not of type Contig")
        self.contig = contig
        self.comparable_samples = {}  #a list of Genomic_position for each comparable sample
        self.bam_file_lists = {} #Sample name, list of bam_file. The sample names are the same than in comparable_samples

    def set_bam_file_lists(self, file_name):
        for lines in open(file_name, 'rU'):
            sample_name = lines.split()[0]
            bam_dir = lines.split()[2].rstrip('\n')
            bam_files = glob.glob(bam_dir+"/*.bam")
            if len(bam_files) == 0:
                raise IOError("No bam file in directory %s\n" % bam_dir)
            self.bam_file_lists[sample_name] = bam_files






    def add_alignement_position(self, sample, gen_position):
        try:
            if gen_position.position <= 0:
                Error.error("Aligned_contig:add_alignement: The gen_position was not set")
            else:

                self.comparable_sample[sample].append(gen_position)
                self.comparable_sample[sample].sort(key=lambda x: x.position) #We keep the list in order for ordered output.
                #TODO: execute sort only when output is aksed?
        except:
            Error.error("Aligned_contig:set_alignement: Object is not good")

#Not the best approach. It would be faster to load the file from the collection if there is a lot of contig.
    def add_sample_to_compare(self, sample_name, gatk_file_name):
        #this is some bad code. To be re-implemented because it does not do what it should.
        #if not sample_name in self.comparable_samples:
        #    self.comparable_samples[sample_name] = []
        #for line in open(gatk_file_name, 'rU'):
        #    new_genomic_position = Genomic_position()
        #    new_genomic_position.fill_from_gatk_line(line)
        raise NotImplementedError




    def select_mutations(self, sample_name, minimum_mutation_ratio=0.0, minimum_depth=0, maximum_mutation_ratio=1.0,
                         maximum_depth=1000000):
        """
        Will remove Genomic_positions. Compare min and max ratio to the second most abundant nucl.
        :param sample_name: The sample name
        :param minimum_mutation_ratio: If no nucl with a mutation rate over this ratio is present, the postion is deletec
        :param minimum_depth: If the depth is insufficient, position is removed.
        :param maximum_depth_ratio: If the depth ratio is greater, position is removed. U
        :param maximum_depth: If the depth id greater, position is removed. Usefull to remove regions where there
        could be an alignement error
        :return: Nothing
        """
        if minimum_depth>maximum_depth:
            raise Error.LogicError("Aligned_contig: select_mutations: Logic error, minimum depth is greater" +
                                   " than maximum depth")
        new_alignement = []
        #sys.stderr.write("Select_mutation for %s in %s\n" % (sample_name, self.contig.id))
        for position in self.comparable_samples[sample_name]:
            try:
                #print(position)
                if maximum_depth >= position.depth >= minimum_depth:
                    pass
                else:
                    continue
                mutation_str = position.get_mutations_string()
                if len(mutation_str.split()) != 4:
                    raise ValueError("Aligned_contig: select_mutations: get_mutation_string is invalid")
                proportion_a = float(mutation_str.split()[0].split(":")[1])
                proportion_c = float(mutation_str.split()[1].split(":")[1])
                proportion_g = float(mutation_str.split()[2].split(":")[1])
                proportion_t = float(mutation_str.split()[3].split(":")[1])
                proportions = [proportion_a, proportion_c, proportion_t, proportion_g]
                proportions.sort(reverse=True)
                second_most_abundant = proportions[1] #Cuz the most abundant is at pos 0.
            #TODO: Add assert to check if the most abundant is the good one?
                if maximum_mutation_ratio >= second_most_abundant >= minimum_mutation_ratio:
                    new_alignement.append(position)
                else:
                    continue
            except:
                raise TypeError("Aligned_contig: select_mutations: There is a problem with one of the objects")
        self.comparable_samples[sample_name] = new_alignement

    def select_mutations_in_genes(self, gene_list=[]):
        """
        This function will select positions that are in a gene list. Default: all genes in Contig.
        There is no modifications to the contig object.
        :param gene_list: A list contining genes name that are already in the contig. If not provided, all genes in the contigs are used.
        :return: Nothing
        """
        #TODO: there is probably a way to optimise this...
        new_alignement=[]
        for gene in self.contig.genes:
            if len(gene_list) > 0:
                if self.contig.genes[gene].name in gene_list:
                    pass
                else:
                    continue
            start = self.contig.genes[gene].start
            end = self.contig.genes[gene].end
            for position in self.alignement:
                if self.alignement[position].position > end:
                    break
                else:
                    if self.alignement[position].position > start:
                        new_alignement.append(self.alignement[position])
                    else:
                        continue
        self.alignement = new_alignement

    def write_out(self, header=False):
        if header:
            head_line = "Contig_id\tSample_name"
            print (head_line)
        for sample in self.comparable_samples:
            line_start = self.contig.id+"\t"+sample+"\t"
            for compared in self.comparable_samples[sample]:
                line = line_start + compared.get_mutations_string()
                print(line)
                line = compared.environment.output_by_sequence()




