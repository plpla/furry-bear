__author__='pier-luc'
__date__ = '2014-06-04'

import sys
import os
from multiprocessing import Process, Queue, Manager
from time import sleep
import copy
import glob

from Modules.Gene import Gene
from Modules.Contig import Contig
from Modules.Option_parser import Option_parser
from Modules.Aligned_contig import Aligned_contig
from Modules.Genomic_position import Genomic_position
from Modules.Genomic_environment import Genomic_environment

from Modules.Utility import isValid


"""
I see this program as multiple function that can be pluged together to obtain the desired result.
Let's see if it can work.
Analysis will probably be an iterative process (execute this program with diverse options each time).
"""

def gatk_reader(file_name, num_of_parser, chunk_queue):
    """
    My job is to read, not to lead!
    Currently losing the last contig of a file...
    """
    #global contig_queue
    contig_chunk = ""
    current_contig_name = ""
    previous_contig_name = ""
    number_of_contig = 0
    for lines in open(file_name, 'rU'):
        if "Locus" in lines:
            continue
        current_contig_name = lines.split("\t")[0].split(":")[0]
        if current_contig_name != previous_contig_name and previous_contig_name != "":
            contig_chunk = contig_chunk.rstrip("\n")   #Remove the last eof. Otherwise, a blank line is inserted and f*** my parser.
            chunk_queue.put(contig_chunk)
            number_of_contig += 1
            sys.stderr.write("Reader %s: %s (# %s) was read and sent to queue.\n" % (os.getpid(), previous_contig_name,
                                                                                     number_of_contig))
            contig_chunk = ""
        contig_chunk += lines
        previous_contig_name = current_contig_name
    #We add the last chunk to the queue... (It should solve the bug)
    contig_chunk = contig_chunk.rstrip("\n")
    chunk_queue.put(contig_chunk)
    number_of_contig += 1
    sys.stderr.write("Reader %s: %s (# %s) was read and sent to queue.\n" % (os.getpid(), current_contig_name,
                                                                             number_of_contig))
    #we add a stop signal for every process in the queue
    for i in range(num_of_parser):
        chunk_queue.put("kill")
    chunk_queue.close()
    sys.stderr.write("Reader %s done\n" % os.getpid())
    #return


def gatk_parser(sample_name, arguments, contig_collection, chunk_queue, result_queue):

    sys.stderr.write("Thread %s : started!\n" % os.getpid())
    min_depth = arguments["d"]
    min_mutation_rate = arguments["r"]
    while True:

        chunk = chunk_queue.get()
        if chunk == "kill":
            sys.stderr.write("Thread %s : poison pill\n" % os.getpid())
            break
        sys.stderr.write("Thread %s : analyzing a chunk of data containing %s lines \n" %
                         (os.getpid(), len(chunk.split("\n"))))
        genomic_position_array = []
        contig_name = ""
        for lines in chunk.split("\n"):
            contig_name = ">"+lines.split("\t")[0].split(":")[0]
            new_genomic_position = Genomic_position()
            new_genomic_position.fill_from_gatk_line(lines)
            genomic_position_array.append(new_genomic_position)
        if contig_name in contig_collection:
            new_aligned_contig = Aligned_contig(contig_collection[contig_name])
            new_aligned_contig.comparable_samples[sample_name] = genomic_position_array
            new_aligned_contig.select_mutations(sample_name=sample_name, minimum_mutation_ratio=min_mutation_rate,
                                                minimum_depth= min_depth, maximum_mutation_ratio= 1.0,
                                                maximum_depth=1000000)

            result_queue.put(new_aligned_contig)
            ####################################################
            #Problem here. The process are not ending.... it takes too much time...
            #######################################
        else:
            sys.stderr.write("Thread %s : %s not added  to list\n" % (os.getpid(), contig_name))
    result_queue.put("kill")


def gatk_pipe_cleaner(result_queue, aligned_contig_dict, num_of_parser):
    sys.stderr.write("Pipe cleaner %s started\n" % os.getpid())
    num_of_kill = 0
    while True:
        new_aligned_contig = result_queue.get()
        if new_aligned_contig == "kill":
            num_of_kill += 1
            sys.stderr.write("Pipe cleaner %s poison %s\n" % (os.getpid(), num_of_kill))
            if num_of_kill == num_of_parser:
                sys.stderr.write("Pipe cleaner %s dying\n" % os.getpid())
                break
            else:
                continue
        #sys.stderr.write("Pipe cleaner: merging %s\n" % new_aligned_contig.contig.id)
        
        if not new_aligned_contig.contig.id in aligned_contig_dict:
            aligned_contig_dict[new_aligned_contig.contig.id] = new_aligned_contig
        else:
            for sample in new_aligned_contig.comparable_samples:
                aligned_contig_dict[new_aligned_contig.contig.id].comparable_samples[sample] = \
                    new_aligned_contig.comparable_samples[sample]


class GATK_handler:
    """
    source: http://ikharn.blogspot.ca/2013/04/producer-consumer-with-python.html       and
    http://stackoverflow.com/questions/15859374/python-producer-consumer-with-process-and-pool
    """
    def __init__(self, num_of_proc):
        self.NUM_OF_PROCESS = num_of_proc   #r# eadonly
        self.manag = Manager()
        self.chunk_queue = Queue(1000)  #TODO: optimize for Colosse
        self.result_queue = Queue()
        self.aligned_contig_dict = self.manag.dict()   #For the pipe cleaner

    def start(self, file_name, sample_name, arguments, contig_collection):
        self.p = Process(target=gatk_reader,  args=(file_name, self.NUM_OF_PROCESS, self.chunk_queue))
        self.p.start()
        self.parsers = [Process(target=gatk_parser, args=(sample_name, arguments, contig_collection, self.chunk_queue,
                                                          self.result_queue))
                        for i in xrange(self.NUM_OF_PROCESS)]
        for proc in self.parsers:
            proc.start()
        self.pipe_cleaner = Process(target= gatk_pipe_cleaner, args=(self.result_queue, self.aligned_contig_dict,
                                                                    self.NUM_OF_PROCESS))
        self.pipe_cleaner.start()

    def check_status(self):
        """
        Function created for debogue.
        TODO: implement as an option?
        """
        while 1:
            for i in self.parsers:
                sys.stderr.write("%s is alive : %s\n" % (i.pid, i.is_alive()))
                sleep(5)
            if not i.is_alive():
                sys.stderr.write("%s is alive : %s dying!\n" % (i.pid, i.is_alive()))
                break

    def join(self):
        sys.stderr.write("MANAGER: %s: joining process\n" % os.getpid())
        self.p.join()
        for i in self.parsers:
            i.join()
        self.pipe_cleaner.join()


def load_genes(file_name):
    gene_collection = {}
    current_gene_id = ""
    for line in open(file_name, 'rU'):
        if line[0] == '>':
            current_gene_id = line.split(' # ')[0]
            #sys.stderr.write(current_gene_id+"\n")
            contig_id = current_gene_id.split("_")[0]
            start = int(line.split(' # ')[1])
            end = int(line.split(' # ')[2])
            strand = int(line.split(" # ")[3])
            if strand == 1:
                strand = '+'
            if strand == -1:
                strand = '-'
            start_type = line.split(' # ')[4].split(';')[2].split('=')[1]
            metadata = line.split(' # ')[4].rstrip('\n')
            gene_collection[current_gene_id] = Gene()
            gene_collection[current_gene_id].set_info(current_gene_id, contig_id, start, end, strand,
                                                      start_type, metadata)
        else:
            gene_collection[current_gene_id].sequence += line.rstrip('\n')
    return gene_collection

def parse_list(file_name):
    item_list = []
    for line in open(file_name, 'rU'):
        item_list.append(line.rstrip('\n'))
    return item_list

def load_contigs_dot_fasta(file_name):
    contig_collection = {}
    header = ""
    sequence = ""
    for lines in open(file_name, 'rU'):
        if lines[0] == ">":
            if not sequence == "":#if sequence is not empty...
                contig_collection[header] = Contig(header, sequence)
                sequence = ""
            header = lines.split()[0]                                           
        else:
            sequence += lines.rstrip('\n')
    return contig_collection

def associate_genes_to_contigs(contig_collection, gene_collection):
    for gene in gene_collection.itervalues():
        try:
            contig_collection[gene.contig_id].genes.append(gene)
        except:
            warning_message = str("Gene %s was not associated to a contig. Most probably because %s does not exist\n" %
                                  (gene.name, gene.contig_id))
            sys.stderr.write(warning_message)
    return contig_collection


def add_genomic_position_to_aligned_contig(arguments, aligned_contig_collection):
    """
    This function add the information of the GATK file based on the desired parameter.
    :param arguments: the program argument dict as obtain from Option_parser
    :param aligned_contig_collection: The aligned contig dict before adding GATK info (contain only contigs)
    :return: A dictionary of Aligned contig filled with Genomic_position selected based on the criterion in arguments
    """
    contig_treated = 0
    #"contig-23:1     2       2.00    2       A:2 C:0 G:0 T:0 N:0 "
    for sample in open(arguments["g"], 'rU'):
        sample_id = sample.split()[0]
        print("Sample in treatment: %s" % sample)
        #Todo: check if GATK file can be openned.
        previous_contig_name = ""
        min_depth = arguments["d"]
        min_mutation_rate = arguments["r"]
        untreated_contig_id =[]

        for gatk_line in open(sample.split('\t')[1], 'rU'):
            if "Locus" in gatk_line:
                continue    #Cuz the header
            current_contig_name = ">"+gatk_line.split("\t")[0].split(":")[0]
            if current_contig_name in untreated_contig_id:
                continue
            if current_contig_name != previous_contig_name and previous_contig_name != "":
                aligned_contig_collection[previous_contig_name].select_mutations(sample_id, min_mutation_rate,
                                                                                 min_depth, 1.0, 1000000)
                contig_treated += 1
                if contig_treated % 1000 == 0:
                    sys.stderr.write("Currently have treated more than %s contigs\n" % contig_treated)
            if not current_contig_name in aligned_contig_collection:
                sys.stderr.write("%s not treated for genomic position." % current_contig_name)
                untreated_contig_id.append(current_contig_name)
                continue
            if not sample_id in aligned_contig_collection[current_contig_name].comparable_samples:
                aligned_contig_collection[current_contig_name].comparable_samples[sample_id] = []
                if len(aligned_contig_collection[current_contig_name].comparable_samples) >=3:
                    print(aligned_contig_collection[current_contig_name].comparable_samples)
            new_genomic_position = Genomic_position()
            new_genomic_position.fill_from_gatk_line(gatk_line)
            aligned_contig_collection[current_contig_name].comparable_samples[sample_id].append(new_genomic_position)
            previous_contig_name = current_contig_name
            #High error probability due to complexity of command...
    #TODO: select samples according to parameters

    #for aligned_contig in aligned_contig_collection:
    #    for genomic_position in aligned_contig.comparable_samples:
     #       aligned_contig.comparable_samples[genomic_position].select_mutations(min_mutation_rate,
      #                                                                           min_depth, 1.0, 1000000)
                                                                                 #Hardcoded... not good!
    return aligned_contig_collection


def get_bam_file_list(file_name, line_number):
    """
    not needed anymore
    """
    raise NotImplementedError("get_bam_file_list: should not be used anymore." +
                              " Store bam_files list in Aligned_contig instead\n")
    #bam_list = []
    #current_line = 0
    #bam_dir = ""
    #for line in open(file_name, 'rU'):
        #print("Line number: %s and current line: %s" % (line_number, current_line))
    #    if current_line == line_number:
    #        bam_dir = line.split()[2].rstrip('\n')
    #        break #to save time...
    #    current_line += 1
    #if bam_dir == "":
    #    raise ValueError("Bam_dir is empty because we could not find the good sample")
    #bam_list = glob.glob(bam_dir+"/*.bam")
    #new_list = []
    #for i in bam_list:
    #    new_list.append(bam_dir+"/"+i)
    #return bam_list


def get_genomic_environment_from_bam(arguments, aligned_contig_collection):
    #1. set contig id in genomic env.
    pass
    """
    for aligned_contig in aligned_contig_collection:
        contig_name = aligned_contig_collection[aligned_contig].contig.id
        aligned_contig.set_bam_file_lists(arguments["g"])
        for sample in aligned_contig_collection[aligned_contig].comparable_samples:
            sys.stderr.write("Is there is a SNP in %s ? Length is %s\n" % (contig_name,
            len(aligned_contig_collection[aligned_contig].comparable_samples[sample])))
            for snp in aligned_contig_collection[aligned_contig].comparable_samples[sample]:
                snp.set_environment_contig(contig_name)


    #2. get env from the bam file... Can we optimize this?

    for aligned_contig in aligned_contig_collection:
        sample_number = 0
        for sample in aligned_contig_collection[aligned_contig].comparable_samples:
            #print("Sample name: %s" % sample)
            #bam_files_list = get_bam_file_list(arguments["g"], sample_number) #That is not optimal...
            #sys.stderr.write("The bam file list: %s\n" % bam_files_list)
            sample_number += 1
            #for position in aligned_contig_collection[aligned_contig].comparable_samples[sample]:
              #  for bam_file in bam_files_list:
              #      if not ".bam" in bam_file:
              ###          sys.stderr.write("The file: %s was not considered because it does not seems to be a bam file\n"
                #                         % bam_file)
                ##        continue
                 #   else:
                  #      print("Populating with file: %s" % bam_file)
                  #      position.environment.populate(bam_file, 30) #30 is min align. score
    return aligned_contig_collection
    #done?
    """


def bam_problem_queue_producer(aligned_contig, problem_queue, num_of_proc):
    num_of_prob =  0
    for contig in aligned_contig:
        for sample in aligned_contig[contig].comparable_samples:
            for position in aligned_contig[contig].comparable_samples[sample]:
                new_problem = [contig, sample, position]
                problem_queue.put(new_problem)
                num_of_prob += 1
                if num_of_prob % 100 == 0:
                    sys.stderr.write("Have created more than %s problems\n" % num_of_prob)
    for proc in range(num_of_proc):
        poison = ["kill"]
        problem_queue.put(poison)
    sys.stderr.write("Bam problem producer (pid: %s) done!\n" % os.getpid())


def bam_problem_queue_solver(problem_queue, result_queue, bam_files_list):
    num_of_prob_solved=0
    while True:
        problem = problem_queue.get()
        if problem[0] == "kill":
            sys.stderr.write("Bam problem solver (pid %s) dying after %s problems\n" % (os.getpid(), num_of_prob_solved))
            break
        aligned_contig = problem[0]
        new_genomic_env = Genomic_environment()
        sys.stderr.write("BAM SOLVER %s: %s for %s at pos %s\n" % (os.getpid(), aligned_contig, problem[1],
                                                                   problem[2].position))
        new_genomic_env.contig = aligned_contig
        new_genomic_env.position = problem[2].position
        sample = problem[1]
        for bam_file in bam_files_list[sample]:
            new_genomic_env.populate(bam_file, 30) #TODO: 30 is min align score. Should be remplace with an argument
        result = [new_genomic_env, sample]
        result_queue.put(result)
        num_of_prob_solved += 1
    result_queue.put(["kill"])



def bam_pipe_cleaner(aligned_contig, result_queue, final_dict, num_of_proc):
    num_of_poison = 0
    while True:
        result = result_queue.get()
        if result[0] == "kill":
            num_of_poison +=1
            sys.stderr.write("Bam pipe cleaner (pid: %s). This is my %s dose of poison\n" % (os.getpid(), num_of_poison))
            if num_of_poison == num_of_proc:
                sys.stderr.write("Bam pipe cleaner (pid: %s) dying\n" % os.getpid())
                break
        continue
        sample = result[1]
        genomic_env = result[0]
        contig_id = aligned_contig[genomic_env.contig].contig.id
        if contig_id in final_dict:
            for position in final_dict[contig_id].comparable_samples[sample]:
                if position.position == genomic_env.position:
                    position.environment = genomic_env
        else:
            #We must create a copy in the final dict
            #unaligned_contig = copy.copy(aligned_contig[contig_id].contig)
            new_aligned_contig = copy.deepcopy(aligned_contig[contig_id])
            final_dict[contig_id] = new_aligned_contig
            for position in final_dict[contig_id].comparable_samples[sample]:
                if position.position == genomic_env.position:
                    position.environment = genomic_env
        #that should do it?


class bam_file_handler:
    def __init__(self, num_of_proc):
        self.NUM_OF_PROCESS =  num_of_proc
        self.problem_queue = Queue(2000)
        self.result_queue = Queue(2000)
        self.manag = Manager()
        self.final_dict = self.manag.dict()

    def start(self, aligned_contig, bam_files_list):
        self.problem_creator = Process(target= bam_problem_queue_producer, args=(aligned_contig, self.problem_queue,
                                                                                 self.NUM_OF_PROCESS))
        self.problem_creator.start()
        sys.stderr.write("Problem creator started")
        self.problem_solvers = [Process(target=bam_problem_queue_solver, args= (self.problem_queue, self.result_queue,
                                                                                bam_files_list))
                                for i in range(self.NUM_OF_PROCESS)]
        for solver in self.problem_solvers:
            solver.start()
        self.pipe_cleaner = Process(target=bam_pipe_cleaner, args=(aligned_contig, self.result_queue, self.final_dict, self.NUM_OF_PROCESS))
        self.pipe_cleaner.start()


    def check_status(self):
        pass

    def join(self):
        self.problem_creator.join()
        for solver in self.problem_solvers:
            solver.join()
        self.pipe_cleaner.join()


def get_bam_file_list(file_name):
    bam_files_list = {}
    for lines in open(file_name, 'rU'):
        sample_name = lines.split()[0]
        bam_files_list[sample_name] = []
        bam_dir = lines.split()[2].rstrip('\n')
        bam_files = glob.glob(bam_dir+"/*.bam")
        if len(bam_files) == 0:
            raise IOError("No bam file in directory %s\n" % bam_dir)
        bam_files_list[sample_name] = bam_files
    return bam_files_list

def parse_gatk_like_file(contig_collection, gff_file):
    
    if not isValid(gff_file):
        raise IOError("File '" + gff_file + "' could not be opened.")
    
    out_dict = {}
    
    for line in open(gff_file, 'rU'):
        line_split = line.split()

        info = line_split[0]
        gene_start = int(line_split[1])
        gene_end = int(line_split[2])
        sens = line_split[3]
        
        if sens == '1':
            sens = '+'
        elif sens == '-1':
            sens = '-'
        else:
            raise ValueError("Strand identifier is not recognisable.")

        contig_id_begin_index = info.index('contig')

        gene_id = info[contig_id_begin_index:]
        contigId = gene_id.split('_')[0]
        
        contig_header = '>' + contigId
        gene_id_fixed = '>' + gene_id
        
        filtered_genes = list()
        
        
        for gene in contig_collection[contig_header].genes:
            if gene.name == gene_id_fixed:
                if gene.start != gene_start or gene.end != gene_end or gene.strand != sens:
                    raise ValueError("Gene '"+ gene_id +"' in gff file does not match in contig collection!")
                else:
                    filtered_genes.append(gene)
        
        out_dict[contig_header] = copy.copy(contig_collection[contig_header])
        out_dict[contig_header].genes = filtered_genes
        
    return out_dict


def run_first_pipeline(arguments):

    #1. load Contigs.fasta
    sys.stderr.write("Loading contigs.fasta\n")
    contig_collection = load_contigs_dot_fasta(arguments["c"])
    sys.stderr.write("Loading done\n")
    #2. Load the genes from prodigal (load genes)
    sys.stderr.write("Loading genes from prodigal file\n")
    gene_collection = load_genes(arguments["p"])
    sys.stderr.write("Loading genes done\n")
    #3 Associate gene to contigs
    sys.stderr.write("Associating genes to contigs\n")
    contig_collection = associate_genes_to_contigs(contig_collection, gene_collection)
    sys.stderr.write("Associating genes to contigs done\n")
    
    contig_collection = parse_gatk_like_file(contig_collection, arguments["gff"])
    
    #4 Create an "Aligned_contig" for each contig since we want to compare alignements
    #4 Load gatk files using a cutoff...

    m = GATK_handler(arguments["n"])
    #Test: if i launch the bam_handler here... will there be a memory problem with the child process.?
    b = bam_file_handler(arguments["n"])
    #####EO test
    for sample in open(arguments["g"], 'rU'):

        sample_id = sample.split()[0]
        print("Sample in treatment: %s\n" % sample_id)
        sys.stderr.write("Loading GATK file: %s\n" % sample.split('\t')[1])
        file_name = sample.split('\t')[1]
        #m = Managerss(2)
        m.start(file_name, sample_id, arguments, contig_collection)
        #m.check_status()   #For debogue
        m.join()
        #We want this stuff in a dict.. not in a list!

    aligned_contig = m.aligned_contig_dict
    del(m) #needed?
    
    sys.stderr.write("Number of grouped contigs %s\n" % len(aligned_contig))
    sys.stderr.write("Loading GATK files done\n")
    sys.stderr.write("Computing context\n")
    bam_files_list = get_bam_file_list(arguments["g"])
    b.start(aligned_contig, bam_files_list)
    b.join()
    aligned_contig = get_genomic_environment_from_bam(arguments, aligned_contig)
    #TODO: add warning if depth between GATK and BAM file is too different
    sys.stderr.write("Computing context done\n")



if __name__ == "__main__":

    parser = Option_parser(sys.argv[1:])
    arguments = parser.get_arguments()

    """
    Fist pipeline. Need: -c , -p , -g
    1. load contigs.fasta
    2. Load genes from Prodigal outfile
    3. Load GATK file using cutoff
    4. Get context from bam files
    """
    run_first_pipeline(arguments)


