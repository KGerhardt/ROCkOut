import sys
import os
import multiprocessing

from rocker_project_manager import project_manager

class aligner:
	def __init__(self, basename = None, tagged_fasta = None, target_database = None, tool = 'BLAST'):
		self.base = basename
		self.input = tagged_fasta
		self.target = target_database
		self.tool = tool
		
	#Module for aligning reads with BLAST. One of these modules per tool; database has to match tool fmt.
	def blast_align(self):
		pass
	
	def mmseqs_align(self):
		pass

#Database building class		
class db_builder:
	def __init__(self, input_files = None, tool = 'mmseqs', output_file = None):
		self.inputs = input_files
		self.tool = tool
		self.output = output_file
		
		self.merged_input = None
		
	def merge_fastas(self):
		for f in self.inputs:
			fh = open(f)
			
			fh.close()
		
	#Module for building a BLAST database from 
	def build_blast_db(self):
		pass
		
	def build_mmseqs_db(self):
		pass
		command = ["mmseqs", "createdb", self.merged_input]

		
def align_reads(project_directory, threads, tool = "mmseqs"):
	rocker = project_manager(directory = project_directory, threads = threads)
	#This gets 
	rocker.parse_project_directory()
	rocker.parse_genomes()
	rocker.parse_tagged_reads()
	
	#Set up database location for the projet
	shared_loc = os.path.normpath(project_directory + "/shared_files") + "/"
	shared_db = os.path.normpath(project_directory + "/shared_files/alignment_databases") + "/"
	shared_genomes = os.path.normpath(project_directory + "/shared_files/combined_genomes") + "/"
	
	if not os.path.exists(shared_loc):
		os.mkdir(shared_loc)
	if not os.path.exists(shared_db):
		os.mkdir(shared_db)
	if not os.path.exists(shared_genomes):
		os.mkdir(shared_genomes)
		
	target_db_name = shared_db + tool + "_combined_database.db"
	
	fasta_paths = []
	to_do = []
	for base in rocker.positive:
		for fasta in rocker.alnomes_pos[base]:
			aln = read_aligner(base_name = rocker.positive[base], input_fasta = fasta)
			aln.prepare_outputs()
			to_do.append(aln)
			count += 1
	for base in rocker.negative:
		for fasta in rocker.alnomes_neg[base]:
			aln = read_aligner(base_name = rocker.negative[base], input_fasta = fasta)
			aln.prepare_outputs()
			to_do.append(aln)
			count += 1
	
	#build DB
	
	#align reads
	pool = multiprocessing.Pool(threads)
	pool.map(run_tagging, to_do)
	pool.close()
	pool.join()
	
align_reads(project_directory = project_directory, threads = threads)