import sys
import os
import multiprocessing
import subprocess

from rocker_project_manager import project_manager

class aligner:
	def __init__(self, base_name = None, base_path = None, tagged_fasta = None, target_database = None, tool = 'diamond'):
		self.base = base_name
		self.path = base_path
		self.input = tagged_fasta
		self.target = target_database
		self.tool = tool
		
		self.read_file = os.path.normpath(self.path + "aligned_reads/" + self.base + "_aligned_reads.txt")
		
	
	def prepare_outputs(self):
		if not os.path.exists(self.path + "aligned_reads/"):
			os.mkdir(self.path + "aligned_reads/")
	
	def align(self):
		self.prepare_outputs()
		
		if self.tool == "diamond":
			'''
			Other opts to consider?
			--evalue (-e)            maximum e-value to report alignments (default=0.001)
			--min-score              minimum bit score to report alignments (overrides e-value setting)
			--id                     minimum identity% to report an alignment
			--query-cover            minimum query cover% to report an alignment
			--subject-cover          minimum subject cover% to report an alignment
			--threads
			'''
			#diamond blastx --db rdl/shared_files/alignment_database/diamond_combined_database.db.dmnd --query XXX
			align_command = ["diamond", "blastx", "--unal", "0", "--outfmt", "6", "--db", self.target, "--query", self.input, "--out", self.read_file]
			subprocess.call(align_command)
			

#Database building class		
class db_builder:
	def __init__(self, tool = 'diamond', input_file = None, output_file = None):
		self.input = input_file
		self.tool = tool
		self.output = output_file
		
		self.merged_input = input_file
		
	def build(self):
		if self.tool == "diamond":
			#This is building a blast protein DB so that diamond can do the aligning.
			makedb = ["diamond", "makedb", "--db", self.output,  "--in", self.input]
			print("Building Diamond database. Log information will follow.")
			subprocess.call(makedb)
			print("")
			print("Done! Alignment will now begin.")
			print("")
			
		if self.tool == "mmseqs":
			pass
		

		
def align_reads(project_directory, threads, tool = "diamond"):
	rocker = project_manager(directory = project_directory, threads = threads)
	#This gets 
	rocker.parse_project_directory()
	#get tagged reads, ready for alignment as QUERIES
	rocker.parse_tagged_reads()
	#Get reference proteins for building the DB and using as alignment TARGETS.
	rocker.parse_targets()
	
	#Set up database location for the projet
	shared_loc = os.path.normpath(project_directory + "/shared_files") + "/"
	shared_db = os.path.normpath(project_directory + "/shared_files/alignment_database") + "/"
	shared_genomes = os.path.normpath(project_directory + "/shared_files/combined_genomes") + "/"
	
	if not os.path.exists(shared_loc):
		os.mkdir(shared_loc)
	if not os.path.exists(shared_db):
		os.mkdir(shared_db)
	if not os.path.exists(shared_genomes):
		os.mkdir(shared_genomes)
		
	target_db_name = shared_db + tool + "_combined_database.db"
	target_fasta_set = shared_genomes + "combined_genomes.fasta.txt"

	#Join the positive files
	combiner = open(target_fasta_set, "w")
	#We want all of them. Align now, filter later.
	for base in rocker.targets:
		for genome in rocker.targets[base]:
			in_fasta = open(genome)
			for line in in_fasta:
				combiner.write(line)
			in_fasta.close()
	
	combiner.close()
	
	#Create the database against which alignment will be done
	db = db_builder(tool, target_fasta_set, target_db_name)
	#Run build appropriate to the tool
	db.build()
		
	fasta_paths = []
	to_do = []
	
	#build out a list of files to align to pass to a parallel function later.
	
	#args...
	#def __init__(self, basename = None, tagged_fasta = None, target_database = None, tool = 'diamond')
	count = 0
	for base in rocker.positive:
		for fasta in rocker.tagged_reads_pos[base]:
			aln = aligner(base_name = base, base_path = os.path.normpath(rocker.positive[base])+"/", tagged_fasta = fasta, target_database = target_db_name, tool = tool)
			aln.prepare_outputs()
			to_do.append(aln)
			count += 1
	for base in rocker.negative:
		for fasta in rocker.tagged_reads_neg[base]:
			aln = aligner(base_name = base, base_path = os.path.normpath(rocker.negative[base])+"/", tagged_fasta = fasta, target_database = target_db_name, tool = tool)
			aln.prepare_outputs()
			to_do.append(aln)
			count += 1
	
	#align reads
	pool = multiprocessing.Pool(threads)
	pool.map(run_align, to_do)
	pool.close()
	pool.join()

def run_align(aligner_object):
	aligner_object.align()
	
	
project_directory = sys.argv[1]
threads = int(sys.argv[2])
	
align_reads(project_directory = project_directory, threads = threads)