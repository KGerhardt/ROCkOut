import sys
import os
import multiprocessing
import subprocess

import argparse

from rocker_project_manager import project_manager
from rocker_progress_tracker import progress_tracker

class aligner:
	def __init__(self, base_name = None, base_path = None, original_reads = None, tagged_fasta = None, target_database = None, tool = 'diamond'):
		self.base = base_name
		self.path = base_path
		self.og = original_reads
		self.input = tagged_fasta
		self.target = target_database
		self.tool = tool
		
		self.final_read_count = 0
		
		self.read_file = os.path.normpath(self.path + "aligned_reads/" + self.base + "_aligned_reads.txt")
		
	
	def prepare_outputs(self):
		if not os.path.exists(self.path + "aligned_reads/"):
			os.mkdir(self.path + "aligned_reads/")
		#outputs go here.
		if not os.path.exists(self.path + "alignment_log/"):
			os.mkdir(self.path + "alignment_log/")
	
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
			diamond_err = open(self.path + "alignment_log/"+self.base + "_DIAMOND_log.txt", "w")
			#diamond_stdout = open(self.path + "alignment_log/"+self.base + "_DIAMOND_stdout.txt", "w")
			
			#Original ruby command structure
			# diamond: '%1$sdiamond %2$s -q "%3$s" -d "%4$s" -a "%5$s.daa" -p %6$d' +
			#' -k 1 --min-score 20 --sensitive && %1$sdiamond view -a "%5$s"' +
			#' -o "%5$s"'},
			
			#-k = --max-target-seqs
			#-p = --threads
			#-a = --daa DIAMOND alignment archive (DAA) file
			
			align_command = ["diamond", "blastx", "--sensitive", "--max-target-seqs", "1",  "--unal", "0", "--outfmt", "6", "--db", self.target, "--query", self.input, "--out", self.read_file]
			
			#print(' '.join(align_command))
			
			subprocess.call(align_command, stderr = diamond_err)
			
			diamond_err.close()
			
			fh = open(self.read_file)
			out = open(self.read_file + "_filtered.txt", "w")
			for line in fh:
				if line.startswith(";"):
					print("Malformed read:", line.strip())
					continue
				segs = line.strip().split("\t")
				if len(segs) < 12:
					print("Malformed read:", line.strip())
					continue
					
				#Only occurs when the two continues do not.
				out.write(line)
			
			out.close()
			fh.close()
			
			os.remove(self.read_file)
			os.rename(self.read_file + "_filtered.txt", self.read_file)
			#diamond_stdout.close()
		
	
	def clean_raw_file(self, file_name, hits):
		fh = open(file_name)
		out = open(file_name + "_filtered.txt", "w")
		defline = fh.readline()
		while defline:
			#As long as there is a defline, seq follows.
			seq = fh.readline()
			read_id = defline.split("_")[1]
			if read_id in hits:
				out.write(defline)
				out.write(seq)
			#Read next line
			defline = fh.readline()
			
		out.close()
		fh.close()
		os.remove(file_name)
		os.rename(file_name + "_filtered.txt", file_name)
		
	def clean_tagged_file(self, file_name, hits):
		fh = open(file_name)
		out = open(file_name + "_filtered.txt", "w")
		defline = fh.readline()
		while defline:
			#As long as there is a defline, seq follows.
			seq = fh.readline()
			#The read ID is pre-tagged here
			read_id = defline.split(";")[0][1:]
			if read_id in hits:
				out.write(defline)
				out.write(seq)
			#Read next line
			defline = fh.readline()
			
		out.close()
		fh.close()
		os.remove(file_name)
		os.rename(file_name + "_filtered.txt", file_name)
	
	def clean_up_your_mess(self):
		hits = {}
		fh = open(self.read_file)
		for line in fh:
			read_id = line.split("\t")[0].split(";")[0]
			#The data in the dict doesn't matter.
			hits[read_id] = 0
		fh.close()
		
		#How many reads remain.
		self.final_read_count = len(hits)
		
		self.clean_raw_file(self.og, hits)
		self.clean_tagged_file(self.input, hits)
		
		
		
		
			
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
			print("Logs recording alignment information can be found in the pos_or_neg/protein/alignment_log/ subdirectories of your project.")
			print("")
			
		if self.tool == "mmseqs":
			pass
		

def align_reads(project_directory, threads, tool = "diamond"):
	rocker = project_manager(directory = project_directory, threads = threads)
	#This gets 
	rocker.parse_project_directory()
	#For cleaning, collect untagged reads
	rocker.parse_raw_reads()
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
	
	#The base name is not a unique ID. Whoops.
	
	for base in rocker.positive:
		for raws, fasta in zip(rocker.read_fastas_pos[base], rocker.tagged_reads_pos[base]):
			unique_id = os.path.basename(fasta)[:-17]
			aln = aligner(base_name = unique_id, base_path = os.path.normpath(rocker.positive[base])+"/", original_reads = raws, tagged_fasta = fasta, target_database = target_db_name, tool = tool)
			aln.prepare_outputs()
			to_do.append(aln)
			count += 1	
	for base in rocker.negative:
		for raws, fasta in zip(rocker.read_fastas_neg[base], rocker.tagged_reads_neg[base]):
			unique_id = os.path.basename(fasta)[:-17]
			aln = aligner(base_name = unique_id, base_path = os.path.normpath(rocker.negative[base])+"/", original_reads = raws, tagged_fasta = fasta, target_database = target_db_name, tool = tool)
			aln.prepare_outputs()
			to_do.append(aln)
			count += 1
	

	#align reads
	prog_bar = progress_tracker(total = len(to_do), message = "Aligning reads. This will take some time.")
	run_data = []
	pool = multiprocessing.Pool(threads)
	for result in pool.imap_unordered(run_align, to_do):
		run_data.append(result)
		prog_bar.update()
	pool.close()
	pool.join()
	
	#Create project index:
	index = open(os.path.normpath(project_directory+ "/ROCkOUT_index.txt"), "w")
	print("protein_name", "reads_included", "group", "in_final_model", sep = "\t", file = index)
	for result in run_data:
		print(*result, "T", sep = "\t", file = index)
		
	index.close()
	
	print("Read alignment done!")

def run_align(aligner_object):
	aligner_object.align()
	#Simulated reads are usually almost entirely off-target, leaving large files that are almost entirely useless.
	#Goes back and removes any reads that aren't in the aligned set.
	aligner_object.clean_up_your_mess()
	if "positive" in aligner_object.path:
		group = "positive"
	else:
		group = "negative"
	
	return (aligner_object.base, aligner_object.final_read_count, group,)
	
	
	
def options():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
			description='''	''')
	
	parser.add_argument('-t', '--threads',  dest = 'threads', default = 1, help =  '')
	parser.add_argument('-o', '--output',  dest = 'out', default = None, help =  'A rocker project directory that you have downloaded genomes for.')

	args, unknown = parser.parse_known_args()
	
	return parser, args
	
def main():
	parser, opts = options()
	
	project_directory = opts.out
	if project_directory is None:
		print("You need to specify a ROCkOut project directory. Quitting.")
		print(parser.print_help())
		quit()
		
	try:
		threads = int(opts.threads)
	except:
		print("Threads must be an integer. Defaulting to 1 thread.")
		threads = 1
		
	multiprocessing.freeze_support()
	
	align_reads(project_directory = project_directory, threads = threads)
	
	
if __name__ == "__main__":
	main()
