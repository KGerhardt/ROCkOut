import sys
import os
import subprocess
import re
import numpy as np
import multiprocessing
import shutil

import argparse

from .rocker_project_manager import project_manager
from .rocker_progress_tracker import progress_tracker

#Class for simulating and tagging reads in ROCker friendly format.
class read_generator:
	#def __init__(self, input_fasta = None, input_gffs = None, positive = True, simulator = "BBMap", coordinates = None, is_pos = False):
	def __init__(self, base_name = None, input_fasta = None, simulator = "BBMap", 
	build_num = 1, coverage_depth = 20.0, snp_rate_per_base = 0.01, insertion_rate = 0.01/19.0, 
	deletion_rate = 0.01/19.0, read_length_min = 80, read_length_mid = 100, read_length_max = 120):
		
		self.base = base_name
		
		self.input_base = os.path.basename(input_fasta)[:-10]
		
		self.sim = simulator
		#File name after all directories and before the final .[extension]
		
		self.build_num = build_num
		
		self.input = input_fasta
		
		#FastQ may not always be generated, but is with the following tools: BBMap
		self.fastq_out = None
		#These two will always be made
		self.fasta_out = None
		
		self.log_file = os.path.normpath(self.base + "/bbmap_log/" + self.input_base + "_log.txt")
		
		self.generation_log = ""
		
		self.cov_depth = coverage_depth
		self.snp_rate = snp_rate_per_base
		#Default behavior for these items is to be 5 indels : 95 SNPs
		self.insrate = insertion_rate
		self.delrate = deletion_rate
			
		self.minlen = read_length_min
		self.mulen =  read_length_mid
		self.maxlen = read_length_max
		
		
	def rename_outputs(self):
		#Names are generated with read length info for model building.
		read_len_label = str(int(round((self.minlen + self.maxlen)/2, 0)))
		#FastQ may not always be generated, but is with the following tools: BBMap
		self.fastq_out = os.path.normpath(self.base + "/raw_reads/" + self.input_base + "_read_length_" + read_len_label + "_raw.fastq")
		#These two will always be made
		self.fasta_out = os.path.normpath(self.base + "/raw_reads/" + self.input_base + "_read_length_" + read_len_label + "_raw.fasta")
		
	def prepare_outputs(self):
		fastas = os.path.normpath(self.base + "/raw_reads") + "/"
		logs = os.path.normpath(self.base + "/bbmap_log") + "/"
		
		if not os.path.exists(fastas):
			os.mkdir(fastas)
		if not os.path.exists(logs):
			os.mkdir(logs)
	
	#Simulate reads using BBmap
	#Args:
	'''
		build_num - needed for threading; bbamp doesn't automatically handle parallelization well. Each parallel process from the same dir needs its own number.
		snp_rate
		insrate
	'''	
	def simulate_reads_bbmap(self):
		build_num = "build="+str(self.build_num) 
		cov_depth = "coverage="+str(self.cov_depth)
		snp_rate = "snprate="+str(self.snp_rate)
		insrate = "insrate="+str(self.insrate)
		delrate = "delrate="+str(self.delrate) 
		minlen = "minlength="+str(self.minlen)
		medlen = "midlength="+str(self.mulen)
		maxlen = "maxlength="+str(self.maxlen)
		
		
		has_midlen = False
		sanity_checker = subprocess.Popen(["randomreads.sh", "--help"], stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
		for line in sanity_checker.stdout:
			line = line.decode()
			if "midlen" in line:
				has_midlen = True
		
		if has_midlen:
			#Set randomreads script to run with the given args.
			command =  ["randomreads.sh", "simplenames=t", "gaussianlength=t", build_num, "ref="+self.input, "out="+self.fastq_out, snp_rate, cov_depth, insrate, delrate, minlen, medlen, maxlen]
		else:
			command =  ["randomreads.sh", "simplenames=t", "gaussianlength=t", build_num, "ref="+self.input, "out="+self.fastq_out, snp_rate, cov_depth, insrate, delrate, minlen, maxlen]

		
		#Run the command and capture the output for logging purposes
		proc = subprocess.Popen(command, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
		for line in proc.stdout:
			#Lines are imported in binary, this just converts to plaintext
			self.generation_log += line.decode()
		
		fh = open(self.log_file, "w")
		fh.write(self.generation_log)
		fh.close()
		
		try:
			proc.stdout.close()
		except:
			pass
			
		try:
			proc.stderr.close()
		except:
			pass
		
	
	#Convert fastq to fasta. Assumes that fastq format is strictly followed, 4 lines exactly.
	def fastq_to_fasta(self):
		if os.path.exists(os.path.normpath(self.fastq_out)):
			line_counter = 1
			fq = open(self.fastq_out)
			fa = open(self.fasta_out, "w")
			
			#We read the FASTA and parse the reads to tag them along the way.
			for line in fq:
				#Defline
				if line_counter % 4 == 1:
					fa.write(">"+line[1:])
				#Seq
				if line_counter % 4 == 2:
					fa.write(line)
				#Other lines are not data to keep.
					
				line_counter += 1
			
			fa.close()
			fq.close()
			#Clean up.
			try:
				os.remove(self.fastq_out)
			except:
				print("Couldn't remove", self.fastq_out)
		else:
			print(self.fastq_out, "could not be generated. ROCkOut cannot continue without this file.")
		
		
	def clean_up(self):
		try:
			shutil.rmtree("ref/genome/"+str(self.build_num))
		except:
			pass
		try:
			shutil.rmtree("ref/index/"+str(self.build_num))
		except:
			print("Couldn't remove ref/index/"+str(self.build_num))

		
def run_generation(read_gen):
	read_gen.rename_outputs()
	read_gen.simulate_reads_bbmap()
	read_gen.fastq_to_fasta()
	read_gen.clean_up()
	return read_gen.base, read_gen.mulen

def parse_read_triplet(rt):
	rt = rt.split(",")
	rt = [int(r) for r in rt]
	rt = tuple(rt)
	return rt
	
def generate_reads(project_directory, threads, short = "90,100,110", 
	standard = "180,200,220", long = "270,300,330", extra_long = "360,400,440",
	cov = 20.0, snp = 0.01, insrate = 0.01/19, delrate = 0.01/19):
	
	rocker = project_manager(directory = project_directory, threads = threads)
	rocker.parse_project_directory()
	rocker.parse_coords()
	rocker.parse_genomes()
	rocker.check_genome_lengths()
	
	to_do = []
	count = 1
	
	try:
		short = parse_read_triplet(short)
	except:
		print("Couldn't interpret short read lengths", short, "setting to default 90, 100, 110")
		short = (90, 100, 110,)
		
	try:
		standard = parse_read_triplet(standard)
	except:
		print("Couldn't interpret med read lengths", standard, "setting to default 180, 200, 220")
		standard = (180, 200, 220,)
		
	try:
		long = parse_read_triplet(long)
	except:
		print("Couldn't interpret long read lengths", standard, "setting to default 270, 300, 330")
		long = (270, 300, 330,)
			
	try:
		extra_long = parse_read_triplet(extra_long)
	except:
		print("Couldn't interpret extra long read lengths", standard, "setting to default 360, 400, 440")
		extra_long = (360, 400, 440,)
	
	sim_lens = [short, standard, long, extra_long]
	print("")
	print("Simulating reads from gaussian distribution with min, mean, max lengths:")
	for rl in sim_lens:
		print("\t",rl)
		
	print("Using parameters:")
	print("Coverage:", cov)
	print("SNP rate:", snp)
	print("Insertion rate:", insrate)
	print("Deletion rate:", delrate)
	print("")
	
	for base in rocker.positive:
		for fasta, length in zip(rocker.genomes_pos[base], rocker.gen_lengths_pos[base]):
			hits = 0
			for length_set in sim_lens:
				#Do not add a too-short genome - the genome must be 1.5x larger than the max read length of the triplet.
				if length > 1.5 * length_set[2]:
					gen = read_generator(base_name = rocker.positive[base], input_fasta = fasta, 
					build_num = count, coverage_depth = cov, snp_rate_per_base = snp, insertion_rate = insrate,  deletion_rate = delrate,
					read_length_min = length_set[0], read_length_mid = length_set[1], read_length_max = length_set[2])
					gen.prepare_outputs()
					to_do.append(gen)
					count += 1
					hits += 1
			if hits == 0:
				print("")
				print("Protein:", base, "genome:", fasta, "was unable to have reads simulated.")
	for base in rocker.negative:
		for fasta, length in zip(rocker.genomes_neg[base], rocker.gen_lengths_neg[base]):
			hits = 0
			for length_set in sim_lens:
				if length > 1.5 * length_set[2]:
					gen = read_generator(base_name = rocker.negative[base], input_fasta = fasta, 
					build_num = count, coverage_depth = cov, snp_rate_per_base = snp, insertion_rate = insrate,  deletion_rate = delrate,
					read_length_min = length_set[0], read_length_mid = length_set[1], read_length_max = length_set[2])
					gen.prepare_outputs()
					to_do.append(gen)
					count += 1
					hits += 1
			if hits == 0:
				print("")
				print("Protein:", base, "genome:", fasta, "was unable to have reads simulated.")
			
	prog_bar = progress_tracker(total = len(to_do), message = "Simulating reads from your proteins...")
	pool = multiprocessing.Pool(threads)
	for result in pool.imap_unordered(run_generation, to_do):
		prog_bar.update()
		#print(result[0], "read length", result[1], "done!")
	pool.close()
	pool.join()
	
	print("Cleaning up.")
	
	try:
		shutil.rmtree("ref")
	except:
		print("Trouble removing simulated refs")
		
	print("Reads simulated!")
	
def options():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
			description='''	''')
	
	parser.add_argument('-t', '--threads',  dest = 'threads', default = 1, help =  '')
	parser.add_argument('-o', '--output',  dest = 'out', default = None, help =  'A rocker project directory that you have downloaded genomes for.')

	args, unknown = parser.parse_known_args()
	
	return parser, args
	
	
def generate(parser, opts):
	project_directory = opts.dir
	if project_directory is None:
		print("You need to specify a ROCkOut project directory. Quitting.")
		print(parser.print_help())
		quit()
		
	try:
		threads = int(opts.threads)
	except:
		print("Threads must be an integer. Defaulting to 1 thread.")
		threads = 1
		
	s, m, l, xl = opts.short, opts.med, opts.long, opts.xl
	
	#parser.add_argument('--coverage',  dest = 'cov', default = "20", help = "Read coverage depth for simulated reads. Default 20") 
	#parser.add_argument('--snprate',  dest = 'snps', default = "0.01", help = "Per base substitution likelihood. Default 0.01") 
	#parser.add_argument('--insertrate',  dest = 'insrate', default = None, help = "Insertion rate. Default 1/19th of snprate.") 
	#parser.add_argument('--delrate',  dest = 'delrate', default = None, help = "Deletion rate. Default 1/19th of snprate.") 
	
	coverage = float(opts.cov)
	snprate = float(opts.snps)
	insrate = opts.insrate
	delrate = opts.delrate
	
	try:
		if insrate is None:
			insrate = snprate/19
		else:
			insrate = float(insrate)
	except:
		print("Insertion rate couldn't be determined from your parameter", insrate)
		print("Defaulting to 1/19th default snp rate (0.01/19).")
		insrate = 0.01/19
		
	try:
		if delrate is None:
			delrate = snprate/19
		else:
			delrate = float(delrate)
	except:
		print("Deletion rate couldn't be determined from your parameter", insrate)
		print("Defaulting to 1/19th snp rate (0.01/19).")	
		delrate = 0.01/19
		
	multiprocessing.freeze_support()
	
	generate_reads(project_directory = project_directory, threads = threads, 
	short = s, standard = m, long = l, extra_long = xl)

'''	
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
		
	s, m, l, xl = opts.short, opts.med, opts.long, opts.xl
	print(s, m, l, xl)
		
	multiprocessing.freeze_support()
	
	generate_reads(project_directory = project_directory, threads = threads, 
	short = s, standard = m, long = l, extra_long = xl)
'''	
	

