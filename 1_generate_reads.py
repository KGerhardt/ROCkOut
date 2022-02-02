import sys
import os
import subprocess
import re
import numpy as np
import multiprocessing

from rocker_project_manager import project_manager

#Class for simulating and tagging reads in ROCker friendly format.
class read_generator:
	#def __init__(self, input_fasta = None, input_gffs = None, positive = True, simulator = "BBMap", coordinates = None, is_pos = False):
	def __init__(self, base_name = None, input_fasta = None, simulator = "BBMap", build_num = 1):
		self.base = base_name
		
		self.input_base = os.path.basename(input_fasta)[:-10]
		
		self.sim = simulator
		#File name after all directories and before the final .[extension]
		
		self.build_num = build_num
		
		self.input = input_fasta
		
		#FastQ may not always be generated, but is with the following tools: BBMap
		self.fastq_out = os.path.normpath(self.base + "/raw_reads/" + self.input_base + "_fastq.txt")
		#These two will always be made
		self.fasta_out = os.path.normpath(self.base + "/raw_reads/" + self.input_base + "_fasta.txt")
		
		self.log_file = os.path.normpath(self.base + "/bbmap_log/" + self.input_base + "_log.txt")
		
		self.generation_log = ""
		
		
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
	def simulate_reads_bbmap(self, cov_depth = 20.0, snp_rate = 0.1, insrate = None, delrate = None, minlen = 95, maxlen = 105):
		#Default behavior for these items is to be 5 indels : 95 SNPs
		if insrate is None:
			insrate = snp_rate/19
		if delrate is None:
			delrate = snp_rate/19
		
		build_num = "build="+str(self.build_num) 
		cov_depth = "coverage="+str(cov_depth)
		snp_rate = "snprate="+str(snp_rate)
		insrate = "insrate="+str(insrate)
		delrate = "delrate="+str(delrate) 
		minlen = "minlength="+str(minlen)
		maxlen = "maxlength="+str(maxlen) 
		
		#Set randomreads script to run with the given args.
		command =  ["randomreads.sh", "simplenames=t", build_num, "ref="+self.input, "out="+self.fastq_out, snp_rate, cov_depth, insrate, delrate, minlen, maxlen]
		
		#Run the command and capture the output for logging purposes
		proc = subprocess.Popen(command, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
		for line in proc.stdout:
			#Lines are imported in binary, this just converts to plaintext
			self.generation_log += line.decode()
		
		fh = open(self.log_file, "w")
		fh.write(self.generation_log)
		fh.close()
		
	
	#Convert fastq to fasta. Assumes that fastq format is strictly followed, 4 lines exactly.
	def fastq_to_fasta(self):
		line_counter = 1
		fq = open(self.fastq_out)
		fa = open(self.fasta_out, "w")
		
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
		os.remove(self.fastq_out)
		
		
project_directory = sys.argv[1]
try:
	threads = int(sys.argv[2])
except:
	print("Couldn't recognize threads as a number. Setting threads to 1.")
	threads = 1

def run_generation(read_gen):
	#read_gen.prepare_outputs()
	read_gen.simulate_reads_bbmap()
	read_gen.fastq_to_fasta()

def generate_reads(project_directory, threads):
	rocker = project_manager(directory = project_directory, threads = threads)
	rocker.parse_project_directory()
	rocker.parse_coords()
	rocker.parse_genomes()
	
	to_do = []
	count = 1
	for base in rocker.positive:
		for fasta in rocker.genomes_pos[base]:
			gen = read_generator(base_name = rocker.positive[base], input_fasta = fasta, build_num = count)
			gen.prepare_outputs()
			to_do.append(gen)
			count += 1
	for base in rocker.negative:
		for fasta in rocker.genomes_neg[base]:
			gen = read_generator(base_name = rocker.negative[base], input_fasta = fasta, build_num = count)
			gen.prepare_outputs()
			to_do.append(gen)
			count += 1
		
	pool = multiprocessing.Pool(threads)
	pool.map(run_generation, to_do)
	pool.close()
	pool.join()
	

generate_reads(project_directory = project_directory, threads = threads)