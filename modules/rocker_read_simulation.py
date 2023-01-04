import sys
import os
import subprocess
import re
import numpy as np
import multiprocessing
import shutil

import argparse

import tarfile

from .rocker_project_manager import project_manager
from .rocker_progress_tracker import progress_tracker

#from rocker_project_manager import project_manager
#from rocker_progress_tracker import progress_tracker

'''
This is the most involved part of rockout

First, generate a diamond database from all proteins for both positive and negative
--max-target-seqs = 0. This produces a report of ALL alignments for each read
We'll sort out the best hit at refine time as the read which has:
	hits a positive target sequence
	has max bit score among those hits

Each genome will:

1 simulate reads covering each genome at [depth]
2 tag reads as originating from a specific genome and as on target or off target
	2b negative is implicit when a genome is classed as negative - this takes care of the target/non-target distinction

3 align reads to the shared database, but keep the aligned results separate on a per-prot basis
4 clean up the original reads and retain only those which aligned to at least one protein
	
'''

def run_sim(read_simulator_item):
	#print("Simulating", read_simulator_item.inf)
	read_simulator_item.simulate_reads()	
	read_simulator_item.fastq_to_fasta()
	read_simulator_item.tag_reads()
	read_simulator_item.align_reads()	
	read_simulator_item.clean_up_your_mess()
	return read_simulator_item.inf
			
class read_manager:
	def __init__(self, 
				threads = 1, 
				dir = None,
				#Read simulator arguments
				coverage_depth = 20.0, 
				snp_rate_per_base = 0.01, 
				insertion_rate = 0.01/19.0, 
				deletion_rate = 0.01/19.0, 
				short = "90,100,110", 
				standard = "180,200,220", 
				long = "270,300,330", 
				extra_long = "360,400,440"):
				
		self.threads = threads
		self.rocker_directory = dir
		self.project = project_manager(directory = self.rocker_directory)
		
		self.project.parse_project_directory()
		self.project.parse_coords()
		self.project.parse_genomes()
		self.project.check_genome_lengths()
		self.project.parse_targets()
		
		#Database stuff
		#raw = self.prep_dir(out_base, "raw_reads/")
		self.shared_loc = self.prep_dir(self.rocker_directory, "shared_files/")
		self.shared_db = self.prep_dir(self.rocker_directory, "shared_files/alignment_database/")
		self.shared_genomes = self.prep_dir(self.rocker_directory, "shared_files/combined_genomes/")
		
		self.multiple_alignment = self.prep_dir(self.rocker_directory, "shared_files/multiple_alignment/")
			
		self.make_db()
		
		#Read sim stuff
		self.short = self.parse_read_triplet(short, "90,100,110")
		self.standard = self.parse_read_triplet(standard, "180,200,220")
		self.long = self.parse_read_triplet(long, "270,300,330")
		self.xl = self.parse_read_triplet(extra_long, "360,400,440")
		
		self.simlens = [self.short, self.standard, self.long, self.xl]
		
		self.sim_arg_template = ' '.join(["randomreads.sh",
								"coverage="+str(coverage_depth),
								"snprate="+str(snp_rate_per_base),
								"insrate="+str(round(insertion_rate, 8)),
								"delrate="+str(round(deletion_rate, 8)),
								"simplenames=t",
								"gaussianlength=t",
								"minlength={min}",
								"midlength={med}",
								"maxlength={max}",
								"build={build_num}",
								"ref={input_genome}",
								"out={output_fastq}"
							])
	
		self.genomes_to_sim = []
		self.genlens = []
		self.coords_files = []
	
		#print(self.project.genomes_pos)
		#print(self.project.gen_lengths_pos)
		#print(self.project.gen_lengths_neg)
		#print(self.project.coords_pos)
	
		for item in self.project.genomes_pos:
			for gen, l in zip(self.project.genomes_pos[item], self.project.gen_lengths_pos[item]):
				self.genomes_to_sim.append(gen)
				self.genlens.append(l)
			for gen in self.project.coords_pos[item]:
				self.coords_files.append(gen)
				
		for item in self.project.genomes_neg:
			for gen, l in zip(self.project.genomes_neg[item], self.project.gen_lengths_neg[item]):
				self.genomes_to_sim.append(gen)
				self.genlens.append(l)
			for gen in self.project.coords_neg[item]:
				self.coords_files.append(gen)		

		self.items_to_sim = None
		self.make_prep()
	
	def make_db(self):	
		target_db_name = self.shared_db + "combined_database.db"
		
		target_fasta_set = self.shared_genomes + "combined_proteins_aa.fasta"
		target_fasta_set_nt = self.shared_genomes + "combined_proteins_nt.fasta"
		
		ma_file = self.multiple_alignment + "complete_multiple_alignment_aa.fasta"
		ma_file_nt = self.multiple_alignment + "complete_multiple_alignment_nt.fasta"
		
		#Join the reference proteins
		combiner = open(target_fasta_set, "w")
		#We want all of them. Align now, filter to pos/neg later.
		for base in self.project.targets:
			for genome in self.project.targets[base]:
				in_fasta = open(genome)
				for line in in_fasta:
					combiner.write(line)
				in_fasta.close()
		
		combiner.close()
		
		#Join the reference proteins
		combiner = open(target_fasta_set_nt, "w")
		#We want all of them. Align now, filter to pos/neg later.
		for base in self.project.targets:
			for genome in self.project.targets[base]:
				genome = genome.replace("_AA.fasta", "_nt.fasta")
				in_fasta = open(genome)
				for line in in_fasta:
					combiner.write(line)
				in_fasta.close()
		
		combiner.close()
		
		#craft a MA
		muscle_ma_command = ["muscle", "-in", target_fasta_set, "-out", ma_file]
		subprocess.run(muscle_ma_command)
		
		muscle_ma_command = ["muscle", "-in", target_fasta_set_nt, "-out", ma_file_nt]
		subprocess.run(muscle_ma_command)		
		
		#These should be nt only
		tree_file = os.path.normpath(self.multiple_alignment + "target_seq_tree.txt")
		tree_log_file = os.path.normpath(self.multiple_alignment + "fasttree_log.txt")
		#fasttree -log ftl.txt eee/shared_files/multiple_alignment/target_seq_MA.txt > ftt.txt
		#muscle_tree_command = ["muscle", "-maketree", "-in", self.ma, "-out", self.tree]
		#subprocess.run(muscle_tree_command)
		
		fasttree_command = ["fasttree", "-nt", "-gtr", "-log", tree_log_file, "-out", tree_file, ma_file_nt]
		subprocess.run(fasttree_command)
				
		makedb = ["diamond", "makedb", "--db", target_db_name,  "--in", target_fasta_set]
		print("Building Diamond database. Log information will follow.")
		subprocess.call(makedb)
		
		self.shared_db = target_db_name
	
	#Function for splitting comma-sep string triplets of numbers eg 90,100,110 as numeric list
	def parse_read_triplet(self, arg, default):
		try:
			rt = rt.split(",")
			rt = [int(r) for r in rt]
			rt = tuple(rt)
		except:
			rt = default
			rt = rt.split(",")
			rt = [int(r) for r in rt]
			rt = tuple(rt)
			
		return rt
	
	#Function for building out directories as part of a structure.
	def prep_dir(self, out_base, dir):
		if not os.path.exists(out_base + dir):
			os.mkdir(out_base + dir)
		return out_base + dir	
	
	#Function for reading a coordinates file to determine where on-target sequences begin and end.
	def extract_from_coords(self, cf):
		coord_starts = []
		coord_ends = []
		fh = open(cf)
		for line in fh:
			segs = line.strip().split()
			parent, prot_id, start, end, strand = segs[0], segs[1], int(segs[2]), int(segs[3]), segs[4]
			#my_coords[prot_id] = (start, end, strand)
			coord_starts.append(start)
			coord_ends.append(end)
		
		fh.close()
			
		#Switch to numpy.
		coord_starts = np.array(coord_starts)
		coord_ends = np.array(coord_ends)
		
		return coord_starts, coord_ends
	
	def make_prep(self):
		num = 1
		self.items_to_sim = []
		
		for prot, l, coord in zip(self.genomes_to_sim, self.genlens, self.coords_files):
			base = prot.split("/genomes/")
			
			cs, ce = self.extract_from_coords(coord)
			
			output = prot.split("/genomes/")
			out_base = output[0] + "/"
			read_name = output[1][:-6] + "_read_len_"
			
			raw = self.prep_dir(out_base, "raw_reads/")
			tagged = self.prep_dir(out_base, "tagged_reads/")
			aln = self.prep_dir(out_base, "aligned_reads/")
			bblog = self.prep_dir(out_base, "bbmap_log/")
			aln_log = self.prep_dir(out_base, "alignment_log/")
			
			for triplet in self.simlens:
				maxlen = triplet[2]
				#print("requested length", maxlen, "seqlen", l)
				if maxlen < l:
					next_item = one_protein(base[0], self.sim_arg_template, num, prot, triplet, self.shared_db, cs, ce)
					self.items_to_sim.append(next_item)
					
				num += 1
						
	def run_simulation(self):
		prog_bar = progress_tracker(total = len(self.items_to_sim), message = "Simulating metagenome and aligning reads. This will take some time.")
	
		pool = multiprocessing.Pool(self.threads)
		for result in pool.imap_unordered(run_sim, self.items_to_sim):
			prog_bar.update()
		pool.close()
		
class one_protein:
	def __init__(self, directory_base, read_sim_template, build_num, input_fasta, sim_length, 
				diamond_database, coord_starts, coord_ends):
				
		self.base = directory_base
		
		self.simulation_index = build_num
		self.simlen = sim_length
		
		self.inf = input_fasta
		
		self.output = self.inf.split("/genomes/")
		self.out_base = self.output[0] + "/"
		self.read_name = self.output[1][:-6] + "_read_len_"
		
		self.this_readlen = self.read_name + str(self.simlen[1])
		
		self.fastq = self.out_base + "raw_reads/"+ self.this_readlen + "_raw.fastq"
		self.fasta = self.out_base + "raw_reads/"+ self.this_readlen + "_raw.fasta"
		self.tagged = self.out_base + "tagged_reads/"+ self.this_readlen + "_tagged.fasta"
		self.aln_reads = self.out_base + "aligned_reads/"+ self.this_readlen + "_aligned_reads.fasta"
		self.aln_err = self.out_base + "alignment_log/" + self.this_readlen + "_DIAMOND_log.txt"
		
		self.template = read_sim_template
		
		#self.log_file = os.path.normpath(self.base + "/bbmap_log/" + self.input_base + "_log.txt")
		
		self.generation_log = ""
		
		self.coord_starts = coord_starts
		self.coord_ends = coord_ends
		
		self.dia = diamond_database
		
		self.final_read_count = 0
	
	def simulate_reads(self):		

		'''
		"minlength={min}",
		"midlength={med}",
		"maxlength={max}",
		"build={build_num}",
		"ref={input_genome}",
		"out={output_fastq}"
		'''
		
		next_gen = self.template.format(min = self.simlen[0],
										med = self.simlen[1],
										max = self.simlen[2],
										build_num = self.simulation_index,
										input_genome = self.inf,
										output_fastq = self.fastq)
										
		
		next_gen = next_gen.split()
		
		#print(next_gen)
		proc = subprocess.Popen(next_gen, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
		for line in proc.stdout:
			#Lines are imported in binary, this just converts to plaintext
			self.generation_log += line.decode()
			
		try:
			proc.stdout.close()
		except:
			pass
			
		try:
			proc.stderr.close()
		except:
			pass
			
		fh = open(self.out_base + "bbmap_log/"+self.this_readlen + "simlog.txt", "w")
		fh.write(self.generation_log)
		fh.close()
			
		try:
			shutil.rmtree("ref/genome/"+str(self.simulation_index))
		except:
			pass
		try:
			shutil.rmtree("ref/index/"+str(self.simulation_index))
		except:
			print("Couldn't remove ref/index/"+str(self.simulation_index))

	def fastq_to_fasta(self):
		if os.path.exists(os.path.normpath(self.fastq)):
			line_counter = 1
			fq = open(self.fastq)
			fa = open(self.fasta, "w")
			
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
				os.remove(self.fastq)
			except:
				print("Couldn't remove", self.fastq)
		else:
			print(self.fastq, "could not be generated. ROCkOut cannot continue without this file.")
		
	def tag_reads(self):
		pos_ct, off_ct, neg_ct = 0, 0, 0
		
		bbmap_match = r'SYN_(\d+)_(.+?(?=_))_(.+?(?=_))_\d_(.)_.+?(?=\.)\._(.+?(?=\$|\s)).+'
		fh = open(self.fasta)
		out = open(self.tagged, 'w')
		in_seq = False
		malformed = False
		for line in fh:
			if line.startswith(">"):
				#Reset
				malformed = False
				
				in_seq = False
				try:
					id, fr, to, comp, genome_id = re.search(bbmap_match, line).groups()
				except:
					print("Could not parse line.")
				
				fr, to = int(fr), int(to)
				mn, mx = min(fr, to), max(fr, to)
				
				#No interference with ROCker's split scheme.
				genome_id = genome_id.replace(';', '_')
				
				#We use the end to figure out if there's an overlap with the gene start to the gene start.
				start_window = np.searchsorted(self.coord_starts, mx, side = 'right')
				#We use the start against the ends to figure out if there's an overlap with the gene end.
				end_window = np.searchsorted(self.coord_ends, mn, side = 'left')
				
				tagged_name = ';'.join([id, str(mn), str(mx), comp, genome_id])
				
				#Read is malformed for one reason or another. Skip and keep going with the others.
				if len(tagged_name) == 0:
					print(tagged_name)
					#Skip writing until the next seq.
					malformed = True
					continue

				if (start_window - 1) == end_window:
					#The read falls inside a target gene window and we should tag it as on-target
					print(">" + tagged_name + ";Target", file = out)
					pos_ct += 1
				else:
					print(">" + tagged_name + ";Non_Target", file = out)
					off_ct += 1
				
			else:
				if not malformed:
					#Write seq.
					out.write(line)
				
		fh.close()
		out.close()
		
		return [pos_ct, off_ct, neg_ct]
		
	def align_reads(self):
		#DIAMOND run log
		diamond_err = open(self.aln_err, "w")
		'''
		Other opts to consider?
		--evalue (-e)            maximum e-value to report alignments (default=0.001)
		--min-score              minimum bit score to report alignments (overrides e-value setting)
		--id                     minimum identity% to report an alignment
		--query-cover            minimum query cover% to report an alignment
		--subject-cover          minimum subject cover% to report an alignment
		--threads
		'''
		
		#Original ruby command structure
		# diamond: '%1$sdiamond %2$s -q "%3$s" -d "%4$s" -a "%5$s.daa" -p %6$d' +
		#' -k 1 --min-score 20 --sensitive && %1$sdiamond view -a "%5$s"' +
		#' -o "%5$s"'},
		
		#-k = --max-target-seqs
		#-p = --threads
		#-a = --daa DIAMOND alignment archive (DAA) file
		
		align_command = ["diamond", 
						"blastx", #align nt to AA
						"--min-score", "20", #min bitscore = 20
						"--sensitive", #try harder to align distant matches
						"--max-target-seqs", "0", #report all discovered alignments for each read
						"--unal", "0", #report no unaligned reads
						"--outfmt", "6", #output as tabular blast
						"--db", self.dia, #use the shared database as the target
						"--query", self.tagged, #align the tagged reads
						"--out", self.aln_reads] #send to the pre-named file

		subprocess.call(align_command, stdout = diamond_err, stderr = subprocess.STDOUT)
		
		diamond_err.close()
		
		fh = open(self.aln_reads)
		out = open(self.aln_reads + "_filtered.txt", "w")
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
		
		try:
			os.remove(self.aln_reads)
		except:
			pass
			
		os.rename(self.aln_reads + "_filtered.txt", self.aln_reads)
		
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
		try:
			os.remove(file_name)
		except:
			pass
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
		try:
			os.remove(file_name)
		except:
			pass
		os.rename(file_name + "_filtered.txt", file_name)
	
	def clean_up_your_mess(self):
		hits = {}
		fh = open(self.aln_reads)
		for line in fh:
			read_id = line.split("\t")[0].split(";")[0]
			#The data in the dict doesn't matter.
			hits[read_id] = 0
		fh.close()
		
		#How many reads remain.
		self.final_read_count = len(hits)
		
		self.clean_raw_file(self.fasta, hits)
		self.clean_tagged_file(self.tagged, hits)
		
def build_project(parser, opts):
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
	'''
	threads = 1, 
	dir = None,
	#Read simulator arguments
	coverage_depth = 20.0, 
	snp_rate_per_base = 0.01, 
	insertion_rate = 0.01/19.0, 
	deletion_rate = 0.01/19.0, 
	short = "90,100,110", 
	standard = "180,200,220", 
	long = "270,300,330", 
	extra_long = "360,400,440"
	'''
	
	mn = read_manager(threads = threads,
					dir = project_directory,
					coverage_depth = coverage,
					snp_rate_per_base = snprate,
					insertion_rate = insrate,
					deletion_rate = delrate,
					short = s,
					standard = m,
					long = l,
					extra_long = xl)
					
	mn.run_simulation()
	
