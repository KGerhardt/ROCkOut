import sys
import os
import subprocess
import re
import numpy as np
import multiprocessing
import shutil

import argparse

import tarfile

from time import sleep

from .rocker_project_manager import project_manager
from .rocker_progress_tracker import progress_tracker
from .reads_labeller import protein_trawler


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
	ok_to_continue = read_simulator_item.fastq_to_fasta()
	if ok_to_continue:
		read_simulator_item.tag_reads()
		read_simulator_item.align_reads()
	
	read_simulator_item.clean_up_your_mess()
	return read_simulator_item.inf
	
def run_tf(target_finder_item):
	target_finder_item.run()
	#print(target_finder_item.probable_hits)
	return target_finder_item.probable_hits
			
class read_manager:
	def __init__(self,
				threads = 1, 
				dir = None,
				#Read simulator arguments
				coverage_depth = 20.0, 
				snp_rate_per_base = 0.01,
				insertion_rate = 0.01/19.0, 
				deletion_rate = 0.01/19.0,
				
				short = None,
				standard = None, 
				long = None, 
				extra_long = None,
				keep_reads = False):
				
		self.threads = threads
		self.rocker_directory = dir
		self.project = project_manager(directory = self.rocker_directory)
		
		self.project.parse_project_directory()
		self.project.parse_coords()
		
		self.project.parse_proteomes()
		self.proteome_jobs = []
		
		self.project.parse_genomes()
		self.project.check_genome_lengths()
		self.project.parse_targets()
		
		self.reference_lengths = {}
		self.reference_prots = []
		self.reference_prots_dict = {}
		
		self.read_labeller = None
		
		#Database stuff
		#raw = self.prep_dir(out_base, "raw_reads/")
		self.shared_loc = self.prep_dir(self.rocker_directory, "shared_files/")
		self.shared_db = self.prep_dir(self.rocker_directory, "shared_files/alignment_database/")
		self.shared_genomes = self.prep_dir(self.rocker_directory, "shared_files/combined_genomes/")
		
		self.multiple_alignment = self.prep_dir(self.rocker_directory, "shared_files/multiple_alignment/")
			
		#self.use_blast = use_blast
			
		self.dia_db = os.path.normpath(self.shared_db + "/combined_database_diamond.db")
		#self.bla_db = os.path.normpath(self.shared_db + "/combined_database_blast.db")
		
		self.make_db()

		#Try with illumina style reads
		self.short = short
		self.standard = standard
		self.long = long
		self.xl = extra_long
		
		self.simlens = [self.short, self.standard, self.long, self.xl]
		
		self.keep = keep_reads
		
		#Old format
		'''
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
		'''
		
		self.sim_arg_template = ' '.join(["randomreads.sh",
								"coverage="+str(coverage_depth),
								"snprate="+str(snp_rate_per_base),
								"insrate="+str(round(insertion_rate, 8)),
								"delrate="+str(round(deletion_rate, 8)),
								"simplenames=t",
								"minlength={min}",
								"maxlength={max}",
								"build={build_num}",
								"ref={input_genome}",
								"out={output_fastq}"
							])
							
	
		self.genomes_to_sim = []
		self.genlens = []
		self.coords_files = []
		self.proteome_files = []
		
		self.proteome_args = []
		
		self.coord_starts = []
		self.coord_ends = []
	
		for item in self.project.genomes_pos:
			for gen, l, coords, proteome in zip(self.project.genomes_pos[item], self.project.gen_lengths_pos[item], self.project.coords_pos[item], self.project.proteomes_pos[item]):
				self.genomes_to_sim.append(gen)
				self.genlens.append(l)
				cs, ce = self.extract_from_coords(coords)
				self.coord_starts.append(cs)
				self.coord_ends.append(ce)
				
				bn = os.path.basename(proteome)
				while bn != os.path.splitext(bn)[0]:
					bn = os.path.splitext(bn)[0]
				output_aln = "/".join([self.project.project_base, "positive", item, "proteome_vs_references", bn + ".protein_alignment.blast.txt"])
				output_coords = "/".join([self.project.project_base, "positive", item, "probable_target_coords", bn + ".probable_targets_coords.txt"])
				
				next_pa = (item, proteome, bn, output_aln, output_coords)	
				self.proteome_args.append(next_pa)
				
		
		for item in self.project.genomes_neg:
			for gen, l, coords, proteome in zip(self.project.genomes_neg[item], self.project.gen_lengths_neg[item], self.project.coords_neg[item], self.project.proteomes_neg[item]):
				self.genomes_to_sim.append(gen)
				self.genlens.append(l)
				cs, ce = self.extract_from_coords(coords)
				self.coord_starts.append(cs)
				self.coord_ends.append(ce)
				
				#Should consider getting only own proteins to exclude from probable targets.
				bn = os.path.basename(proteome)
				while bn != os.path.splitext(bn)[0]:
					bn = os.path.splitext(bn)[0]
				output_aln = "/".join([self.project.project_base, "negative", item, "proteome_vs_references", bn + ".protein_alignment.blast.txt"])
				output_coords = "/".join([self.project.project_base, "negative", item, "probable_target_coords", bn + ".coords.txt"])
				
				next_pa = (item, proteome, bn, output_aln, output_coords)	
				self.proteome_args.append(next_pa)
				
		self.per_ID_coords = []
		#self.find_probable_targets()
				
		self.items_to_sim = []
		#self.make_prep()
		
	def make_db(self):	
		target_db_name_dia = self.dia_db
		#target_db_name_bla = self.bla_db
		
		target_fasta_set = os.path.normpath(self.shared_genomes + "/combined_proteins_aa.fasta")
		target_fasta_set_nt = os.path.normpath(self.shared_genomes + "/combined_proteins_nt.fasta")
		
		ma_file = os.path.normpath(self.multiple_alignment + "/complete_multiple_alignment_aa.fasta")
		ma_file_nt = os.path.normpath(self.multiple_alignment + "/complete_multiple_alignment_nt.fasta")
		
		#Join the reference proteins
		cur_seqlen = 0
		cur_id = ""
		
		combiner = open(target_fasta_set, "w")
		#We want all of them. Align now, filter to pos/neg later.
		for base in self.project.targets:
			for genome in self.project.targets[base]:
				in_fasta = open(genome)
				for line in in_fasta:
					if line.startswith(">"):
						if cur_seqlen > 0:
							prot = cur_id.split("__")[3]
							self.reference_prots.append(prot)
							self.reference_lengths[cur_id] = cur_seqlen
						cur_id = line.strip().split()[0][1:]
						cur_seqlen = 0
					else:
						cur_seqlen += len(line.strip())
						
					combiner.write(line)
				in_fasta.close()
				
				if cur_seqlen > 0:
					prot = cur_id.split("__")[3]
					self.reference_prots.append(prot)
					self.reference_lengths[cur_id] = cur_seqlen
		
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
		muscle_ma_command_aa = ["muscle", "-align", target_fasta_set, "-output", ma_file]
		muscle_ma_command_nt = ["muscle", "-align", target_fasta_set_nt, "-output", ma_file_nt]	
		
		#These should be nt only
		tree_file = os.path.normpath(self.multiple_alignment + "/target_seq_tree.txt")
		tree_log_file = os.path.normpath(self.multiple_alignment + "/fasttree_log.txt")

		
		fasttree_command = ["fasttree", "-nt", "-gtr", "-log", tree_log_file, "-out", tree_file, ma_file_nt]
		
		#Subprocess calls for MA+tree
		
		
		subprocess.call(muscle_ma_command_aa)
		subprocess.call(muscle_ma_command_nt)	
		subprocess.call(fasttree_command)
				
				
		makedb = ["diamond", "makedb", "--db", target_db_name_dia,  "--in", target_fasta_set]
		print("Building Diamond database. Log information will follow.")
		process = subprocess.call(makedb)
		
		#try:
			#makeblastdb -in <reference.fa> -dbtype nucl -parse_seqids -out <database_name> -title "Database title"
			#makedb = ["makeblastdb", "-in", target_fasta_set, "-parse_seqids", "-dbtype", "prot", "-out", target_db_name_bla]
			#print(" ".join(makedb))
			#print("Building BLAST database for positive targets. Log information will follows.")
			#process = subprocess.call(makedb)
		#except:
			#print("Couldn't make BLAST database of positive targets!")
		
		#if self.use_blast:
			#self.shared_db = target_db_name_bla
		#else:
		
		self.shared_db = target_db_name_dia
	
	#Function for building out directories as part of a structure.
	def prep_dir(self, out_base, dir):
		dirpath = os.path.normpath(out_base + "/" + dir)
		if not os.path.exists(dirpath):
			os.mkdir(dirpath)
		return dirpath
	
	#Function for reading a coordinates file to determine where on-target sequences begin and end.
	def extract_from_coords(self, cf):
		coord_starts = []
		coord_ends = []

		fh = open(cf)
		for line in fh:
			segs = line.strip().split()
			parent_uniprot_id, prot_id, tran_id, start, end, strand = segs[0], segs[1], segs[2], int(segs[3]), int(segs[4]), segs[5],
			coord_starts.append(start)
			coord_ends.append(end)
			
			#real parent since the name in the file seems to be wrong for now... fix in 0_download
			parent_uniprot_id = os.path.basename(os.path.dirname(os.path.dirname(cf)))
			
			#We probably need this to happen on a per-file basis with a third return
			self.reference_prots.append(prot_id)
			parent_genome = prot_id.split(".")[0]
			
			if parent_genome not in self.reference_prots_dict:
				self.reference_prots_dict[parent_genome] = {}
			
			#print(cf, parent_genome, prot_id)
			self.reference_prots_dict[parent_genome][parent_uniprot_id] = (start, end, prot_id,)
		
		fh.close()
			
		#Switch to numpy.
		#coord_starts = np.array(coord_starts)
		#coord_ends = np.array(coord_ends)
		
		return coord_starts, coord_ends
		
	def find_probable_targets(self):
		self.reference_prots = set(self.reference_prots) #Render this searchable.
	
		for next_pa in self.proteome_args:
			item, proteome, bn, output_aln, output_coords = next_pa[0], next_pa[1], next_pa[2], next_pa[3], next_pa[4]
			#next_pa = (item, proteome, bn, output_aln, output_coords)
			next_set = probable_target_finder(item, proteome, self.dia_db, output_aln, output_coords, self.reference_lengths, self.reference_prots,)
			self.proteome_jobs.append(next_set)
		
		for directory in self.project.positive:
			dir1 = "/".join([self.project.project_base, "positive", directory, "probable_target_coords"])
			dir1 = os.path.normpath(dir1)
			dir2 = "/".join([self.project.project_base, "positive", directory, "proteome_vs_references"])
			dir2 = os.path.normpath(dir2)
			if not os.path.exists(dir1):
				os.mkdir(dir1)
			if not os.path.exists(dir2):
				os.mkdir(dir2)
			
		for directory in self.project.negative:
			dir1 = "/".join([self.project.project_base, "negative", directory, "probable_target_coords"])
			dir1 = os.path.normpath(dir1)
			dir2 = "/".join([self.project.project_base, "negative", directory, "proteome_vs_references"])
			dir2 = os.path.normpath(dir2)
			if not os.path.exists(dir1):
				os.mkdir(dir1)
			if not os.path.exists(dir2):
				os.mkdir(dir2)
		
	def run_target_finder(self):
		prog_bar = progress_tracker(total = len(self.proteome_jobs), message = "Finding unlabeled probable targets.")
		#print("")
		#print("Finding unlabeled probable targets. This will take a bit.")
		#print("")
		final_results = {}
		
	
		pool = multiprocessing.Pool(self.threads)
		for result in pool.imap(run_tf, self.proteome_jobs):
			if result is not None:
				#self.probable_hits[self.reference_genome][qid] = self.protein_metadata[qid]
				for reference_genome in result:
					if reference_genome not in final_results:
						final_results[reference_genome] = {}
					for qid in result[reference_genome]:
						final_results[reference_genome][qid] = result[reference_genome][qid]
						
				#self.per_ID_coords.append(result)
			prog_bar.update() #Cannot update a progress bar using map?
		pool.close()
		
		#Collect per ID coords to match regular coords.
		self.per_ID_coords = final_results
		
		self.read_labeller = protein_trawler(self.rocker_directory)
		self.read_labeller.load_labels()
		
	def make_prep(self):
		num = 1 #Arbitrary read simulation index for bbtools.
		
		
		for prot, l, cs, ce in zip(self.genomes_to_sim, self.genlens, self.coord_starts, self.coord_ends):
			base = prot.split("/genomes/")
						
			output = prot.split("/genomes/")
			out_base = output[0] + "/"
			read_name = output[1][:-6] + "_read_len_"
			
			raw = self.prep_dir(out_base, "raw_reads/")
			tagged = self.prep_dir(out_base, "tagged_reads/")
			aln = self.prep_dir(out_base, "aligned_reads/")
			bblog = self.prep_dir(out_base, "bbmap_log/")
			aln_log = self.prep_dir(out_base, "alignment_log/")
			#prob_tgt = self.prep_dir(out_base, "probable_target_coords/")
			
			for triplet in self.simlens:
				#maxlen = triplet[2]
				maxlen = triplet[1]
				if maxlen < l: #if the requested read length is less than the genome's length, otw skip
					#seen_files.append(input_protein)
					next_item = one_protein(base[0], self.sim_arg_template, num, prot, 
					triplet, self.shared_db, cs, ce, 
					self.reference_prots_dict, #Real coords
					self.per_ID_coords, self.read_labeller, self.keep,) #Homology coords
					self.items_to_sim.append(next_item)
					
				num += 1
				
	def run_simulation(self):
		#if self.use_blast:
		#	prog_bar = progress_tracker(total = len(self.items_to_sim), message = "Simulating metagenome and aligning reads using BLASTx. This will take a long time.")
		#else:
		
		prog_bar = progress_tracker(total = len(self.items_to_sim), message = "Simulating metagenome and aligning reads using DIAMOND. This will take some time.")
		
		
		pool = multiprocessing.Pool(self.threads)
		for result in pool.imap_unordered(run_sim, self.items_to_sim):
			prog_bar.update()
		pool.close()
		pool.join()
		
		sleep(1.5)
		
	
class probable_target_finder:
	#probable_target_finder(proteome, self.shared_db, output_aln, output_coords, self.use_blast)
	def __init__(self, parent, proteome_file, diamond_database, out_aln, out_coords, reference_lengths, ref_proteins):
		self.parent = parent
			
		#self.in_coords = reference_coords
		self.proteome = proteome_file
		
		self.reference_genome = os.path.basename(proteome_file)
		while self.reference_genome != os.path.splitext(self.reference_genome)[0]:
			self.reference_genome = os.path.splitext(self.reference_genome)[0]
		
		self.db = diamond_database
		self.out_aln = out_aln
		self.out_coords = out_coords

		self.aligned_prots = None
		self.probable_target_coords = None
		
		self.reflens = reference_lengths
		self.ref_prots = ref_proteins
		
		self.protein_metadata = {}
		self.probable_hits = {}
		
	def read_protein_locs(self):
		#self.protein_metadata = {}
		fh = open(self.proteome)
		for line in fh:
			if line.startswith(">"):
				relevant = line.strip().split()[0]
				relevant = relevant.split(";")
				id = relevant[0][1:]
				start = int(relevant[1])
				end = int(relevant[2])
				
				strand = relevant[3]
				
				self.protein_metadata[id] = [start, end, strand]
			else:
				pass
	
	def align_prots(self):
		align_command = ["diamond",
				"blastp", #protein aln
				"--sensitive", #try harder to align distant matches
				"--max-target-seqs", "0", #report all discovered alignments for each read
				"--unal", "0", #report no unaligned reads
				"--outfmt", "6", 
							"qseqid", 
							"sseqid", 
							"pident", 
							"length", 
							"mismatch", 
							"gapopen", 
							"qstart",
							"qend",
							"sstart", 
							"send", 
							"evalue", 
							"bitscore", 
							"qlen", 
							"slen", #output as tabular blast with seqlen
				"--db", self.db, #the reference protein db
				"--query", self.proteome, #full proteome
				"--out", self.out_aln] #This output alignment file
		
		
		process = subprocess.call(align_command, stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
		
	def parse_alignments_to_coords(self):
		
		outwriter = open(self.out_coords, "w")
		print("origin_genome", "query_parent", "query_id", "qstart", "qend", "qstrand", "aligns_to_target", "pct_aln_to_tgt", "pct_ID_to_tgt", sep = "\t", file = outwriter)
		passing_targets = []
		fh = open(self.out_aln)
		for line in fh:
			#print(line)
			segs = line.strip().split("\t")
			query = segs[0].split(";")
			qid = query[0]
			if qid in self.ref_prots:
				label = "maps_to_reference_protein"
			else:
				label = "discovered_match_through_homology"
				
			if qid not in self.ref_prots: #We don't need to check already referenced items			
				target = segs[1]
				tgt_length = self.reflens[target]
				pid = float(segs[2])
				aln_len = int(segs[3])
				pct_aln = round((aln_len/tgt_length) * 100, 2)
								
				if self.reference_genome not in self.probable_hits:
					self.probable_hits[self.reference_genome] = {}
				
				if pct_aln >= 70.0 and pid >= 50.0: #70% alignment length and 50 pct. ID
					self.probable_hits[self.reference_genome][qid] = self.protein_metadata[qid]
					#print(self.reference_genome, self.parent, qid, *self.protein_metadata[qid], target, pct_aln, pid, label, sep = "\t", file = outwriter)
			
		outwriter.close()
		fh.close()
		
		
		
	def run(self):
		try:
			self.read_protein_locs()
			self.align_prots()
			self.parse_alignments_to_coords()
			if len(self.probable_hits) == 0:
				self.probable_hits = None
		except:
			if os.path.exists(self.out_aln):
				os.remove(self.out_aln)
			if os.path.exists(self.out_coords):
				os.remove(self.out_coords)
				
			self.probable_hits = None #This essentially is a result of a GFF being empty, and we want to skip this.
	
def check_process_started(file_path):
	filesize = 0
	hang_count = 0
	hanging = False
	while not hanging:
		if filesize > 0:
			break
		if os.path.exists(file_path): 
			output_size = os.stat(file_path)
			output_size = output_size.st_size
			
		else:#not finding the file is equivalent to a hang.
			output_size = filesize 
			
		if output_size > filesize:
			filesize = output_size
		else:
			hang_count += 1
			if hang_count > 6:
				hanging = True
			
		sleep(30)
		
	return hanging

	
class one_protein:
	def __init__(self, directory_base, read_sim_template, build_num, input_fasta, sim_length, 
				alignment_database, coord_starts, coord_ends, target_coords, probable_coords, 
				read_labeller, keep_reads):
				
		self.base = directory_base
		
		self.keep = keep_reads
		
		self.labeller = read_labeller
		
		self.simulation_index = build_num
		self.simlen = sim_length
		
		self.inf = input_fasta
		
		self.output = self.inf.split("/genomes/")
		self.out_base = self.output[0] + "/"
		self.ref_genome = self.output[1][:-6]
		self.read_name = self.output[1][:-6] + "_read_len_"
		
		#Average of the readlens, curtailed to an int
		self.this_readlen = self.read_name + str(int((self.simlen[0] + self.simlen[1]) / 2))
		
		self.fastq = self.out_base + "raw_reads/"+ self.this_readlen + "_raw.fastq"
		self.fasta = self.out_base + "raw_reads/"+ self.this_readlen + "_raw.fasta"
		self.tagged = self.out_base + "tagged_reads/"+ self.this_readlen + "_tagged.fasta"
		self.aln_reads = self.out_base + "aligned_reads/"+ self.this_readlen + "_aligned_reads.blast.txt"
		
		#self.use_blast = use_blast
		
		#if self.use_blast:
		#	self.aln_err = self.out_base + "alignment_log/" + self.this_readlen + ".BLAST_log.txt"
		#else:
		
		self.aln_err = self.out_base + "alignment_log/" + self.this_readlen + ".DIAMOND_log.txt"
		
		
		self.template = read_sim_template
		
		#self.log_file = os.path.normpath(self.base + "/bbmap_log/" + self.input_base + "_log.txt")

		
		self.generation_log = ""
		
		self.reference_tgts = target_coords
		self.homol_tgts = probable_coords
		
		self.own_starts = []
		self.own_ends = []
		
		self.foreign_starts = []
		self.foreign_ends = []
		self.foreign_names = []
		
		self.homology_starts = []
		self.homology_ends = []
		self.homology_names = []
		
		self.collect_starts()
		
		#self.coord_starts = coord_starts
		#self.coord_ends = coord_ends
		
		#self.probable_starts = prob_starts
		#self.probable_ends = prob_ends
		
		self.db = alignment_database

		self.final_read_count = 0
	
	def collect_starts(self):
		own_genome = self.ref_genome
		parent_uniprot_id = os.path.basename(os.path.dirname(os.path.dirname(self.inf)))
		
		if own_genome in self.reference_tgts:
			for protein_id in self.reference_tgts[own_genome]:
				tup = self.reference_tgts[own_genome][protein_id]
				start, end = tup[0], tup[1]
				if protein_id == parent_uniprot_id:
					self.own_starts.append(start)
					self.own_ends.append(end)
				else:
					ft = tup[2]
					self.foreign_starts.append(start)
					self.foreign_ends.append(end)
					self.foreign_names.append(ft)
				
		if own_genome in self.homol_tgts:
			for protein_id in self.reference_tgts[own_genome]:
				tup = self.reference_tgts[own_genome][protein_id]
				start, end = tup[0], tup[1]
				
				self.homology_starts.append(start)
				self.homology_ends.append(end)
				self.homology_names.append(protein_id)
				
	def simulate_reads(self):
		#print("Simulating", self.inf, "read_length", int((self.simlen[1] + self.simlen[0])/2))
		#print("")
	
		next_gen = self.template.format(min = self.simlen[0],
										max = self.simlen[1],
										build_num = self.simulation_index,
										input_genome = self.inf,
										output_fastq = self.fastq)

			
		#print(next_gen)
		'''
		next_gen = next_gen.split()
		
		read_size = self.simlen[1] + 1
		gensize = 0
		if os.path.exists(self.inf):
			with open(self.inf) as fh:
				for line in fh:
					if not line.startswith(">"):
						gensize += len(line.strip())
						
		if read_size < gensize:
		'''
		
		next_gen = next_gen.split()
		
		proc = subprocess.Popen(next_gen, stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
		sleep(5)
		hanging = check_process_started(self.fastq)
		
		if hanging:
			print(self.aln_reads, "was hanging and can't continue.")
			print("")
			try:
				proc.terminate()
				proc.kill()
			except:
				pass
			
			try:
				os.remove(os.path.normpath(self.fastq))
			except:
				pass
			
		else:
			proc.wait()
		
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
		
		fh = open(self.out_base + "bbmap_log/"+self.this_readlen + ".simlog.txt", "w")
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
		has_fasta = True
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
			has_fasta = False
			print(self.fastq, "could not be generated.")
			print("The reads may have failed to simulate or some other issue occurred.")
			print("This file will be skipped, but ROCkOut will continue.")
			print("")
		
		return has_fasta
		
	def calculate_read_overlap(self, read_start, read_end, gene_start, gene_end):
		read_positions = np.arange(read_start, read_end, dtype = np.int32)
		gene_positions = np.arange(gene_start, gene_end, dtype = np.int32)
		
		#I'm not entirely sure this math is right...
		overlap_size = len(np.intersect1d(read_positions, gene_positions)) + 1
		
		return overlap_size
		
	def tag_reads(self):
		pos_ct, prob_ct, off_ct = 0, 0, 0
		
		previously_tagged = {}
		
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
				#print(line)
				
				#read_id, fr, to, comp, genome_id = re.search(bbmap_match, line).groups()
				
				try:
					read_id, fr, to, comp, genome_id = re.search(bbmap_match, line).groups()
				except:
					print("Could not parse line.")
				
				fr, to = int(fr), int(to)
				#mn, mx = min(fr, to), max(fr, to)
				mn, mx = min([fr, to]), max([fr, to])
				
				#No interference with ROCker's split scheme.
				genome_id = genome_id.replace(';', '_')
				
				#label_read(self, genome, read_start, read_end, read_ID)
				printable_name = self.labeller.tag_read(self.ref_genome, mn, mx, read_id, comp, genome_id)
				
				#print(printable_name)
				print(printable_name, file = out)
				
			else:
				if not malformed:
					#Write seq.
					out.write(line)
				
		fh.close()
		out.close()
		
		return [pos_ct, prob_ct, off_ct]
		
	def align_reads(self):
		#if self.use_blast:
		#	self.align_reads_blast()
		#else:
		self.align_reads_diamond()
		
	def align_reads_diamond(self):
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
						"--outfmt", "6", 
							"qseqid", 
							"sseqid", 
							"pident", 
							"length", 
							"mismatch", 
							"gapopen", 
							"qstart",
							"qend",
							"sstart", 
							"send", 
							"evalue", 
							"bitscore", 
							"qlen", 
							"slen", 
						#"--outfmt", "6", #output as tabular blast
						"--db", self.db, #use the shared database as the target
						"--query", self.tagged, #align the tagged reads
						"--out", self.aln_reads,
						'--threads', '1'] #send to the pre-named file
		
		#align_command = " ".join(align_command)
		#print(align_command)
		#print("")
		
		#os.system(align_command)
		#process = subprocess.call(align_command, stdout = diamond_err, stderr = subprocess.STDOUT)
		process = subprocess.Popen(align_command, stdout = diamond_err, stderr = subprocess.STDOUT)
		sleep(5)
		hanging = check_process_started(self.aln_reads)
		
		if hanging:
			print(self.aln_reads, "was hanging and can't continue.")
			process.terminate()
		else:
			process.wait()
		
		diamond_err.close()
		
		#print("Filtering")
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
				
			#Extract the additional features here
			read_id_and_info = segs[0]
			read_id_and_info = read_id_and_info.split(";")
			#overlap_bp = int(read_id_and_info[3])
			overlap_bp = read_id_and_info[3]
			percent_of_read_overlapping_the_gene = read_id_and_info[4]
			
			qlen = int(segs[12])
			
			#percent_of_read_overlapping_the_gene = 100*overlap_bp/qlen
			#percent_of_read_overlapping_the_gene = round(percent_of_read_overlapping_the_gene, 2)
				
			'''	
			read_id_and_info = ";".join([read_id_and_info[0],
										read_id_and_info[1],
										read_id_and_info[2],
										read_id_and_info[4],
										read_id_and_info[5],
										read_id_and_info[6]])
			'''							
			read_id_and_info = ";".join(read_id_and_info)
			
			#Only occurs when the two continues do not.
			#out.write()
			print(read_id_and_info, *segs[1:], overlap_bp, percent_of_read_overlapping_the_gene, sep = "\t", file = out)
		
		out.close()
		fh.close()
		
		try:
			os.remove(self.aln_reads)
		except:
			pass
			
		os.rename(self.aln_reads + "_filtered.txt", self.aln_reads)
		
	#No longer used
	def align_reads_blast(self):
		#DIAMOND run log
		err = open(self.aln_err, "w")
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
		
		align_command = ["blastx", #align nt to AA
						"-db", self.db,#use the shared database as the target
						#"--min-score", "20", #min bitscore = 20
						#"--sensitive", #try harder to align distant matches
						#"-max_target_seqs", "1000", #report all discovered alignments for each read
						#"--unal", "0", #report no unaligned reads
						"-outfmt", "'6", 
							"qseqid", 
							"sseqid", 
							"pident", 
							"length", 
							"mismatch", 
							"gapopen", 
							"qstart",
							"qend",
							"sstart", 
							"send", 
							"evalue", 
							"bitscore", 
							"qlen", 
							"slen'", #output as tabular blast
						"-query", self.tagged, #align the tagged reads
						"-out", self.aln_reads] #send to the pre-named file
						
		
		process = subprocess.call(align_command, stdout = err, stderr = subprocess.STDOUT)
		
		
		err.close()
		
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
				
			#Extract the additional features here
			read_id_and_info = segs[0]
			read_id_and_info = read_id_and_info.split(";")
			overlap_bp = int(read_id_and_info[3])
			qlen = int(segs[12])
			percent_of_read_overlapping_the_gene = round(overlap_bp/qlen, 2)
				
			read_id_and_info = ";".join([read_id_and_info[0],
										read_id_and_info[1],
										read_id_and_info[2],
										read_id_and_info[4],
										read_id_and_info[5],
										read_id_and_info[6]])
										
			#Only occurs when the two continues do not.
			#out.write()
			print(read_id_and_info, *segs[1:], overlap_bp, percent_of_read_overlapping_the_gene, sep = "\t", file = out)
		
		out.close()
		fh.close()
		
		try:
			os.remove(self.aln_reads)
		except:
			pass
			
		os.rename(self.aln_reads + "_filtered.txt", self.aln_reads)
		
	def clean_raw_file(self, file_name, hits):
		if os.path.exists(file_name):
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
		if os.path.exists(file_name):
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
		if os.path.exists(self.aln_reads):
			fh = open(self.aln_reads)
			for line in fh:
				read_id = line.split("\t")[0].split(";")[0]
				#The data in the dict doesn't matter.
				hits[read_id] = 0
			fh.close()
		
		#How many reads remain.
		self.final_read_count = len(hits)
		
		if self.final_read_count == 0:
			print("")
			print(self.base, "had no reads successfully align for alignment length", int((self.simlen[0]+self.simlen[1])/2))
			print("This file and its generated reads will be removed.")
			print("")
			
			if os.path.exists(self.aln_reads):
				os.remove(self.aln_reads)
			if os.path.exists(self.fasta):
				os.remove(self.fasta)
			if os.path.exists(self.tagged):
				os.remove(self.tagged)
		else:
			if not self.keep:
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
		
	#s, m, l, xl = opts.short, opts.med, opts.long, opts.xl
	
	sl = opts.sl # short-lower
	su = opts.su # short-upper
	
	ml = opts.ml # short-lower
	mu = opts.mu # short-upper
	
	ll = opts.ll
	lu = opts.lu
	
	xll = opts.xll
	xlu = opts.xlu
	
	keep_reads = opts.keep_all
	
	coverage = float(opts.cov)
	snprate = float(opts.snps)
	insrate = opts.insrate
	delrate = opts.delrate
	
	#do_blast = opts.use_blast
	dia_sens = opts.dia_sens
	
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

	mn = read_manager(threads = threads,
					dir = project_directory,
					coverage_depth = coverage,
					snp_rate_per_base = snprate,
					insertion_rate = insrate,
					deletion_rate = delrate,
					short = [sl, su],
					standard = [ml, mu],
					long = [ll, lu],
					extra_long = [xll, xlu],
					keep_reads = keep_reads)
	
	
	mn.find_probable_targets()
	mn.run_target_finder()
	
	
	mn.make_prep()	
	mn.run_simulation()
	
		#Check to see if the ref/genomes and ref/index dirs are empty. Remove them if they are
	if os.path.exists("ref/"):
		if os.path.exists("ref/genome/"):
			if len(os.listdir("ref/genome/")) == 0:
				shutil.rmtree("ref/genome/")
		if os.path.exists("ref/index"):
			if len(os.listdir("ref/index/")) == 0:
				shutil.rmtree("ref/index/")
		if len(os.listdir("ref/")) == 0:
			shutil.rmtree("ref/")

