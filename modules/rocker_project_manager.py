import sys
import os

def listdir_mac_purge(directory):
	files = os.listdir(directory)
	purged = []
	for f in files:
		if f != ".DS_Store":
			purged.append(f)
			
	return purged

#We'll use this template multiple times
class project_manager:
	def __init__(self, directory = None, threads = 1):
		self.project_base = directory
		self.threads = threads
		
		#Protein folders
		self.positive = {}
		self.negative = {}
		
		#Coordinate files for target gene
		self.coords_pos = {}
		self.coords_neg = {}
		
		#Proteome AA FASTA files
		self.proteomes_pos = {}
		self.proteomes_neg = {}
		self.prob_tgts_pos = {}
		self.prob_tgts_neg = {}
		
		#GFF files for target genomes
		self.gffs_pos = {}
		self.gffs_neg = {}
		
		#fastas of the original genomes
		self.genomes_pos = {}
		self.genomes_neg = {}
		self.gen_lengths_pos = None
		self.gen_lengths_neg = None
		
		#fastq reads generated by randomreads.sh for the genomes
		self.read_fastqs_pos = {}
		self.read_fastqs_neg = {}
		
		#fasta translations of the fastqs
		self.read_fastas_pos = {}
		self.read_fastas_neg = {}
		
		#renamed/tagged reads for the fasta translations
		self.tagged_reads_pos = {}
		self.tagged_reads_neg = {}
		
		#Blast/diamond alignments of the tagged reads
		self.alignments_pos = {}
		self.alignments_neg = {}
		
		self.targets = {}
		self.targets_nt = {}
		self.active_targets = None
		
		self.mult_aln_base = None
		self.mult_aln_files = {}
		
		self.outputs_base = None
		self.rocker_filter = None
		self.targets_blast = None
		self.targets_dia = None
		
	#List out proteins in a ROCker directory
	def parse_project_directory(self):
		if self.project_base is None:
			print("ROCker requires an active porject directory with at least a positive set.")
			return None
		else:
			pos = os.path.normpath(self.project_base + "/positive")
			neg = os.path.normpath(self.project_base + "/negative")
			if os.path.exists(pos):
				prots = listdir_mac_purge(pos)
				for p in prots:
					self.positive[p] = os.path.normpath(pos + "/" + p)
			if os.path.exists(neg):
				prots = listdir_mac_purge(neg)
				for p in prots:
					self.negative[p] = os.path.normpath(neg + "/" + p)
		
	#Collect FASTA-format nucleotide sequences for genomes containing the target proteins
	def parse_genomes(self):
		relevant_name = "genomes"
		for base in self.positive:
			relevant_directory = os.path.normpath(self.positive[base] + "/" + relevant_name) + "/"
			if os.path.exists(relevant_directory):
				self.genomes_pos[base] = [relevant_directory + c for c in listdir_mac_purge(relevant_directory)]
		for base in self.negative:
			relevant_directory = os.path.normpath(self.negative[base] + "/" + relevant_name) + "/"
			if os.path.exists(relevant_directory):
				self.genomes_neg[base] = [relevant_directory + c for c in listdir_mac_purge(relevant_directory)]
		
	def check_genome_lengths(self):
		self.gen_lengths_pos = {}
		for base in self.genomes_pos:
			self.gen_lengths_pos[base] = []
			for genome in self.genomes_pos[base]:
				self.gen_lengths_pos[base].append(self.get_gen_len(genome))
				
		self.gen_lengths_neg = {}
		for base in self.genomes_neg:
			self.gen_lengths_neg[base] = []
			for genome in self.genomes_neg[base]:
				self.gen_lengths_neg[base].append(self.get_gen_len(genome))
		
	def get_gen_len(self, file):
		tot_len = 0
		fh = open(file)
		for line in fh:
			if line.startswith(">"):
				continue
			else:
				tot_len += len(line.strip())
		fh.close()
		return tot_len
		
	#Locate GFFs for possible use in selecting target proteins without use of the ROCker downloader. Not currently used.
	def parse_gffs(self):
		relevant_name = "gffs"
		for base in self.positive:
			relevant_directory = os.path.normpath(self.positive[base] + "/" + relevant_name) + "/"
			if os.path.exists(relevant_directory):
				self.gffs_pos[base] = [relevant_directory + c for c in listdir_mac_purge(relevant_directory)]
		for base in self.negative:
			relevant_directory = os.path.normpath(self.negative[base] + "/" + relevant_name) + "/"
			if os.path.exists(relevant_directory):
				self.gffs_neg[base] = [relevant_directory + c for c in listdir_mac_purge(relevant_directory)]
		
	#Collect coordinate files for each identified protein; protein ID and location within genome files per the above parse_genomes()
	def parse_coords(self):
		relevant_name = "coords"
		for base in self.positive:
			relevant_directory = os.path.normpath(self.positive[base] + "/" + relevant_name) + "/"
			if os.path.exists(relevant_directory):
				self.coords_pos[base] = [relevant_directory + c for c in listdir_mac_purge(relevant_directory)]
		for base in self.negative:
			relevant_directory = os.path.normpath(self.negative[base] + "/" + relevant_name) + "/"
			if os.path.exists(relevant_directory):
				self.coords_neg[base] = [relevant_directory + c for c in listdir_mac_purge(relevant_directory)]
	
	#collect proteome files for detection of probable target sequences
	def parse_proteomes(self):
		relevant_name = "proteome"
		for base in self.positive:
			relevant_directory = os.path.normpath(self.positive[base] + "/" + relevant_name) + "/"
			if os.path.exists(relevant_directory):
				self.proteomes_pos[base] = [relevant_directory + c for c in listdir_mac_purge(relevant_directory)]
		for base in self.negative:
			relevant_directory = os.path.normpath(self.negative[base] + "/" + relevant_name) + "/"
			if os.path.exists(relevant_directory):
				self.proteomes_neg[base] = [relevant_directory + c for c in listdir_mac_purge(relevant_directory)]
		
	def parse_probable_targets(self):
		relevant_name = "probable_target_coords"
		for base in self.positive:
			relevant_directory = os.path.normpath(self.positive[base] + "/" + relevant_name) + "/"
			if os.path.exists(relevant_directory):
				self.prob_tgts_pos[base] = [relevant_directory + c for c in listdir_mac_purge(relevant_directory)]
		for base in self.negative:
			relevant_directory = os.path.normpath(self.negative[base] + "/" + relevant_name) + "/"
			if os.path.exists(relevant_directory):
				self.prob_tgts_neg[base] = [relevant_directory + c for c in listdir_mac_purge(relevant_directory)]
		
	#Collect raw, FASTQ-format reads as generated by BBTools
	def parse_raw_reads(self):
		relevant_name = "raw_reads"
		for base in self.positive:
			relevant_directory = os.path.normpath(self.positive[base] + "/" + relevant_name) + "/"
			if os.path.exists(relevant_directory):
				self.read_fastas_pos[base] = [relevant_directory + c for c in listdir_mac_purge(relevant_directory)]
		for base in self.negative:
			relevant_directory = os.path.normpath(self.negative[base] + "/" + relevant_name) + "/"
			if os.path.exists(relevant_directory):
				self.read_fastas_neg[base] = [relevant_directory + c for c in listdir_mac_purge(relevant_directory)]

	#Collect ROCker-tagged reads that have annotations used by rocker in the deflines for each seq.
	def parse_tagged_reads(self):
		relevant_name = "tagged_reads"
		for base in self.positive:
			relevant_directory = os.path.normpath(self.positive[base] + "/" + relevant_name) + "/"
			if os.path.exists(relevant_directory):
				self.tagged_reads_pos[base] = [relevant_directory + c for c in listdir_mac_purge(relevant_directory)]
		for base in self.negative:
			relevant_directory = os.path.normpath(self.negative[base] + "/" + relevant_name) + "/"
			if os.path.exists(relevant_directory):
				self.tagged_reads_neg[base] = [relevant_directory + c for c in listdir_mac_purge(relevant_directory)]				
		
	#Collect valid taget proteins for use in creating a MA of prots and for aligning short reads
	def parse_targets(self):
		relevant_name = "target_protein"
		for base in self.positive:
			relevant_directory = os.path.normpath(self.positive[base] + "/" + relevant_name) + "/"
			if os.path.exists(relevant_directory):
				files = listdir_mac_purge(relevant_directory)
				for c in files:
					if c.endswith("AA.fasta"):
						if base not in self.targets:
							self.targets[base] = []
						self.targets[base].append(relevant_directory + c)
					else:
						if base not in self.targets_nt:
							self.targets_nt[base] = []
						self.targets_nt[base].append(relevant_directory + c)
		
		for base in self.negative:
			relevant_directory = os.path.normpath(self.negative[base] + "/" + relevant_name) + "/"
			if os.path.exists(relevant_directory):
				files = listdir_mac_purge(relevant_directory)
				for c in files:
					if c.endswith("AA.fasta"):
						if base not in self.targets:
							self.targets[base] = []
						self.targets[base].append(relevant_directory + c)
					else:
						if base not in self.targets_nt:
							self.targets_nt[base] = []
						self.targets_nt[base].append(relevant_directory + c)
						
		if self.active_targets is None:
			self.active_targets = self.targets
		
	#Takes a list of protein names from the positive directory and sets the active targets to this set. Used for model refinement. 
	def set_targets(self, protein_base_list):
		self.active_targets = {}
		for p in protein_base_list:
			self.active_targets[base] = self.targets[base]
			
	#Collect aligned reads for use in visualization/model refinement and most-discriminant cutoff identification.
	def parse_aligns(self):
		relevant_name = "aligned_reads"
		for base in self.positive:
			relevant_directory = os.path.normpath(self.positive[base] + "/" + relevant_name) + "/"
			if os.path.exists(relevant_directory):
				self.alignments_pos[base] = [relevant_directory + c for c in listdir_mac_purge(relevant_directory)]
		for base in self.negative:
			relevant_directory = os.path.normpath(self.negative[base] + "/" + relevant_name) + "/"
			if os.path.exists(relevant_directory):
				self.alignments_neg[base] = [relevant_directory + c for c in listdir_mac_purge(relevant_directory)]		
				
	def parse_multiple_alignment(self):
		self.mult_aln_base = os.path.normpath(self.project_base + "/shared_files/multiple_alignment")
		if os.path.exists(self.mult_aln_base):
			files = os.listdir(self.mult_aln_base)
			self.mult_aln_files["aln_nt"]  = os.path.normpath(self.mult_aln_base + "/complete_multiple_alignment_nt.fasta")
			self.mult_aln_files["aln_aa"]  = os.path.normpath(self.mult_aln_base + "/complete_multiple_alignment_aa.fasta")
			self.mult_aln_files["log"]  = os.path.normpath(self.mult_aln_base + "/fasttree_log.txt")
			self.mult_aln_files["tree"] = os.path.normpath(self.mult_aln_base + "/target_seq_tree.txt")
		else:
			print("Multiple alignment directory not found!")
	
	def collect_final_info(self):
		self.outputs_base = os.path.normpath(self.project_base + "/final_outputs/")
		filter_dir = os.path.normpath(self.outputs_base + "/model")
		db_dir = os.path.normpath(self.outputs_base + "/database")
		
		if os.path.exists(filter_dir):
			filter = os.path.normpath(filter_dir + "/ROCkOut_Filter.txt")
			if os.path.exists(filter):
				self.rocker_filter = filter
			
		if os.path.exists(db_dir):
			dia = os.path.normpath(db_dir + "/positive_proteins_diamond_db.dmnd")
			bla = os.path.normpath(db_dir + "/positive_proteins_blast_database")
			
			#print("DIAMOND", dia)
			
			if os.path.exists(dia):
				self.targets_dia = dia
			if os.path.exists(bla+".phr"): #the blast db is actually 6 files with different extensions
				self.targets_blast = bla
		

	