import sys
import os

import re
import numpy as np
import pandas as pd

class protein_trawler:
	def __init__(self, project_directory, splits = 5, train_fraction = 0.4):
		self.base = os.path.normpath(project_directory)
		self.positives = None
		self.negatives = None
		
		#Information on each protein, organized by genome.
		self.proteins_by_genome = {}
		self.starts_by_genome = {}
		self.ends_by_genome = {}
		self.pos_neg_prob_by_genome = {}
		self.annot_by_genome = {}
		
		#Dict tracking whether each homolog recovered maps to a currently positive protein sequence
		self.foreign_matches = {}
		
		self.read_files = None
		self.raw_read_files = None
		self.valid_targets = None
		self.valid_origin_proteins = None
		self.targets_for_writeout = []
		
		self.datasets = None
		
		self.splits = splits
		self.train_fraction = train_fraction

		self.seeds = {}
		self.train_indices = None
		self.test_indices = None
		
	####LOAD COORDINATE LABELS	
	
	#Peruse a subdirectory, positive or negative, and return all coordinates, 
	#probable coordinates, and annotations for all proteins	
	def trawl(self, is_positive = True):
		#positive/ID/coords
		#positive/ID/probable_target_coords
		#positive/ID/proteome
		if is_positive:
			posneg = "positive"
			#final_label = "Positive"
		else:
			posneg = "negative"
			#final_label = "Negative"
		
		final_label = "Target"
			
		coordinate_dict = {}
		foreign_matches = {}
		
		trunk = os.path.normpath(self.base + "/"+posneg)
				
		if os.path.exists(trunk):
			uniprot_ids = os.listdir(trunk)
			for id in uniprot_ids:
				coords = os.path.normpath(trunk+"/"+id+"/coords")
				probs = os.path.normpath(trunk+"/"+id+"/probable_target_coords")
				labels = os.path.normpath(trunk+"/"+id+"/proteome")
				target_prots = os.path.normpath(trunk+"/"+id+"/target_protein")
				
				coord_files = [os.path.normpath(coords+"/"+f) for f in os.listdir(coords)]
				probs_files = [os.path.normpath(probs+"/"+f) for f in os.listdir(probs)]
				label_files = [os.path.normpath(labels+"/"+f) for f in os.listdir(labels)]
				target_prots = [os.path.normpath(target_prots+"/"+f) for f in os.listdir(target_prots)]
				
				target_prots = [t for t in target_prots if "_AA.fasta" in t]
				
				coord_files.sort()
				probs_files.sort()
				label_files.sort()
				target_prots.sort()
				
				if is_positive:
					for file in target_prots:
						fh = open(file)
						#print(file)
						for line in fh:
							self.targets_for_writeout.append(line)
							if line.startswith(">"):
								next_item = line.strip().split()[0][1:]
								origin_protein = next_item.split("__")[-1]
								self.valid_origin_proteins.append(origin_protein)
								self.valid_targets.append(next_item)
						fh.close()	
								
				for c, p, f in zip(coord_files, probs_files, label_files):
					genome = os.path.basename(c)						
					#genome = genome.split(".coords.txt")[0]
					genome = genome[:-11]
					
					if genome not in coordinate_dict:
						#protein ID, start, stop, label_type, annotation
						coordinate_dict[genome] = {}
						foreign_matches[genome] = {}
					
					proteome_labels = {}
					with open(f) as fh:
						for line in fh:
							if line.startswith(">"):
								annot = line.strip()[1:]
								protein_id = annot.split(";")[0]
								proteome_labels[protein_id] = annot
					with open(c) as fh:
						for line in fh:
							#segs = line.strip().split("\t")
							segs = line.strip().split()
							prot = segs[1]
								
							loc = [int(segs[3]), int(segs[4])]
							start = min(loc)
							end = max(loc)
							annot = proteome_labels[prot]
							
							info = [start, end, final_label, annot]
							if prot not in coordinate_dict[genome]:
								coordinate_dict[genome][prot] = info
							else:
								#Overwrite a negative or homology target for this protein
								if final_label == "Positive":
									coordinate_dict[genome][prot] = info
					
					probables = {}
					with open(p) as fh:
						#has header
						fh.readline()
						#origin_genome	query_parent	[query_id]	[qstart]	[qend]	qstrand	aligns_to_target	pct_aln_to_tgt	pct_ID_to_tgt
						for line in fh:
							segs = line.strip().split("\t")
							prot = segs[2]
							target = segs[6]
							target = target.split("__")
							target = target[-1]
							
							if prot not in foreign_matches[genome]:
								foreign_matches[genome][prot] = []
							if target not in foreign_matches[genome][prot]:
								foreign_matches[genome][prot].append(target)
							
							if prot not in probables:
								loc = [int(segs[3]), int(segs[4])]
								start = min(loc)
								end = max(loc)
								annot = proteome_labels[prot]
								info = [start, end, "Homology_Target", annot]
								probables[prot] = info
					
					#We don't want to overwrite a negative or positive here
					for prot in probables:
						if prot not in coordinate_dict[genome]:
							coordinate_dict[genome][prot] = probables[prot]
		to_remove = []
		for genome in foreign_matches:
			if len(foreign_matches[genome]) == 0:
				to_remove.append(genome)
				
		for genome in to_remove:
			discard = foreign_matches.pop(genome)
						
						
		return coordinate_dict, foreign_matches
	
	#Trawl positive and negative; combine results
	def load_labels(self):
		self.valid_targets = []
		self.valid_origin_proteins = []
		
		positives, foreign_keys_pos = self.trawl(is_positive = True) #must exist
		negatives, foreign_keys_neg = self.trawl(is_positive = False) #may exist
		
		self.valid_targets = set(self.valid_targets)
		self.valid_origin_proteins = set(self.valid_origin_proteins)
		self.targets_for_writeout = ''.join(self.targets_for_writeout)
		
		positive_proteins = []
		negative_proteins = []
		
		#Combine results into one list of prots per genome
		all_genomes = []
		for genome in positives:
			all_genomes.append(genome)
			for prot in positives[genome]:
				if positives[genome][prot][2] != ";Homology_Target":
					positive_proteins.append(prot)
					
		for genome in negatives:
			all_genomes.append(genome)
			for prot in negatives[genome]:
				if negatives[genome][prot][2] != ";Homology_Target":
					negative_proteins.append(prot)
		
		positive_proteins = set(positive_proteins)
		negative_proteins = set(negative_proteins)
		
		for genome in foreign_keys_pos:
			for protein in foreign_keys_pos[genome]:
				list_of_alignment_targets = set(foreign_keys_pos[genome][protein])
				
				hits_a_positive_target = list_of_alignment_targets.intersection(positive_proteins)
				if len(hits_a_positive_target) > 0:
					hits_a_positive_target = True
				else:
					hits_a_positive_target = False
				
				foreign_keys_pos[genome][protein] = hits_a_positive_target
				
		for genome in foreign_keys_neg:
			for protein in foreign_keys_neg[genome]:
				list_of_alignment_targets = set(foreign_keys_neg[genome][protein])
				
				hits_a_positive_target = list_of_alignment_targets.intersection(positive_proteins)
				if len(hits_a_positive_target) > 0:
					hits_a_positive_target = True
				else:
					hits_a_positive_target = False
					
				foreign_keys_neg[genome][protein] = hits_a_positive_target

		self.foreign_matches = {}
		#Working as expected
		for genome in foreign_keys_pos:
			for protein in foreign_keys_pos[genome]:
				self.foreign_matches[protein] = foreign_keys_pos[genome][protein]
				
		for genome in foreign_keys_neg:
			for protein in foreign_keys_neg[genome]:
				self.foreign_matches[protein] = foreign_keys_neg[genome][protein]
		
		all_genomes = list(set(all_genomes)) #all unique genome observations
		
		self.proteins_by_genome = {}
		self.starts_by_genome = {}
		self.ends_by_genome = {}
		self.pos_neg_prob_by_genome = {}
		self.annot_by_genome = {}
		
		for genome in all_genomes:
			self.proteins_by_genome[genome] = []
			self.starts_by_genome[genome] = []
			self.ends_by_genome[genome] = []
			self.pos_neg_prob_by_genome[genome] = []
			self.annot_by_genome[genome] = []
			
			if genome in positives:
				for prot in positives[genome]:
					self.proteins_by_genome[genome].append(prot)
					self.starts_by_genome[genome].append(positives[genome][prot][0])
					self.ends_by_genome[genome].append(positives[genome][prot][1])
					self.pos_neg_prob_by_genome[genome].append(positives[genome][prot][2])
					self.annot_by_genome[genome].append(positives[genome][prot][3])
			if genome in negatives:
				for prot in negatives[genome]:
					if prot not in self.proteins_by_genome[genome]: #Don't overwrite data from a positive
						self.proteins_by_genome[genome].append(prot)
						self.starts_by_genome[genome].append(negatives[genome][prot][0])
						self.ends_by_genome[genome].append(negatives[genome][prot][1])
						self.pos_neg_prob_by_genome[genome].append(negatives[genome][prot][2])
						self.annot_by_genome[genome].append(negatives[genome][prot][3])
		
		self.positives = positive_proteins
		self.negatives = negative_proteins
		
		return None
	

	####READ SIMULATION INFO
	
	#Calculate the number of bases by which a read overlaps its origin gene
	def calculate_read_overlap(self, read_start, read_end, gene_start, gene_end):
		read_positions = np.arange(read_start, read_end, dtype = np.int32)
		gene_positions = np.arange(gene_start, gene_end, dtype = np.int32)
		
		#I'm not entirely sure this math is right...
		overlap_size = len(np.intersect1d(read_positions, gene_positions)) + 1
		
		return overlap_size
			
	#Locate the origin protein of a read, if any, and label it as such
	def tag_read(self, genome, read_start, read_end, read_ID, comp, genome_id):
		useful_target = False
		for gene_start, gene_end, protein_name, prot_label, annot in zip(self.starts_by_genome[genome],
														self.ends_by_genome[genome],
														self.proteins_by_genome[genome],
														self.pos_neg_prob_by_genome[genome],
														self.annot_by_genome[genome]):
			
			origin_protein = ""
			annotation = ""
			final_label = "Non_Target"
														
			if read_end >= gene_start and read_start <= gene_end:
				useful_target = True
				overlap_bp = self.calculate_read_overlap(read_start, read_end, gene_start, gene_end)
				qlen = read_end-read_start + 1
				pct_overlap = 100*overlap_bp/qlen
				pct_overlap = round(pct_overlap, 2)
				
				origin_protein = protein_name
				annotation = annot
				
				tagged_name = ';'.join([read_ID, str(read_start), str(read_end), str(overlap_bp), str(pct_overlap), comp, genome_id, "origin_protein="+protein_name, prot_label])	

				break #If we find it, we label it and finish	
				
		if not useful_target:
			tagged_name = ';'.join([read_ID, str(read_start), str(read_end), "0", "0.0", comp, genome_id, "origin_protein=NA", "Non_Target"])	
			
		printable_name = ">" + tagged_name
			
		return printable_name


	#LOCATE ALL UNIQUE ALIGNMENT FILES
	def locate_unique_read_sets(self):
		#Find all the alignment files
		reads = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(self.base)) for f in fn if "aligned_reads" in f and "_read_len_" in f]
		
		self.read_files = {}
		for r in reads:
			#Looks like this:
			#BX842653_read_len_100_aligned_reads.blast.txt
			base_file = os.path.basename(r)
			base_file = base_file.split("_read_len_")
			genome = base_file[0]
			read_length = base_file[1].split("_aligned_reads.blast.txt")[0]
			read_length = int(read_length) #Transform from string
			
			if read_length not in self.read_files:
				self.read_files[read_length] = {}
			#Two copies of a genome simulated at the same read length by BBTools will have
			#Exactly the same simulated reads, meaning any addtl. copies
			#of any genome beyond the first are always 100% redundant and don't need looked at.
			if genome not in self.read_files[read_length]: 
				self.read_files[read_length][genome] = r
				
		return None
		
	#LOCATE ALL UNIQUE TAGGED READ FILES
	def locate_unique_tagged_reads(self):
		#Find all the alignment files
		reads = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser(self.base)) for f in fn if "_tagged.fasta" in f and "_read_len_" in f]
		
		self.raw_read_files = {}
		for r in reads:
			#Looks like this:
			#BX842653_read_len_100_tagged.fasta
			base_file = os.path.basename(r)
			base_file = base_file.split("_read_len_")
			genome = base_file[0]
			read_length = base_file[1].split("_tagged.fasta")[0]
			read_length = int(read_length) #Transform from string
			
			if read_length not in self.raw_read_files:
				self.raw_read_files[read_length] = {}
			#Two copies of a genome simulated at the same read length by BBTools will have
			#Exactly the same simulated reads, meaning any addtl. copies
			#of any genome beyond the first are always 100% redundant and don't need looked at.
			if genome not in self.raw_read_files[read_length]: 
				self.raw_read_files[read_length][genome] = r
				
		return None
		
		
	#REFINEMENT	
		
	#SINGLE BEST HIT READ BY BITSCORE
	def besthit_reads(self, df):
		#print(self.valid_targets)
		#Filter to valid targets
		df = df[df['target'].isin(self.valid_targets)]
		df = df.reset_index(drop=True)
		
		#Collect max bitscore by read
		idx = df.groupby(['read_id'])["bitscore"].transform(max) == df["bitscore"]
		df = df[idx]
		df = df.reset_index(drop = True)
			
		#collect a random read among cases where max bitscore is tied
		idx = df.groupby(["read_id"])['aln_len'] #random integer column - it's effectively discarded
		random_sels = []
		for item in idx:
			rows = item[1].index
			if len(rows) > 2:
				rand = np.random.choice(rows, size = 1)[0]
				random_sels.append(rand)
			else:
				random_sels.append(rows[0])

		random_sels = np.array(random_sels, dtype = np.int32)
		random_sels = np.sort(random_sels)
		
		df = df.iloc[random_sels]
		df = df.reset_index(drop=True)
		
		return df
	
	#LOAD AND BESTHIT ONE READ FILE, FILTERING TO ONLY VALID TARGETS
	def load_read_file(self, read_file, origin_genome):
		#An alignment record has these fields
		#read_id = 496;1125999;1126099;0;0.0;+;ENA|CP002663|CP002663.1;origin_protein=NA;Non_Target	
		#target = UNIPROT__E3D1H3_NEIM7__E3D1H3_NEIM7__CP001561.3252	
		#pct ID = 72.2	
		#length = 18	
		#mismatch = 5	
		#gap = 0	
		#qstart = 47	
		#qend = 100
		#sstart = 22	
		#send = 39	
		#evalue = 1.19e-01	
		#bitscore = 20.8	
		#qlen = 101
		#slen = 391	
		#overlap bp = 0	
		#pct overlap = 0.0
		
		df = pd.read_csv(read_file, sep = "\t", 
						usecols = [0, 1, 2, 3, 8, 9, 11, 12, 15],
						names = ["read_id", "target", "pct_id", "aln_len", "alignment_min_pos", "alignment_max_pos", "bitscore", "query_length", "pct_overlap"])
		
		df = self.besthit_reads(df)
		#Convert to amino acid lengths
		#df['query_length'] = df['query_length']/3
		#Calculate percent alignment with conversion of readlen from AA to nt
		df['pct_aln'] = round((3*(df['aln_len']/df['query_length']))*100, 2)
		
		true_labels = []
		homolog_annotations = []
		
		corrected_mins = np.minimum(df["alignment_min_pos"], df["alignment_max_pos"])
		corrected_maxes = np.maximum(df["alignment_min_pos"], df["alignment_max_pos"])
		
		df["alignment_min_pos"] = corrected_mins
		df["alignment_max_pos"] = corrected_maxes
		
		corrected_mins = None
		corrected_maxes = None

		#Find truth value
		for label in df['read_id']:
			#Default values
			truth = "Non_Target"
			annotation = ""
			
			segs = label.split(";")
			group = segs[-1]
			#Is positive or is negative
			if group == "Target":
				protein = segs[-2]
				protein = protein[15:]
				if protein in self.positives:
					truth = "Positive"
				if protein in self.negatives:
					truth = "Negative"
				
			
			#Is a homology target
			if group == "Homology_Target":
				protein = segs[-2]
				protein = protein[15:]
				hits_to = self.foreign_matches[protein]
				#Only homologs to positive sequences actually count
				if hits_to in self.positives:
					truth = "Homology_Target"
					annotation = "Homolog Annotation:" + self.annot_by_genome[origin_genome][protein]
				else:
					truth = "Non_Target"
					
			true_labels.append(truth)
			homolog_annotations.append(annotation)
			
		df['classifier'] = true_labels
		df['annotations'] = homolog_annotations
				
		return df
		
	def load_raws_file(self, raw_file, passing_IDs):
		results = {}
		cur_seq = []
		cur_id = None
		with open(raw_file) as fh:
			for line in fh:
				if line.startswith(">"):
					if cur_id is not None:
						cur_seq = ''.join(cur_seq)
						results[cur_id] = cur_seq
						
					cur_id = line.strip()
					cur_seq = []
					#Removes the initial > from the line.
					if cur_id[1:] not in passing_IDs:
						cur_id = None
				else:
					if cur_id is not None:
						cur_seq.append(line.strip())
		
		#Final iteration
		if cur_id is not None:
			cur_seq = ''.join(cur_seq)
			results[cur_id] = cur_seq			

		return results
		
	#PROCESS AN ENTIRE PROJECT'S WORTH OF READS TO BESTHITS AND PREPARE TRAIN/TEST SPLITS
	def collect_and_label_reads(self):
		self.datasets = {}
		for read_length in self.read_files:
			if read_length not in self.datasets:
				self.datasets[read_length] = []
			for genome in self.read_files[read_length]:
				next_set_of_reads = self.load_read_file(self.read_files[read_length][genome], genome)
				self.datasets[read_length].append(next_set_of_reads)
				
		for read_length in self.datasets:
			self.datasets[read_length] = pd.concat(self.datasets[read_length])
			#RESET INDEX OR IT REPEATS INDEX NUMBERS MANY TIMES
			self.datasets[read_length] = self.datasets[read_length].reset_index(drop = True)
	
	def collect_reads_and_raws(self, read_length, genome):
		next_set_of_alignments = self.load_read_file(self.read_files[read_length][genome], genome)
		
		acceptable_reads = set(next_set_of_alignments['read_id'])
		
		next_set_of_raws = self.load_raws_file(raw_file = self.raw_read_files[read_length][genome], 
												passing_IDs = acceptable_reads)
			
		acceptable_reads = None
		refactored_IDs = {}
		id_to_tag = dict(zip(next_set_of_alignments['read_id'], next_set_of_alignments['classifier']))
		for id in next_set_of_raws:
			with_addendum = id+";"+id_to_tag[id[1:]]
			refactored_IDs[with_addendum] = next_set_of_raws[id]
			next_set_of_raws[id] = None
			
			
		
		#return next_set_of_alignments, next_set_of_raws
		return next_set_of_alignments, refactored_IDs
	
	def collect_train_test_splits(self, splits = 5, train_fraction = 0.4, seed = None):
		self.train_indices = {}
		self.test_indices = {}
		self.seeds = {}
		for rl in self.datasets:
		
			if rl not in self.train_indices:
				self.train_indices[rl] = {}
				self.test_indices[rl] = {}
				self.seeds[rl] = {}
				
			for i in range(0, splits):
				seed = np.random.randint(0, 1024)
				self.seeds[rl][i] = seed
				
				train = self.datasets[rl].groupby('classifier').apply(pd.DataFrame.sample, frac=self.train_fraction, replace=False, random_state = seed).reset_index(level='classifier', drop=True)

				self.train_indices[rl][i] = train.index
				
				mask = np.full(len(self.datasets[rl]), fill_value = True, dtype = bool)
				mask[train.index] = False
								
				self.test_indices[rl][i] = pd.Index(np.where(mask)[0], dtype = np.int64)
		
	def prepare_for_cross_validation(self, splits = 5, train_fraction = 0.4, seed = None):
		print("Collecting metadata")
		self.load_labels()
		print("Loading alignments")
		self.locate_unique_read_sets()		
		self.collect_and_label_reads()
		print("Preparing train/test splits")
		self.collect_train_test_splits(splits = 5, train_fraction = 0.4, seed = None)

		print("Reads loaded!")
		
	def extract_prep(self):
		print("Collecting metadata")
		self.load_labels()
		print("Loading alignments")
		self.locate_unique_read_sets()
		self.locate_unique_tagged_reads()
		
			
	def get_training_data(self):
		for rl in self.datasets:
			for split_index in self.train_indices[rl]:
				train = self.datasets[rl].iloc[self.train_indices[rl][split_index]]
				yield train, rl, split_index
	
	#We actually just want one dataset per rl
	def get_testing_data(self):
		rls_by_test_index = {}
		for rl in self.datasets:
			for split_index in self.test_indices[rl]:
				if split_index not in rls_by_test_index:
					rls_by_test_index[split_index] = []
				rls_by_test_index[split_index].append(rl)
		
		for split_index in rls_by_test_index:
			next_test = []
			for rl in rls_by_test_index[split_index]:
				next_test.append(self.datasets[rl].iloc[self.test_indices[rl][split_index]])
				
			next_test = pd.concat(next_test)
			
		#test = self.datasets[rl].iloc[self.test_indices[rl][split_index]]
			yield next_test, split_index
	
				
#dir = sys.argv[1]
#mn = protein_trawler(dir)
#mn.load_labels()

#Load alignments and classify according to current pos/neg sets.
#mn.load_valid_targets()
#mn.locate_unique_read_sets()
#mn.collect_and_label_reads()

#mn.collect_train_test_splits()

#for train, test in mn.get_tests_and_trains():
#	print(train)
#	print(test)
