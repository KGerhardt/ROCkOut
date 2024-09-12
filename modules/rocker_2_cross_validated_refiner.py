import sys
import os
import shutil

import subprocess

try:
	from .rocker_project_manager import project_manager
	from .pplacer.rocker_phylomap_build import phylomap_build 
	from .visualizer import rockout_visualizer
	from .reads_labeller import protein_trawler
except:
	from rocker_project_manager import project_manager
	from pplacer.rocker_phylomap_build import phylomap_build 
	from visualizer import rockout_visualizer
	from reads_labeller import protein_trawler

	
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None

from math import ceil

#Class for cross-validating rockout models on the data and producing a final weighted model
class cross_validate_refiner:
	def __init__(self, directory, skip_pplace = False, cutoff_bias = "balanced", show_all_reads=False, verbose = False, method = "sequence_outgroups"):
		self.base = os.path.normpath(directory)
		self.final_outputs = os.path.normpath(directory+"/final_outputs")
		
		self.labeller = protein_trawler(self.base, splits = 5, train_fraction = 0.6)
		self.manager = project_manager(self.base)
		self.viz = rockout_visualizer()
		
		#Repos for train/test labels
		self.MA_bins = {}
		
		#Find misc. files including multiple alignments
		self.manager.parse_project_directory()
		self.manager.parse_aligns()
		self.train_method = method

		#Multiple alignment file for finding protein -> MA position offsets
		self.ma_aa_file = None
		
		#Phylogenetic placement files
		self.skip_pplace = skip_pplace
		self.pplace_pos = None
		self.pplace_neg = None
		
		#Values that ensure these are replaced every time
		self.min_aln = 10000.0
		self.max_aln = -1.0
		
		self.min_pct = 10000.0
		self.max_pct = -1.0
		
		self.min_bitscore = 10000.0
		self.max_bitscore = -1.0
		
		#Binning information
		self.protein_to_ma_offsets = None
		
		#Determine the bin boundaries
		
		self.multiple_align_size = None		
		self.aln_x_axis = None
		
		self.bitscore_y_axis = None
		self.pct_id_y_axis = None
		
		self.cutoff_bias = cutoff_bias
		
		#self.window_size = 20
		self.window_size = 20
		#self.aln_window_size = 10.0
		self.aln_window_size = 10.0
		#How fine-grain will the y-axis be? 0.5 bitscore windows by default
		#Note: it's entirely possible to calculate the ROC without binning, 
		#but gains are negligible and it's slower.
		#self.bitscore_resolution = 0.5
		self.bitscore_resolution = 1.0
		self.percent_ID_resolution = 0.5
		self.percent_alignment_resolution = 2.5
		
		self.model_F1_scores = {}
		self.combined_bit = None
		self.combined_id = None
		self.combined_aln = None
		self.min_readlen = 10000
		self.max_readlen = -1
		
		#Verbosity
		self.loud = verbose
		self.show_all_reads = show_all_reads
		
	def parse_project_directory(self):
		#Find misc. files including multiple alignments
		self.manager.parse_project_directory()
		self.manager.parse_aligns()
		self.manager.parse_targets()
		self.manager.parse_multiple_alignment()
		self.ma_aa_file = self.manager.mult_aln_files['aln_aa']

	#Find the appropriate files for the pplacer code to build a phylogenetic reference package
	def prep_for_pplace(self):
		genomes_for_pplacer_pos = []
		for protein_ID in self.manager.positive: #Only positive sequences, please
			for genome in self.manager.targets_nt[protein_ID]:
				genomes_for_pplacer_pos.append(genome)
				
			genomes_for_pplacer_pos = list(set(genomes_for_pplacer_pos))
		
		genomes_for_pplacer_neg = []	
		for protein_ID in self.manager.negative: #Only positive sequences, please
			for genome in self.manager.targets_nt[protein_ID]:
				genomes_for_pplacer_neg.append(genome)
					
			genomes_for_pplacer_neg = list(set(genomes_for_pplacer_neg))
			
		if len(genomes_for_pplacer_pos) == 0:
			genomes_for_pplacer_pos = None
			
		if len(genomes_for_pplacer_neg) == 0:
			genomes_for_pplacer_neg = None
			
		self.pplace_pos = genomes_for_pplacer_pos
		self.pplace_neg = genomes_for_pplacer_neg
		
		
		phylomap_build(pos = genomes_for_pplacer_pos,
						neg = genomes_for_pplacer_neg,
						output = self.final_outputs)
		
		return None
	
	#Load the project's multiple alignment of positive sequences to place reads
	def load_multiple_alignment_and_find_offsets(self):
		current_seq = ""
		current_prot = ""
		self.offsets = {}
		fh = open(self.manager.mult_aln_files['aln_aa'])
		for line in fh:
			if line.startswith(">"):
				if len(current_seq) > 0:
					self.offsets[current_prot] = current_seq
					current_seq = ""
					
				current_prot = line.strip().split()[0][1:]
			else:
				current_seq += line.strip()
				
		fh.close()
		#Final iteration
		if len(current_seq) > 0:
			self.offsets[current_prot] = current_seq
			self.multiple_align_size = len(current_seq) #All seqlens are identical with multiple alignment.
		
		for p in self.offsets:
			offset_list = []
			offset = 0
			for character in self.offsets[p]:
				if character == "-":
					offset += 1
				else:
					offset_list.append(offset)
				
			offset_list = np.array(offset_list, dtype = np.int32)
			self.offsets[p] = offset_list

	#Causes the protein_trawler class to find all of the project's read alignments,
	#current positive, negative, and homolog to positive sequences
	#Filter to positive targets + besthit by bitscore
	#Labels reads as positive/negative/homolog/off-target, collects annotations for homologs
	def collect_reads_and_metadata(self):
		self.labeller.prepare_for_cross_validation(splits = 5, train_fraction = 0.4, seed = None, method = self.train_method)
		for rl in self.labeller.datasets:
			minimum_aln = np.min(self.labeller.datasets[rl]["pct_aln"])
			max_aln = np.max(self.labeller.datasets[rl]["pct_aln"])
			max_pct = np.max(self.labeller.datasets[rl]["pct_id"])
			min_pct = np.min(self.labeller.datasets[rl]["pct_id"])
			min_bit = np.min(self.labeller.datasets[rl]["bitscore"])
			max_bit = np.max(self.labeller.datasets[rl]["bitscore"])
			
			if minimum_aln < self.min_aln:
				self.min_aln = minimum_aln
			if max_aln > self.max_aln:
				self.max_aln = max_aln
			if min_pct < self.min_pct:
				self.min_pct = min_pct
			if max_pct > self.max_pct:
				self.max_pct = max_pct
			if min_bit < self.min_bitscore:
				self.min_bitscore = min_bit
			if max_bit > self.max_bitscore:
				self.max_bitscore = max_bit
		
		#We regard > 100% aln to inherently be a miss during filtering
		if self.max_aln > 100.0:
			self.max_aln = 100.0
		
		aln_steps = ceil((100.0-self.min_aln)/self.percent_alignment_resolution)
		id_steps = ceil((self.max_pct-self.min_pct)/self.percent_ID_resolution)
		
		self.bitscore_y_axis = np.arange(self.min_bitscore, self.max_bitscore, self.bitscore_resolution)
		#100.0 is a special number for both of these axes, 
		#so we switch to linspace to ensure that the max value is observed in the breaks
		self.pct_id_y_axis = np.round(np.linspace(self.min_pct, self.max_pct, num = id_steps), 2)
		self.aln_x_axis = np.round(np.linspace(self.min_aln, 100, num = aln_steps), 2)
	
	#Functions for generating 2D histograms of readmapping results by positive/confounder
	def prep_bins_bitscore_vs_pos_in_MA(self):
		per_position_data = {"Target" : {}, "Confounder": {}}
		for vbin in self.bitscore_y_axis:
			#The range is always 0:max_pos
			per_position_data["Target"][vbin] = np.zeros(self.multiple_align_size, dtype = np.int32)
			per_position_data["Confounder"][vbin] = np.zeros(self.multiple_align_size, dtype = np.int32)
		
		return per_position_data
		
	def prep_bins_pct_id_vs_pos_in_MA(self):
		per_position_data = {"Target" : {}, "Confounder": {}}
		for vbin in self.pct_id_y_axis:
			#The range is always 0:max_pos
			per_position_data["Target"][vbin] = np.zeros(self.multiple_align_size, dtype = np.int32)
			per_position_data["Confounder"][vbin] = np.zeros(self.multiple_align_size, dtype = np.int32)
		
		return per_position_data
			
	def prep_bins_pct_id_vs_pct_aln(self):
		per_position_data = {"Target":{}, "Confounder":{}}
		for id in self.pct_id_y_axis:
			per_position_data["Target"][id] = np.zeros(len(self.aln_x_axis), dtype = np.int32)
			per_position_data["Confounder"][id] = np.zeros(len(self.aln_x_axis), dtype = np.int32)
			
		return per_position_data
	
	#Find the nearest value in an array
	def find_nearest(self, array, value):
		#array = np.asarray(array)
		idx = (np.abs(array - value)).argmin()
		#return array[idx].astype(np.int32)
		return np.int32(idx)
	
	#Convert a read's positional information to each of the 2D histograms' bins
	def bin_read(self, bs, id, aln, targ, s, e):
		#Coerce reads to bins
		bs_index = self.find_nearest(self.bitscore_y_axis, bs)
		id_index = self.find_nearest(self.pct_id_y_axis, id)
		aln_index = self.find_nearest(self.aln_x_axis, aln)
		
		falls_into_bitscore = self.bitscore_y_axis[bs_index]
		falls_into_id = self.pct_id_y_axis[id_index]
		falls_into_aln = self.aln_x_axis[aln_index]
		
		#Map read local alignment on target protein to that protein's pos in MA
		lowhi = np.arange(s-1, e-1)
		fills_bins = self.offsets[targ][lowhi] + lowhi
		
		median_mapping_pos = np.median(fills_bins)
		
		#return falls_into_bitscore, bs_index, falls_into_id, id_index, falls_into_aln, aln_index, fills_bins
		#return falls_into_bitscore, falls_into_id, aln_index, fills_bins, median_mapping_pos
		return falls_into_bitscore, falls_into_id, aln_index, fills_bins, median_mapping_pos
	
	#Instead of repeating the read labelling process for each train/test split, just do it once up front.
	def pre_label_reads(self):
		for rl in self.labeller.datasets:
			self.MA_bins[rl] = {}
			bit_bins = []
			id_bins = []
			aln_indices = []
			MA_medians = []
			
			for read_ID, s, e, bs, targ, id, aln in zip(self.labeller.datasets[rl]["read_id"],
												self.labeller.datasets[rl]["alignment_min_pos"], 
												self.labeller.datasets[rl]["alignment_max_pos"], 
												self.labeller.datasets[rl]["bitscore"], 
												self.labeller.datasets[rl]["target"],
												self.labeller.datasets[rl]["pct_id"],
												self.labeller.datasets[rl]["pct_aln"]):
												
				bit_bin, id_bin, aln_index, MA_bins, MA_pos = self.bin_read(bs, id, aln, targ, s, e)
				bit_bins.append(bit_bin)
				id_bins.append(id_bin)
				aln_indices.append(aln_index)
				MA_medians.append(MA_pos)
				
				self.MA_bins[rl][read_ID] = MA_bins
			
			#Add to the original dataframes
			self.labeller.datasets[rl]["bitscore_bin"] = bit_bins
			self.labeller.datasets[rl]["id_bin"] = id_bins
			self.labeller.datasets[rl]["aln_indices"] = aln_indices
			self.labeller.datasets[rl]["MA_median_mapping_pos"] = MA_medians
				
	#Convert bins from dicts of vertical_bin:array into matrices ofvert x horiz
	#Descending order: row 0 is max vertical axis
	def convert_bins_to_dfs(self, vert_bins, per_position_data):
		#Need to join the disparate tabs into 2D arrs for computation.
		collected_target_data = []
		collected_confounder_data = []
		
		desc_bitscores = sorted(vert_bins, reverse=True)
		
		for vbin in desc_bitscores:
			collected_target_data.append(per_position_data["Target"][vbin])
			per_position_data["Target"][vbin] = None
			collected_confounder_data.append(per_position_data["Confounder"][vbin])
			per_position_data["Confounder"][vbin] = None
		
		collected_target_data = np.vstack(collected_target_data)
		collected_confounder_data = np.vstack(collected_confounder_data)
		
		#Column-wise cumulative sum of bases at each pos. in the protein, split by tgt vs confounder
		collected_target_data = np.cumsum(collected_target_data, axis = 0)
		collected_confounder_data = np.cumsum(collected_confounder_data, axis = 0)
		
		return collected_target_data, collected_confounder_data, desc_bitscores
	
	#No longer used. Relic code.
	#Non-Youden based method of calculating pct ID vs. pct aln cutoffs
	def calculate_ID_aln_cutoffs(self, id_aln_tgt, id_aln_con, desc_ids):
		performance = (id_aln_tgt/id_aln_con)
		
		max_ratios = np.amax(performance, axis = 0)
		
		#nan means zero target and zero con, so we default to a super conservative cutoff of all reads fail.
		nothing_allowed = np.where(np.isnan(max_ratios))
		#max ratio < 1.0 indicates that there is no place where more targets are observed than confounders. Best bet is all reads fail.
		bad_ratios = np.where(max_ratios < 1.0)
		
		#INFINITE HANDLING START
		#Infinite means only targets were observed with zero cons. 
		#In these cols, we need to determine the cut where all of the targets have been seen.
		#Select those columns with an infinite ratio max
		infs = np.where(np.isinf(max_ratios))
		#Subset data to those columns 
		inf_tgt = id_aln_tgt[:, infs]
		inf_ratios = performance[infs]
		
		
		counts = np.sum(inf_ratios == np.inf, axis = 1)
		
		#print(counts)
		
		counts = counts/2
		counts = np.median(counts)
		best_guess = int(counts)
		best_guess = desc_ids[best_guess]
		
		#Find maximum observed target read count by column
		#Find mininum still infinite ratio by column
		id_cuts = np.argmax(inf_tgt, axis = 0)
		id_cuts = id_cuts[0]
		
		infs = infs[0].tolist()
		
		#INFINITE HANDLING END
		
		#Find the ID cutoffs that best separate pos from not-pos except...
		final_id_cutoffs = np.argmax(performance, axis = 0).tolist()
		final_id_cutoffs = [desc_ids[c] for c in final_id_cutoffs]
		final_id_cutoffs = np.array(final_id_cutoffs,dtype = np.float64)
		
		#If no data seen, we are as conservative as possible and assume all reads fail here
		#Usually happens when certain %aln values don't appear in shorter/longer simulated reads
		final_id_cutoffs[nothing_allowed] = 101.0
		#If the best thing to do for classification is to assume all reads fail, then they do.
		final_id_cutoffs[bad_ratios] = 101.0
		
		#Add in the inf values
		#for column, value in zip(infs, tf):
		for column in infs:
			final_id_cutoffs[column] = best_guess
			
		return final_id_cutoffs
	
	#Youden-based method of caculating pct ID vs. pct aln cutoffs
	def calculate_youden_idaln(self, tgt, con, desc_ids):
		#Okay, data's collected, so we select windows and calc ROC from those.
		half_aln_window = int(self.aln_window_size/self.percent_alignment_resolution)
		#We treat reads with >100% aln differently.
		#acceptable_aln = self.aln_x_axis[self.aln_x_axis <= 100.0]
		acceptable_aln = self.aln_x_axis
		max_aln = len(acceptable_aln)
		final_cutoffs = []
		
		for window_midpoint in range(0, len(acceptable_aln)):
			#Get sliding window start, end indices; truncate at edges.
			window_start = window_midpoint-half_aln_window
			if window_start < 0:
				window_start = 0
			window_end = window_midpoint + half_aln_window
			if window_end > max_aln:
				window_end = max_aln
												
			current_window_tgt = tgt[:, np.arange(window_start, window_end)]
			current_window_con = con[:, np.arange(window_start, window_end)]
			
			cutoff_vector = self.calculate_youden(current_window_tgt, current_window_con, self.aln_window_size)
			if self.cutoff_bias == "balanced":
				cutoff_index = cutoff_vector[0]
				
			#75th percentile
			if self.cutoff_bias == "favor_false_negatives":
				cutoff_index = cutoff_vector[1]
				
			#max
			if self.cutoff_bias == "strongly_favor_false_negatives":
				cutoff_index = cutoff_vector[2]
			
			#25th percentile
			if self.cutoff_bias == "favor_false_positives":
				cutoff_index = cutoff_vector[3]
			
			#minimum
			if self.cutoff_bias == "strongly_favor_false_positives":
				cutoff_index = cutoff_vector[4]
			
			#The ID values are expected to drop from 100, which often results in a spike from zero on the left. 
			#This fixes that - it's aesthetic, mostly
			cutoff_value = desc_ids[cutoff_index]
			if cutoff_value < 25.0:
				cutoff_value = 101.0
			
			final_cutoffs.append(cutoff_value)

			
		#return final_cutoffs, cutoff_indices
		return final_cutoffs
	
	#Function for calculating the Youden index of a pos/confounder pair of readmapping reuslts
	def calculate_youden(self, current_window_tgt, current_window_con, expected_window_size):	
		#Select the maximum depth of coverage for the current window at each bitscore. Should always be increasing down the matrix.
		max_by_bitscore_tgt = np.amax(current_window_tgt, axis = 1)
		max_by_bitscore_con = np.amax(current_window_con, axis = 1)
		
		#max_tgt_bitscore = np.min(np.where(np.sum(current_window_tgt, axis = 1) > 0)[0])
		#print(max_tgt_bitscore)
		
		#Cumulative sums means that the max is always found at the final position in the array.
		tgt_max = max_by_bitscore_tgt[max_by_bitscore_tgt.shape[0]-1]
		con_max = max_by_bitscore_con[max_by_bitscore_con.shape[0]-1]
		
		#We care about maximizing the Youden index, per the original ROCker paper
		#Youden = sensitivity + specificity - 1
		
		#Within the context of cumulative sums, at a bitscore in descending order:
		#TP = max_by_bitscore_tgt
		#FP = max_by_bitscore_con
		#FN = tgt_max - max_by_bitscore_tgt
		#TN = con_max - max_by_bitscore_con
		
		#these are the two from the above that we don't already have.
		#fn = np.subtract(tgt_max, max_by_bitscore_tgt)	
		tn = np.subtract(con_max, max_by_bitscore_con)
		
		total_obs = tgt_max + con_max
		accuracy = np.divide(np.add(max_by_bitscore_tgt, tn), total_obs)
		
		#sensitivity = TP / (TP + FN)
		#This simplifies.
		#TP + FN = max_by_bitscore_tgt + tgt_max - max_by_bitscore_tgt = tgt_max
		sensitivity = np.divide(max_by_bitscore_tgt, tgt_max)
		
		#Specificity = TN / (FP + TN)
		#This simplifies
		#FP + TN = max_by_bitscore_con + con_max - max_by_bitscore_con = con_max
		specificity = np.divide(tn, con_max)
		
		Youden_by_bs = sensitivity + specificity - 1
		
		#Argmax returns the first occurrence is there is a tie
		#this equivalent to the highest bitscore among ties, usually seen where Youden = 1
		cutoff = np.argmax(Youden_by_bs)
		
		if np.isnan(Youden_by_bs[cutoff]):
			#Natural return is 0 when there's no conf, but we actually want to choose
			#the min bitscore in this case, which is the last index.
			#cutoff = max_by_bitscore_tgt.shape[0]-1
			cutoff = max_by_bitscore_tgt.shape[0]-1
			cut_vec = [cutoff] * 5
			#print("beep", cutoff)
			#cutoff = 0
		else:
			#Find the median value of tied matches.
			max_youden = Youden_by_bs[cutoff]
			ties = np.where(Youden_by_bs == max_youden)[0]
			
			balanced = int(np.median(ties))
			fn = int(np.quantile(ties, 0.25))
			fn_plus = int(np.min(ties))
			fp = int(np.quantile(ties, 0.75))
			fp_plus = int(np.max(ties))
			
			#Check to see if a cutoff is too high/ excludes positive reads without gain, 
			#drop cutoffs to max level appropriate to save these.
			#Should be descending
			cut_vec = [balanced, fn, fn_plus, fp, fp_plus]
			
			#Handle edge effects
			#if True:
			if False:
				if current_window_con.shape[1] < expected_window_size:
					#print("window", current_window_con.shape[1])
					consums = (np.sum(current_window_con, axis = 1) == 0)
					tgtsums = np.sum(current_window_tgt, axis = 1)
					tgtsums = (tgtsums < np.max(tgtsums))
					lowest_no_confounder = np.argmin(consums)
					
					lowest_max_targets = np.argmin(tgtsums)
					#print(lowest_no_confounder, lowest_max_targets)
					replacement = min([lowest_no_confounder, lowest_max_targets])
					
					#print(lowest_no_confounder, lowest)
					cut_vec = [min([v, replacement]) for v in cut_vec]
					

		return cut_vec
	
	#Use Youden index to find the best bitscore, pct ID cutoffs over the MA of reads
	def calculate_pos_in_MA_cutoffs(self, tgt, con, desc_array):
		#Ceiling of the number of full windows.
		half_window = int(self.window_size/2)
		#Okay, data's collected, so we select windows and calc ROC from those.
		final_cutoffs = []
		cutoff_indices = []
		for window_midpoint in range(0, self.multiple_align_size):
			#Get sliding window start, end indices; truncate at edges.
			window_start = window_midpoint-half_window
			if window_start < 0:
				window_start = 0
			window_end = window_midpoint + half_window
			if window_end > self.multiple_align_size:
				window_end = self.multiple_align_size
			
			#Bitscore
			#Select columns matching the window from cum sums
			current_window_tgt = tgt[:, np.arange(window_start, window_end)]
			current_window_con = con[:, np.arange(window_start, window_end)]
							
			#In-group cutoffs
			#cutoff, accuracy, sensitivity, specificity = self.calculate_youden(current_window_tgt, current_window_con)
			#cutoff_index = self.calculate_youden(current_window_tgt, current_window_con)
			cutoff_vector = self.calculate_youden(current_window_tgt, current_window_con, self.window_size)
			if self.cutoff_bias == "balanced":
				cutoff_index = cutoff_vector[0]
				
			#75th percentile
			if self.cutoff_bias == "favor_false_negatives":
				cutoff_index = cutoff_vector[1]
				
			#max
			if self.cutoff_bias == "strongly_favor_false_negatives":
				cutoff_index = cutoff_vector[2]
			
			#25th percentile
			if self.cutoff_bias == "favor_false_positives":
				cutoff_index = cutoff_vector[3]
			
			#minimum
			if self.cutoff_bias == "strongly_favor_false_positives":
				cutoff_index = cutoff_vector[4]
				
			cutoff_value = desc_array[cutoff_index]
			final_cutoffs.append(cutoff_value)
			cutoff_indices.append(cutoff_index)
			
		#return final_cutoffs, cutoff_indices
		#return final_cutoffs
		return final_cutoffs
	
	#Applies to one training dataset and calculates 3 ROCkOut cutoff sets
	def calculate_curves_for_one_set(self, dataframe, readlength):
		#These reference metadata values gathered earilier
		per_position_data_bitscore = self.prep_bins_bitscore_vs_pos_in_MA()
		per_position_data_id = self.prep_bins_pct_id_vs_pos_in_MA()
		id_aln = self.prep_bins_pct_id_vs_pct_aln()
		
		aln_mids = []
										
		#Iterate through the reads and fill the appropriate row/columns
		for read_id, classifier, bs, targ, id, aln, bit_bin, id_bin, aln_index, MA_pos in zip(dataframe["read_id"],
											dataframe["classifier"], 
											#dataframe["alignment_min_pos"], 
											#dataframe["alignment_max_pos"], 
											dataframe["bitscore"], 
											dataframe["target"],
											dataframe["pct_id"],
											dataframe["pct_aln"],
											dataframe["bitscore_bin"],
											dataframe["id_bin"],
											dataframe["aln_indices"],
											dataframe["MA_median_mapping_pos"]):
			
			MA_bins = self.MA_bins[readlength][read_id]
			#bit_bin, id_bin, aln_index, MA_bins, MA_pos = self.bin_read(bs, id, aln, targ, s, e)
			#aln_mids.append(MA_pos)
			
			med_bin = int(np.median(MA_bins))
			
			if classifier == "Positive":
				#per_position_data_bitscore["Target"][bit_bin][MA_bins] += 1
				per_position_data_bitscore["Target"][bit_bin][med_bin] += 1
				per_position_data_id["Target"][id_bin][MA_bins] += 1
				#per_position_data_id["Target"][id_bin][med_bin] += 1
				id_aln["Target"][id_bin][aln_index] += 1
			else:
				#per_position_data_bitscore["Confounder"][bit_bin][MA_bins] += 1
				per_position_data_bitscore["Confounder"][bit_bin][med_bin] += 1
				per_position_data_id["Confounder"][id_bin][MA_bins] += 1
				#per_position_data_id["Confounder"][id_bin][med_bin] += 1
				id_aln["Confounder"][id_bin][aln_index] += 1
				
		id_aln_tgt, id_aln_con, desc_ids = self.convert_bins_to_dfs(self.pct_id_y_axis, id_aln)
		
		#idaln_cutoffs = self.calculate_ID_aln_cutoffs(id_aln_tgt, id_aln_con, desc_ids)
		idaln_cutoffs = self.calculate_youden_idaln(id_aln_tgt, id_aln_con, np.flip(self.pct_id_y_axis))

		collected_target_data_bitscore, collected_confounder_data_bitscore, desc_bitscores = self.convert_bins_to_dfs(self.bitscore_y_axis, per_position_data_bitscore)
		collected_target_data_id, collected_confounder_data_id, desc_ids = self.convert_bins_to_dfs(self.pct_id_y_axis, per_position_data_id)
		
		bitscore_cutoffs = self.calculate_pos_in_MA_cutoffs(collected_target_data_bitscore,
															collected_confounder_data_bitscore,
															desc_bitscores)
		id_position_cutoffs = self.calculate_pos_in_MA_cutoffs(collected_target_data_id,
															collected_confounder_data_id,
															desc_ids)
		

		return idaln_cutoffs, id_position_cutoffs, bitscore_cutoffs
		#return idaln_cutoffs, id_position_cutoffs, bitscore_cutoffs, bs_plot, id_plot, aln_plot
	
	#Used for loading the filter matrices from file outputs;
	#Not really used in this code.
	def import_filter_from_file(self, file):
		per_rl_x = {}
		per_rl_y = {}
		allx = []
		
		with open(file) as fh:
			header = fh.readline()
			for line in fh:
				segs = line.strip().split("\t")
				rl = int(segs[0])
				if rl not in per_rl_x:
					per_rl_x[rl] = []
					per_rl_y[rl] = []
				
				xaxis = float(segs[1])
				yaxis = float(segs[2])
				per_rl_x[rl].append(xaxis)
				per_rl_y[rl].append(yaxis)
				
				allx.append(xaxis)
				
		allx = list(set(allx))
		allx.sort()
				
		rls = list(per_rl_x.keys())
		min_rl = min(rls)
		max_rl = max(rls)
		self.min_readlen = min_rl
		self.max_readlen = max_rl
		z_shape = max_rl - min_rl + 1
		x_shape = len(allx)
		
		x_index = 0
		y_index = 0
		x_to_index = {}
		for x in allx:
			x_to_index[x] = x_index
			x_index += 1
		
		matrix = np.zeros(shape = (z_shape, x_shape), dtype = np.float64)
		for rl in per_rl_x:
			zi = rl - min_rl
			for x, y in zip(per_rl_x[rl], per_rl_y[rl]):
				xi = x_to_index[x]
				matrix[zi, xi] = y
				
		matrix = self.interpolate_matrix(matrix)
		allx = np.array(allx, dtype = np.float64)
		
		return matrix, allx
	
	#Convert the calculated cutoffs into a filter matrix for filtering test data
	def cutoffs_to_filter(self, dataset, xaxis, yaxis):
		min_rl, max_rl = 10000, -1
		for rl in dataset:
			if rl < min_rl:
				min_rl = rl
			if rl > max_rl:
				max_rl = rl
		
		#Find the size of the final matrix
		z_shape = max_rl - min_rl + 1
		
		per_rl_x = {}
		per_rl_y = {}
		#Transform the data for convenience
		for rl in dataset:
			if rl not in per_rl_x:
				per_rl_x[rl] = xaxis
				per_rl_y[rl] = dataset[rl] 	
		
		x_shape = xaxis.shape[0]
		
		x_index = 0
		y_index = 0
		x_to_index = {}
		for x in xaxis:
			x_to_index[x] = x_index
			x_index += 1
		
		#Create matrix; fill the rows which have actually calculated values
		matrix = np.zeros(shape = (z_shape, x_shape), dtype = np.float64)
		for rl in per_rl_x:
			zi = rl - min_rl
			for x, y in zip(per_rl_x[rl], per_rl_y[rl]):
				xi = x_to_index[x]
				matrix[zi, xi] = y
		
		#Interpolate between adjacent filled rows
		matrix = self.interpolate_matrix(matrix)
		
		return matrix, min_rl, max_rl
	
	#Function that linearly interpolates between cutoffs of adjacent calculated models
	#e.g. read lengths of 100, 150, 200 - interpolate over each position in X for the 49 empty values between
	#100 and 150 and 150 and 200. Essentially draw straight lines parallel to the read length axis joining
	#adjacent read-length models' cutoffs.
	def interpolate_matrix(self, cutoff_matrix):
		rowsums = np.sum(cutoff_matrix, axis = 1)
		calculated_rows = np.sort(np.where(rowsums > 0)[0])

		for i in range(0, len(calculated_rows)-1):
			start = calculated_rows[i]
			end = calculated_rows[i+1]
			
			step_sizes = (cutoff_matrix[end] - cutoff_matrix[start]) / (end-start)
			
			for i in range(1, (end-start)):
				cutoff_matrix[start + i] = cutoff_matrix[start] + (step_sizes * i)
			
		return cutoff_matrix
	
	def fill_confmat(self, true_labels, assignments, presumptive_false_negs, presumptive_true_negs):
		confmat = {'true_positive':0,
				'false_positive':0,
				'false_negative':presumptive_false_negs,
				'true_negative':presumptive_true_negs}
				
		for truth, assign in zip(true_labels, assignments):
			if truth == "Positive":
				if assign:
					confmat['true_positive'] += 1
				else:
					confmat['false_negative'] += 1
			else:
				if assign:
					confmat['false_positive'] += 1
				else:
					confmat['true_negative'] += 1
		
		return confmat
	
	def confmat_to_stats(self, confmat):
		totals = confmat["true_positive"] + confmat["false_positive"] + confmat["true_negative"] +confmat["false_negative"]
		accuracy = (confmat["true_positive"]+confmat["true_negative"]) / totals
		#sensitivity = true_positive / (true_positive + false_negative)
		sensitivity = confmat["true_positive"] / (confmat["true_positive"] + confmat["false_negative"])
		#Specificity = TN / (false_positive + true_negative)
		specificity = confmat["true_negative"] / (confmat["false_positive"] + confmat["true_negative"])
	
		total_size = confmat['true_negative']+confmat['true_positive']+confmat["false_positive"]+confmat["false_negative"]
		fpr = confmat["false_positive"] / (confmat["false_positive"]+confmat["true_negative"])
		fnr = confmat["false_negative"] / (confmat["false_negative"] + confmat["true_positive"])
		f1 = (2*confmat['true_positive']) / (2*confmat['true_positive'] + confmat['false_positive']+confmat['false_negative'])
		
		return accuracy, sensitivity, specificity, fpr, fnr, f1
	
	def filter_testing_data(self, test, bs, id, ia, minimum_rl, maximum_rl, test_idx, return_confmat = False, test_out_path = ""):
		#Presumptive failures for pct aln > 100 pct or (shouldn't be possible) pct ID > 100
		too_much_aln = test.loc[(test['pct_aln'] > 100.0) | (test['pct_id'] > 100.0)]
		
		presumptive_false_negs = len(too_much_aln[too_much_aln['classifier'] == "Positive"])
		presumptive_true_negs = len(too_much_aln[too_much_aln['classifier'] != "Positive"])
		
		#Remaining reads
		test = test.loc[(test['pct_aln'] <= 100.0) & (test['pct_id'] <= 100.0)]
		
		percent_alignment_matches = np.searchsorted(self.aln_x_axis, test["pct_aln"], side = 'left')

		bs_passes = []
		idpos_passes = []
		idaln_passes = []
		
		
		#for real_label, start, end, alignment_target, bitscore, evalue, pct_aln, pct_id in zip(test["classifier"],
		for real_label, start, end, targ, bitscore, pct_aln_index, pct_aln, pct_id, qlen, median_mapping_pos in zip(test["classifier"],
																			test["alignment_min_pos"],
																			test["alignment_max_pos"],
																			test["target"],
																			test["bitscore"],
																			percent_alignment_matches,
																			test["pct_aln"],
																			test["pct_id"],
																			test["query_length"],
																			test["MA_median_mapping_pos"]):
			
			#Query lengths in the table are in AA length; we need nt length
			
			#Can't go outside model bounds
			if qlen <  minimum_rl:
				qlen = minimum_rl
			if qlen > maximum_rl:
				qlen = maximum_rl
				
			#print(pct_aln)
			#print(pct_aln, self.aln_x_axis[pct_aln_index], self.aln_x_axis[pct_aln_index-1], max(self.aln_x_axis))
			
			#Read length to filter matrix row position
			readlen_index = qlen - minimum_rl
			readlen_index = int(readlen_index)
			
			if readlen_index > bs.shape[0]:
				readlen_index = bs.shape[0] - 1
				
			#Collect the location in the middle of the alignment
			#read_mid = int(((start+end-1)/2))
			#Place it into multiple alignment of targets
			#median_mapping_pos = self.offsets[targ][read_mid] + read_mid
			
			median_mapping_pos = int(median_mapping_pos)
				
			bitscore_cutoff = bs[readlen_index, median_mapping_pos]
			id_pos_cutoff = id[readlen_index, median_mapping_pos]
			id_aln_cutoff = ia[readlen_index, pct_aln_index]
			
			passes_bitscore = bitscore >= bitscore_cutoff
			passes_idpos = pct_id >= id_pos_cutoff
			passes_id_aln = pct_id >= id_aln_cutoff
			
			bs_passes.append(passes_bitscore)
			idpos_passes.append(passes_idpos)
			idaln_passes.append(passes_id_aln)
			
		bs_passes = np.array(bs_passes, dtype = np.bool_)
		idpos_passes = np.array(idpos_passes, dtype = np.bool_)
		idaln_passes = np.array(idaln_passes, dtype = np.bool_)
		vote = (np.sum([bs_passes, idpos_passes, idaln_passes], axis = 0) > 1)
		label_is_positive = np.array(test["classifier"] == "Positive", dtype = np.bool_)
		
		false_positives = np.where(np.logical_and(np.logical_not(label_is_positive), vote))[0]
		false_negatives = np.where(np.logical_and(label_is_positive, np.logical_not(vote)))[0]
		
		#print(label_is_positive)
		
		test_writeout = test
		test_writeout['passes_bitscore_pos_filter'] = bs_passes
		test_writeout['passes_pct_id_pos_filter'] = idpos_passes
		test_writeout['passes_pct_id_pct_aln_filter'] = idaln_passes
		test_writeout['passes_ensemble_filter'] = vote
		
		test_writeout.to_csv(os.path.normpath(test_out_path+ "/testing_data_"+str(test_idx)+".txt"), sep = "\t", index=False)
		
		bs_perf = self.fill_confmat(test["classifier"], bs_passes, presumptive_false_negs, presumptive_true_negs)
		idpos_perf = self.fill_confmat(test["classifier"], idpos_passes, presumptive_false_negs, presumptive_true_negs)
		idaln_perf = self.fill_confmat(test["classifier"], idaln_passes, presumptive_false_negs, presumptive_true_negs)
		vote_perf = self.fill_confmat(test["classifier"], vote, presumptive_false_negs, presumptive_true_negs)
		
		#print(test_idx, readlen, votes)
		#print(test_idx, votes)
		res = self.confmat_to_stats(vote_perf)	
		
		
		self.model_F1_scores[test_idx] = res[5]
		
		#Final reporting
		if return_confmat:
			return bs_perf, idpos_perf, idaln_perf, vote_perf
		else:
			if self.loud:
				print("Presumptively filtered for pct. alignment > 100")
				print("False negs:", presumptive_false_negs, "True negs:", presumptive_true_negs)
				print("")
				print("bitscore only", bs_perf)
				print("idpos only", idpos_perf)
				print("idaln only", idaln_perf)
				print("ensemble", vote_perf)
				
				print("FPR", res[3]*100)
				print("FNR", res[4] * 100)
				print("F1 score (model weight)", res[5])
				print("")
	
	def craft_cutoff_plots(self, dataframe, bitscore_cutoffs, id_position_cutoffs, idaln_cutoffs):
		#Remove too highs.
		dataframe = dataframe.loc[(dataframe['pct_aln'] <= 100.0) & (dataframe['pct_id'] <= 100.0)]
		if not self.show_all_reads:
			if len(dataframe.index) > 12000:
				dataframe = self.viz.sample_to_reasonable(dataframe)
	
		bs_plot = self.viz.visualize_reads_and_roc_curves(df = dataframe,
														x = 'MA_median_mapping_pos', 
														y = 'bitscore', 
														x_curve = np.arange(0, self.multiple_align_size),
														y_curve = bitscore_cutoffs,
														xname = "Position in Multiple Alignment", 
														yname = "Bitscore")
		id_plot = self.viz.visualize_reads_and_roc_curves(df = dataframe, 
														x = 'MA_median_mapping_pos', 
														y = 'pct_id', 
														x_curve = np.arange(0, self.multiple_align_size), 
														y_curve = id_position_cutoffs,
														xname = "Position in Multiple Alignment", 
														yname = "Percent Identity to Reference")
		aln_plot = self.viz.visualize_reads_and_roc_curves(df = dataframe,
														x = 'pct_aln', 
														y = 'pct_id', 
														x_curve = self.aln_x_axis, 
														y_curve = idaln_cutoffs,
														xname = "Percent Alignment", 
														yname = "Percent Identity to Reference")
														
		return bs_plot, id_plot, aln_plot
	
	def output_cv_models(self, idaln, idpos, bitpos, model_base):
		pass
	
	def process_CV(self):
		#Output tracking
		if not os.path.exists(os.path.normpath(self.base + "/final_outputs/train_test")):
			os.mkdir(os.path.normpath(self.base + "/final_outputs/train_test"))
		if not os.path.exists(os.path.normpath(self.base + "/final_outputs/train_test/training_data")):
			os.mkdir(os.path.normpath(self.base + "/final_outputs/train_test/training_data"))
		if not os.path.exists(os.path.normpath(self.base + "/final_outputs/train_test/testing_data_and_results")):
			os.mkdir(os.path.normpath(self.base + "/final_outputs/train_test/testing_data_and_results"))
		if not os.path.exists(os.path.normpath(self.base + "/final_outputs/train_test/models")):
			os.mkdir(os.path.normpath(self.base + "/final_outputs/train_test/models"))
	
		training_bitscore_cutoffs = {}
		training_id_pos_cutoffs = {}
		training_id_aln_cutoffs = {}
		
		train_out_base = os.path.normpath(os.path.normpath(self.base + "/final_outputs/train_test/training_data"))
		test_out_base = os.path.normpath(os.path.normpath(self.base + "/final_outputs/train_test/testing_data_and_results"))
		model_base = os.path.normpath(os.path.normpath(self.base + "/final_outputs/train_test/models"))
		
		#Iterates over all readlengths for one CV id out of 5
		for training_dataset, read_length, train_index in self.labeller.get_training_data():
		
			training_dataset.to_csv(os.path.normpath(train_out_base+ "/training_data_"+str(train_index)+".txt"),
									sep = "\t")
		
			if train_index not in training_bitscore_cutoffs:
				training_bitscore_cutoffs[train_index] = {}
				training_id_pos_cutoffs[train_index] = {}
				training_id_aln_cutoffs[train_index] = {}
				
			idaln, idpos, bitpos = self.calculate_curves_for_one_set(training_dataset, read_length)
				
			self.output_cv_models(idaln, idpos, bitpos, model_base)
				
			#We write out the calculated filters because they're useful project-end data
			training_bitscore_cutoffs[train_index][read_length] = bitpos
			training_id_pos_cutoffs[train_index][read_length] = idpos
			training_id_aln_cutoffs[train_index][read_length] = idaln
		
		ma_axis = np.arange(0, self.multiple_align_size)
		bs_filters = {}
		id_filters = {}
		idaln_filters = {}
		readlen_range = {}
		for train_index in training_bitscore_cutoffs:
			bitscore_filter_matrix, min_rl, max_rl = self.cutoffs_to_filter(dataset = training_bitscore_cutoffs[train_index],
																			xaxis = ma_axis,
																			yaxis = self.pct_id_y_axis)
			id_pos_filter_matrix, min_rl, max_rl = self.cutoffs_to_filter(training_id_pos_cutoffs[train_index],
																			xaxis = ma_axis,
																			yaxis = self.pct_id_y_axis)
			id_aln_filter_matrix, min_rl, max_rl = self.cutoffs_to_filter(training_id_aln_cutoffs[train_index],
																			xaxis = self.aln_x_axis,
																			yaxis = self.pct_id_y_axis)
														
			bs_filters[train_index] = bitscore_filter_matrix
			id_filters[train_index] = id_pos_filter_matrix
			idaln_filters[train_index] = id_aln_filter_matrix
			readlen_range[train_index] =  (min_rl, max_rl)
			
			
						
		#for testing_dataset, read_length, test_index in self.labeller.get_testing_data():
		for testing_dataset, test_index in self.labeller.get_testing_data():
			bs = bs_filters[test_index]
			id = id_filters[test_index]
			ia = idaln_filters[test_index]
			min_readlen = readlen_range[test_index][0]
			max_readlen = readlen_range[test_index][1]
			self.filter_testing_data(testing_dataset, bs, id, ia, min_readlen, max_readlen, test_index, test_out_path = test_out_base)
		
		final_model_bitscore_cutoffs = {}
		final_model_id_pos_cutoffs = {}
		final_model_id_aln_cutoffs = {}
		#We could use any of the model collections for the indices - they're shared.
		total_weight = 0
		for model_index in training_bitscore_cutoffs:
			current_model_weight = self.model_F1_scores[model_index]
			total_weight += current_model_weight
			#Still deciding on how to read id_aln model...
			for read_length in training_bitscore_cutoffs[model_index]:
				if read_length not in final_model_bitscore_cutoffs:
					final_model_bitscore_cutoffs[read_length] = np.multiply(current_model_weight , training_bitscore_cutoffs[model_index][read_length])
					final_model_id_pos_cutoffs[read_length] = np.multiply(current_model_weight , training_id_pos_cutoffs[model_index][read_length])
					#final_model_id_aln_cutoffs[read_length] = np.multiply(current_model_weight , training_id_aln_cutoffs[model_index][read_length])
					final_model_id_aln_cutoffs[read_length] = training_id_aln_cutoffs[model_index][read_length]
				else:
					final_model_bitscore_cutoffs[read_length] += np.multiply(current_model_weight , training_bitscore_cutoffs[model_index][read_length])
					final_model_id_pos_cutoffs[read_length] += np.multiply(current_model_weight , training_id_pos_cutoffs[model_index][read_length])
					#final_model_id_aln_cutoffs[read_length] += np.multiply(current_model_weight , training_id_aln_cutoffs[model_index][read_length])
					final_model_id_aln_cutoffs[read_length] += np.minimum(final_model_id_aln_cutoffs[read_length], training_id_aln_cutoffs[model_index][read_length])
		
		for read_length in final_model_bitscore_cutoffs:
			final_model_bitscore_cutoffs[read_length] = final_model_bitscore_cutoffs[read_length]/total_weight
			final_model_id_pos_cutoffs[read_length] = final_model_id_pos_cutoffs[read_length]/total_weight
			final_model_id_aln_cutoffs[read_length] = final_model_id_aln_cutoffs[read_length]/total_weight
		
		bs, min_rl, max_rl = self.cutoffs_to_filter(dataset = final_model_bitscore_cutoffs,
																		xaxis = ma_axis,
																		yaxis = self.pct_id_y_axis)
		id, min_rl, max_rl = self.cutoffs_to_filter(final_model_id_pos_cutoffs,
																		xaxis = ma_axis,
																		yaxis = self.pct_id_y_axis)
		ia, min_rl, max_rl = self.cutoffs_to_filter(final_model_id_aln_cutoffs,
																		xaxis = self.aln_x_axis,
																		yaxis = self.pct_id_y_axis)
																		
		#Set up
		for rl in self.labeller.datasets:
			readlens = self.labeller.datasets[rl]['alignment_max_pos'] - self.labeller.datasets[rl]['alignment_min_pos']
			min_readlen, max_readlen = min(readlens), max(readlens)
			min_readlen, max_readlen = min_readlen + 1, max_readlen + 1
			if min_readlen < self.min_readlen:
				self.min_readlen = min_readlen
			if max_readlen > self.max_readlen:
				self.max_readlen = max_readlen
			
		#complete_data = []
		overall_performance = {'true_positive':0,
								'false_positive':0,
								'false_negative':0,
								'true_negative':0}
		
		if self.loud:
			print("")
			print("Complete model performance")
			print("")
			
			header = ["read_length", "FPR", "FNR", "F1_score"]
			print(*header, sep = "\t")
			for rl in self.labeller.datasets:
				#complete_data.append(self.labeller.datasets[rl])
				bs_perf, idpos_perf, idaln_perf, vote_perf = self.filter_testing_data(
									self.labeller.datasets[rl], 
									bs,
									id,
									ia, 
									self.min_readlen, 
									self.max_readlen, 
									-1, 
									return_confmat = True,
									test_out_path = test_out_base)
				for classification in vote_perf:
					overall_performance[classification] += vote_perf[classification]
				
				res = self.confmat_to_stats(vote_perf)
				print(rl, "%.4f" % (res[3]*100), "%.4f" % (res[4]*100), "%.4f" % res[5], sep = "\t")
			
			res = self.confmat_to_stats(overall_performance)
			print("overall", "%.4f" % (res[3]*100), "%.4f" % (res[4]*100), "%.4f" % res[5], sep = "\t")
			
		self.combined_bit = final_model_bitscore_cutoffs
		self.combined_id = final_model_id_pos_cutoffs
		self.combined_aln = final_model_id_aln_cutoffs
		
	def check_exists_or_make_dir(self, path):
		path = os.path.normpath(self.base+"/"+path)
		if not os.path.exists(path):
			os.makedirs(path, exist_ok = True)
		
	def create_output_directories(self):
		dirs = ["final_outputs",
				"final_outputs/figures",
				"final_outputs/model",
				"final_outputs/reads",
				"final_outputs/train_test",
				"final_outputs/phylogenetic_placement",
				"final_outputs/database"]
		
		for d in dirs:
			self.check_exists_or_make_dir(d)
	
	def output_MA_and_DB(self):
		positive_targets_file = os.path.normpath(self.base + "/final_outputs/database/positive_proteins_aa.fasta")
		fh = open(positive_targets_file, "w")
		fh.write(self.labeller.targets_for_writeout)
		fh.close()
		
		positive_targets_dmnd = os.path.normpath(self.base + "/final_outputs/database/positive_proteins_diamond_db.dmnd")
		
		try:
			makedb = ["diamond", "makedb", "--db", positive_targets_dmnd,  "--in", positive_targets_file]
			print("Building Diamond database for positive targets. Log information will follow.")
			subprocess.call(makedb)
		except:
			print("Couldn't make DIAMOND database of positive targets!")
			
		complete_MA_to_cp = os.path.normpath(self.base + "/final_outputs/model/complete_multiple_alignment_aa.fasta")
		shutil.copy(self.manager.mult_aln_files['aln_aa'], complete_MA_to_cp) 
		
	def output_final_models(self):
		#self.combined_bit = final_model_bitscore_cutoffs
		#self.combined_id = final_model_id_pos_cutoffs
		#self.combined_aln = final_model_id_aln_cutoffs
		ma_positions = np.arange(0, self.multiple_align_size)
		
		with open(os.path.normpath(self.base+"/final_outputs/model/bitscore_vs_MA_pos.txt"), "w") as fh:
			print("read_length", "position_in_MA", "bitscore", sep = "\t", file = fh)
			for readlen in self.combined_bit:
				for x, y in zip(ma_positions, self.combined_bit[readlen]):
					print(readlen, x, y, sep = "\t", file = fh)
					
		with open(os.path.normpath(self.base+"/final_outputs/model/pct_id_vs_MA_pos.txt"), "w") as fh:
			print("read_length", "position_in_MA", "percent_id", sep = "\t", file = fh)
			for readlen in self.combined_id:
				for x, y in zip(ma_positions, self.combined_id[readlen]):
					print(readlen, x, y, sep = "\t", file = fh)
		
		with open(os.path.normpath(self.base+"/final_outputs/model/pct_id_vs_pct_aln.txt"), "w") as fh:
			print("read_length", "percent_aln", "percent_id", sep = "\t", file = fh)
			for readlen in self.combined_aln:
				for x, y in zip(self.aln_x_axis, self.combined_aln[readlen]):
					print(readlen, x, y, sep = "\t", file = fh)
		
		
	def output_labeled_reads(self):
		for read_length in self.labeller.datasets:
			outpath = os.path.normpath(self.base+"/final_outputs/reads/complete_reads_read_length_"+str(read_length)+".txt")
			self.labeller.datasets[read_length].to_csv(outpath, sep = "\t", header = True, index = False)
	
	#Obsolete
	def output_train_test_indices(self):
		for rl in self.labeller.train_indices:
			for sample in self.labeller.train_indices[rl]:
				filename = os.path.normpath(self.base+"/final_outputs/test_train/train_indices_CV_split_"+str(sample)+"_read_length_"+str(rl)+".txt")
				with open(filename, "w") as fh:
					for index in self.labeller.train_indices[rl][sample]:
						print(index, file = fh)
		for rl in self.labeller.test_indices:
			for sample in self.labeller.test_indices[rl]:
				filename = os.path.normpath(self.base+"/final_outputs/test_train/test_indices_CV_split_"+str(sample)+"_read_length_"+str(rl)+".txt")
				with open(filename, "w") as fh:
					for index in self.labeller.test_indices[rl][sample]:
						print(index, file = fh)
		
	def output_plots(self):
		for readlength in self.labeller.datasets:
			bplot, iplot, aplot = self.craft_cutoff_plots(dataframe = self.labeller.datasets[readlength], 
									bitscore_cutoffs = self.combined_bit[readlength], 
									id_position_cutoffs = self.combined_id[readlength], 
									idaln_cutoffs = self.combined_aln[readlength])
			bplot.write_html(os.path.normpath('.'.join([self.base+"/final_outputs/figures/bitscore_vs_MA_pos", "read_length_" + str(readlength), "html"])))
			iplot.write_html(os.path.normpath('.'.join([self.base+"/final_outputs/figures/pct_ID_vs_MA_pos", "read_length_" + str(readlength), "html"])))
			aplot.write_html(os.path.normpath('.'.join([self.base+"/final_outputs/figures/pct_ID_vs_pct_aln", "read_length_" + str(readlength), "html"])))
		
	def run(self):
		self.parse_project_directory()
		if not self.skip_pplace:
			self.prep_for_pplace()
		else:
			print("Not generating a new pplacer phylogeny this time.")
			
		self.load_multiple_alignment_and_find_offsets()
		
		#Create pplacer phylogeny
		self.collect_reads_and_metadata()
		self.pre_label_reads()
		self.process_CV()
		
		self.create_output_directories()
		
		self.output_MA_and_DB()
		
		self.output_final_models()
		
		self.output_labeled_reads()
		
		#self.output_train_test_indices()
		
		print("Creating plots of the models...")
		self.output_plots()
	
		print("ROCkOut Refine is complete!")
		

def build_rockout_model(parser, opts):
	#multiprocessing.freeze_support()

	dirname = opts.dir
	
	if dirname is None:
		sys.exit("ROCkOut needs a directory name to build a model!")
	
	cutoff_tiebreaker = opts.cutoff
	
	cutoff_dict = {"balanced":"balanced",
					"fn":"favor_false_negatives",
					"fn+":"strongly_favor_false_negatives",
					"fp":"favor_false_positives",
					"fp+":"strongly_favor_false_positives"}
					
	if cutoff_tiebreaker not in cutoff_dict:
		print("Your selected cutoff", cutoff_tiebreaker, "was not recognized.")
		sys.exit("Please choose one of: fn+, fn, balanced, fp, fp+")
	else:
		cutoff_bias = cutoff_dict[cutoff_tiebreaker]
		
	skip_pplacer = opts.skip_pplacer
	
	all_reads = opts.big_viz
	
	filter_meth = opts.cv_meth
	if filter_meth not in ['sequence_outgroups', 'subsample']:
		print("Filter method", filter_meth, "not recognized.")
		print("Defaulting to sequence_outgroups")
		filter_meth = "sequence_outgroups"
	
	quiet = opts.quiet
	verbose = (not quiet)
				
	mn = cross_validate_refiner(directory = dirname, 
								skip_pplace = skip_pplacer, 
								cutoff_bias = cutoff_bias,
								show_all_reads = all_reads,
								verbose = verbose,
								method = filter_meth)
	mn.run()
