import numpy as np
import sys
import pandas as pd
import os

import plotly.graph_objects as go

try:
	from .rocker_project_manager import project_manager
	from .extract_read_labels import output_reads
except:
	from rocker_project_manager import project_manager
	from extract_read_labels import output_reads	


#Warnings in this code are annoying and pointless
pd.options.mode.chained_assignment = None
np.seterr(all="ignore")
	
class cross_validator:
	def __init__(self, dir):
		self.prjdir = dir
		self.manager = None
		self.data = None
		self.basename = os.path.basename(os.path.normpath(dir))
		
		self.outdir = os.path.normpath(self.prjdir + "/final_outputs/cross_validation/")
		if not os.path.exists(self.outdir):
			os.makedirs(self.outdir, exist_ok = True)
		
		self.valid_targets = []		
		self.valid_origin_proteins = []
		
		self.filter_model = {}
		self.observed_readlens = []
		
		self.offsets = {}
		self.window_size = 20
		self.resolution = 0.5
		
		self.subsample_fraction = 0.4
		self.num_subsamps = 5
		self.subsamp_indices = {}
		
		self.models = []		
		self.senspec = {}
		
		self.testing_data_acc_sens_spec = {}
		self.testing_data_asp_e3 = {}
		self.testing_data_asp_e10 = {}
		self.testing_data_asp_e20 = {}
		self.testing_data_asp_e30 = {}
		self.testing_data_asp_staramr = {}
		self.test_index = 1
		
		
	
	def ready_inputs(self):
		self.data, self.manager = output_reads(self.prjdir, external = False)
		dataframes = {}
		for tup in self.data:
			source = os.path.basename(tup[0])
			source = source.split("_read_len_")
			read_len = source[1]
			source = source[0]
			read_len = read_len.split("_aligned_reads.")[0]
			read_len = int(read_len)
			if read_len not in dataframes:
				dataframes[read_len] = []
			
			next_alignments = tup[1]
			dataframes[read_len].extend(next_alignments)
		
		for rl in dataframes:
			dataframes[rl] = pd.DataFrame(dataframes[rl])
			#[read_name, target, label, start, end, bitscore, pct_id, ]
			dataframes[rl].columns = ["read_id","target", "classifier", "start", "end", "bitscore", "evalue", "pct_id", "pct_aln"]
		
		self.data = dataframes
		
	def prepare_offsets(self):
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
			#print(self.offsets[p])
			for character in self.offsets[p]:
				if character == "-":
					offset += 1
				else:
					offset_list.append(offset)
				
			offset_list = np.array(offset_list, dtype = np.int32)
			self.offsets[p] = offset_list

	def find_valid_targets(self):
		self.valid_targets = []
		self.valid_origin_proteins = []
		
		#self.manager.parse_coords()
		for prot in self.manager.targets:
			for file in self.manager.targets[prot]:				
				if '/positive/' in file:
					fh = open(file)
					for line in fh:
						#self.targets_for_writeout += line
						if line.startswith(">"):
							next_item = line.strip().split()[0][1:]
							origin_protein = next_item.split("__")[-1]
							#origin_genome = origin_genome.split(".")[0]
							self.valid_origin_proteins.append(origin_protein)
							self.valid_targets.append(next_item)
					fh.close()
				
		self.valid_targets = set(self.valid_targets)
		#self.invalid_targets = set(self.invalid_targets)
		
		self.valid_origin_proteins = set(self.valid_origin_proteins)
		#self.invalid_origin_genomes = set(self.invalid_origin_genomes)
		
	def find_filter(self):
		if os.path.exists(self.prjdir):
			filter_path = os.path.normpath(self.prjdir + "/final_outputs/model/ROCkOut_Filter.txt")
			#print(os.path.exists(filter_path))
			if os.path.exists(filter_path):
				#Only set this to not None if the file can be found.
				self.filter_file = filter_path
			else:
				print("Project appears incomplete.")
		else:
			print("Project not found.")
			
	def get_filter(self):
		self.filter_model = {}
		self.observed_readlens = []
		fh = open(self.filter_file)
		#skip header
		fh.readline()
		for line in fh:
			segs = line.strip().split()
			read_len, midpt, bitscore = int(segs[0]), float(segs[1]), float(segs[2])
			if midpt not in self.filter_model:
				self.filter_model[midpt] = [[], []]
			self.filter_model[midpt][0].append(read_len)
			self.filter_model[midpt][1].append(bitscore)
			self.observed_readlens.append(read_len)
				
		fh.close()
		
		for midpt in self.filter_model:
			sort_order = np.argsort(self.filter_model[midpt][0])
			self.filter_model[midpt][0] = np.array(self.filter_model[midpt][0], dtype = np.int32)[sort_order]
			self.filter_model[midpt][1] = np.array(self.filter_model[midpt][1], dtype = np.float_)[sort_order]
			
		self.observed_readlens = np.unique(np.array(self.observed_readlens, dtype = np.int32))
		
	def interpolate(self):
		midpoints = np.unique(np.array(list(self.filter_model.keys()), dtype = np.int32))
		min_rl = min(self.observed_readlens)
		max_rl = max(self.observed_readlens)
		max_midpt = max(midpoints)
		
		self.filter_matrix = np.zeros(shape = (max_rl-min_rl, max_midpt+1), dtype = np.float_)
		
		for m in midpoints:
			readlens = self.filter_model[m][0]
			this_min_rl = min(readlens)
			bitscores = self.filter_model[m][1]
			dist = np.max(readlens)-np.min(readlens)
			next_addition = np.zeros(dist, dtype = np.float_)
			
			for i in range(0, len(readlens)-1):
				interp = np.linspace(bitscores[i], bitscores[i+1], readlens[i+1]-readlens[i], endpoint = True, dtype = np.float_)
				self.filter_matrix[(readlens[i]-this_min_rl):(readlens[i+1]-this_min_rl), m] += interp
			
		mat_to_dict = {}
		row = 0
		for i in range(min_rl, max_rl):
			mat_to_dict[i] = self.filter_matrix[row]
			row += 1
			
		self.filter_matrix = mat_to_dict
		
	#Maybe going through the raws and collecting the readlengths first would be good?
	def run_filter_blast_tabular(self, raw, aln, passing, failing, raw_filt):
		mapping_reads = []
		fh = open(aln)
		for line in fh:
			segs = line.strip().split("\t")
			alignment_target = segs[1]
			if alignment_target not in self.offsets:
				continue
			else:
				mapping_reads.append(segs[0])
		fh.close()
		
		mapping_reads = set(mapping_reads)
		
		cur_sl = 0
		cur_read = ""
		is_valid = False
		read_lengths = {}
		fh = open(raw)
		for line in fh:
			if line.startswith(">"):
				segs = line.strip().split()
				read_name = segs[0][1:]
				if len(cur_read) > 0:
					read_lengths[cur_read] = cur_sl
				
				if read_name not in mapping_reads:
					cur_read = ""
					is_valid = False
				else:
					cur_read = read_name
					is_valid = True
					
				cur_sl = 0
			else:
				if is_valid:
					line = line.strip()
					cur_sl += len(line)
			
		fh.close()
	
		#Final iter
		if len(cur_read) > 0:
			read_lengths[cur_read] = cur_sl
		
		passing_reads = []
		passing = open(passing, "w")
		failing = open(failing, "w")
		fh = open(aln)
		for line in fh:
			segs = line.strip().split("\t")
			alignment_target = segs[1]
			if alignment_target not in self.offsets:
				continue
			else:
				read_name = segs[0]
				
				readlength = read_lengths[read_name]
				
				if readlength not in self.filter_matrix:
					continue
				
				these_offsets = self.offsets[alignment_target]
				
				alignment_range = [int(segs[8]), int(segs[9])]
				
				#We need to get the offsets here.
				
				start = min(alignment_range)
				end = max(alignment_range)
				
				lowhi = np.arange(start-1, end-1)
				ma_corrected_positions = these_offsets[lowhi] + lowhi
				
				bitscore = float(segs[11])
				
				cutoffs = self.filter_matrix[readlength][ma_corrected_positions]
				
				diffs = np.subtract(bitscore, cutoffs)
				
				decider = np.mean(diffs)
				
				
				
				if decider > 0:
					passing.write(line)
					passing_reads.append(read_name)
				else:
					failing.write(line)
			
		fh.close()
		passing.close()
		failing.close()
		
		#And clean the raws
		passing_reads = set(passing_reads)
		current_read = None
		current_defline = None
		current_seq = ""
		
		out = open(raw_filt, "w")
		fh = open(raw)
		for line in fh:
			if line.startswith(">"):
				if current_read is not None:
					if current_read in passing_reads:
						out.write(current_defline)
						out.write(current_seq)
				current_defline = line
				current_read = current_defline.strip().split()[0][1:]
				current_seq = ''
			else:
				current_seq += line
		
		fh.close()
		if current_read is not None:
			if current_read in passing_reads:
				out.write(current_defline)
				out.write(current_seq)
		
		out.close()
	
	def confmat_dict(self, ncol):
		tp = np.zeros(ncol, dtype = np.int32)
		tn = np.zeros(ncol, dtype = np.int32)
		fp = np.zeros(ncol, dtype = np.int32)
		fn = np.zeros(ncol, dtype = np.int32)
		confmat = {"tp":tp, "tn":tn, "fp":fp, "fn":fn}
		return confmat
		
	def confmat_to_stats(self, confmat):
		totals = confmat["tp"] + confmat["fp"] + confmat["tn"] +confmat["fn"]
		accuracy = (confmat["tp"]+confmat["tn"]) / totals
		#sensitivity = TP / (TP + FN)
		sensitivity = confmat["tp"] / (confmat["tp"] + confmat["fn"])
		#Specificity = TN / (FP + TN)
		specificity = confmat["tn"] / (confmat["fp"] + confmat["tn"])
		
		return accuracy, sensitivity, specificity
	
	def prepare_subsamples(self):
	
		for rl in self.data:
			self.testing_data_acc_sens_spec[rl] = {}
			self.testing_data_asp_e3[rl] = {}
			self.testing_data_asp_e10[rl] = {}
			self.testing_data_asp_e20[rl] = {}
			self.testing_data_asp_e30[rl] = {}
			self.testing_data_asp_staramr[rl] = {}
		
			self.test_index = 1
			for i in range(0, self.num_subsamps):
				print("Read length", rl, "subsample", i+1, "of", self.num_subsamps)
				#Proportionally sample the original data based on the classifier labels: Positive, Negative, Non_Target
				train = self.data[rl].groupby('classifier').apply(pd.DataFrame.sample, frac=self.subsample_fraction).reset_index(level='classifier', drop=True)
				
				#output train data
				train_out = os.path.normpath(self.outdir + "/" + self.basename+"_read_len_"+str(rl)+"_train_group_"+str(self.test_index)+".txt")
				train.to_csv(train_out, sep = "\t", index = False)
				
				cutoffs = self.calculate_roc_curves(train, rl)

				#Select all the indices from the original dataframe NOT in the training list
				mask = np.full(len(self.data[rl]), fill_value = True, dtype = bool)
				mask[train.index] = False
				test = self.data[rl][mask]
				
				
				test_out = os.path.normpath(self.outdir + "/" +self.basename+"_read_len_"+str(rl)+"_test_group_"+str(self.test_index)+".txt")
				#Write test data out
				test.to_csv(test_out, sep = "\t", index = False)
				self.test_index += 1
				
				#Confusion matrix rows
				#tp = np.zeros(cutoffs.shape[0], dtype = np.int32)
				#tn = np.zeros(cutoffs.shape[0], dtype = np.int32)
				#fp = np.zeros(cutoffs.shape[0], dtype = np.int32)
				#fn = np.zeros(cutoffs.shape[0], dtype = np.int32)
				
				#evals
				evals = [1*10**-3, 1*10**-10, 1*10**-20, 1*10**-30]
				
				aln_cut = 50.0
				id_cut = 90.0
				
				rocker_conf = self.confmat_dict(cutoffs.shape[0])
				eval_10e3_conf = self.confmat_dict(cutoffs.shape[0])
				eval_10e10_conf = self.confmat_dict(cutoffs.shape[0])
				eval_10e20_conf = self.confmat_dict(cutoffs.shape[0])
				eval_10e30_conf = self.confmat_dict(cutoffs.shape[0])
				staramr_conf = self.confmat_dict(cutoffs.shape[0])
				
				#Iterate over rows in the testing data and filter
				for real_label, start, end, alignment_target, bitscore, evalue, pct_aln, pct_id in zip(test["classifier"],
								test["start"],
								test["end"],
								test["target"],
								test["bitscore"],
								test["evalue"],
								test["pct_aln"],
								test["pct_id"]):
					
					
					lowhi = np.arange(start-1, end-1)
					these_offsets = self.offsets[alignment_target]
					ma_corrected_positions = these_offsets[lowhi] + lowhi
					
					#These would be filtered on the basis of not aligning to a target, so their offsets don't reall get considerd.
					if alignment_target not in self.valid_targets:
						if real_label == "Positive":
							fn[ma_corrected_positions] += 1
						else:
							tn[ma_corrected_positions] += 1
					
					else:
						model_cutoffs = cutoffs[ma_corrected_positions]
						
						diffs = np.subtract(bitscore, cutoffs)
						
						decider = np.mean(diffs)
						
						if real_label == "Positive":
							if evalue < evals[0]:
								eval_10e3_conf["tp"][ma_corrected_positions] += 1
							else:
								eval_10e3_conf["fn"][ma_corrected_positions] += 1
								
							if evalue < evals[1]:
								eval_10e10_conf["tp"][ma_corrected_positions] += 1
							else:
								eval_10e10_conf["fn"][ma_corrected_positions] += 1
								
							if evalue < evals[2]:
								eval_10e20_conf["tp"][ma_corrected_positions] += 1
							else:
								eval_10e20_conf["fn"][ma_corrected_positions] += 1
								
							if evalue < evals[3]:
								eval_10e30_conf["tp"][ma_corrected_positions] += 1
							else:
								eval_10e30_conf["fn"][ma_corrected_positions] += 1
								
							if pct_aln > aln_cut and pct_id > id_cut:
								staramr_conf["tp"][ma_corrected_positions] += 1
							else:
								staramr_conf["fn"][ma_corrected_positions] += 1
							
							if decider > 0:
								rocker_conf["tp"][ma_corrected_positions] += 1
							else:
								rocker_conf["fn"][ma_corrected_positions] += 1
						else:
							if evalue < evals[0]:
								eval_10e3_conf["fp"][ma_corrected_positions] += 1
							else:
								eval_10e3_conf["tn"][ma_corrected_positions] += 1
								
							if evalue < evals[1]:
								eval_10e10_conf["fp"][ma_corrected_positions] += 1
							else:
								eval_10e10_conf["tn"][ma_corrected_positions] += 1
								
							if evalue < evals[2]:
								eval_10e20_conf["fp"][ma_corrected_positions] += 1
							else:
								eval_10e20_conf["tn"][ma_corrected_positions] += 1
								
							if evalue < evals[3]:
								eval_10e30_conf["fp"][ma_corrected_positions] += 1
							else:
								eval_10e30_conf["tn"][ma_corrected_positions] += 1
								
							if pct_aln > aln_cut and pct_id > id_cut:
								staramr_conf["fp"][ma_corrected_positions] += 1
							else:
								staramr_conf["tn"][ma_corrected_positions] += 1
						
							if decider > 0:
								rocker_conf["fp"][ma_corrected_positions] += 1
							else:
								rocker_conf["tn"][ma_corrected_positions] += 1

				
				accuracy, sensitivity, specificity = self.confmat_to_stats(rocker_conf)
				e3acc, e3sens, e3spec = self.confmat_to_stats(eval_10e3_conf)
				e10acc, e10sens, e10spec = self.confmat_to_stats(eval_10e10_conf)
				e20acc, e20sens, e20spec = self.confmat_to_stats(eval_10e20_conf)
				e30acc, e30sens, e30spec = self.confmat_to_stats(eval_10e30_conf)
				staracc, starsens, starspec = self.confmat_to_stats(staramr_conf)
				
				#Record the results.
				self.testing_data_acc_sens_spec[rl][i] = {"Train_indices":train.index.tolist(), "Accuracy":accuracy, "Sensitivity":sensitivity, "Specificity":specificity}
				self.testing_data_asp_e3[rl][i] = {"Train_indices":train.index.tolist(), "Accuracy":e3acc, "Sensitivity":e3sens, "Specificity":e3spec}
				self.testing_data_asp_e10[rl][i] = {"Train_indices":train.index.tolist(), "Accuracy":e10acc, "Sensitivity":e10sens, "Specificity":e10spec}
				self.testing_data_asp_e20[rl][i] = {"Train_indices":train.index.tolist(), "Accuracy":e20acc, "Sensitivity":e20sens, "Specificity":e20spec}
				self.testing_data_asp_e30[rl][i] = {"Train_indices":train.index.tolist(), "Accuracy":e30acc, "Sensitivity":e30sens, "Specificity":e30spec}
				self.testing_data_asp_staramr[rl][i] = {"Train_indices":train.index.tolist(), "Accuracy":staracc, "Sensitivity":starsens, "Specificity":starspec}
			
		
			
	def calculate_youden(self, current_window_tgt, current_window_con):		
		#Select the maximum depth of coverage for the current window at each bitscore. Should always be increasing down the matrix.
		max_by_bitscore_tgt = np.amax(current_window_tgt, axis = 1)
		max_by_bitscore_con = np.amax(current_window_con, axis = 1)
		
		#The numbers of reads falling into each bitscore window, descending.
		#If I flip these matrices rows, find youden cuts, then flip the cut indices, does that work?
		#print(current_window_tgt)
		#print(current_window_con)
		#print("")
		
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
			cutoff = max_by_bitscore_tgt.shape[0]-1
			
		return cutoff, accuracy, sensitivity, specificity
		
	def calculate_roc_curves(self, training_data, rl):		
		#Ceiling of the number of full windows.
		half_window = int(self.window_size/2)
		
		#print(training_data["classifier"], training_data["pct_id"])
		
		def find_nearest(array, value):
			array = np.asarray(array)
			idx = (np.abs(array - value)).argmin()
			return array[idx]
			
		self.models = []
		#self.aln_len_models = []
		#self.pct_id_models = []
		
		self.senspec = {}
		
		min_bitscore = min(training_data["bitscore"])
		max_bitscore = max(training_data["bitscore"])
		
		#Determine the bin boundaries
		vert_bins = np.arange(min_bitscore, max_bitscore, self.resolution)
		
		#Set up data repo - rows corresp. to each vertical bin and are the length of the protein, 
		#Also divided by target vs. confounder.
		per_position_data = {"Target" : {}, "Confounder": {}}
		for vbin in vert_bins:
			#The range is always 0:max_pos
			#per_position_data["Target"][vbin] = np.zeros(self.max_position, dtype = np.int32)
			per_position_data["Target"][vbin] = np.zeros(self.multiple_align_size, dtype = np.int32)
			#per_position_data["Confounder"][vbin] = np.zeros(self.max_position, dtype = np.int32)
			per_position_data["Confounder"][vbin] = np.zeros(self.multiple_align_size, dtype = np.int32)
		
		#Iterate through the reads and fill the appropriate row/columns
		for classifier, s, e, bs, targ in zip(training_data["classifier"], training_data["start"], training_data["end"], training_data["bitscore"], training_data["target"]):
			
			#TODO adjust the position of the read placements
			falls_into = find_nearest(vert_bins, bs)
			
			lowhi = np.arange(s-1, e-1)
			fills_bins = self.offsets[targ][lowhi] + lowhi
			#print(targ)
			#print(s, e)
			#print(fills_bins.shape)
			#fills_bins = fills_bins[np.arange(s, e)]

			#print("")
			
			if classifier == "Positive":
				#per_position_data["Target"][falls_into][s:e] += 1
				per_position_data["Target"][falls_into][fills_bins] += 1
			else:
				#per_position_data["Confounder"][falls_into][s:e] += 1
				per_position_data["Confounder"][falls_into][fills_bins] += 1
		
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
		
		#Okay, data's collected, so we select windows and calc ROC from those.
		#for window_midpoint in range(0, self.max_position):
		for window_midpoint in range(0, self.multiple_align_size):
			#Get sliding window start, end indices; truncate at edges.
			window_start = window_midpoint-half_window
			if window_start < 0:
				window_start = 0
			window_end = window_midpoint + half_window
			if window_end > self.multiple_align_size:
				window_end = self.multiple_align_size
			
			#Select columns matching the window from cum sums
			current_window_tgt = collected_target_data[:, np.arange(window_start, window_end)]
			current_window_con = collected_confounder_data[:, np.arange(window_start, window_end)]
			
			#In-group cutoffs
			cutoff, accuracy, sensitivity, specificity = self.calculate_youden(current_window_tgt, current_window_con)
			
			#print(current_window_tgt)
			#print(current_window_con)
			#print(cutoff)
			#print(accuracy, sensitivity, specificity)
			#print("")
			
			if rl not in self.senspec:
				self.senspec[rl] = []
				
			self.senspec[rl].append((window_midpoint, accuracy[cutoff], sensitivity[cutoff], specificity[cutoff],))

			cutoff_bitscore = desc_bitscores[cutoff]
			
			data = (rl, window_midpoint, cutoff_bitscore,)
			
			self.models.append(data)
				
		self.models = pd.DataFrame(self.models)
		self.models.columns = ["read_length", "window_midpt", "bitscore_cutoff"]
		
		#For the CV module we just need one read length's values.
		return self.models["bitscore_cutoff"]
		
	def asp_plot(self, asp, name):
		for rl in asp:
			next_figure = go.Figure()
			MA_size = asp[rl][0]["Accuracy"].shape[0]
			midpoints = np.arange(0, MA_size, dtype = np.int32)
			for subsample_index in asp[rl]:
				next_figure.add_trace(go.Scatter(x = midpoints, 
							y = asp[rl][subsample_index]["Accuracy"],
							marker = dict(color = "red"),
							legendgroup = "Subsample " + str(subsample_index),
							name = "Accuracy")
							)
									
				next_figure.add_trace(go.Scatter(x = midpoints, 
							y = asp[rl][subsample_index]["Sensitivity"],
							marker = dict(color = "blue"),
							legendgroup = "Subsample " + str(subsample_index),
							name = "Sensitivity")
							)
				
				next_figure.add_trace(go.Scatter(x = midpoints, 
							y = asp[rl][subsample_index]["Specificity"],
							marker = dict(color = "yellow"),
							legendgroup = "Subsample " + str(subsample_index),
							name = "Specificity")
							)
				
			next_figure.update_xaxes(title_text='Position in Multiple Alignment')
			next_figure.update_yaxes(title_text='Accuracy, Sensitvity, and Specificity')
			next_figure.update_layout(hovermode="x unified")
			next_figure.update_layout(yaxis_range=[0,1])
			next_figure.write_html(os.path.normpath(self.outdir + "/" + name +"_acc_sens_spec_"+str(rl)+".html"))
	
	def plot_all_asps(self):
		self.asp_plot(self.testing_data_acc_sens_spec, "rocker")
		self.asp_plot(self.testing_data_asp_e3, "evalue_10neg3")
		self.asp_plot(self.testing_data_asp_e10, "evalue_10neg10")
		self.asp_plot(self.testing_data_asp_e20, "evalue_10neg20")
		self.asp_plot(self.testing_data_asp_e30, "evalue_10neg30")
		self.asp_plot(self.testing_data_asp_staramr, "staramr")
	
prj = sys.argv[1]			
mn = cross_validator(prj)
mn.ready_inputs() #Load data, set project manager. Has to go first

mn.prepare_offsets()
mn.find_valid_targets()
mn.find_filter()
mn.get_filter()
mn.interpolate()

mn.prepare_subsamples()

mn.plot_all_asps()
