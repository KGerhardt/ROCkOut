import sys
import os
import pandas as pd
import numpy as np

import plotly.express as px
import plotly.graph_objs as go

from .rocker_project_manager import project_manager
from .pplacer.rocker_phylomap_build import phylomap_build 

import subprocess

def import_rocker_project(project_directory):
	manager = project_manager(project_directory)
	manager.parse_project_directory()
	manager.parse_aligns()
	return manager
	
class plot_data:
	def __init__(self):
		self.project_dir = None
		self.index = None
		self.manager = None
		
		self.multiple_alignment = None
		
		self.description = None
		
		self.offsets = None
		self.multiple_align_size = 0
		self.valid_targets = None
		self.invalid_targets = None
		self.targets_for_writeout = None
		
		self.valid_origin_proteins = None

		self.positive_targets_file = None
		self.positive_targets_dmnd = None
		#Blastdb makes 6 files, so the prefix is enough
		self.positive_targets_blast = None
		
		self.read_files = None
		self.loaded_data = None
		self.successful_load = None
		
		self.readlens = None
		self.current_read_len = None
		
		self.active_prot_dict = {}
		self.active_parents = None
		self.active_proteins = None
		self.max_position = None
		
		self.senspec = None

		
		#How many AA will be considered at a time when calculating ROC - this a sliding window
		self.window_size = 20
		#How fine-grain will the y-axis be? 0.5 bitscore windows by default
		#Note: it's entirely possible to calculate the ROC without binning, but gains are negligible and it's slower.
		self.resolution = 0.5
		self.models = None
		self.aln_len_models = None
		self.pct_id_models = None
		
		self.figs = {}
	
	def set_project(self, project_directory):
		#try:
		if True:
			print("Parsing:", project_directory, "...")
			self.project_dir = project_directory
			
			out_path_base = os.path.normpath(self.project_dir + "/final_outputs")
			if not os.path.exists(out_path_base):
				os.mkdir(out_path_base)
			
			self.positive_targets_file = os.path.normpath(self.project_dir + "/final_outputs/database/positive_proteins_aa.fasta")
			self.positive_targets_dmnd = os.path.normpath(self.project_dir + "/final_outputs/database/positive_proteins_diamond_db.dmnd")
			#Blastdb makes 6 files, so the prefix is enough
			self.positive_targets_blast = os.path.normpath(self.project_dir + "/final_outputs/database/positive_proteins_blast_database")
			
			self.index = os.path.normpath(project_directory + "/ROCkOUT_index.txt")
			self.manager = project_manager(project_directory)
			self.manager.parse_project_directory()
			self.manager.parse_targets()
			self.manager.parse_aligns()
			self.manager.parse_multiple_alignment()
			
			
			self.prepare_offsets() #load a MA and determine position by position offset for each target protein.
			self.find_valid_targets()

			self.load_reads()
			
			#self.load_from_index()
			
			self.successful_load = True
			print("Project loaded!")
			
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
				
			#We need to hit both the main
			phylomap_build(pos = genomes_for_pplacer_pos,
							neg = genomes_for_pplacer_neg,
							output = out_path_base)	
		'''
		except:
			print("Couldn't parse directory!")
			self.successful_load = False
		'''
		
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
		self.invalid_targets = []
		self.targets_for_writeout = ""
		
		self.valid_origin_proteins = []
		#self.invalid_origin_genomes = []
		
		#copy to outputs
		
		
		self.manager.parse_coords()
		for prot in self.manager.targets:
			for file in self.manager.targets[prot]:				
				if '/positive/' in file:
					fh = open(file)
					for line in fh:
						self.targets_for_writeout += line
						if line.startswith(">"):
							next_item = line.strip().split()[0][1:]
							origin_protein = next_item.split("__")[-1]
							#origin_genome = origin_genome.split(".")[0]
							self.valid_origin_proteins.append(origin_protein)
							self.valid_targets.append(next_item)
					fh.close()
				else: #is negative
					fh = open(file)
					for line in fh:
						self.targets_for_writeout += line
						if line.startswith(">"):
							next_item = line.strip().split()[0][1:]
					fh.close()
				
		self.valid_targets = set(self.valid_targets)
		#self.invalid_targets = set(self.invalid_targets)
		
		self.valid_origin_proteins = set(self.valid_origin_proteins)
		#self.invalid_origin_genomes = set(self.invalid_origin_genomes)
		
	def load_reads(self):
		self.loaded_data = {}
		
		read_best_hits = {}
		
		for alignments in self.manager.alignments_pos:
			for file in self.manager.alignments_pos[alignments]:
				#File names contain info and always have the same seps. 
				#To avoid OS issues, we just get the name after the final '/' in the path with os.basename
				relevant = os.path.basename(file)
				components = relevant.split("_read_len_")
				#Specific protein name
				title = components[0]
				#read length
				avg_read_length = components[1].split("_aligned_reads.blast.txt")[0]
				#Initialize this list.
				
				if avg_read_length not in self.loaded_data:
					self.loaded_data[avg_read_length] = []
					read_best_hits[avg_read_length] = {}
					
				origin_genome = os.path.basename(file)
				origin_genome = origin_genome.split("_read_len_")[0]
				fh = open(file)
				for line in fh:
					#This is blast data.
					segs = line.strip().split("\t")
					#target, off-target, or negative.
					read_classifier = segs[0].split(";")
					read_id = ';'.join(read_classifier[:-1])
					read_classifier = read_classifier[len(read_classifier)-1]
					
					alignment_target = segs[1]
					
					if alignment_target not in self.valid_targets:
						continue
						
					#Positive foreign targets can be skipped - these are either:
					#(1) Negative foreign targets which will be seen as negatives in the negatives directories
					#(2) Positive foreign targets which will be seen as targets in the positives directories
					if read_classifier == "Foreign_Target":
						continue
						
					alignment_range = [int(segs[8]), int(segs[9])]
					
					start = min(alignment_range)
					end = max(alignment_range)
					
					mid_finder = np.arange(start-1, end-1)
					
					offset_locs = mid_finder + self.offsets[alignment_target][mid_finder]
					midpt = np.median(offset_locs)
					
					#start, end = start + offset_locs[start], end + offset_locs[end]
					#midpt = np.median(self.offsets[alignment_target][mid_finder])
					
					#midpt = int((start+end)/2)
					
					bitscore = float(segs[11])
					
					pct_id = float(segs[2])
					alignment_length = int(segs[3])
					
					double_key = alignments + ", " + title
					
					data = (read_id, 
							avg_read_length, 
							alignments, 
							title, 
							read_classifier, 
							alignment_target, 
							midpt, 
							start, 
							end, 
							bitscore,
							pct_id,
							alignment_length,
							double_key,)
										
					if read_id not in read_best_hits[avg_read_length]:
						read_best_hits[avg_read_length][read_id] = data
					else:
						#read ID has already been found, so the read either aligns multiple times or comes from multiple copies of the same genome
						existing_bitscore = read_best_hits[avg_read_length][read_id][9]
						existing_label = read_best_hits[avg_read_length][read_id][4]
						
						#if true, these are almost certainly the same read from the same genome 
						#under different UniProt IDs - favor the one labeled as target
						if np.isclose(bitscore, existing_bitscore):
							#replace a non-target with a target
							if read_classifier == "Target" and existing_label != "Target":
								read_best_hits[avg_read_length][read_id] = data
						else:
							if bitscore > existing_bitscore: #Displace a clearly lower bitscore either way
								read_best_hits[avg_read_length][read_id] = data
						
					#self.loaded_data[avg_read_length].append(data)
				fh.close()
				
		for alignments in self.manager.alignments_neg:
			for file in self.manager.alignments_neg[alignments]:
				#File names contain info and always have the same seps. 
				#To avoid OS issues, we just get the name after the final '/' in the path with os.basename
				relevant = os.path.basename(file)
				components = relevant.split("_read_len_")
				#Specific protein name
				title = components[0]
				#read length
				avg_read_length = components[1].split("_aligned_reads.blast.txt")[0]
				
				if avg_read_length not in self.loaded_data:
					self.loaded_data[avg_read_length] = []
					
				fh = open(file)
				for line in fh:
					#This is blast data.
					segs = line.strip().split("\t")
					#target, off-target, or negative.
					read_classifier = segs[0].split(";")
					read_id = ';'.join(read_classifier[:-1])
					read_classifier = read_classifier[len(read_classifier)-1]

					#read_id = segs[0]
					#read_classifier = read_classifier[-1]
					
					#read_classifier = "Negative" # We do not need to check this for negatives
					alignment_target = segs[1]
			
					if alignment_target not in self.valid_targets:
						continue	

					if read_classifier == "Foreign_Target":
						origin_protein = segs[-2]
						if origin_protein in self.valid_origin_proteins:
							read_classifier = "Target"
						else:
							read_classifier = "Negative"
					else:
						read_classifier = "Negative"
			
					alignment_range = [int(segs[8]), int(segs[9])]
					start = min(alignment_range)
					end = max(alignment_range)
					
					#print(alignment_range)
					mid_finder = np.arange(start-1, end-1)
					
					offset_locs = mid_finder + self.offsets[alignment_target][mid_finder]
					midpt = np.median(offset_locs)
					
					bitscore = float(segs[11])
					pct_id = float(segs[2])
					alignment_length = int(segs[3])
					
					double_key = alignments + ", " + title
					
					data = (read_id, 
							avg_read_length, 
							alignments, 
							title, 
							read_classifier, 
							alignment_target, 
							midpt, 
							start, 
							end, 
							bitscore,
							pct_id,
							alignment_length,
							double_key,)
					
					#data = (avg_read_length, alignments, title, read_classifier, alignment_target, midpt, start, end, bitscore, )
					
					if read_id not in read_best_hits[avg_read_length]:
						read_best_hits[avg_read_length][read_id] = data
					else:	
						existing_bitscore = read_best_hits[avg_read_length][read_id][9]
						existing_label = read_best_hits[avg_read_length][read_id][4]
						
						#if true, these are almost certainly the same read from the same genome 
						#under different UniProt IDs - favor the one labeled as target
						if np.isclose(bitscore, existing_bitscore):
							#replace a non-target with a target
							if read_classifier == "Target" and existing_label != "Target":
								read_best_hits[avg_read_length][read_id] = data
						else:
							if bitscore > existing_bitscore: #Displace a clearly lower bitscore either way
								read_best_hits[avg_read_length][read_id] = data
					
					#self.loaded_data[avg_read_length].append(data)
				fh.close()

		self.active_parents = set()
		self.active_proteins = set()
			
		self.max_position = 0

		for rl in read_best_hits:
			self.loaded_data[rl] = []
			for read_id in read_best_hits[rl]:
				#print(read_best_hits[rl][read_id])
				self.loaded_data[rl].append(read_best_hits[rl][read_id])
			read_best_hits[rl] = None
			
		
		for rl in self.loaded_data:
			self.loaded_data[rl] = pd.DataFrame(self.loaded_data[rl])
			
			self.loaded_data[rl].columns = ["read_id", "read_length", "parent", "protein", 
											"classifier", "target", "midpoint", "start", "end", 
											"bitscore", "pct_id", "aln_len", "UniProt ID, protein name"]
			
			#print(self.loaded_data[rl])
			#This should probably be reimplemented up above.
			#clean the dataframes to best hit by read - this takes care of repeat genome simulations
			self.loaded_data[rl] = self.loaded_data[rl].loc[self.loaded_data[rl].reset_index().groupby(['read_id'])['bitscore'].idxmax()]
			
			#Add max as needed.
			self.max_position = max(self.max_position, max(self.loaded_data[rl]["end"]))
			self.active_parents = self.active_parents.union(self.loaded_data[rl]["parent"])
			self.active_proteins = self.active_proteins.union(self.loaded_data[rl]["protein"]) 
			
		#print(self.loaded_data)
			
		self.readlens = list(self.loaded_data.keys())
		self.current_read_len = min(self.readlens)
			
		self.active_parents = list(self.active_parents)
		self.active_proteins = list(self.active_proteins)
		self.active_prot_dict = {}
		for prot in self.active_proteins:
			self.active_prot_dict[prot] = True
		
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
		
	def calculate_roc_curves(self):		
		#Ceiling of the number of full windows.
		half_window = int(self.window_size/2)
		
		def find_nearest(array, value):
			array = np.asarray(array)
			idx = (np.abs(array - value)).argmin()
			return array[idx]
			
		self.models = []
		self.aln_len_models = []
		self.pct_id_models = []
		
		self.senspec = {}

		#Separate the data into groups of mean read lengths.
		#read_lengths = list(set(current_data["read_length"]))
		
		#From meeting:
		#We're going to repeat this model building for alignment length and percent identity 
		#and use the three models as an ensemble for the filter step.
		
		for rl in self.loaded_data:
			#Since each model is calculated at one read length and we interpolate read lengths between them,
			#We select the data associated with each readlength and get the model for it.
			one_length = self.loaded_data[rl].copy()
			#print(one_length)
			
			#Select out active prots
			one_length = one_length.loc[one_length['protein'].isin(self.active_proteins)]
			
			min_bitscore = min(one_length["bitscore"])
			max_bitscore = max(one_length["bitscore"])
			
			#min_aln_len = min(one_length[""])
			#This could be repeated for aln length, pct ID for more models...
			
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
			for classifier, s, e, bs, targ in zip(one_length["classifier"], one_length["start"], one_length["end"], one_length["bitscore"], one_length["target"]):
				
				#TODO adjust the position of the read placements
				falls_into = find_nearest(vert_bins, bs)
				
				lowhi = np.arange(s-1, e-1)
				fills_bins = self.offsets[targ][lowhi] + lowhi
				#print(targ)
				#print(s, e)
				#print(fills_bins.shape)
				#fills_bins = fills_bins[np.arange(s, e)]

				#print("")
				
				if classifier == "Target":
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
				
				if rl not in self.senspec:
					self.senspec[rl] = []
					
				self.senspec[rl].append((window_midpoint, accuracy[cutoff], sensitivity[cutoff], specificity[cutoff],))

				cutoff_bitscore = desc_bitscores[cutoff]
				
				data = (rl, window_midpoint, cutoff_bitscore,)
				
				self.models.append(data)
				
		self.models = pd.DataFrame(self.models)
		self.models.columns = ["read_length", "window_midpt", "bitscore_cutoff"]
		
	def craft_plot(self):
		if self.loaded_data is not None:
			
			for rl in self.loaded_data:
				current_data = self.loaded_data[rl].copy()
				#We do have the parent involved...
				#print(current_data)
				#Not needed.
				#current_data = current_data.loc[current_data['read_length']==self.current_read_len]
				current_data = current_data.loc[current_data['protein'].isin(self.active_proteins)]
				
				fig = px.scatter(current_data, x="midpoint", y ="bitscore", color = 'classifier', hover_data = ['read_id', 'target'], 
				color_discrete_sequence = ["blue", "aqua", "darkorange", "green"], 
				category_orders={"classifier": ["Target", "Probable_Target", "Non_Target", "Negative"]}, 
				symbol = "UniProt ID, protein name")
				
				fig.update_layout(scene = dict(
							xaxis_title='Position in Protein',	
							yaxis_title='Bitscore'))
							
				fig_title = "scatter_2d_read_len_" + str(rl)
				#Add to list with title.
				self.figs[fig_title] = fig
		
	def craft_reads_plot(self):
		if self.loaded_data is not None:
			for rl in self.loaded_data:
				current_data = self.loaded_data[rl].copy()
				
				current_data = current_data.loc[current_data['protein'].isin(self.active_proteins)]
				#hit_type_to_color = {"Target":"blue", "Non_Target":"darkorange", "Negative":"green"}
				hit_type_to_color = {"Target":"blue", "Probable_Target":"aqua", "Non_Target":"darkorange", "Negative":"green"}
				hit_colors = [hit_type_to_color[t] for t in current_data["classifier"]]
				current_data.insert(0, "color", hit_colors)
				
				# Horizontal line shape
				shapes=[dict(
						type='line',
						x0 = current_data['start'].iloc[i],
						y0 = current_data["bitscore"].iloc[i],
						x1 = current_data['end'].iloc[i],
						y1 = current_data["bitscore"].iloc[i],
						line = dict(
							color = current_data["color"].iloc[i],
							width = 0.5
						)
					) for i in range(len(current_data))]
					
				layout = go.Layout(
					shapes = shapes,
					scene = dict(xaxis_title='Position in Protein',
								yaxis_title='Bitscore'))
								
				# Plot the chart
				fig = go.Figure(go.Scatter(), layout)
				
				fig_title = "reads_plot_read_len_" + str(rl)
				#Add to list with title.
				self.figs[fig_title] = fig
				
	def craft_3d_plot(self):
		if self.loaded_data is not None:
			#Pop the readlens together into one DF
			current_data = pd.concat(list(self.loaded_data.values()))
			current_data = current_data.loc[current_data['protein'].isin(self.active_proteins)]
			
			fig = px.scatter_3d(current_data, x="midpoint", y ="bitscore", z="read_length", color = 'classifier',
				color_discrete_sequence = ["blue", "aqua", "darkorange", "green"], 
				category_orders={"classifier": ["Target", "Probable_Target", "Non_Target", "Negative"]})
			
			fig.update_layout(scene = dict(
						xaxis_title='Position in Protein',
						yaxis_title='Bitscore',
						zaxis_title='Read Length'))
						
			fig_title = "scatter_3d"
			#Add to list with title.
			self.figs[fig_title] = fig
								
	def craft_roc_plot(self):
		if self.models is None:
			if self.loaded_data is not None:
				self.calculate_roc_curves()
		
		if self.models is not None:
			current_data = self.models.copy()
			fig = px.line_3d(current_data, x="window_midpt", y="bitscore_cutoff", z="read_length", line_group="read_length")
			fig.update_layout(scene = dict(
						xaxis_title='Position in Protein',
						yaxis_title='Read Length',
						zaxis_title='Bitscore'))
			
			fig_title = "roc_plot_3d"
			#Add to list with title.
			#self.figs[fig_title] = fig
			if 'scatter_3d' in self.figs:
				self.figs['scatter_3d'] = go.Figure(data=self.figs['scatter_3d'].data + fig.data)
			
			subs = []
			for rl in self.loaded_data:
				sub = current_data[current_data['read_length'] == rl]
				fig = px.line(sub, x="window_midpt", y= "bitscore_cutoff")
				fig.update_layout(scene = dict(
						xaxis_title='Position in Protein',
						yaxis_title='Bitscore'))
						
				fig_title = "roc_plot_"+str(rl)
				
				#out_path_base = os.path.normpath(self.project_dir + "/final_outputs/figures")
				#fig.write_html(out_path_base+"/"+fig_title+".html")
				
				self.figs[fig_title] = fig
				
				
				corresp = "scatter_2d_read_len_"+str(rl)
				if corresp in self.figs:
					self.figs[corresp] = go.Figure(data=self.figs[corresp].data + fig.data)
	
	def craft_acc_sens_spec(self):
		if self.senspec is not None:
			for rl in self.senspec:
				one_data = pd.DataFrame(self.senspec[rl])
				one_data.columns = ["midpoint", "acc", "sens", "spec"]
				one_data = one_data.melt(id_vars = ["midpoint"], var_name = 'Measure')
				
				fig = px.line(one_data, x="midpoint", y ="value", color = 'Measure')
				fig.update_layout(scene = dict(
							xaxis_title='Position in Protein',
							yaxis_title='Metric value'))
				fig.update_layout(yaxis_range=[0,1])
							
				fig_title = "classifier_plot_read_len_" + str(rl)
				#Add to list with title.
				self.figs[fig_title] = fig
	
	def output_plots(self):
		#Create output folders if needed
		out_path_base = os.path.normpath(self.project_dir + "/final_outputs/figures")

		if not os.path.exists(out_path_base):
			os.mkdir(out_path_base)
		out_path_base = os.path.normpath(out_path_base)
	
		config = {
			'scrollZoom': True,
				'toImageButtonOptions': {
				'format': 'svg', # one of png, svg, jpeg, webp
				'filename': 'custom_image',
				'height': 1080,
				'width': 1920,
				'scale': 1 # Multiply title/legend/axis/canvas sizes by this factor
		  }
		}	
	
		#Save interactive plots.
		for plot_title in self.figs:
			#We don't need the separate roc plots without the read overlays.
			if not plot_title.startswith("roc_plot_"):
				print("Printing", plot_title)
				self.figs[plot_title].write_html(os.path.normpath(out_path_base+"/"+plot_title+".html"), config = config)

	def output_models(self):
		out_path_base = os.path.normpath(self.project_dir + "/final_outputs")
		if not os.path.exists(out_path_base):
			os.mkdir(out_path_base)
			
		dbout = os.path.normpath(self.project_dir + "/final_outputs/database/")
		if not os.path.exists(dbout):
			os.mkdir(dbout)
			
		out_path_base += "/model"
		out_path_base = os.path.normpath(out_path_base)
		if not os.path.exists(out_path_base):
			os.mkdir(out_path_base)

		output_senspec = os.path.normpath(out_path_base + "/accuracy_sensitivity_and_specificity.txt")
		sp = open(output_senspec, "w")
		print("read_length", "window_midpt", "accuracy", "sensitivity", "specificity", sep = "\t", file = sp)
		for rl in self.senspec:
			for row in self.senspec[rl]:
				print(rl, row[0], row[1], row[2], row[3], sep = "\t", file = sp)
		sp.close()
		
		
		fh = open(self.positive_targets_file, "w")
		fh.write(self.targets_for_writeout)
		fh.close()
		
		try:
			makedb = ["diamond", "makedb", "--db", self.positive_targets_dmnd,  "--in", self.positive_targets_file]
			print("Building Diamond database for positive targets. Log information will follow.")
			subprocess.call(makedb)
		except:
			print("Couldn't make DIAMOND database of positive targets!")
		
		try:
			#makeblastdb -in <reference.fa> -dbtype nucl -parse_seqids -out <database_name> -title "Database title"
			makedb = ["makeblastdb", "-in", self.positive_targets_file, "-parse_seqids", "-dbtype", "prot", "-out", self.positive_targets_blast]
			print(" ".join(makedb))
			print("Building BLAST database for positive targets. Log information will follow.")
			subprocess.call(makedb)
		except:
			print("Couldn't make BLAST database of positive targets!")
		
		output_filter = os.path.normpath(out_path_base + "/ROCkOut_Filter.txt")
		self.models.to_csv(output_filter, sep = "\t", index = False)
		
		
		reads_path = os.path.normpath(self.project_dir + "/final_outputs/reads/")
		if not os.path.exists(reads_path):
			os.mkdir(reads_path)
		
		for rl in self.loaded_data:
			this_output = os.path.normpath(reads_path + "/read_len_" + str(rl)+"_final_reads.tsv")
			self.loaded_data[rl].to_csv(this_output, sep = "\t")
			
		#Copy the multiple alignment to the model directory
		ma_path = os.path.normpath(self.project_dir + "/shared_files/multiple_alignment/complete_multiple_alignment_aa.fasta" )
		out_ma = os.path.normpath(self.project_dir + "/final_outputs/model/complete_multiple_alignment_aa.fasta")
		fh = open(ma_path)
		out = open(out_ma, "w")
		for line in fh:
			out.write(line)
		out.close()
		fh.close()
		
	def update_index(self):
		pass
		
	def load_from_index(self):
		record = {}
		fh = open(self.index)
		fh.readline()
		for line in fh:
			segs = line.strip().split("\t")
			parent, protein_name, readlen, num_reads, group, included = segs[0], segs[1], int(segs[2]), int(segs[3]), segs[4], segs[5]=="T"
			#print(parent, protein_name, readlen, num_reads, group, included)
			if not included:
				print(parent, protein_name, readlen, num_reads, group, included)
				self.active_proteins.pop(self.active_proteins.index(protein_name))
				
		fh.close()
			
	def just_make_plots(self):
		self.craft_plot()
		self.craft_reads_plot()
		self.craft_3d_plot()
		self.craft_roc_plot()
		
	def just_roc(self):
		self.calculate_roc_curves()
		self.craft_roc_plot()
			
	def run_plotter(self, directory):
		self.set_project(directory)
		self.calculate_roc_curves()
		self.craft_plot()
		self.craft_reads_plot()
		self.craft_3d_plot()
		self.craft_roc_plot()
		self.craft_acc_sens_spec()



def non_interactive(parser, opts):
	mn = plot_data()
	dir = opts.dir
	mn.run_plotter(dir)
	mn.output_plots()
	mn.output_models()
	