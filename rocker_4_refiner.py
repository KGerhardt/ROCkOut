import sys
import os
import pandas as pd
import numpy as np

import plotly.express as px
import plotly.graph_objs as go

from .rocker_project_manager import project_manager
#from rocker_project_manager import project_manager

def import_rocker_project(project_directory):
	manager = project_manager(project_directory)
	manager.parse_project_directory()
	manager.parse_aligns()
	return manager
	
	
class plot_data:
	def __init__(self):
		self.project_dir = None
		self.manager = None
		
		self.description = None
		
		self.read_files = None
		self.loaded_data = None
		self.successful_load = None
		
		self.readlens = None
		self.current_read_len = None
		
		self.active_prot_dict = {}
		self.active_parents = None
		self.active_proteins = None
		self.max_position = None
		
		#How many AA will be considered at a time when calculating ROC - this a sliding window
		self.window_size = 20
		#How fine-grain will the y-axis be? 0.5 bitscore windows by default
		#Note: it's entirely possible to calculate the ROC without binning, but gains are negligible and it's slower.
		self.resolution = 0.5
		self.models = None
		
		self.figs = {}
	
	def set_project(self, project_directory):
		try:
			print("Parsing:", project_directory, "...")
			self.project_dir = project_directory
			self.manager = project_manager(project_directory)
			self.manager.parse_project_directory()
			self.manager.parse_aligns()
			self.load_reads()
			self.successful_load = True
			print("Project loaded!")
		except:
			print("Couldn't parse directory!")
			self.successful_load = False
			
	def load_reads(self):
		self.loaded_data = {}
		for alignments in self.manager.alignments_pos:
			for file in self.manager.alignments_pos[alignments]:
				#File names contain info and always have the same seps. 
				#To avoid OS issues, we just get the name after the final '/' in the path with os.basename
				relevant = os.path.basename(file)
				components = relevant.split("_read_length_")
				#Specific protein name
				title = components[0]
				#read length
				avg_read_length = components[1].split("_aligned_reads.txt")[0]
				#Initialize this list.
				if avg_read_length not in self.loaded_data:
					self.loaded_data[avg_read_length] = []
				fh = open(file)
				for line in fh:
					#This is blast data.
					segs = line.strip().split("\t")
					#target, off-target, or negative.
					read_classifier = segs[0].split(";")
					read_classifier = read_classifier[len(read_classifier)-1]
					alignment_target = segs[1]
					alignment_range = [int(segs[8]), int(segs[9])]
					start = min(alignment_range)
					end = max(alignment_range)
					midpt = int((start+end)/2)
					
					bitscore = float(segs[11])
					
					data = (avg_read_length, alignments, title, read_classifier, alignment_target, midpt, start, end, bitscore, )
					
					self.loaded_data[avg_read_length].append(data)
				fh.close()
				
		for alignments in self.manager.alignments_neg:
				
			for file in self.manager.alignments_neg[alignments]:
				#File names contain info and always have the same seps. 
				#To avoid OS issues, we just get the name after the final '/' in the path with os.basename
				relevant = os.path.basename(file)
				components = relevant.split("_read_length_")
				#Specific protein name
				title = components[0]
				#read length
				avg_read_length = components[1].split("_aligned_reads.txt")[0]
				if avg_read_length not in self.loaded_data:
					self.loaded_data[avg_read_length] = []
				fh = open(file)
				for line in fh:
					#This is blast data.
					segs = line.strip().split("\t")
					#target, off-target, or negative.
					read_classifier = segs[0].split(";")
					read_classifier = read_classifier[len(read_classifier)-1]
					alignment_target = segs[1]
					alignment_range = [int(segs[8]), int(segs[9])]
					start = min(alignment_range)
					end = max(alignment_range)
					midpt = int((start+end)/2)
					
					bitscore = float(segs[11])
					
					data = (avg_read_length, alignments, title, read_classifier, alignment_target, midpt, start, end, bitscore, )
					
					self.loaded_data[avg_read_length].append(data)
				fh.close()

		self.active_parents = set()
		self.active_proteins = set()
			
		self.max_position = 0
			
		for rl in self.loaded_data:
			self.loaded_data[rl] = pd.DataFrame(self.loaded_data[rl])
			self.loaded_data[rl].columns = ["read_length", "parent", "protein", "classifier", "target", "midpoint", "start", "end", "bitscore"]
			#Add max as needed.
			self.max_position = max(self.max_position, max(self.loaded_data[rl]["end"]))
			self.active_parents = self.active_parents.union(self.loaded_data[rl]["parent"])
			self.active_proteins = self.active_proteins.union(self.loaded_data[rl]["protein"]) 

		self.readlens = list(self.loaded_data.keys())
		self.current_read_len = min(self.readlens)
			
		self.active_parents = list(self.active_parents)
		self.active_proteins = list(self.active_proteins)
		self.active_prot_dict = {}
		for prot in self.active_proteins:
			self.active_prot_dict[prot] = True
		
	def calculate_roc_curves(self):		
		#Ceiling of the number of full windows.
		half_window = int(self.window_size/2)
		
		def find_nearest(array, value):
			array = np.asarray(array)
			idx = (np.abs(array - value)).argmin()
			return array[idx]
			
		self.models = []
		
		#Separate the data into groups of mean read lengths.
		#read_lengths = list(set(current_data["read_length"]))
		for rl in self.loaded_data:
			#Since each model is calculated at one read length and we interpolate read lengths between them,
			#We select the data associated with each readlength and get the model for it.
			one_length = self.loaded_data[rl].copy()
			
			#Select out active prots
			one_length = one_length.loc[one_length['protein'].isin(self.active_proteins)]
			
			min_bitscore = min(one_length["bitscore"])
			max_bitscore = max(one_length["bitscore"])
			
			#Determine the bin boundaries
			vert_bins = np.arange(min_bitscore, max_bitscore, self.resolution)
			
			#Set up data repo - rows corresp. to each vertical bin and are the length of the protein, 
			#Also divided by target vs. confounder.
			per_position_data = {"Target" : {}, "Confounder": {}}
			for vbin in vert_bins:
				#The range is always 0:max_pos
				per_position_data["Target"][vbin] = np.zeros(self.max_position, dtype = np.int32)
				per_position_data["Confounder"][vbin] = np.zeros(self.max_position, dtype = np.int32)
			
			#Iterate through the reads and fill the appropriate row/columns
			for classifier, s, e, bs in zip(one_length["classifier"], one_length["start"], one_length["end"],one_length["bitscore"]):
				falls_into = find_nearest(vert_bins, bs)
				
				if classifier == "Target":
					per_position_data["Target"][falls_into][s:e] += 1
				else:
					per_position_data["Confounder"][falls_into][s:e] += 1
			
			
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
			for window_midpoint in range(0, self.max_position):
				#Get sliding window start, end indices; truncate at edges.
				window_start = window_midpoint-half_window
				if window_start < 0:
					window_start = 0
				window_end = window_midpoint + half_window
				if window_end > self.max_position:
					window_end = self.max_position
				
				#Select columns matching the window from cum sums
				current_window_tgt = collected_target_data[:, np.arange(window_start, window_end)]
				current_window_con = collected_confounder_data[:, np.arange(window_start, window_end)]
				
				#Select the maximum depth of coverage for the current window at each bitscore. Should always be increasing down the matrix.
				max_by_bitscore_tgt = np.amax(current_window_tgt, axis = 1)
				max_by_bitscore_con = np.amax(current_window_con, axis = 1)
				
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
				
				cutoff_bitscore = desc_bitscores[cutoff]
				#print(rl, window_start, window_end, cutoff, Youden_by_bs[cutoff], desc_bitscores[cutoff])
				
				data = (rl, window_midpoint, cutoff_bitscore,)
				
				self.models.append(data)
				
		self.models = pd.DataFrame(self.models)
		self.models.columns = ["read_length", "window_midpt", "bitscore_cutoff"]
					
	def craft_plot(self):
		if self.loaded_data is not None:
			
			for rl in self.loaded_data:
				current_data = self.loaded_data[rl].copy()
				#Not needed.
				#current_data = current_data.loc[current_data['read_length']==self.current_read_len]
				current_data = current_data.loc[current_data['protein'].isin(self.active_proteins)]
				
				fig = px.scatter(current_data, x="midpoint", y ="bitscore", color = 'classifier', symbol = 'protein')
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
				hit_type_to_color = {"Target":"blue", "Non_Target":"orange", "Negative":"green"}
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
			
			fig = px.scatter_3d(current_data, x="midpoint", y ="bitscore", z="read_length", color = 'classifier')
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
	
	def output_plots(self):
		#Create output folders if needed
		out_path_base = os.path.normpath(self.project_dir + "/final_outputs")
		if not os.path.exists(out_path_base):
			os.mkdir(oout_path_base)
		out_path_base += "/figures"
		if not os.path.exists(out_path_base):
			os.mkdir(out_path_base)
		#Save interactive plots.
		for plot_title in self.figs:
			print("Printing", plot_title)
			self.figs[plot_title].write_html(out_path_base+"/"+plot_title+".html")

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
	

'''
mn = plot_data()
mn.run_plotter("rockout_out")	
'''
