import sys
import argparse
import os
import numpy as np

import plotly.graph_objects as go
import pandas as pd

#from modules.rocker_project_manager import project_manager

class rocker_filterer:
	def __init__(self, project_directory, option, raws = None, reads = None):
		self.proj_dir = project_directory
		self.plot = option
		self.filter_file = None
		self.ma_file = None
		self.offsets = None
		self.filter_model = None
		
		self.observed_readlens = None
		
		self.filter_matrix = None
		
		self.raw_reads = raws
		self.reads_to_filter = reads
		
		self.run_filterer()
	
	def find_filter(self):
		if os.path.exists(self.proj_dir):
			filter_path = os.path.normpath(self.proj_dir + "/final_outputs/model/ROCkOut_Filter.txt")
			#print(filter_path)
			if os.path.exists(filter_path):
				#Only set this to not None if the file can be found.
				self.filter_file = filter_path
			else:
				print("Project appears incomplete.")
		else:
			print("Project not found.")
		
	def find_ma(self):
		if os.path.exists(self.proj_dir):
			ma_path = os.path.normpath(self.proj_dir + "/shared_files/multiple_alignment/complete_multiple_alignment_aa.fasta" )
			if os.path.exists(ma_path):
				#Only set this to not None if the file can be found.
				self.ma_file = ma_path
				self.get_offsets()
			else:
				print("Project appears incomplete.")
		else:
			print("Project not found.")
		
	def get_offsets(self):
		current_seq = ""
		current_prot = ""
		self.offsets = {}
		fh = open(self.ma_file)
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
		
		#We should convert the matrix into a dict of readlength:cutoffs
			
	def plot_filter(self):
		as_pandas = pd.DataFrame(self.filter_matrix)
		fig = go.Figure(data=[go.Surface(z=as_pandas.values)])
		fig.update_layout(scene = dict(
							xaxis_title='Position in Protein',	
							yaxis_title='Read Length',
							zaxis_title='Bitscore Cutoff'))
		fig.write_html(os.path.normpath("ROCkOut_surface.html"))

	#This now needs to take the multiple alignment info.	
	def run_filter_blast_tabular(self, output_pass = None, output_fail = None):
		
		
		if self.raw_reads is None:
			self.raw_reads = [None] * len(self.reads_to_filter)
	
		for raw_reads, reads_file in zip(self.raw_reads, self.reads_to_filter):
			if output_pass is None:
				output_pass = os.path.basename(reads_file)
				while output_pass != os.path.splitext(output_pass)[0]:
					output_pass = os.path.splitext(output_pass)[0]
				output_pass += "_rockout_filter_passing_reads.fasta"
				
				output_fail = output_pass
				output_fail = output_fail.replace("_passing_", "_failing_")
			
			
			
			has_readlen = False
			mapping_reads = []
			fh = open(reads_file)
			for line in fh:
				segs = line.strip().split("\t")
				if len(segs) > 12:
					has_readlen = True
					break
				alignment_target = segs[1]
				if alignment_target not in self.offsets:
					continue
				else:
					mapping_reads.append(segs[0])
			fh.close()
			
			#If has readlen, this is an empty set
			mapping_reads = set(mapping_reads)
			
			if not has_readlen:
				cur_sl = 0
				cur_read = ""
				is_valid = False
				read_lengths = {}
				fh = open(raw_reads)
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
			
			passing = open(output_pass, "w")
			failing = open(output_fail, "w")
			fh = open(reads_file)
			for line in fh:
				segs = line.strip().split("\t")
				alignment_target = segs[1]
				if alignment_target not in self.offsets:
					continue
				else:
					read_name = segs[0]
					
					if has_readlen:
						readlength = int(segs[12])
					else:
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
					
					#readlength = end-start #This won't work for this setup
					
					cutoffs = self.filter_matrix[readlength][ma_corrected_positions]
					
					diffs = np.subtract(bitscore, cutoffs)
					
					decider = np.mean(diffs)
					if decider > 0:
						passing.write(line)
					else:
						failing.write(line)
				
			fh.close()
			passing.close()
			failing.close()
			
	def run_filterer(self):
		
		self.find_filter()
		self.find_ma()
		
		if self.filter_file is None:
			print("ROCkOut filter file not found! Quitting")
			sys.exit()
		if self.ma_file is None:
			print("ROCkOut multiple alignment file not found! Quitting")
			sys.exit()
			
		self.get_filter()
		self.interpolate()
	
		if self.plot:
			self.plot_filter()
		else:
			self.run_filter_blast_tabular()
			
			
#inf = sys.argv[1]
proj_dir = sys.argv[1]
r = [sys.argv[2]]

if len(sys.argv) > 3:
	raw = [sys.argv[3]]
else:
	raw = None
	
#print(r)
#print(raw)

'''
if len(sys.argv) > 2:
	plot_or_filter = sys.argv[2]
	plot_filter = sys.argv[2] == "2"
'''
plot_filter = False

#r = ["test_redux/positive/P02979/aligned_reads/V01278_read_len_110_aligned_reads.fasta"]
#raw = ["test_redux/positive/P02979/tagged_reads/V01278_read_len_110_tagged.fasta"]

mn = rocker_filterer(proj_dir, option = plot_filter, raws = raw, reads = r)
#mn.run_filterer()