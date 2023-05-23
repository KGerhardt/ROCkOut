import sys
import argparse
import os
import numpy as np

import plotly.graph_objects as go
import pandas as pd

#from modules.rocker_project_manager import project_manager

class rocker_filterer:
	def __init__(self, project_directory, filter_dir):
		self.proj_dir = project_directory
		self.filt_dir = filter_dir
		
		self.alns = self.filt_dir + "/alignments"
		self.orig = self.filt_dir + "/original_reads"
		self.passing = self.filt_dir + "/ROCkOut_passing_alignments"
		self.failing = self.filt_dir + "/ROCkOut_failing_alignments"
		self.fastas = self.filt_dir + "/ROCkOut_passing_read_fastas"
		
		#self.plot = option
		self.filter_file = None
		self.ma_file = None
		self.offsets = None
		self.filter_model = None
		
		self.observed_readlens = None
		
		self.filter_matrix = None
		
		self.raw_reads = None
		self.reads_to_filter = None
		self.filtered_passing = None
		self.passing_alns = None
		self.failing_alns = None
	
	def find_filter(self):
		if os.path.exists(self.proj_dir):
			filter_path = os.path.normpath(self.proj_dir + "/final_outputs/model/ROCkOut_Filter.txt")
			#print(filter_path)
			#print(os.path.exists(filter_path))
			if os.path.exists(filter_path):
				#Only set this to not None if the file can be found.
				self.filter_file = filter_path
			else:
				print("Project appears incomplete.")
		else:
			print("Project not found.")
		
	def find_ma(self):
		if os.path.exists(self.proj_dir):
			ma_path = os.path.normpath(self.proj_dir + "/final_outputs/model/complete_multiple_alignment_aa.fasta" )
			#print(ma_path)
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
		
	def plot_filter(self):
		as_pandas = pd.DataFrame(self.filter_matrix)
		fig = go.Figure(data=[go.Surface(z=as_pandas.values)])
		fig.update_layout(scene = dict(
							xaxis_title='Position in Protein',	
							yaxis_title='Read Length',
							zaxis_title='Bitscore Cutoff'))
		fig.write_html(os.path.normpath("ROCkOut_surface.html"))

	#This now needs to take the multiple alignment info.	
	def parse_filter_dir(self):
		if not os.path.exists(self.passing):
			os.mkdir(self.passing)
		if not os.path.exists(self.failing):
			os.mkdir(self.failing)
		if not os.path.exists(self.fastas):
			os.mkdir(self.fastas)
			
		self.raw_reads = []
		self.reads_to_filter = []
		self.filtered_passing = []
		self.passing_alns = []
		self.failing_alns = []
		
		raws = os.listdir(self.orig)
		raws.sort()
		alignments = os.listdir(self.alns)
		alignments.sort()
		
		if not len(raws) == len(alignments):
			print("Mismatch in length of reads to filter and their original FASTA counterparts. Exiting.")
		else:
			for raw, aln in zip(raws, alignments):
				basename = aln.replace("_ROCkOut_alignments.txt", "")
				alnpath = os.path.normpath(self.alns + "/" + aln)
				rawpath = os.path.normpath(self.orig + "/" + raw)
				filt_raw = os.path.normpath(self.fastas + "/" + basename + "_filtered.fasta")
				filt_pass = os.path.normpath(self.passing + "/" + basename + "ROCkOut_passing.txt")
				filt_fail = os.path.normpath(self.failing + "/" + basename + "ROCkOut_failing.txt")
				
				self.raw_reads.append(rawpath)
				self.reads_to_filter.append(alnpath)
				self.passing_alns.append(filt_pass)
				self.failing_alns.append(filt_fail)
				self.filtered_passing.append(filt_raw)
				
	def run_filter_blast_tabular(self):

		for raw, aln, passing, failing, raw_filt in zip(self.raw_reads, self.reads_to_filter, self.passing_alns, self.failing_alns, self.filtered_passing):
			print("Filtering", raw)
			
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
					
					#readlength = end-start #This won't work for this setup
					
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
			
	def run_filterer(self):
		print("Loading filter resources...")
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
	
		#self.plot_filter()
		self.parse_filter_dir()
		
		print("Filtering reads...")
		self.run_filter_blast_tabular()
	
def do_filter(parser, opts):
	project_dir = opts.dir
	filter_dir = opts.filter_dir
	
	if project_dir is None or filter_dir is None:
		print("ROCkOut needs both a project directory and a filter directory to filter reads.")
		parser.print_help()
		sys.exit()
	else:
		mn = rocker_filterer(project_dir, filter_dir)
		mn.run_filterer()
		