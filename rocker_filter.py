import sys
import argparse
import os
import numpy as np

import plotly.graph_objects as go
import pandas as pd

class rocker_filterer:
	def __init__(self, filter_file):
		self.filter_file = filter_file
		self.filter_model = None
		self.observed_readlens = None
		self.filter_matrix = None
		self.get_filter()
		self.interpolate()
		
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
			
	def plot_filter(self):
		as_pandas = pd.DataFrame(self.filter_matrix)
		fig = go.Figure(data=[go.Surface(z=as_pandas.values)])
		fig.update_layout(scene = dict(
							xaxis_title='Position in Protein',	
							yaxis_title='Read Length',
							zaxis_title='Bitscore Cutoff'))
		fig.write_html(os.path.normpath("ROCkOut_surface.html"))

		
	def run_filter_blast_tabular(self, reads_file, output = None):
		fh = open(reads_file)
		for line in fh:
			segs = line.strip().split("\t")
			alignment_range = [int(segs[8]), int(segs[9])]
			start = min(alignment_range)
			end = max(alignment_range)
			bitscore = float(segs[11])
			readlength = end-start
			
			diffs = self.filter_matrix[readlength, start:end]
			diffs = np.subtract(bitscore, diffs)
			
			decider = np.mean(diffs)
			if decider > 0:
				print(*segs)
			
		fh.close()
			

inf = sys.argv[1]
mn = rocker_filterer(inf)
#mn.plot_filter()
mn.run_filter_blast_tabular(reads_file = "eee/positive/Q00014/aligned_reads/M64090_read_length_100_aligned_reads.txt")
#print(mn.filter_model)