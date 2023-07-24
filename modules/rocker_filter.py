import sys
import argparse
import os
import numpy as np

import plotly.graph_objects as go
import pandas as pd

#Suppress warning for making a new column in a dataframe
pd.options.mode.chained_assignment = None
#from modules.rocker_project_manager import project_manager

class rocker_filterer:
	def __init__(self, project_directory, filter_dir = ""):
		self.proj_dir = project_directory
		self.filt_dir = filter_dir
		self.model_dir = os.path.normpath(self.proj_dir + "/final_outputs/model/")

		
		self.alns = self.filt_dir + "/alignments"
		self.orig = self.filt_dir + "/original_reads"
		
		self.passing = self.filt_dir + "/ROCkOut_passing_alignments"
		self.failing = self.filt_dir + "/ROCkOut_failing_alignments"
		self.fastas = self.filt_dir + "/ROCkOut_passing_read_fastas"
		
		#self.plot = option
		self.filter_file = None
		self.id_filter_file = None
		self.idaln_file = None
		
		self.bit_positions = None
		self.id_positions = None
		self.aln_positions = None
		
		self.min_readlen = None
		self.max_readlen = None
		
		self.ma_file = None
		self.offsets = None
		self.filter_model = None
		
		self.observed_readlens = None
		
		self.filter_matrix = None
		self.idpos_filtmat = None
		self.idaln_filtmat = None
		
		self.raw_reads = None
		self.reads_to_filter = None
		self.filtered_passing = None
		self.passing_alns = None
		self.failing_alns = None
	
	def dir_prep(self):
		out_base_p = os.path.normpath(self.filt_dir + "/ROCkOut_passing_alignments")
		out_base_f = os.path.normpath(self.filt_dir + "/ROCkOut_failing_alignments")
		if not os.path.exists(out_base_p):
			os.mkdir(out_base_p)
		if not os.path.exists(out_base_f):
			os.mkdir(out_base_f)
		
	def find_filters(self):
		if os.path.exists(self.proj_dir):
			if os.path.exists(self.model_dir):
				asps = os.listdir(self.model_dir)
				asps = [os.path.normpath(self.model_dir + "/" +a) for a in asps if "accuracy_sensitivity_and_specificity" in a]
				for a in asps:
					if a.endswith("_id_aln.txt"):
						self.idaln_file = a
					else:
						if a.endswith("_id.txt"):
							self.id_filter_file = a
						else:
							self.filter_file = a
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
		
	def load_filters(self):
		self.filter_matrix, self.bit_positions = self.import_filter(mn.filter_file)
		self.idpos_filtmat, self.id_positions =self.import_filter(mn.id_filter_file)
		self.idaln_filtmat, self.aln_positions = self.import_filter(mn.idaln_file)
		
	def import_filter(self, file):
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
		
		#print(allx)
		
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
		
		matrix = np.zeros(shape = (z_shape, x_shape), dtype = np.float_)
		for rl in per_rl_x:
			#z = readlen - min_rl
			zi = rl - min_rl
			for x, y in zip(per_rl_x[rl], per_rl_y[rl]):
				xi = x_to_index[x]
				matrix[zi, xi] = y
				
		matrix = self.interpolate_matrix(matrix)
		allx = np.array(allx, dtype = np.float_)
		
		return matrix, allx
		
	def interpolate_matrix(self, cutoff_matrix):
		rowsums = np.sum(cutoff_matrix, axis = 1)
		calculated_rows = np.sort(np.where(rowsums > 0)[0])

		for i in range(0, len(calculated_rows)-1):
			start = calculated_rows[i]
			end = calculated_rows[i+1]
			
			#print(end-start)
			step_sizes = (cutoff_matrix[end] - cutoff_matrix[start]) / (end-start)
			
			for i in range(1, (end-start)):
				cutoff_matrix[start + i] = cutoff_matrix[start] + (step_sizes * i)
			
		return cutoff_matrix
	
	'''
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
	'''
	
	def plot_filter(self):
		as_pandas = pd.DataFrame(self.filter_matrix)
		fig = go.Figure(data=[go.Surface(z=as_pandas.values)])
		fig.update_layout(scene = dict(
							xaxis_title='Position in Protein',	
							yaxis_title='Read Length',
							zaxis_title='Bitscore Cutoff'))
		fig.write_html(os.path.normpath("ROCkOut_surface.html"))

		
	def load_reads_and_besthit(self, reads_file):
		df = pd.read_csv(reads_file, sep = "\t", header=None, 
						usecols = [0, 1, 2, 3, 6, 7, 8, 9, 11, 12, 13])
		df.columns = ["read", "target", "pctid", "alnlen", "qst", "qnd", "sst", "snd", "bs", "qlen", "slen"]

		#Collect max bitscore by read
		idx = df.groupby(['read'])['bs'].transform(max) == df['bs']
		df = df[idx]
		df = df.reset_index(drop = True)
			
		#collect a random read among cases where max bitscore is tied
		idx = df.groupby(["read"])['sst'] #random integer column - it's effectively discarded
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
		
		#Calculate percent alignment for use in filtering
		df["pct_aln"] = np.round(100 * df["alnlen"] / (df["qlen"]/3), 2)

		return df
		
	def filter_reads(self, read_df):
		percent_alignment_matches = np.searchsorted(self.aln_positions, read_df["pct_aln"], side = 'left')
		percent_alignment_matches -= 1
		multiple_alignment_positions = []
		
		for tgt, aln1, aln2 in zip(read_df["target"], read_df["sst"], read_df["snd"]):
			read_mid = int(((aln1+aln2-1)/2))
			ma_midpoint = self.offsets[tgt][read_mid] + read_mid
			multiple_alignment_positions.append(ma_midpoint)
			
		multiple_alignment_positions = np.array(multiple_alignment_positions, dtype = np.int32)
		
		passing_indices = []
		
		index = 0
		for read, ma_pos, pct_aln_index, pct_id, bitscore, readlen in zip(read_df["read"],
																	multiple_alignment_positions, 
																	percent_alignment_matches,
																	read_df["pctid"],
																	read_df["bs"],
																	read_df["qlen"]):
			#We take the best guess we can 
			if readlen < self.min_readlen:
				readlen = self.min_readlen
			if readlen > self.max_readlen:
				readlen = self.max_readlen
				
			readlen_index = readlen - self.min_readlen
			
			bitscore_cutoff = self.filter_matrix[readlen_index, ma_pos]
			id_pos_cutoff = self.idpos_filtmat[readlen_index, ma_pos]
			id_aln_cutoff = self.idaln_filtmat[readlen_index, pct_aln_index]
			
			passbs = (bitscore >= bitscore_cutoff)
			passid = (pct_id >= id_pos_cutoff)
			passidaln = (pct_id >= id_aln_cutoff)
				
			vote = sum([passbs, passid, passidaln])
			passing_indices.append(vote)
				
			index += 1
		
		vote_count = np.array(passing_indices, dtype = np.int32)
		read_df["passed_filters_out_of_3"] = vote_count

		passing_reads = read_df[read_df["passed_filters_out_of_3"] >= 2]
		failing_reads = read_df[read_df["passed_filters_out_of_3"] < 2]
		
		if True:
			outs = {"tp":0, "fp":0, "fn":0, "tn":0}
			for read, vote in zip(read_df["read"], vote_count):
				segs = read.split(";")
				label = segs[-1]
				
				if vote >= 2:
					if label == "Positive":
						outs["tp"] += 1
					else:
						outs["fp"] += 1
				else:
					if label == "Positive":
						outs["fn"] += 1
					else:
						outs["tn"] += 1
			
			
			fpr = outs["fp"] / (outs["fp"]+outs["tn"])
			fnr = outs["fn"] / (outs["fn"] + outs["tp"])
			print(outs)
			print("FPR", fpr)
			print("FNR", fnr)
		
			
		
		return passing_reads, failing_reads


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
		
proj = sys.argv[1]
reads = sys.argv[2]
fd = sys.argv[3]
	
mn = rocker_filterer(proj, filter_dir = fd)
mn.dir_prep()
mn.find_ma()
mn.get_offsets()
mn.find_filters()
mn.load_filters()

loaded_reads = mn.load_reads_and_besthit(reads)
p, f = mn.filter_reads(loaded_reads)

outp = os.path.normpath(fd+"/ROCkOut_passing_alignments/passing_reads.blast.txt")
outf = os.path.normpath(fd+"/ROCkOut_failing_alignments/failing_reads.blast.txt")
p.to_csv(outp, sep = "\t", index = False, header = None)
f.to_csv(outf, sep = "\t", index = False, header = None)
