import sys
import argparse
import os
import numpy as np

import plotly.express as px
import plotly.graph_objects as go
import pandas as pd

#Suppress warning for making a new column in a dataframe
pd.options.mode.chained_assignment = None
#from modules.rocker_project_manager import project_manager

#Change check sim to false for final build.
class rocker_filterer:
	def __init__(self, project_directory, filter_dir = "", reads = [], basenames = [], raws = [], write_outputs = True, check_sim = True):
		self.proj_dir = project_directory
		self.filt_dir = filter_dir
		self.model_dir = os.path.normpath(self.proj_dir + "/final_outputs/model/")
		self.reads = reads
		self.basenames = basenames
		self.raws = raws
		
		self.write_outputs = write_outputs
		self.sim_check = check_sim
		
		self.alns = self.filt_dir + "/alignments"
		self.orig = self.filt_dir + "/original_reads"
		
		self.passing = os.path.normpath(self.filt_dir + "/ROCkOut_passing_alignments")
		self.failing = os.path.normpath(self.filt_dir + "/ROCkOut_failing_alignments")
		self.fastas = os.path.normpath(self.filt_dir + "/ROCkOut_passing_read_fastas")
		self.failing_fastas = os.path.normpath(self.filt_dir + "/ROCkOut_failing_read_fastas")
		self.viz = os.path.normpath(self.filt_dir + "/visualizations")
		
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
		if not os.path.exists(self.passing):
			os.mkdir(self.passing)
		if not os.path.exists(self.failing):
			os.mkdir(self.failing)
		if not os.path.exists(self.fastas):
			os.mkdir(self.fastas)
		if not os.path.exists(self.failing_fastas):
			os.mkdir(self.failing_fastas)
			
		if not os.path.exists(self.viz):
			os.mkdir(self.viz)
			
	def find_filters(self):
		if os.path.exists(self.proj_dir):
			if os.path.exists(self.model_dir):
				models = os.listdir(self.model_dir)
				models = [os.path.normpath(self.model_dir + "/" +m) for m in models]
				for m in models:
					if m.endswith("pct_id_vs_MA_pos.txt"):
						self.id_filter_file = m
					if m.endswith("bitscore_vs_MA_pos.txt"):
						self.filter_file = m
					if m.endswith("pct_id_vs_pct_aln.txt"):
						self.idaln_file = m
							
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
		self.filter_matrix, self.bit_positions = self.import_filter(self.filter_file)
		self.idpos_filtmat, self.id_positions =self.import_filter(self.id_filter_file)
		self.idaln_filtmat, self.aln_positions = self.import_filter(self.idaln_file)
		
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

	def load_reads(self, reads_file):
		df = pd.read_csv(reads_file, sep = "\t", header=None, 
						usecols = [0, 1, 2, 3, 6, 7, 8, 9, 11, 12, 13])
						
		colnames = ["read", "target", "pctid", "alnlen", "qst", "qnd", "sst", "snd", "bs", "qlen", "slen"]
		df.columns = colnames

		
		valid_offsets = set(list(self.offsets.keys()))
				
		df = df.loc[df["target"].isin(valid_offsets)]
		df = df.reset_index(drop = True)
		

		#Calculate percent alignment for use in filtering
		df["pct_aln"] = np.round(100 * df["alnlen"] / (df["qlen"]/3), 2)

		return df
	
	def besthit_reads(self, df):
		#Collect max bitscore by read
		idx = df.groupby(['read'])['bs'].transform("max") == df['bs']
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
		
		return df
		
	def label_sim_reads(self, read_names, assignments):
		try:
			outs = {"tp":0, "fp":0, "fn":0, "tn":0}

			for read, guessed_positive in zip(read_names, assignments):		
				segs = read.split(";")
				label = segs[-1]
				#print(label, guessed_positive)
				
				if guessed_positive:
					if label == "Positive":
						outs["tp"] += 1
					else:
						outs["fp"] += 1
				else:
					if label == "Positive":
						outs["fn"] += 1
					else:
						outs["tn"] += 1
				
				
				total_size = outs['tn']+outs['tp']+outs["fp"]+outs["fn"]
				fpr = outs["fp"] / (outs["fp"]+outs["tn"])
				fnr = outs["fn"] / (outs["fn"] + outs["tp"])
				
				print("Read count:", total_size)
				print(outs)
				print("Acc.", "%.4f" % ((outs['tn']+outs['tp'])/total_size), "%")
				print("FPR", "%.4f" % (fpr*100), "%")
				print("FNR", "%.4f" % (fnr*100), "%")
				print("F1 ", "%.4f" % (2*outs['tp'] / (2*outs['tp'] + outs['fp']+outs['fn'])))
		except:
			pass

		return None
		
	def filter_raws(self, selection, raw_file, output_file):
		out = open(output_file, "w")
		next_defline = ""
		next_seqid = ""
		next_seq = []
		with open(raw_file) as fh:
			for line in fh:
				#print(line)
				if line.startswith(">"):
					if len(next_seq) > 0:
						#print(next_seqid)
						if next_seqid in selection:
							next_seq = "".join(next_seq)
							out.write(next_defline)
							out.write(next_seq)
							
					next_seq = []
					next_defline = line
					next_seqid = line.strip()
					next_seqid = next_seqid[1:]
					next_seqid = next_seqid.split()[0]
					
				else:
					next_seq.append(line)
		
		#Final iteration
		if len(next_seq) > 0:
			if next_seqid in selection:
				next_seq = "".join(next_seq)
				out.write(next_defline)
				out.write(next_seq)
		
		out.close()
		
	def filter_reads(self, read_df, name_prefix):
		percent_alignment_matches = np.searchsorted(self.aln_positions, read_df["pct_aln"], side = 'left')
		percent_alignment_matches -= 1
		multiple_alignment_positions = []
		
		for tgt, aln1, aln2 in zip(read_df["target"], read_df["sst"], read_df["snd"]):
			read_mid = int(((aln1+aln2-1)/2))
			ma_midpoint = self.offsets[tgt][read_mid] + read_mid
			multiple_alignment_positions.append(ma_midpoint)
			
		multiple_alignment_positions = np.array(multiple_alignment_positions, dtype = np.int32)
		
		passing_indices = []
		passes_bitscore = []
		passes_id = []
		passes_idaln = []
		
		#print("idaln", self.idaln_file)
		#print("bitscore", self.filter_file)
		#print("pct", self.id_filter_file)
		
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
			
			#print(pct_aln, id_aln_cutoff)
			
			passbs = (bitscore >= bitscore_cutoff)
			passid = (pct_id >= id_pos_cutoff)
			passidaln = (pct_id >= id_aln_cutoff)
				
			passes_bitscore.append(passbs)
			passes_id.append(passid)
			passes_idaln.append(passidaln)
				
			vote = sum([passbs, passid, passidaln])
			passing_indices.append(vote)
				
			index += 1
		
		vote_count = np.array(passing_indices, dtype = np.int32)
		read_df["passes_bitscore_v_pos"] = passes_bitscore
		read_df["passes_pct_id_v_pos"] = passes_id
		read_df["passes_pct_id_v_pct_aln"] = passes_idaln

		read_df["passed_filters_out_of_3"] = vote_count
		read_df["passes_ensemble"] = (vote_count >= 2)
		
		#Viz stuff
		if True:
			average_rl = int(np.mean(read_df['qlen']))
			average_rl_index = average_rl - self.min_readlen
			nearest_bitpos =  self.filter_matrix[average_rl_index, :]
			nearest_idpos = self.idpos_filtmat[average_rl_index, :]
			nearest_idaln = self.idaln_filtmat[average_rl_index, :]
			read_df['ma_pos'] = multiple_alignment_positions
			
			
			cat = []
			for guessed_positive in read_df["passes_ensemble"]:
				
				if guessed_positive:
					cat.append("Positive")
				else:
					cat.append("Negative")
			
			read_df['label'] = cat 
			
			
			self.plot_results(read_df,
								nearest_bitpos, 
								nearest_idpos, 
								nearest_idaln,
								name_prefix)
		
		passing_reads = read_df[read_df["passed_filters_out_of_3"] >= 2]
		failing_reads = read_df[read_df["passed_filters_out_of_3"] < 2]
		
		#if self.sim_check:
		if False:
			print("Bitscore vs. position in MA")
			self.label_sim_reads(read_df["read"], read_df["passes_bitscore_v_pos"])

			print("Pct. ID vs. position in MA")
			self.label_sim_reads(read_df["read"], read_df["passes_pct_id_v_pos"])

			print("Pct. ID vs. Pct. aln")
			self.label_sim_reads(read_df["read"], read_df["passes_pct_id_v_pct_aln"])

			print("Ensemble model vote")
			self.label_sim_reads(read_df["read"], read_df["passes_ensemble"])
			
		return passing_reads, failing_reads

	def plot_results(self, df, bp, ip, ia, prefix):
		bsl = px.line(x = self.bit_positions, y = bp)
		pil = px.line(x = self.id_positions, y = ip)
		ial = px.line(x = self.aln_positions, y = ia)

		bs = px.scatter(df, x = 'ma_pos', y = 'bs', color = "label", hover_data=["read"])
		bs = go.Figure(data=bs.data + bsl.data)
		bs.write_html(os.path.normpath(self.viz + "/" + prefix + "_bitscore_filter_plot.html"))
		
		id = px.scatter(df, x = 'ma_pos', y = 'pctid', color = "label", hover_data=["read"])
		id = go.Figure(data=id.data + pil.data)
		id.write_html(os.path.normpath(self.viz + "/" + prefix + "_pct_id_filter_plot.html"))
		
		idaln = px.scatter(df, x = 'pct_aln', y = 'pctid', color = "label", hover_data=["read"])
		idaln = go.Figure(data=idaln.data + ial.data)
		idaln.write_html(os.path.normpath(self.viz + "/" + prefix + "_id_aln_filter_plot.html"))
		
		
	
	def output_files(self, df, alignments_out, raws_in, raws_out):
		df.to_csv(alignments_out, sep = "\t", index = False, header = None)
		acceptable_names = set(df["read"])
		self.filter_raws(acceptable_names, raws_in, raws_out)
		return None
	
	def load_resources(self):
		print("Loading filter resources...")
		self.dir_prep()
		self.find_ma()
		self.get_offsets()
		self.find_filters()
		self.load_filters()
		
		if self.filter_file is None:
			print("ROCkOut filter file not found! Quitting")
			sys.exit()
		if self.ma_file is None:
			print("ROCkOut multiple alignment file not found! Quitting")
			sys.exit()
			
	def run_filterer(self, r, b, raw):
		loaded_reads = self.load_reads(r)
		loaded_reads = self.besthit_reads(loaded_reads)
		p, f = self.filter_reads(loaded_reads, get_bn(r))

		outp = os.path.normpath(self.passing + "/" +b+".ROCkOut_passing.txt")
		outf = os.path.normpath(self.failing+ "/" +b+".ROCkOut_failing.txt")
		pass_raw = os.path.normpath(self.fastas+"/"+b+".filtered.fasta")
		fail_raw = os.path.normpath(self.failing_fastas+"/"+b+".filtered.fasta")
		
		#print(p)
		#self.label_sim_reads(p[""], p[""])
		
		self.output_files(p, outp, raw, pass_raw)
		self.output_files(f, outf, raw, fail_raw)			
			
		return p, f
		
	def run_filterer_internal(self, loaded_reads):
		loaded_reads.columns = ["read", "target", "pctid", "alnlen", "qst", "qnd", "sst", "snd", "bs", "qlen", "slen"]
		#Calculate percent alignment for use in filtering
		loaded_reads["pct_aln"] = np.round(100 * loaded_reads["alnlen"] / (loaded_reads["qlen"]/3), 2)
		loaded_reads = self.besthit_reads(loaded_reads)
		
		print(loaded_reads)
			
def get_bn(file):
	b = os.path.basename(file)
	while b != os.path.splitext(b)[0]:
		b = os.path.splitext(b)[0]
	return b
	
def do_filter(parser, opts):
	project_dir = opts.dir
	filter_dir = opts.filter_dir
	
	if filter_dir is None or project_dir is None:
		print("Can't proceed without a project directory and filter directory.")
		parser.print_help()
		sys.exit()
		
	reads = []
	basenames = []
	raws = []
	for f in os.listdir(os.path.normpath(filter_dir+"/alignments")):
		reads.append(os.path.normpath(filter_dir+"/alignments/"+f))
		basenames.append(get_bn(f))
	for f in os.listdir(os.path.normpath(filter_dir+"/original_reads")):
		raws.append(os.path.normpath(filter_dir+"/original_reads/"+f))
		
	reads.sort()
	basenames.sort()
	raws.sort()
		
	if project_dir is None or filter_dir is None:
		print("ROCkOut needs both a project directory and a filter directory to filter reads.")
		parser.print_help()
		sys.exit()
	else:
		mn = rocker_filterer(project_dir, filter_dir, reads, basenames, raws, check_sim = False)
		mn.load_resources()
		print("Filtering reads")

		for r, b, raw in zip(reads, basenames, raws):
			passing_reads, failing_reads = mn.run_filterer(r, b, raw)
			
def internal_filter(bit, bitx, id, idx, idaln, idalnx, reads, filter_dir):
	mn = rocker_filterer(project_dir = None, filter_dir = filter_dir, reads = [], basenames = [], raws = [], write_outputs = False, check_sim = False)
	#we need to pass resources instead
	#mn.load_resources()
	mn.filter_matrix = mn.interpolate_matrix(bit)
	mn.idpos_filtmat = mn.interpolate_matrix(id)
	mn.idaln_filtmat = mn.interpolate_matrix(idaln)
	mn.bit_positions = np.array(bitx, dtype = np.float_)
	mn.id_positions = np.array(idx, dtype = np.float_)
	mn.aln_positions = np.array(idalnx, dtype = np.float_)
		
	#allx = np.array(allx, dtype = np.float_)
	
	for read in reads:
		passing_reads, failing_reads = mn.run_filterer_internal(read)
		yield passing_reads, failing_reads
		
	