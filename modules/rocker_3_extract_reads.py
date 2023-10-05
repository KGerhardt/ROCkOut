import sys
import os 

import subprocess

try:
	from .rocker_project_manager import project_manager
	from .reads_labeller import protein_trawler
except:
	from rocker_project_manager import project_manager
	from reads_labeller import protein_trawler
	
class read_extractor:
	def __init__(self, directory, aln_prefix, raw_prefix):
		self.base = os.path.normpath(directory)
		self.labeller = protein_trawler(self.base, splits = 5, train_fraction = 0.6)
		self.manager = project_manager(self.base)
		self.aln_pre = aln_prefix
		self.raw_pre = raw_prefix
		
		
	def parse_project_directory(self):
		#Find misc. files including multiple alignments
		self.manager.parse_project_directory()
		self.manager.parse_aligns()
		self.manager.parse_targets()
		self.manager.parse_multiple_alignment()
		self.ma_aa_file = self.manager.mult_aln_files['aln_aa']
		
	def load_and_output_reads(self):
		self.labeller.extract_prep()
		
		for read_length in self.labeller.read_files:
			if self.aln_pre is not None:
				aln_dat = [self.aln_pre, str(read_length), "alignments.txt"]
				out_aln = open(os.path.normpath('_'.join(aln_dat)), "w")
			if self.raw_pre is not None:
				raw_dat = [self.raw_pre, str(read_length), "raw_reads.fasta"]
				out_raw = open(os.path.normpath('_'.join(raw_dat)), "w")
			
			for genome in self.labeller.read_files[read_length]:
				alns, raws = self.labeller.collect_reads_and_raws(read_length, genome)
				if self.aln_pre is not None:
					alns.to_csv(path_or_buf = out_aln, sep = "\t", index = False, header = False)
				if self.raw_pre is not None:
					for id in raws:
						print(id, file = out_raw)
						print(raws[id], file = out_raw)
			
			if self.aln_pre is not None:
				out_aln.close()
			if self.raw_pre is not None:
				out_raw.close()
			
	def run(self):
		#self.parse_project_directory()
		self.load_and_output_reads()
			
			
def extract_reads(parser, opts):
	#multiprocessing.freeze_support()

	dirname = opts.dir
	if dirname is None:
		sys.exit("ROCkOut needs a directory name to build a model!")
	
	aln_prefix = opts.aln_base
	raw_prefix = opts.raw_base
	
	if aln_prefix is None and raw_prefix is None:
		print("Choose at least one filename prefix to supply using either or both of:")
		print("--alignments_prefix")
		print("--raw_reads_prefix")
		sys.exit()
				
	mn = read_extractor(directory = dirname, 
						aln_prefix = aln_prefix,
						raw_prefix = raw_prefix)
	mn.run()

