import subprocess
import sys
import os
import tempfile

import argparse

from rocker_project_manager import project_manager

class MA_manager:
	def __init__(self, dir):
		self.in_dir = dir
		self.roc_proj = None
		self.target_paths = None
		self.combined_seqs = None
		self.ma = None
		self.tree = None
		
	def access_files(self):
		self.roc_proj = project_manager(directory = self.in_dir)
		self.roc_proj.parse_project_directory()
		self.roc_proj.parse_targets()
		self.target_paths = self.roc_proj.targets
		
	def merge_seqs(self):
		if not os.path.exists(os.path.normpath(self.in_dir + "/shared_files")):
			os.mkdir(os.path.normpath(self.in_dir + "/shared_files"))
		if not os.path.exists(os.path.normpath(self.in_dir + "/shared_files/multiple_alignment")):
			os.mkdir(os.path.normpath(self.in_dir + "/shared_files/multiple_alignment"))
		self.combined_seqs = open(os.path.normpath(self.in_dir + "/shared_files/multiple_alignment/combined_target_seqs.txt"), "w")
		
		for prot in self.target_paths:
			for file in self.target_paths[prot]:
				fh = open(file)
				for line in fh:
					if line.startswith(">"):
						line = line.strip()[1:]
						#Give ID clearly for later.
						line = ">"+prot+"\t"+line+"\n"
					#Write either seqline or updated defline.
					self.combined_seqs.write(line)
				fh.close()
		
		self.combined_seqs.close()
		self.combined_seqs = os.path.normpath(self.in_dir + "/shared_files/multiple_alignment/combined_target_seqs.txt")
		
	def create_ma(self):
		self.ma = os.path.normpath(self.in_dir + "/shared_files/multiple_alignment/target_seq_MA.txt")
		muscle_ma_command = ["muscle", "-in", self.combined_seqs, "-out", self.ma]
		subprocess.run(muscle_ma_command)
		
	def create_tree(self):
		self.tree = os.path.normpath(self.in_dir + "/shared_files/multiple_alignment/target_seq_tree.txt")
		muscle_tree_command = ["muscle", "-maketree", "-in", self.ma, "-out", self.tree]
		subprocess.run(muscle_tree_command)

directory = sys.argv[1]

manager = MA_manager(directory)
manager.access_files()
manager.merge_seqs()
manager.create_ma()
manager.create_tree()
