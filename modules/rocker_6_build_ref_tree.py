import sys
import os
import subprocess

from .rocker_project_manager import project_manager
from .rocker_progress_tracker import progress_tracker

class pplacer_tree_builder:
	def __init__(self, dir):
		self.projdir = dir
		self.basename = os.path.basename(self.projdir)
		
		self.proj_man = None
		
		self.locus = "locus_placeholder"
		self.refpkg = os.path.normpath(self.projdir+"/shared_files/pplace/"+self.basename+".refpkg")
		self.mult_aln = None
		self.fasttree_log = None
		self.fasttree_tree = None
		
		self.jplace = None
		
	def prep(self):
		pplace_path = os.path.normpath(self.projdir+"/shared_files/pplace/")
		if not os.path.exists(pplace_path):
			os.mkdir(pplace_path)
		
		
		
	def load_dir(self):
		self.proj_man = project_manager(directory = self.projdir)
		self.proj_man.parse_project_directory()
		self.proj_man.parse_tagged_reads()
		self.proj_man.parse_mult_aln()
		
		if self.proj_man.mult_aln_files is not None:
			self.mult_aln = self.proj_man.mult_aln_files["aln"]
			self.fasttree_log = self.proj_man.mult_aln_files["log"]
			self.fasttree_tree = self.proj_man.mult_aln_files["tree"]
		else:
			print("Didn't find any multiple alignment files!")
		
		
	def build_refpkg(self):
		taxit_arg = "taxit create -l {locus} -P {packname} --aln-fasta {aln} --tree-stats {log} --tree-file {tree}"
		taxit_arg = taxit_arg.format(locus=self.locus, packname=self.refpkg, aln=self.mult_aln, log=self.fasttree_log, tree=self.fasttree_tree)
		taxit_arg = taxit_arg.split()
		
		if not os.path.exists(self.refpkg):
			subprocess.call(taxit_arg)
		
		
	def prep_pplace(self):
		self.jplace = os.path.normpath(self.refpkg + "/" + os.path.basename(self.mult_aln)[:-5] + "jplace")
		pplace_arg = "pplacer -c {refpak} {mult_aln} --no-pre-mask -o {outjp}"
		
		#So what we need to do is to cut the reads to their aligned range on the protein and then make a MA from the cut results.
		
		for base in self.proj_man.tagged_reads_pos:
			pplace_arg = pplace_arg.format(refpak = self.refpkg, mult_aln = self.mult_aln, outjp = self.jplace)
			
			print(pplace_arg)
		
		pplace_arg = pplace_arg.split()
		#subprocess.call(pplace_arg)
		
		
def make_reftree(parser, opts):
	dir = opts.dir
	mn = pplacer_tree_builder(dir = dir)
	mn.prep()
	
	mn.load_dir()
	mn.build_refpkg()
	mn.prep_pplace()
