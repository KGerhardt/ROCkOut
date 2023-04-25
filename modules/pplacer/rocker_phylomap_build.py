import os
import sys

class phyl_builder:
	def __init__(self, genomes, outpath):
		self.gens = genomes
		self.outpath = outpath
		
		self.is_nt = True

		self.make_dir("phylogenetic_placement")
		self.make_dir("phylogenetic_placement/source")
		
		self.combined_genomes = os.path.normpath(outpath + "/phylogenetic_placement/source/combined_genomes.fasta")
		self.cg_writer = None
		
		self.mult_aln = os.path.normpath(outpath + "/phylogenetic_placement/source/combined_genomes_multiple_alignment.fasta")
		#self.mult_aln = ma
		self.fasttree_log = os.path.normpath(outpath + "/phylogenetic_placement/source/combined_genomes_fasttree_log.txt")
		#self.fasttree_log = tree_log
		self.fasttree_tree = os.path.normpath(outpath + "/phylogenetic_placement/source/combined_genomes_fasttree_tree.txt")
		#self.fasttree_tree = tree
		
		self.make_dir("phylogenetic_placement/reference_package")
		
		self.refpkg = os.path.normpath(outpath + "/phylogenetic_placement/reference_package/phylomap_ref.pkg")
		if os.path.exists(self.refpkg):
			from shutil import rmtree
			rmtree(self.refpkg)
			
			
		self.locus =  "placeholder"
		
	def make_dir(self, dir):
		to_make = os.path.normpath(self.outpath + "/" + dir)
		if not os.path.exists(to_make):
			os.makedirs(to_make, exist_ok = True)
		
	def combine_fasta(self, file):
		try:
			fh = open(file)
			for line in fh:
				self.cg_writer.write(line)
			fh.close()
		except:
			print("Couldn't read file", file)
			print("Skipping.")
			
	def make_ma(self):
		print("Creating positive sequences multiple alignment")
		#This might be a version thing
		muscle_ma_command = ["muscle", "-in", self.combined_genomes, "-out", self.mult_aln]
		muscle_ma_command = " ".join(muscle_ma_command)
		#print(muscle_ma_command)
		os.system(muscle_ma_command)
		
	def make_tree(self):
		print("Creating positive sequences phylogenetic tree")

		fasttree_command = ["fasttree", "-nt", "-gtr", "-log", self.fasttree_log, "-out", self.fasttree_tree, self.mult_aln]
		fasttree_command = " ".join(fasttree_command)
		os.system(fasttree_command)
			
	def make_refpak(self):
		print("Finalizing positive sequence phylogenetic placement information")
		
		if self.is_nt:
			taxit_arg = "taxit create -l {locus} -P {packname} --aln-fasta {aln} --tree-stats {log} --tree-file {tree}"
			taxit_arg = taxit_arg.format(locus=self.locus, packname=self.refpkg, aln=self.mult_aln, log=self.fasttree_log, tree=self.fasttree_tree)
			os.system(taxit_arg)
		else:
			print("Make proteins work?")
		
	def run(self):
		
		self.cg_writer = open(self.combined_genomes, "w")
		for g in self.gens:
			self.combine_fasta(g)
		self.cg_writer.close()
		
		self.make_ma()
		
		self.make_tree()
		self.make_refpak()		

def phylomap_build(genomes = None, output = None, pakname = "rocker_reference_seqs", quiet = True):			
	mn = phyl_builder(genomes = genomes,
					outpath = output)		
	mn.run()
		