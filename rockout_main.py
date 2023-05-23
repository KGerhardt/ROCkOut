import argparse
import sys
import subprocess

def options(action):
	parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
			description='''	''')

	#This is universal and I want it first
	parser.add_argument('-d', '--directory',  dest = 'dir', default = None, help =  'A rocker project directory that you wish to create or modify')
	
	if action == "download":
		parser.description = "ROCkOut download"
		parser.add_argument('-p', '--positive',  dest = 'pos', default = None, help =  'File containing positive uniprot IDs, 1 per line')
		parser.add_argument('-n', '--negative',  dest = 'neg', default = None, help =  'File containing negative uniprot IDs, 1 per line')
		
	if action == "add":
		parser.add_argument('--id',  dest = 'id', default = None, help =  'Uniprot ID to add')
		parser.add_argument('--type',  dest = 'pos_neg', default = None, help =  "Group for this ID. One of 'positive', 'negative'")
	
	if action == "build":
		parser.description = "ROCkOut build. Read lengths intend to hit a mean readlength +/- 10%"
		'''
		parser.add_argument('--short',  dest = 'short', default = "90,100,110", 
		help =  'Comma-sep triplet of very short read lengths, default: 90,100,110')
		parser.add_argument('--med',  dest = 'med', default = "180,200,220", 
		help =  'Comma-sep triplet of moderate read lengths, default: 180,200,220')
		parser.add_argument('--long',  dest = 'long', default = "270,300,330", 
		help =  'Comma-sep triplet of long read lengths, default: 270,300,330')
		parser.add_argument('--xl',  dest = 'xl', default = "360,400,440", 
		help =  'Comma-sep triplet of very long read lengths, default: 360,400,440')
		'''
		parser.add_argument('--short-lower',  dest = 'sl', type = int, default = 90,
		help =  'Lower-bound on simulated read length for the short set, default: 90')
		parser.add_argument('--short-upper',  dest = 'su', type = int, default = 110,
		help =  'Upper-bound on simulated read length for the short set, default: 110')
		
		parser.add_argument('--med-lower',  dest = 'ml', type = int, default = 135,
		help =  'Lower-bound on simulated read length for the medium set, default: 135')
		parser.add_argument('--med-upper',  dest = 'mu', type = int, default = 165,
		help =  'Upper-bound on simulated read length for the medium set, default: 165')
		
		parser.add_argument('--long-lower',  dest = 'll', type = int, default = 225,
		help =  'Lower-bound on simulated read length for the long set, default: 225')
		parser.add_argument('--long-upper',  dest = 'lu', type = int, default = 275,
		help =  'Upper-bound on simulated read length for the long set, default: 275')
		
		parser.add_argument('--xl-lower',  dest = 'xll', type = int, default = 270,
		help =  'Lower-bound on simulated read length for the extra long set, default: 270')
		parser.add_argument('--xl-upper',  dest = 'xlu', type = int, default = 330,
		help =  'Upper-bound on simulated read length for the extra long set, default: 330')

		
		parser.add_argument('--coverage',  dest = 'cov', default = 20.0, help = "Read coverage depth for simulated reads. Default 20") 
		parser.add_argument('--snprate',  dest = 'snps', default = 0.01, help = "Per base substitution likelihood. Default 0.01") 
		parser.add_argument('--insertrate',  dest = 'insrate', default = None, help = "Insertion rate. Default 1/19th of snprate.") 
		parser.add_argument('--delrate',  dest = 'delrate', default = None, help = "Deletion rate. Default 1/19th of snprate.") 
	
		parser.add_argument('--sensitivity', dest = 'dia_sens', default = 1, type = int, help = "DIAMOND sensitivity level integer from 1 to 4. Default 1. 1 = sensitive, 2 = more sensitive, 3 = very sensitive, 4 = ultra sensitive.")
		parser.add_argument('--use_blast', dest = 'use_blast', action = 'store_true', help = "Use BLASTx for read alignment instead of DIAMOND.")
		
	
	'''
	#This needs more thinking.
	if action == "add-genome":
		parser.description = "ROCkOut add-genome"
		parser.add_argument('--genome',  dest = 'gen', default = None, help = "Genome to add.") 
		parser.add_argument('--protein',  dest = 'prot', default = None, help = "If this is a target sequence") 
	'''
		
	if action == "refine":
		parser.description = "ROCkOut refine GUI mode"
		
	if action == "refine-non-interactive":
		parser.description = "ROCkOut refine non-interactive mode"
	
	if action == "align":
		parser.description = "ROCkOut align. Align one or several reads to the positive protein set of a ROCkOut project. Requires downloading a complete ROCkOut model or running ROCkOut through the refine module on your project, first."
		parser.add_argument('-i', '--input_reads',  dest = 'reads', default = None, help =  'A comma-sep list of paths to nucleotide FASTA format reads to align to a ROCkOut project. ')
		parser.add_argument('--reads_dir',  dest = 'reads_dir', default = None, help =  'A directory containing any number of FASTA format reads to align.')
		
		parser.add_argument('-f','--filter_directory', dest = 'filter_dir', default = None, help = 'Output directory for read alignments and filtered read outputs.')
		parser.add_argument('--use_blast', dest = 'use_blast', action = 'store_true', help = "Use BLASTx for read alignment instead of DIAMOND.")
		
	if action == "filter":
		parser.description = "ROCkOut filter: give a ROCkOut project directory and a reads directory created with ROCkOut align to filter reads."
		parser.add_argument('-f','--filter_directory', dest = 'filter_dir', default = None, help = 'A directory containing two subdirectories: "alignments", and "original_reads." These will be used by ROCkOut to filter the alignments')
		
	
	if action == "place":
		parser.description = "ROCkOut multiple alignment"
		parser.add_argument('-f','--filter_directory', dest = 'filter_dir', default = None, help = 'A directory containing five subdirectories as created by ROCkOut filter.')
		parser.add_argument('-p','--placement_count', dest = 'placements', default = 1, help = 'Max number of phylogenetic placements to keep per read. Default 1 (best single placement).')
		parser.add_argument('-g','--placement_group', dest = 'placement_target', default = "positive", help = 'Use (default) "positive", "negative", or "both" reference proteins as the phylogenetic target group. "Negative" and "both" depend on both positive and negative groups being present in the project.')
		
	#These are universal and I want them last
	parser.add_argument('-t', '--threads',  dest = 'threads', default = 1, help =  'Num threads to use for parallel processing')
	parser.add_argument('-q', '--quiet',  dest = 'quiet', action = 'store_true', help =  'Suppress progress reports')
		
	args, unknown = parser.parse_known_args()
	return parser, args
	
def install_check():
	all_ok = True
	try:
		import requests
	except:
		all_ok = False
		print("Python module requests not installed.")
	try:
		import numpy
	except:
		all_ok = False
		print("Python module numpy not installed.")
	try:
		import pandas
	except:
		all_ok = False
		print("Python module pandas not installed.")
	try:
		import plotly
	except:
		all_ok = False
		print("Python module plotly not installed.")
	try:
		import dash
	except:
		all_ok = False
		print("Python module Dash not installed.")
	#ext tools
	try:
		arg = subprocess.run(["muscle", "-version"], capture_output=True)
		version = arg.stdout.decode().strip()
		if version != "MUSCLE v3.8.31 by Robert C. Edgar":
			all_ok = False
			print("Muscle installed, but version is wrong.")
			print("Your version is:", version)
			print("You need MUSCLE v3.8.31.")
		
	except:
		all_ok = False
		print("Muscle not installed. Only muscle version=3.8.31 is acceptable.")
		
	try:
		arg = subprocess.run(["diamond", "--version"], capture_output=True)
	except:
		all_ok = False
		print("DIAMOND not installed.")
	
	try:
		arg = subprocess.run(["randomreads.sh", "--version"], capture_output=True)
		#version = arg.stdout.decode().strip()
		version = arg.stderr.decode().split("\n")[1]

		if version != "BBMap version 38.93":
			all_ok = False
			print("BBMap is installed, but version is wrong.")
			print("Your version is:", version)
			print("You need BBMap version 38.93.")
		
	except:
		all_ok = False
		print("BBMap not installed. Only BBMap version=38.93 is acceptable.")
	
	return all_ok
	
def run_download(parser, opts):
	from modules.rocker_0_download_from_uniprot import download
	download(parser, opts)
	
def run_build(parser, opts):
	from modules.rocker_1_read_simulation import build_project
	build_project(parser, opts)
	
def run_refine(parser, opts):
	from modules.rocker_3_dash import run_rocker_dash
	#run_rocker_dash(parser, opts)
	run_rocker_dash()
	
def run_refine_non_interactive(parser, opts):
	from modules.rocker_4_refiner import non_interactive
	non_interactive(parser, opts)
	
#Align reads to a project's proteins
def run_align(parser, opts):
	from modules.rocker_align_to_refs import align_to_refs
	align_to_refs(parser, opts)
	
def run_filter(parser, opts):
	from modules.rocker_filter import do_filter
	do_filter(parser, opts)
	
def run_pplacer(parser, opts):
	from modules.pplacer.rocker_phylomap_place import phylomap_place
	phylomap_place(parser, opts)

def main():
	valid_actions = ['download', 'build', 'refine', 'align', 'filter', 'place']

	if len(sys.argv) < 2:
		print("ROCkOut needs to be given an action! One of:", valid_actions)
		sys.exit("You can also do rocker install-check to see if all of the appropriate dependencies are installed.")
	
	action = sys.argv[1]
	
	if action == "install-check":
		all_good = install_check()
		if all_good:
			sys.exit("All dependencies are installed correctly.")
		else:
			sys.exit("At least one dependency is not installed properly. Not all behaviors of ROCkOut will be available.")
	
	
	if action not in valid_actions:
		print("Action '" + str(action) + "' not recognized.")
		sys.exit("The action must be one of:", valid_actions)
		
	parser, opts = options(action)
	
	if action == "download":
		run_download(parser, opts)
		
	if action == "build":
		run_build(parser, opts)
		
	#if action == "refine":
	#	run_refine(parser, opts)
		
	if action == "refine":
		run_refine_non_interactive(parser, opts)
		
	if action =="align":
		run_align(parser, opts)
		
	if action =="filter":
		run_filter(parser, opts)
		
	if action == "place":
		run_pplacer(parser, opts)
		
	
if __name__ == "__main__":
	main()	
	
	
	
	
	
	