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
		parser.description = "ROCkOut build"
		parser.add_argument('--short',  dest = 'short', default = "90,100,110", 
		help =  'Comma-sep triplet of very short read lengths, default: 90,100,110')
		parser.add_argument('--med',  dest = 'med', default = "180,200,220", 
		help =  'Comma-sep triplet of moderate read lengths, default: 180,200,220')
		parser.add_argument('--long',  dest = 'long', default = "270,300,330", 
		help =  'Comma-sep triplet of long read lengths, default: 270,300,330')
		parser.add_argument('--xl',  dest = 'xl', default = "360,400,440", 
		help =  'Comma-sep triplet of very long read lengths, default: 360,400,440')
		
		parser.add_argument('--coverage',  dest = 'cov', default = 20.0, help = "Read coverage depth for simulated reads. Default 20") 
		parser.add_argument('--snprate',  dest = 'snps', default = 0.01, help = "Per base substitution likelihood. Default 0.01") 
		parser.add_argument('--insertrate',  dest = 'insrate', default = None, help = "Insertion rate. Default 1/19th of snprate.") 
		parser.add_argument('--delrate',  dest = 'delrate', default = None, help = "Deletion rate. Default 1/19th of snprate.") 
	
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
		
	
	if action == "align-and-filter":
		parser.description = "ROCkOut align and filter"
		
	if action == "filter":
		parser.description = "ROCkOut filter"
		
	if action == "mult-aln":
		parser.description = "ROCkOut multiple alignment"
		
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
	
def run_ma(parser, opts):
	from modules.rocker_0B_create_MA import make_ma
	make_ma(parser, opts)
	
def run_add(parser, opts):
	pass
	
def run_build(parser, opts):
	from modules.rocker_1_generate_reads import generate
	from modules.rocker_2_tag_reads import tag_reads
	from modules.rocker_3_align_to_refs import align_to_refs
	
	generate(parser, opts)
	tag_reads(parser, opts)
	align_to_refs(parser, opts)
	
	
def run_refine(parser, opts):
	from modules.rocker_4_refiner import non_interactive
	non_interactive(parser, opts)
	
def run_refine_non_interactive(parser, opts):
	from modules.rocker_4_refiner import non_interactive
	non_interactive(parser, opts)
	
def run_align_and_filter(parser, opts):
	pass
	
def run_filter(parser, opts):
	pass
	
def run_pplace_build(parser, opts):
	from modules.rocker_6_build_ref_tree import make_reftree
	make_reftree(parser, opts)


def main():
	if len(sys.argv) < 2:
		print("ROCkOut needs to be given an action! One of  'download', 'add', 'build', 'refine', 'refine-ni', 'align-and-filter', 'filter'")
		sys.exit("You can also do rocker install-check to see if all of the appropriate dependencies are installed.")
	
	action = sys.argv[1]
	
	if action == "install-check":
		all_good = install_check()
		if all_good:
			sys.exit("All dependencies are installed correctly.")
		else:
			sys.exit("At least one dependency is not installed properly. Not all behaviors of ROCkOut will be available.")
	
	valid_actions = ['download', 'mult-aln', 'add', 'build', 'refine', 'refine-non-interactive', 'align', 'filter', 'pplace-prep']
	if action not in valid_actions:
		print("Action '" + str(action) + "' not recognized.")
		sys.exit("The action must be one of 'download', 'mult-aln', 'add', 'build', 'refine', 'refine-non-interactive', 'align-and-filter', 'filter', 'pplace-prep'")
		
	parser, opts = options(action)
	
	if action == "download":
		run_download(parser, opts)
		
	if action == 'mult-aln':
		run_ma(parser, opts)
		
	if action == "build":
		run_build(parser, opts)
		
	if action == "refine":
		pass
		#run_refine(parser, opts)
		
	if action == "refine-non-interactive":
		run_refine_non_interactive(parser, opts)
		
	if action == "pplace-prep":
		run_pplace_build(parser, opts)
		
	
if __name__ == "__main__":
	main()	
	
	
	
	
	
	