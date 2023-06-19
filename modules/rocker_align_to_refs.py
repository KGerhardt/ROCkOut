import sys
import os
import subprocess

from .rocker_project_manager import project_manager

def import_rocker_project(project_directory):
	manager = project_manager(project_directory)
	manager.parse_project_directory()
	manager.collect_final_info()
	#What we need here is the final outputs/database contents, returning the blast DB base and the diamond db
	return manager
	
class rocker_aligner:
	def __init__(self, project_directory, use_blast = False, thds = 1, inpaths = None, filter_prep = None):
		self.prjdir = project_directory
		self.mn = import_rocker_project(self.prjdir)
		
		self.diamond_db = self.mn.targets_dia
		self.blast_db = self.mn.targets_blast
		
		self.reads_to_align = inpaths
		self.output_paths = None
		self.filtered_reads = None
		
		self.commands = []
		
		self.command_log = ""
		self.stdout = ""
		self.stderr = ""
		
		self.ready_for_filter = filter_prep
		
		self.use_blast = use_blast
		self.threads = thds
	
	def craft_pairs(self):
		self.reads_to_align = self.reads_to_align.strip().split(",")

		if not os.path.exists(self.ready_for_filter):
			os.makedirs(self.ready_for_filter, exist_ok = True)
		alns = os.path.normpath(self.ready_for_filter + "/alignments")
		orig = os.path.normpath(self.ready_for_filter + "/original_reads/")
		if not os.path.exists(alns):
			os.makedirs(alns, exist_ok = True)
		if not os.path.exists(orig):
			os.makedirs(orig, exist_ok = True)
	
		self.filtered_reads = []
		for read in self.reads_to_align:
			basename = os.path.basename(read)
			alignment = os.path.normpath(self.ready_for_filter + "/alignments/"+basename+ "_ROCkOut_alignments.txt")
			filtered_raw = os.path.normpath(self.ready_for_filter + "/original_reads/"+basename+ "_raw_reads_.txt")
			
			self.filtered_reads.append((read, alignment, filtered_raw,))
			
			bla_fmt = "'6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'"
			dia_fmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
			if self.use_blast:
				arg = "blastx -db {database} -query {input} -out {output} -outfmt {bla_fmt} -num_threads {threads}"
				arg = arg.format(database = self.blast_db, input = read, output = alignment, bla_fmt = bla_fmt, threads = self.threads)
				#Maybe needs a besthit filter here, too?
			else:
				arg = "diamond blastx --very-sensitive --unal 0 --db {database} --query {input} --out {output} --outfmt {dia_fmt} --threads {threads}"
				arg = arg.format(database = self.diamond_db, input = read, output = alignment, dia_fmt = dia_fmt, threads = self.threads)
			
			self.commands.append(arg)
					

	def clean_reads(self, r, a, c):
		reads_to_keep = {}
		with open(a) as fh:
			for line in fh:
				if line.startswith("#"):
					continue
				segs = line.strip().split("\t")
				read_name = segs[0]
				score = segs[11]
				if read_name not in reads_to_keep:
					reads_to_keep[read_name] = score
					
		current_read = None
		current_defline = None
		current_seq = ""
		out = open(c, "w")
		with open(r) as fh:
			for line in fh:
				if line.startswith(">"):
					if current_read is not None:
						if current_read in reads_to_keep:
							out.write(current_defline)
							out.write(current_seq)
					current_defline = line
					current_read = current_defline.strip().split()[0][1:]
					current_seq = ''
				else:
					current_seq += line
		
		if current_read is not None:
			if current_read in reads_to_keep:
				out.write(current_defline)
				out.write(current_seq)
		
		out.close()
			
	def align_reads(self):
		for alignment_arg in self.commands:
			print("Executing alignment:")
			print(alignment_arg)
			self.command_log += alignment_arg + "\n"
			alignment_arg = alignment_arg.split()
			proc = subprocess.run(alignment_arg, stdout=subprocess.PIPE, stderr = subprocess.PIPE)
			self.stdout += proc.stdout.decode() + "\n"
			self.stderr += proc.stderr.decode() + "\n"
			#print(alignment_arg)
			#print(proc.stdout.decode())
			#print(proc.stderr.decode())
			print("")
			
		if self.filtered_reads is not None:
			for tup in self.filtered_reads:
				raw = tup[0]
				aln = tup[1]
				clean = tup[2]
				self.clean_reads(raw, aln, clean)
				
	
def align_to_refs(parser, opts):
	dir = opts.dir
	reads = opts.reads
	reads_dir = opts.reads_dir
	#output = opts.alignments
	threads = opts.threads
	blast = opts.use_blast
	
	if dir is None:
		parser.print_help()
		print("")
		print("I need a ROCkOut project directory!")
		sys.exit()
		
	if reads is None and reads_dir is None:
		parser.print_help()
		print("")
		print("I need at least one input read file or a read directory to filter!")
		sys.exit()
		
	if reads is not None and reads_dir is not None:
		parser.print_help()
		print("")
		print("Please supply only one kind of read input: paths with --reads or a directory with --reads_dir, not both")
		sys.exit()
		
	if reads is not None:
		if os.path.isdir(reads):
			parser.print_help()
			print("")
			print("You supplied a directory to --reads. Did you mean to use --reads_dir?")
			sys.exit()
	
	if reads_dir is not None:
		if os.path.isfile(reads_dir):
			parser.print_help()
			print("")
			print("You supplied a file path to --reads_dir. Did you mean to use --reads?")
			sys.exit()
		else:
			reads = ",".join([os.path.normpath(reads_dir+"/"+r) for r in os.listdir(reads_dir)])
		
	filter_prep = opts.filter_dir
		
	if filter_prep is None:
		print("Please specify an output directory! This directory will be created for you.")
		parser.print_help()
		sys.exit()
	
	mn = rocker_aligner(
	project_directory = dir,
	use_blast = blast,
	thds = threads, 
	inpaths = reads,
	filter_prep = filter_prep)
	
	mn.craft_pairs()
	mn.align_reads()
	
	