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
	def __init__(self, project_directory, use_blast = False, thds = 1, inpaths = None, outpaths = None, filter_prep = None):
		self.prjdir = project_directory
		self.mn = import_rocker_project(self.prjdir)
		
		self.diamond_db = self.mn.targets_dia
		self.blast_db = self.mn.targets_blast
		
		self.reads_to_align = inpaths
		self.output_paths = outpaths
		self.filtered_reads = None
		
		self.commands = []
		
		self.command_log = ""
		self.stdout = ""
		self.stderr = ""
		
		self.ready_for_filter = filter_prep
		
		self.use_blast = use_blast
		self.threads = thds
	
	def craft_pairs(self):
		if self.reads_to_align is not None:
			
			self.reads_to_align = self.reads_to_align.strip().split(",")
			if self.ready_for_filter:
			
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
					
					if self.use_blast:
						arg = "blastx -db {database} -query {input} -out {output} -outfmt 6 -num_threads {threads}"
						arg = arg.format(database = self.blast_db, input = read, output = alignment, threads = self.threads)
						#Maybe needs a besthit filter here, too?
					else:
						arg = "diamond blastx --very-sensitive --unal 0 --db {database} --query {input} --out {output} --outfmt 6 --threads {threads}"
						arg = arg.format(database = self.diamond_db, input = read, output = alignment, threads = self.threads)
					
					self.commands.append(arg)
					
			else:
			
				if self.output_paths is not None:
					self.output_paths = self.output_paths.strip().split(",")
				else:
					self.output_paths = [os.path.normpath(read_file + "_ROCkOut_alignments.txt") for read_file in self.reads_to_align]
				
				if len(self.reads_to_align) == len(self.output_paths):
					for r, o in zip(self.reads_to_align, self.output_paths):
						if self.use_blast:
							arg = "blastx -db {database} -query {input} -out {output} -outfmt 6 -num_threads {threads}"
							arg = arg.format(database = self.blast_db, input = r, output = o, threads = self.threads)
							#Maybe needs a besthit filter here, too?
						else:
							arg = "diamond blastx --very-sensitive --unal 0 --db {database} --query {input} --out {output} --outfmt 6 --threads {threads}"
							arg = arg.format(database = self.diamond_db, input = r, output = o, threads = self.threads)
						
						self.commands.append(arg)
				else:
					print("I can't work with a different number of reads and outputs.")
		else:
			print("I need reads to align.")
	
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
			self.command_log += alignment_arg + "\n"
			alignment_arg = alignment_arg.split()
			proc = subprocess.run(alignment_arg, stdout=subprocess.PIPE, stderr = subprocess.PIPE)
			self.stdout += proc.stdout.decode() + "\n"
			self.stderr += proc.stderr.decode() + "\n"
			#print(alignment_arg)
			#print(proc.stdout.decode())
			#print(proc.stderr.decode())
			
		if self.filtered_reads is not None:
			for tup in self.filtered_reads:
				raw = tup[0]
				aln = tup[1]
				clean = tup[2]
				self.clean_reads(raw, aln, clean)
				
			
		
	
def align_to_refs(parser, opts):
	dir = opts.dir
	reads = opts.reads
	output = opts.alignments
	threads = opts.threads
	blast = opts.use_blast
	
	if dir is None:
		parser.print_help()
		sys.exit()
		
	if reads is None:
		parser.print_help()
		sys.exit()
		
	filter_prep = opts.filter_dir
	
	mn = rocker_aligner(
	project_directory = dir,
	use_blast = blast,
	thds = threads, 
	inpaths = reads, 
	outpaths = output,
	filter_prep = filter_prep)
	
	mn.craft_pairs()
	mn.align_reads()
	
	