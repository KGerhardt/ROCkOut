import sys
import os
import subprocess
import re
import numpy as np
import multiprocessing

from rocker_project_manager import project_manager
from rocker_progress_tracker import progress_tracker

#Class for simulating and tagging reads in ROCker friendly format.
class read_tagger:
	#def __init__(self, input_fasta = None, input_gffs = None, positive = True, simulator = "BBMap", coordinates = None, is_pos = False):
	def __init__(self, base_name = None, input_fasta = None, coordinates = None):
		self.base = base_name
		
		self.input_base = os.path.basename(input_fasta)[:-10]
		
		self.input = input_fasta
		
		self.tags = os.path.normpath(self.base + "/tagged_reads/" + self.input_base + "_tagged_fasta.txt")
		
		self.coords = coordinates
		self.coord_starts = None
		self.coord_ends = None
		
	def prepare_outputs(self):
		tagged = os.path.normpath(self.base + "/tagged_reads") + "/"
		if not os.path.exists(tagged):
			os.mkdir(tagged)

	#coords file to start-end coordinates lists.
	def extract_from_coords(self):
		#my_coords = {}
		self.coord_starts = []
		self.coord_ends = []
		fh = open(self.coords)
		for line in fh:
			segs = line.strip().split()
			parent, prot_id, start, end, strand = segs[0], segs[1], int(segs[2]), int(segs[3]), segs[4]
			#my_coords[prot_id] = (start, end, strand)
			self.coord_starts.append(start)
			self.coord_ends.append(end)
		
		fh.close()
			
		#Switch to numpy.
		self.coord_starts = np.array(self.coord_starts)
		self.coord_ends = np.array(self.coord_ends)
		
	def tag_bbmap_reads(self):
		#A bbmap regex matcher with groups defined for the info we want to extract from each read.
		
		self.extract_from_coords()
		
		is_negative = ("/negative/" in self.tags)
		
		bbmap_match = r'SYN_(\d+)_(.+?(?=_))_(.+?(?=_))_\d_(.)_.+?(?=\.)\._(.+?(?=\$|\s)).+'
		fh = open(self.input)
		out = open(self.tags, 'w')
		in_seq = False
		malformed = False
		for line in fh:
			if line.startswith(">"):
				#Reset
				malformed = False
				
				in_seq = False
				id, fr, to, comp, genome_id = re.search(bbmap_match, line).groups()
				fr, to = int(fr), int(to)
				mn, mx = min(fr, to), max(fr, to)
				
				#No interference with ROCker's split scheme.
				genome_id = genome_id.replace(';', '_')
				
				#We use the end to figure out if there's an overlap with the gene start to the gene start.
				start_window = np.searchsorted(self.coord_starts, mx, side = 'right')
				#We use the start against the ends to figure out if there's an overlap with the gene end.
				end_window = np.searchsorted(self.coord_ends, mn, side = 'left')
				
				tagged_name = ';'.join([id, str(mn), str(mx), comp, genome_id])
				
				#Read is malformed for one reason or another. Skip and keep going with the others.
				if len(tagged_name) == 0:
					print(tagged_name)
					#Skip writing until the next seq.
					malformed = True
					continue
					
				if is_negative:
					print(">" + tagged_name + ";Negative", file = out)
				else:
					if (start_window - 1) == end_window:
						#The read falls inside a target gene window and we should tag it as on-target
						print(">" + tagged_name + ";Target", file = out)
					else:
						print(">" + tagged_name + ";Non_Target", file = out)
				
			else:
				if not malformed:
					#Write seq.
					out.write(line)
				
		fh.close()
		out.close()

		
project_directory = sys.argv[1]
try:
	threads = int(sys.argv[2])
except:
	print("Couldn't recognize threads as a number. Setting threads to 1.")
	threads = 1

def run_tagging(read_tag):
	#read_gen.prepare_outputs()
	read_tag.tag_bbmap_reads()
	return read_tag.base

def generate_reads(project_directory, threads):
	rocker = project_manager(directory = project_directory, threads = threads)
	rocker.parse_project_directory()
	
	#Coords currently expects 1 reads file/coords file, which isn't true anymore.
	rocker.parse_coords()
	rocker.parse_raw_reads()
	
	
	pos_reads = []
	for base in rocker.coords_pos:
		cut = len(rocker.positive[base])
		coords = {}
		for c in rocker.coords_pos[base]:
			coord_key = c[cut+8:-11]
			coords[coord_key] = c
		for r in rocker.read_fastas_pos[base]:
			read_key = r[cut+11:-10].split("_read_length_")[0]
			pos_reads.append((rocker.positive[base], r, coords[read_key],))
	
	neg_reads = []
	for base in rocker.coords_neg:
		cut = len(rocker.negative[base])
		coords = {}
		for c in rocker.coords_neg[base]:
			coord_key = c[cut+8:-11]
			coords[coord_key] = c
		for r in rocker.read_fastas_neg[base]:
			read_key = r[cut+11:-10].split("_read_length_")[0]
			neg_reads.append((rocker.negative[base], r, coords[read_key],))
		
	to_do = []
	
	for tup in pos_reads:
		gen = read_tagger(base_name = tup[0], input_fasta = tup[1], coordinates = tup[2])
		gen.prepare_outputs()
		to_do.append(gen)
		
	for tup in neg_reads:
		gen = read_tagger(base_name = tup[0], input_fasta = tup[1], coordinates = tup[2])
		gen.prepare_outputs()
		to_do.append(gen)
	

	prog_bar = progress_tracker(total = len(to_do), message = "Tagging reads...")
	pool = multiprocessing.Pool(threads)
	for result in pool.imap_unordered(run_tagging, to_do):
		prog_bar.update()
	pool.close()
	pool.join()
	print("Reads tagged!")
	

generate_reads(project_directory = project_directory, threads = threads)