import sys
import os

class pplacer_operator:
	def __init__(self, rocker_dir, reads_dir, threads = 1, max_placements = 1, target_group = "positive"):
		self.rocker = os.path.normpath(rocker_dir)
		self.refpak = None
		self.ref_ma = None
		self.group = target_group
		
		self.reads_base = os.path.normpath(reads_dir)

		self.commands = None
		
		self.threads = threads
		self.max_place = max_placements
		
	def make_dir(self, dir):
		to_make = os.path.normpath(dir)
		if not os.path.exists(to_make):
			os.makedirs(to_make, exist_ok = True)
		
	#Find taxtastic refpak for placement
	def locate_refpak(self):
		pak_loc = os.path.normpath(self.rocker+"/final_outputs/phylogenetic_placement/{ext}/reference_package")
		pak_loc = pak_loc.format(ext = self.group)
		rp = os.listdir(pak_loc)
		rp = os.path.normpath(pak_loc + "/" + rp[0])
		if os.path.exists(rp):
			self.refpak = rp
			
	def locate_ref_ma(self):
		ma = os.path.normpath(self.rocker+"/final_outputs/phylogenetic_placement/{ext}/source/combined_genomes_multiple_alignment.fasta")
		ma = ma.format(ext = self.group)
		if os.path.exists(ma):
			self.ref_ma = ma
		
	def prepare_commands(self):
		#ma_reads
		reads_loc = os.path.normpath(self.reads_base + "/ROCkOut_passing_read_fastas")

		mafft_dir = os.path.normpath(self.reads_base +"/read_multiple_alignments")
		jplace_dir = os.path.normpath(self.reads_base +"/pplacer_jplace")
		tsv_dir = os.path.normpath(self.reads_base +"/pplacer_tsv")
		
		self.make_dir(mafft_dir)
		self.make_dir(jplace_dir)
		self.make_dir(tsv_dir)
		
		self.commands = []
		
		for read in os.listdir(reads_loc):
			read_basename = os.path.basename(read)
			read_basename = read_basename.replace(".filtered.fasta", "")
			read_path = os.path.normpath(reads_loc + "/" + read)
			ma_output = os.path.normpath(mafft_dir + "/"+read_basename+"_multiple_alignment.fasta")
			jplace_out = os.path.normpath(jplace_dir + "/"+read_basename+".jplace")
			tsv_out = os.path.normpath(tsv_dir + "/"+read_basename+".jplace.tsv")
			
			#Multiply align reads
			mafft_comm = "mafft --quiet --thread {thds} --add {combined_reads} --reorder {orig_ma} > {out}"
			mafft_comm = mafft_comm.format(thds = self.threads, combined_reads = read_path, orig_ma =self.ref_ma, out =ma_output)
			
			#Place the MA with pplacer
			pplace_comm = "pplacer --discard-nonoverlapped -c {refpkg} -o {out} {reads_aln} -j {thds} --keep-at-most {place_count}"
			pplace_comm = pplace_comm.format(refpkg = self.refpak, out = jplace_out, reads_aln = ma_output, thds = self.threads, place_count = self.max_place)
			
			#Convert the jplace to human-readable
			guppy_comm = "guppy to_csv {jplace} > {tsv}_temp.csv"
			guppy_comm = guppy_comm.format(jplace = jplace_out, tsv = tsv_out)
			
			comm_group = (mafft_comm, pplace_comm, guppy_comm, tsv_out,)
			self.commands.append(comm_group)
				
	def run_commands(self):
		for command in self.commands:
			mafft = command[0]
			place = command[1]
			guppy = command[2]
			tsv = command[3]
			os.system(mafft)
			os.system(place)
			os.system(guppy)
			swap = open(tsv, "w")
			fh = open(tsv+"_temp.csv")
			for line in fh:
				line = line.replace(",", "\t")
				swap.write(line)
			swap.close()
			fh.close()
			os.remove(tsv+"_temp.csv")
	
	def run_place(self):
		self.locate_refpak()
		self.locate_ref_ma()
		self.prepare_commands()
		self.run_commands()

def check_reads_dir(directory):
	all_good = True
	if not os.path.exists(os.path.normpath(directory+"/ROCkOut_passing_read_fastas")):
		all_good = False
		
	return all_good
		
def phylomap_place(parser, opts):
	rocker_dir = opts.dir
	reads_dir = opts.filter_dir
	choice = opts.placement_target
	
	valid_placement_choices = ["positive", "negative", "both"]
	if choice not in valid_placement_choices:
		print("Placement group", choice, "not found.")
		print("Must be one of", *valid_placement_choices)
		parser.print_help()
		sys.exit()
	
	if rocker_dir is None or reads_dir is None:
		print("I need both a ROCkOut project dirctory and a reads directory to place reads!")
		parser.print_help()
		sys.exit()
		
	ok_to_go = check_reads_dir(reads_dir)
	if not ok_to_go:
		print("You need to run ROCkOut filter on your reads directory first!")
		print("If you really want to phylogenetically place unfiltered reads,\ncreate a subdirectory called 'ROCkOut_passing_read_fastas' in the\nreads directory you're passing and place nt FASTA reads there.")
		sys.exit()
	
	threads = opts.threads
	place_ct = opts.placements
	
	# rocker_dir, reads_dir, threads = 1, max_placements = 1
	mn = pplacer_operator(rocker_dir = rocker_dir,
						reads_dir = reads_dir,
						threads = threads,
						max_placements = place_ct,
						target_group = choice)
	mn.run_place()
	