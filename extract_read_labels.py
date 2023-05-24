import sys
import os

from modules.rocker_project_manager import project_manager

prjdir = sys.argv[1] # provide the rocker project directory.

def collect_group(manager, base, pos = True):
	collection = {}
	final_return = None
	if pos:
		for g in manager.genomes_pos[base]:
			genome_base = os.path.basename(g)
			genome_base = genome_base.split(".fasta")[0]
			collection[genome_base] = {}
			collection[genome_base]["target_proteins"] = []
			for target in manager.targets_nt[base]:
				with open(target) as fh:
					for line in fh:
						if line.startswith(">"):
							protein_id = line.strip().split()[0]
							#print(protein_id)
							protein_id = protein_id.split("__")[-1]
							collection[genome_base]["target_proteins"].append(protein_id)
			
			collection[genome_base]["target_proteins"] = set(collection[genome_base]["target_proteins"])
			
			collection[genome_base]["protein_info"] = []
			
			#print(base, targets)
			for target in manager.gffs_pos[base]:
				
				with open(target) as fh:
					for line in fh:
						if line.startswith("#"):
							continue
							
						segs = line.strip().split()
						tag = segs[2]
						if tag == "CDS":
							start = int(segs[3])
							stop = int(segs[4])
							strand = (segs[6] == "+")
							annot = segs[8]
							protein_name = annot.split(";")[0]
							protein_name = protein_name[3:]
							#print(collection[genome_base]["target_proteins"], protein_name)
							if protein_name in collection[genome_base]["target_proteins"]:
								next_group = (start, stop, strand, protein_name,)
								collection[genome_base]["protein_info"].append(next_group)
	else:
		for g in manager.genomes_neg[base]:
			genome_base = os.path.basename(g)
			genome_base = genome_base.split(".fasta")[0]
			collection[genome_base] = {}
			collection[genome_base]["target_proteins"] = []
			for target in manager.targets_nt[base]:
				with open(target) as fh:
					for line in fh:
						if line.startswith(">"):
							protein_id = line.strip().split()[0]
							#print(protein_id)
							protein_id = protein_id.split("__")[-1]
							collection[genome_base]["target_proteins"].append(protein_id)
			
			collection[genome_base]["target_proteins"] = set(collection[genome_base]["target_proteins"])
			
			collection[genome_base]["protein_info"] = []
			
			#print(base, targets)
			for target in manager.gffs_neg[base]:
				
				with open(target) as fh:
					for line in fh:
						if line.startswith("#"):
							continue
							
						segs = line.strip().split()
						tag = segs[2]
						if tag == "CDS":
							start = int(segs[3])
							stop = int(segs[4])
							strand = (segs[6] == "+")
							annot = segs[8]
							protein_name = annot.split(";")[0]
							protein_name = protein_name[3:]
							#print(collection[genome_base]["target_proteins"], protein_name)
							if protein_name in collection[genome_base]["target_proteins"]:
								next_group = (start, stop, strand, protein_name,)
								collection[genome_base]["protein_info"].append(next_group)
		
	final_return = {}
	for genome in collection:
		final_return[genome] = collection[genome]["protein_info"]
		
	return final_return

def trawl(pd):
	'''
	self.project_base = directory
	self.threads = threads
	
	#Protein folders
	self.positive = {}
	self.negative = {}
	
	#Coordinate files for target gene
	self.coords_pos = {}
	self.coords_neg = {}
	
	#Proteome AA FASTA files
	self.proteomes_pos = {}
	self.proteomes_neg = {}
	
	#GFF files for target genomes
	self.gffs_pos = {}
	self.gffs_neg = {}
	
	#fastas of the original genomes
	self.genomes_pos = {}
	self.genomes_neg = {}
	self.gen_lengths_pos = None
	self.gen_lengths_neg = None
	
	#fastq reads generated by randomreads.sh for the genomes
	self.read_fastqs_pos = {}
	self.read_fastqs_neg = {}
	
	#fasta translations of the fastqs
	self.read_fastas_pos = {}
	self.read_fastas_neg = {}
	
	#renamed/tagged reads for the fasta translations
	self.tagged_reads_pos = {}
	self.tagged_reads_neg = {}
	
	#Blast/diamond alignments of the tagged reads
	self.alignments_pos = {}
	self.alignments_neg = {}
	
	self.targets = {}
	self.targets_nt = {}
	self.active_targets = None
	
	self.mult_aln_base = None
	self.mult_aln_files = {}
	
	self.outputs_base = None
	self.rocker_filter = None
	self.targets_blast = None
	self.targets_dia = None
	'''
	
	mn = project_manager(directory = pd)
	mn.parse_project_directory()
	mn.parse_genomes()
	mn.parse_gffs()
	mn.parse_coords()
	mn.parse_tagged_reads()
	mn.parse_targets()
	
	#We want to collect a list of positive and negative info for each target and negative genome, 
	#knowing that simmed reads at a length are identical for each genome in each grp
	read_label_reference = {}
	

	for p in mn.positive:
		this_protein = collect_group(mn, p, pos = True)
		#print(p, this_protein)
		for genome in this_protein:
			if genome not in read_label_reference:
				read_label_reference[genome] = {"pos":[], "neg":[]}
			for target_protein in this_protein[genome]:
				#Starts and stops for true positives
				read_label_reference[genome]["pos"].append((target_protein[0], target_protein[1], target_protein[3],))
		
	for p in mn.negative:
		this_protein = collect_group(mn, p, pos = False)
		for genome in this_protein:
			if genome not in read_label_reference:
				read_label_reference[genome] = {"pos":[], "neg":[]}
			for target_protein in this_protein[genome]:
				#Starts and stops for true negatives
				read_label_reference[genome]["neg"].append((target_protein[0], target_protein[1], target_protein[3],))
	
	return read_label_reference, mn
		
read_labeller, manager = trawl(prjdir)

def one_read(f, p, n):
	with open(f) as fh:
		for line in fh:
			if line.startswith(">"):
				isp = False
				segs = line.strip().split(";")
				start = int(segs[1])
				end = int(segs[2])
				segs[-1] = "Non_Target"
				for start_stop_name in p:
					if (start >= start_stop_name[0] and start < start_stop_name[1]) or (end > start_stop_name[0] and end <= start_stop_name[1]):
						segs[-1] = "source_protein="+start_stop_name[2]+";Positive"
				for start_stop_name in n:
					if (start >= start_stop_name[0] and start < start_stop_name[1]) or (end > start_stop_name[0] and end <= start_stop_name[1]):
						segs[-1] = "source_protein="+start_stop_name[2]+";Negative"
				segs = ";".join(segs)
				print(segs)
			else:
				print(line.strip())
				
def process_reads(labels, mn):
	seen_genomes = {}
	for p in mn.positive:
		for read in mn.tagged_reads_pos[p]:
			genome = os.path.basename(read)
			read_len = int(genome.split("_read_len_")[1].split("_tagged.fasta")[0])
			if read_len not in seen_genomes:
				seen_genomes[read_len] = []
			if genome not in seen_genomes[read_len]: #Each genome only needs seen once
				seen_genomes[read_len].append(genome)
				genome = genome.split("_read_len_")[0]
				#print(genome, read_len)
				positive_set = labels[genome]["pos"]
				negative_set = labels[genome]["neg"]
				one_read(read, positive_set, negative_set)
		
	for p in mn.negative:
		for read in mn.tagged_reads_neg[p]:
			genome = os.path.basename(read)
			read_len = int(genome.split("_read_len_")[1].split("_tagged.fasta")[0])
			if read_len not in seen_genomes:
				seen_genomes[read_len] = []
			if genome not in seen_genomes[read_len]: #Each genome only needs seen once
				seen_genomes[read_len].append(genome)
				genome = genome.split("_read_len_")[0]
				print(genome, read_len)
				positive_set = labels[genome]["pos"]
				negative_set = labels[genome]["neg"]
	
process_reads(read_labeller, manager)