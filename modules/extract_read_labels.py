import sys
import os

try:
	from .rocker_project_manager import project_manager
except:
	from rocker_project_manager import project_manager

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
	mn = project_manager(directory = pd)
	mn.parse_project_directory()
	mn.parse_genomes()
	mn.parse_gffs()
	mn.parse_coords()
	mn.parse_tagged_reads()
	mn.parse_targets()
	mn.parse_aligns()
	mn.parse_probable_targets()
	mn.parse_multiple_alignment()
	
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
		
#This deals with FASTAs and is designed to output labelled reads for a user.
def one_read_ext(f, p, n):
	with open(f) as fh:
		for line in fh:
			if line.startswith(">"):
				segs = line.strip().split(";")
				start = int(segs[1])
				end = int(segs[2])
				segs[-1] = "Non_Target"
				for start_stop_name in p:
					if (start >= start_stop_name[0] and start < start_stop_name[1]) or (end > start_stop_name[0] and end <= start_stop_name[1]):
						segs[-1] = "source_protein="+start_stop_name[2]+";Positive"
				if label == "Non_Target": #Prevent overwrite of a positive label.
					for start_stop_name in n:
						if (start >= start_stop_name[0] and start < start_stop_name[1]) or (end > start_stop_name[0] and end <= start_stop_name[1]):
							segs[-1] = "source_protein="+start_stop_name[2]+";Negative"
				segs = ";".join(segs)
				print(segs)
			else:
				print(line.strip())

#This deals with alignments and is designed for internal ROCkOut usage.
#This also gets the lengths of reads because we need those.
def one_read_int(f, p, n, do_bh = True, valid_tgts = None, readlen_file = None):
	read_results = {}
	read_scores = {}
	#pct = 0
	#nct = 0
	with open(f) as fh:
		for line in fh:
			segs = line.strip().split("\t")
			
			target = segs[1]
			if valid_tgts is not None:
				if tgt not in valid_tgts:
					continue #Skip alignments to invalid target seqs.
			
			read_name = segs[0]
			read_name = read_name.split(";")
			read_ID = read_name[0]
			
			#These needed to be alignments, not read loc
			start = int(read_name[1])
			end = int(read_name[2])
			
			#This is an alignment local to the thing.
			alignment_range = [int(segs[8]), int(segs[9])]
			
			alnstart = min(alignment_range)
			alnend = max(alignment_range)
			
			read_name = ";".join(read_name[:-1]) #Exclude the label.
			label = "Non_Target"
			
			bitscore = float(segs[11])
			evalue = float(segs[10])
			pct_id = float(segs[2])
			
			overlap_pct = float(segs[14])
			
			aln_len = int(segs[3])
			qlen = int(segs[12]) / 3 #bp to AA lengths
			
			#If the qlen is specified in read aln fmt from earlier.
			if len(segs) > 12:
				pct_aln = round(aln_len/qlen * 100, 2)
			else:
				pct_aln = -1.0
			
			if read_name not in read_scores:
				read_scores[read_name] = bitscore #First appearance of a read is always the current winning best-hit assignment.
			else:
				if read_scores[read_name] >= bitscore:
					#Just skip a read that's already got a better hit in this dataset.
					continue
				else:
					read_scores[read_name] = bitscore #Set the new bitscore max and process the read as the new winner
			
			for start_stop_name in p:
				if (start >= start_stop_name[0] and start < start_stop_name[1]) or (end > start_stop_name[0] and end <= start_stop_name[1]):
					label = "Positive"
					#pct += 1
					#segs[-1] = "source_protein="+start_stop_name[2]+";Positive"
			if label == "Non_Target":
				#nct += 1
				for start_stop_name in n:
					if (start >= start_stop_name[0] and start < start_stop_name[1]) or (end > start_stop_name[0] and end <= start_stop_name[1]):
						label = "Negative"
					#segs[-1] = "source_protein="+start_stop_name[2]+";Negative"

			read_results[read_name] = [read_name, target, label, alnstart, alnend, bitscore, evalue, pct_id, pct_aln, overlap_pct]
	
	'''
	#I don't think importing the read lengths is necessarily worth it for the filter setup. Using the sim readlen is OK I think.
	cur_sl = 0
	cur_read = ""
	is_valid = False
	read_lengths = {}
	with open(readlen_file) as fh:
		for line in fh:
			if line.startswith(">"):
				segs = line.strip().split()
				read_name = segs[0][1:]
				read_name = read_name.split(";")[0]
				
				if len(cur_read) > 0:
					#read_lengths[cur_read] = cur_sl
					read_results[read_name].append(cur_sl)
				if read_name not in read_results:
					cur_read = ""
					is_valid = False
				else:
					cur_read = read_name
					is_valid = True
					
				cur_sl = 0
			else:
				if is_valid:
					line = line.strip()
					cur_sl += len(line)
	'''
	#print(f, pct, nct)
	finalized = []
	for r in read_results:
		finalized.append(read_results[r])
	#for read in read_results:
	#	print(read, *read_results[read])
		
	return finalized
					
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
				positive_set = labels[genome]["pos"]
				negative_set = labels[genome]["neg"]
				one_read_ext(read, positive_set, negative_set)

	for p in mn.negative:
		for read in mn.tagged_reads_neg[p]:
			genome = os.path.basename(read)
			read_len = int(genome.split("_read_len_")[1].split("_tagged.fasta")[0])
			if read_len not in seen_genomes:
				seen_genomes[read_len] = []
			if genome not in seen_genomes[read_len]: #Each genome only needs seen once, or its reads are effectively duplicates.
				seen_genomes[read_len].append(genome)
				genome = genome.split("_read_len_")[0]
				positive_set = labels[genome]["pos"]
				negative_set = labels[genome]["neg"]
				one_read_ext(read, positive_set, negative_set)

def process_alignments(labels, mn):
	results = []
	seen_genomes = {}
	for p in mn.positive:
		for read in mn.alignments_pos[p]:
			genome = os.path.basename(read)
			read_len = int(genome.split("_read_len_")[1].split("_aligned_reads.blast.txt")[0])
			if read_len not in seen_genomes:
				seen_genomes[read_len] = []
			if genome not in seen_genomes[read_len]: #Each genome only needs seen once
				seen_genomes[read_len].append(genome)
				genome = genome.split("_read_len_")[0]
				positive_set = labels[genome]["pos"]
				negative_set = labels[genome]["neg"]
				next_group = one_read_int(read, positive_set, negative_set)
				results.append((read, next_group))
		
	for p in mn.negative:
		for read in mn.alignments_neg[p]:
			genome = os.path.basename(read)
			read_len = int(genome.split("_read_len_")[1].split("_aligned_reads.blast.txt")[0])
			if read_len not in seen_genomes:
				seen_genomes[read_len] = []
			if genome not in seen_genomes[read_len]: #Each genome only needs seen once, or its reads are effectively duplicates.
				seen_genomes[read_len].append(genome)
				genome = genome.split("_read_len_")[0]
				positive_set = labels[genome]["pos"]
				negative_set = labels[genome]["neg"]
				next_group = one_read_int(read, positive_set, negative_set)
				results.append((read, next_group))
	
	return results
		
#This will come from main.
def options():
	pass
					
def output_reads(prjdir, external = True):
	read_labeller, manager = trawl(prjdir)
	if external:
		process_reads(read_labeller, manager)
		prep = None
	else:
		prep = process_alignments(read_labeller, manager)

	return prep, manager
	
#dir = sys.argv[1]
#output_reads(dir, False)