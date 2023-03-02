import sys
import os
import multiprocessing
#Pip install needed
import requests
import shutil

import argparse

'''
This script is a revision on the original ROCker downloading functionality.
'''

def reverse_complement(seq):
	translate = {"A":"T", "T":"A", "C":"G", "G":"C"}
	#reverse
	seq = seq[::-1]
	#complement
	seq = ''.join([translate[nt] for nt in seq])
	return seq

#This was a sanity checking function
def print_comparison(nt_seq, aa_seq):
	for i in range(0, len(aa_seq)):
		codon = nt_seq[(i*3):((i+1)*3)]
		print(codon, aa_seq[i])

class download_manager:
	def __init__(self, prot = None, outdir = None, my_index = 1, max_index = 1, positive = False):
		self.output = os.path.normpath(outdir) + "/"
		self.prot = prot
		self.index = my_index
		self.total = max_index
		self.is_pos = positive

	def check_annotation_information(self):
		#prot = arg[0]
		#output_dir = arg[1]
		url_base = "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch/uniprotkb/"
		url_base += self.prot
		url_base += "/annot"
		
		print("Downloading", url_base, "(Protein " + str(self.index), "of", str(self.total)+")")
		#Download the annotation
		try:
			annot_data = requests.get(url_base)
			if not os.path.exists(self.output + self.prot):
				os.mkdir(self.output + self.prot)
		except:
			#If the data cannot be accessed, break
			pass
			
		#Extract the text in binary; decode to string
		annot_data = annot_data.content.decode()
				
		#Save the annotation
		annot_name = os.path.normpath(self.output + self.prot + "/initial_annotation") + "/"
		if not os.path.exists(annot_name):
			os.mkdir(annot_name)
		
		fh = open(annot_name + self.prot +"_uniprot_annotation.txt", "w")
		fh.write(annot_data)
		fh.close()
		#Split the annotation into single lines
		annot_data = annot_data.split('\n')
		
		#primary_target_name = ">UNIPROT__"
		primary_target_name = [">UNIPROT__", []]

		gff_search_name = ""
		
		#We'll have to get the proteins listed under this annotation; this stores those prots
		my_translations = {}
		#Look into these
		additionals = []
		#Parse the file
		for line in annot_data:
			if line.startswith("ID"):
				#Get the primary name for this sequence. Matches this format:
				#ID   ERMC1_STAAU             Reviewed;         244 AA.
			
				segs = line.strip().split()
				primary_target_name[0] += segs[1]
				primary_target_name[0] += "__"
				primary_target_name[0] += self.prot
				primary_target_name[0] += "__"
				
			if line.startswith("DE"):
				segs = line.strip().split()
				if segs[1] == "RecName:":
					#Python is seeing something in this format from segs:
					#['DE', 'RecName:', 'Full=rRNA', 'adenine', 'N-6-methyltransferase;']
					#We need to keep an extra copy of the label for searching GFFs later
					
					#Remove first 2 items
					usable = ' '.join(segs[2:])
					#Remove "Full=..." and final semicolon
					usable = usable[5:(len(usable)-1)]
					
					gff_search_name = usable
					
					usable = "(" + usable + ")"
					#primary_target_name += " "
					#primary_target_name += usable
					primary_target_name[1].append(usable)
				if segs[1] == "AltName:":
					#Same as above, but we don't need to separately keep the name
					#Matches this format:
					#DE   AltName: Full=Erythromycin resistance protein;
				
					usable = ' '.join(segs[2:])
					usable = usable[5:(len(usable)-1)]
					usable = "(" + usable + ")"
					#primary_target_name += " "
					#primary_target_name += usable
					primary_target_name[1].append(usable)
					
				if segs[1].startswith("EC="):
					#Same as alt name, modified to match this format:
					#DE            EC=2.1.1.184;
					#There may be spaces in the EC...; section.
					#combine
					usable = ' '.join(segs[1:])
					usable = usable[3:(len(usable) - 1)]
					usable = "(" + usable + ")"
					#primary_target_name += " "
					#primary_target_name += usable
					primary_target_name[1].append(usable)
		
			if line.startswith("DR   EMBL;") and line.endswith("Genomic_DNA."):
				segs = line.split("; ")
				further_invest = segs[1]
				additionals.append(further_invest)
				protein_name = segs[2]
				protein_strand = segs[3]
				my_translations[protein_name] = protein_strand
		
		additionals = list(set(additionals))
		
		primary_target_name[1] = ' '.join(primary_target_name[1])
		
		protein_listings = []
		
		#Create per-protein directories to organize sub-outputs
		gffs_name = os.path.normpath(self.output + self.prot + "/gffs") + "/"
		if not os.path.exists(gffs_name):
			os.mkdir(gffs_name)
			
		coords_name = os.path.normpath(self.output + self.prot + "/coords") + "/"
		if not os.path.exists(coords_name):
			os.mkdir(coords_name)
			
		proteins_name = os.path.normpath(self.output + self.prot + "/proteome") + "/"
		if not os.path.exists(proteins_name):
			os.mkdir(proteins_name)
			
		genomes_name = os.path.normpath(self.output + self.prot + "/genomes") + "/"
		if not os.path.exists(genomes_name):
			os.mkdir(genomes_name)
			
		print("    Collecting additional gene information for protein " + str(self.index), "of", str(self.total)+".", str(len(additionals)), "total sources to peruse.")

		winning_translation = ""
		winning_translation_info = None
		
		#wts = []
		wt_info = []
		#tags = []
		
		wts = {}
		#Nucleotide seqs
		nts = {}
		
		#Process identified proteins by downloading the relevant GFF and parsing them for useful results.
		for p in additionals:
			one_file = self.download_gffs(p, gff_search_name)
			formatted_results = []
			successful_tags = []
			gen_tags = {}
			
			if len(one_file) > 0: #Don't make a coords file for an empty gff
				#Record coordinates for use later.
				fh = open(coords_name + p +"_coords.txt", "w")
				for item_set in one_file:
					#print(item_set[0:4])
					#item_set[2] is trans ID
					if item_set[0] in additionals or item_set[2] in my_translations:
						#print(item_set)
						tag = item_set[2]
						wts[tag] = {}
						
						if p not in gen_tags:
							gen_tags[p] = [tag]
						else:
							gen_tags[p].append(tag)
						
						this_trans = item_set[len(item_set)-1]
						
						#wts[tag]["seq_prot"] = this_trans
						
						to_add = list(item_set[:-1]) #all but the actual sequence
						to_add.append(p)
						to_add.append(gff_search_name)
						
						#print(to_add)
						wts[tag]['info'] = to_add
						wts[tag]['seq'] = this_trans
						#wt_info.append(to_add)
						#tags.append(to_add[5])
						
						'''
						if len(this_trans) > len(winning_translation):
							winning_translation = this_trans
							winning_translation_info = list(item_set[:-1])
							winning_translation_info.append(p)
							winning_translation_info.append(gff_search_name)
						'''
						
						no_trans = item_set[:-1]
						print(*no_trans, file = fh)
						successful_tags.append(p)
						
				fh.close()
			
			successful_tags = list(set(successful_tags))
			
			
			#Get the FNA from the sequences
			for protein in successful_tags:
				try:
					sequence = requests.get("https://www.ebi.ac.uk/Tools/dbfetch/dbfetch/embl/"+protein+"/fasta")
					sequence = sequence.content.decode()
			
					fh = open(genomes_name + protein +".fasta", "w")
					fh.write(sequence)
					fh.close()
					
					sequence = sequence.split("\n")
					to_combine = []
					for line in sequence:
						if not line.startswith(">"):
							to_combine.append(line)
					sequence = None
					to_combine = ''.join(to_combine)
					
					tt = None
					for tag in gen_tags[protein]:
						tt = tag
						if tag in wts:
							parent_name = wts[tag]['info'][0]
							gene_name = wts[tag]['info'][1]
							prot_name = wts[tag]['info'][2]
							#python will zero index, but the genome indices are 1 indexed
							start = int(wts[tag]['info'][3])-1
							#python selection will exclude the final index; leaving this as-is works
							end = int(wts[tag]['info'][4])
							strand = wts[tag]['info'][5]
							
							#infl = [parent_name, gene_name, prot_name, start, end, strand]
							
							nt_seq = to_combine[start:end]
							if strand == "-":
								nt_seq = reverse_complement(nt_seq)
								
							nts[tag] = nt_seq
							
				except:
					print("Could not download", "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch/embl/"+protein+"/fasta for protein", self.prot)
					print(tt, "was the failing tag")
					
			to_combine = None

		#Implicitly skipped when wts is empty
		if len(wts) > 0:
			target_dir_name = os.path.normpath(self.output + self.prot + "/target_protein") + "/"
			
			if not os.path.exists(target_dir_name):
				os.mkdir(target_dir_name)
				
			aa_output_name = target_dir_name + self.prot + "_target_protein_AA.fasta"
			
			#This isn't filling correctly - It's probably a tag issue
			nt_output_name = target_dir_name + self.prot + "_target_protein_nt.fasta"

			aaout = open(aa_output_name, "w")
			ntout = open(nt_output_name, "w")
			
			for tag in wts:
				current_protein = wts[tag]['info'][1]
				
				print_name = primary_target_name[0] + current_protein + " " + primary_target_name[1]
				#print_name = primary_target_name
				#print_name = print_name.format(this_protein = current_protein)
				trans = wts[tag]['seq']
				
				reformatted = [trans[i:i+60] for i in range(0, len(trans), 60)]
				reformatted = "\n".join(reformatted)

				print(print_name, file = aaout)
				print(reformatted, file = aaout)
				
				nt_seq = nts[tag]
				reformatted = [nt_seq[i:i+60] for i in range(0, len(nt_seq), 60)]
				reformatted = "\n".join(reformatted)
				print(print_name, file = ntout)
				print(reformatted, file = ntout)
			
			aaout.close()
			ntout.close()
			
		print("    Protein", str(self.index), "complete!")
			
		return None
		
	def download_gffs(self, prot, target):
		#self.output + self.prot
		url_base = "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch/embl/"
		url_base += prot
		url_base += "/gff3"
				
		gff_data = requests.get(url_base).content.decode()
		
		output_gff = self.output + self.prot + "/gffs/" + prot + ".gff3"
		
		coordinate_info = []
		if not gff_data.startswith("ERROR"):
			fh = open(output_gff, "w")
			fh.write(gff_data)
			fh.close()
			
			gff_data = gff_data.split("\n")
			
			output_protein = self.output + self.prot + "/proteome/" + prot + ".proteome.fasta"
			po = open(output_protein, "w")
			for line in gff_data:
				segs = line.strip().split("\t")
				if len(segs) > 7:
					annot = segs[8]
					if "translation=" in annot: #This is a line containing a protein record
						#parent = segs[0]
						start, end = segs[3], segs[4]
						strand = segs[6]
						annot = annot.split(";")
						id = annot[0][3:]
						trans_position = len(annot)-1
						
						trans = annot[trans_position][12:]
						reformatted = [trans[i:i+60] for i in range(0, len(trans), 60)]
						reformatted = "\n".join(reformatted)
						
						annot = annot[1:(trans_position-1)]
						annot = ";".join(annot)
						
						to_write = ">"+id+";"+start+";"+end+";"+strand + " " + annot +  "\n" + reformatted
						print(to_write, file = po)
					
			po.close()
					
			coordinate_info = self.gff_to_coords(gff_data, target)
		
		return coordinate_info
		
	def check_and_extract_primary_target_seq(self):
		target_seq = None
		
		#If seq not present, then none and we handle that up above; else the sequence was found and we want it.
		return target_seq
		
	#Assumes a downloaded GFF3 file that has been processed by download_gffs so that it is a list of strings with each item in the list corresponding to a line in the original GFF file.	
	def gff_to_coords(self, gff_data, target):
		results = []
		
		for line in gff_data:
			#Header lines. Skip em.
			if line.startswith("#"):
				continue
			segs = line.strip().split("\t")
			#Skip line if the record is incomplete.
			if len(segs) < 9:
				continue
				
			#Protein sequence for the target gene or none if it's not in the GFF or is a negative sequence that we don't need.
			target_sequence = None
								
			#Set these to none; they'll be filled or the none value will be used as a check.
			prot_id, tran_id = None, None
			trans = None
			gene_id = None
			#The annotation chunk is a semicolon separated line of unknown length.	
			annotation = segs[8].split(";")
			for item in annotation:
				
				if item.startswith("ID="):
					gene_id = item[3:]
			
				#The labels seem very backwards to me, but they match the original code.
				if item.startswith("protein_id="):
					#Tran ID is always the first thing and always starts with 'protein_id='. This extracts the segment after 'ID='
					tran_id = item[11:]
				if item.startswith("db_xref=UniProtKB"):
					prot_id = item.split(":")[1]
					
				if item.startswith("translation="):
					trans = item[12:]
					
			coords = [int(segs[3]), int(segs[4])]
					
			#Regardless of strand we go min->max as from->to
			from_coord = str(min(coords))
			to_coord   = str(max(coords))
			
			strand = segs[6]
			
			#Nothing to do here.
			if prot_id is None and tran_id is None:
				continue
			
			if prot_id is None:
				prot_id = "None"
			if tran_id is None:
				tran_id = "None"
			if gene_id is None:
				gene_id = "None"
			
			results.append((prot_id, gene_id, tran_id, from_coord, to_coord, strand, trans,))
				
		return results
		
def do_download(dl_man_obj):
	dl_man_obj.check_annotation_information()
		
class uniprot_downloader:
	def __init__(self, positives = None, negatives = None, threads = 1, outdir = "rocker_downloads"):
		self.outdir = os.path.normpath(outdir) + "/"
		self.posdir = os.path.normpath(outdir + "/" + "positive") + "/"
		self.negdir = os.path.normpath(outdir + "/" + "negative") + "/"
		self.positive_set = positives
		self.negative_set = negatives
		self.sets_from_files()
		
		self.worklist = []
		
		
		self.target_protein_files = None
		self.gffs = None
		self.coordinate_files = None
		
		try:
			threads = int(threads)
		except:
			print("Cannot parse threads. Setting to default of 1 thread.")
			threads = 1
			
		self.threads = threads
			
	def sets_from_files(self):
		try:
			prot_set = []
			fh = open(self.positive_set)
			for line in fh:
				prot_set.append(line.strip())
			fh.close()
			self.positive_set = prot_set
		except:
			print("ROCker needs a positive prot_set. Cannot continue.")
			return None
		
		try:
			if self.negative_set is not None:
				prot_set = []
				fh = open(self.negative_set)
				for line in fh:
					prot_set.append(line.strip())
				fh.close()
				self.negative_set = prot_set
		except:
			pass
			
	def prepare_directories(self):
		if not os.path.exists(self.outdir):
			os.mkdir(self.outdir)			
		if len(self.positive_set) > 0:
			if not os.path.exists(self.posdir):
				os.mkdir(self.posdir)
		if self.negative_set is not None:
			if len(self.negative_set) > 0:
				if not os.path.exists(self.negdir):
					os.mkdir(self.negdir)

	def prepare_downloaders(self):
		count = 1
		total = len(self.positive_set)
		if self.negative_set is not None:
			total += len(self.negative_set)
			
		for prot in self.positive_set:
			dl = download_manager(prot, self.posdir, my_index = count, max_index = total, positive = True)
			self.worklist.append(dl)
			count += 1
		if self.negative_set is not None:
			for prot in self.negative_set:
				dl = download_manager(prot, self.negdir, my_index = count, max_index = total, positive = False)
				self.worklist.append(dl)
				count += 1
		
	def execute_downloads(self):
		pool = multiprocessing.Pool(self.threads)
		pool.imap_unordered(do_download, self.worklist)
		#for result in pool.imap_unordered(do_download, self.worklist):
			#target_proteins
			#gffs
			
		pool.close()
		pool.join()
	
	#Look through outputs and see if anything is empty.
	def check_results(self):
		#log file
		if not os.path.exists(os.path.normpath(self.outdir+"/shared_files")):
			os.mkdir(os.path.normpath(self.outdir+"/shared_files"))
		if not os.path.exists(os.path.normpath(self.outdir+"/shared_files/downloads/")):
			os.mkdir(os.path.normpath(self.outdir+"/shared_files/downloads/"))
			
		errlog = open(os.path.normpath(self.outdir+"/shared_files/downloads/download_errlog.txt"), "w")
		
		for pos in self.positive_set:
			path = os.path.normpath(self.posdir + pos + "/genomes")
			try:
				genome_files = os.listdir(path)
			except:
				genome_files = []
			if len(genome_files) < 1:
				print(pos, "positive", "failed", sep = "\t", file = errlog)
				print("Uniprot ID", pos, "was empty!")
				print("This ID may have been changed or may be missing required data on Uniprot.")
				print("Removing this positive ID from the dataset.")
				try:
					shutil.rmtree(os.path.normpath(self.posdir + pos))
				except:
					pass
			else:
				print(pos, "positive", "succeeded", sep = "\t", file = errlog)
		
		if self.negative_set is not None:			
			for neg in self.negative_set:
				path = os.path.normpath(self.negdir + neg + "/genomes")
				try:
					genome_files = os.listdir(path)
				except:
					genome_files = []
				if len(genome_files) < 1:
					print(neg, "negative", "failed", sep = "\t", file = errlog)
					print("Uniprot ID", neg, "was empty!")
					print("This ID may have been changed or may be missing required data on Uniprot.")
					print("Removing this negative ID from the dataset.")
					try:
						shutil.rmtree(os.path.normpath(self.negdir + neg))
					except:
						pass
				else:
					print(neg, "negative", "succeeded", sep = "\t", file = errlog)
			
		errlog.close()

def options():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
			description='''	''')

	parser.add_argument('-p', '--positive',  dest = 'pos', default = None, help =  '')
	parser.add_argument('-n', '--negative',  dest = 'neg', default = None, help =  '')
	parser.add_argument('-t', '--threads',  dest = 'threads', default = 1, help =  '')
	parser.add_argument('-o', '--output',  dest = 'out', default = "rockout_project", help =  '')

	args, unknown = parser.parse_known_args()
	
	return parser, args

def download(parser, opts):
	multiprocessing.freeze_support()
	pos_list = opts.pos
	neg_list = opts.neg
	try:
		threads = int(opts.threads)
	except:
		print("Threads has to be an integer. Defaulting to 1 thread")
		threads = 1

	dirname = opts.dir
	
	if dirname is None:
		sys.exit("ROCkOut needs a directory name to place downloads in!")
	
	if pos_list is None:
		print("Rocker needs a positive list!")
		print(parser.print_help())
		quit()
	
	print("Beginning ROCkOut!")

	dl = uniprot_downloader(pos_list, neg_list, threads, dirname)
	dl.prepare_directories()
	dl.prepare_downloaders()
	dl.execute_downloads()
	dl.check_results()
	



