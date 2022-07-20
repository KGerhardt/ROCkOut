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
		
		primary_target_name = ">UNIPROT:"
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
				primary_target_name += segs[1]
				primary_target_name += " "
				primary_target_name += self.prot
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
					primary_target_name += " "
					primary_target_name += usable
				if segs[1] == "AltName:":
					#Same as above, but we don't need to separately keep the name
					#Matches this format:
					#DE   AltName: Full=Erythromycin resistance protein;
				
					usable = ' '.join(segs[2:])
					usable = usable[5:(len(usable)-1)]
					usable = "(" + usable + ")"
					primary_target_name += " "
					primary_target_name += usable
					
				if segs[1].startswith("EC="):
					#Same as alt name, modified to match this format:
					#DE            EC=2.1.1.184;
					#There may be spaces in the EC...; section.
					#combine
					usable = ' '.join(segs[1:])
					usable = usable[3:(len(usable) - 1)]
					usable = "(" + usable + ")"
					primary_target_name += " "
					primary_target_name += usable
		
			if line.startswith("DR   EMBL;") and line.endswith("Genomic_DNA."):
				segs = line.split("; ")
				further_invest = segs[1]
				additionals.append(further_invest)
				protein_name = segs[2]
				protein_strand = segs[3]
				my_translations[protein_name] = protein_strand
		
		additionals = list(set(additionals))
		
		protein_listings = []
		
		#Create per-protein directories to organize sub-outputs
		gffs_name = os.path.normpath(self.output + self.prot + "/gffs") + "/"
		if not os.path.exists(gffs_name):
			os.mkdir(gffs_name)
			
		coords_name = os.path.normpath(self.output + self.prot + "/coords") + "/"
		if not os.path.exists(coords_name):
			os.mkdir(coords_name)
			
		genomes_name = os.path.normpath(self.output + self.prot + "/genomes") + "/"
		if not os.path.exists(genomes_name):
			os.mkdir(genomes_name)
			
		print("    Collecting additional gene information for protein " + str(self.index), "of", str(self.total)+".", str(len(additionals)), "total sources to peruse.")

		winning_translation = ""
		
		#Process identified proteins by downloading the relevant GFF and parsing them for useful results.
		for p in additionals:
			
			one_file = self.download_gffs(p, gff_search_name)
			formatted_results = []
			successful_tags = []
			
			#Record coordinates for use later.
			fh = open(coords_name + p +"_coords.txt", "w")
			for item_set in one_file:
				if item_set[0] in additionals or item_set[1] in my_translations:
					this_trans = item_set[len(item_set)-1]
					if len(this_trans) > len(winning_translation):
						winning_translation = this_trans
					no_trans = item_set[:-1]
					print(*no_trans, file = fh)
					successful_tags.append(p)
					
			fh.close()
			
			#Get the FNA from the sequences
			for tag in successful_tags:
				try:
					sequence = requests.get("https://www.ebi.ac.uk/Tools/dbfetch/dbfetch/embl/"+tag+"/fasta")
					sequence = sequence.content.decode()
					fh = open(genomes_name + tag +"_fasta.txt", "w")
					fh.write(sequence)
					fh.close()
				except:
					print("Could not download", "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch/embl/"+tag+"/fasta for protein", self.prot)
		
		print("    Protein", str(self.index), "complete!")

		if len(winning_translation) == 0:
			winning_translation = None	
		#We only want the positive prots here
		if self.is_pos:
			#When the sequence is positive and we have a winning answer, we print the AA seq to a file.
			if winning_translation is not None:
				reformatted = [winning_translation[i:i+60] for i in range(0, len(winning_translation), 60)]
				reformatted = "\n".join(reformatted)
				target_dir_name = os.path.normpath(self.output + self.prot + "/target_protein") + "/"
				if not os.path.exists(target_dir_name):
					os.mkdir(target_dir_name)
					
				output_name = target_dir_name + self.prot + "_target_protein_AA.fasta.txt"
				out = open(output_name, "w")
				print(primary_target_name, file = out)
				print(reformatted, file = out)
				out.close()
			
		return None
		
	def download_gffs(self, prot, target):
		#self.output + self.prot
		url_base = "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch/embl/"
		url_base += prot
		url_base += "/gff3"
				
		gff_data = requests.get(url_base).content.decode()
		
		output_gff = self.output + self.prot + "/gffs/" + prot + ".gff3"
		
		fh = open(output_gff, "w")
		fh.write(gff_data)
		fh.close()
		
		gff_data = gff_data.split("\n")
		
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
			#The annotation chunk is a semicolon separated line of unknown length.	
			annotation = segs[8].split(";")
			for item in annotation:
				if self.is_pos:
					pass
					#print(item)
			
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
			
			results.append((prot_id, tran_id, from_coord, to_coord, strand, trans))
				
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
		if len(self.negative_set) > 0:
			if not os.path.exists(self.negdir):
				os.mkdir(self.negdir)

	def prepare_downloaders(self):
		count = 1
		total = len(self.positive_set) + len(self.negative_set)
		for prot in self.positive_set:
			dl = download_manager(prot, self.posdir, my_index = count, max_index = total, positive = True)
			self.worklist.append(dl)
			count += 1
		for prot in self.negative_set:
			dl = download_manager(prot, self.negdir, my_index = count, max_index = total, positive = False)
			self.worklist.append(dl)
			count += 1
		
	def execute_downloads(self):
		pool = multiprocessing.Pool(self.threads)
		pool.imap_unordered(do_download, self.worklist)
		pool.close()
		pool.join()
	
	#Look through outputs and see if anything is empty.
	def check_results(self):
		for pos in self.positive_set:
			path = os.path.normpath(self.posdir + pos + "/genomes")
			try:
				genome_files = os.listdir(path)
			except:
				genome_files = []
			if len(genome_files) < 1:
				print("Uniprot ID", pos, "was empty!")
				print("This ID may have been changed or may be missing required data on Uniprot.")
				print("Removing this positive ID from the dataset.")
				try:
					shutil.rmtree(os.path.normpath(self.posdir + pos))
				except:
					pass
		for neg in self.negative_set:
			path = os.path.normpath(self.negdir + neg + "/genomes")
			try:
				genome_files = os.listdir(path)
			except:
				genome_files = []
			if len(genome_files) < 1:
				print("Uniprot ID", neg, "was empty!")
				print("This ID may have been changed or may be missing required data on Uniprot.")
				print("Removing this negative ID from the dataset.")
				try:
					shutil.rmtree(os.path.normpath(self.negdir + neg))
				except:
					pass
	


	
def options():
	parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
			description='''	''')

	parser.add_argument('-p', '--positive',  dest = 'pos', default = None, help =  '')
	parser.add_argument('-n', '--negative',  dest = 'neg', default = None, help =  '')
	parser.add_argument('-t', '--threads',  dest = 'threads', default = 1, help =  '')
	parser.add_argument('-o', '--output',  dest = 'out', default = "rockout_project", help =  '')

	args, unknown = parser.parse_known_args()
	
	return parser, args

def main():
	multiprocessing.freeze_support()
	parser, opts = options()
	pos_list = opts.pos
	neg_list = opts.neg
	try:
		threads = int(opts.threads)
	except:
		print("Threads has to be an integer. Defaulting to 1 thread")
		threads = 1

	dirname = opts.out
	
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

	

if __name__ == "__main__":
	main()



