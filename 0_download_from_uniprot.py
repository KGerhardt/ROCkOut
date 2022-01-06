import sys
import os
import multiprocessing
#Pip install needed
import requests

'''
This script is a revision on the original ROCker downloading functionality.
'''

class download_manager:
	def __init__(self, prot = None, outdir = None, my_index = 1, max_index = 1):
		self.output = outdir
		self.prot = prot
		self.index = my_index
		self.total = max_index

	def check_annotation_information(self):
		#prot = arg[0]
		#output_dir = arg[1]
		url_base = "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch/uniprotkb/"
		url_base += self.prot
		url_base += "/annot"
		
		print("Downloading", url_base, "(Protein " + str(self.index), "of", str(self.total)+", step 1 of 3)")
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
		
		'''
		#We need to collect the on-target AA seqs for only the positive prots
		
		#Top lines of the annot file
		ID   ERMA_STAAN              Reviewed;         243 AA.
		AC   P0A0H2; P06699;
		DT   01-JAN-1988, integrated into UniProtKB/Swiss-Prot.
		DT   01-JAN-1988, sequence version 1.
		DT   29-SEP-2021, entry version 97.
		DE   RecName: Full=rRNA adenine N-6-methyltransferase;
		DE            EC=2.1.1.184;
		DE   AltName: Full=Erythromycin resistance protein;
		DE   AltName: Full=Macrolide-lincosamide-streptogramin B resistance protein;
			

		#Name to match 
		>UNIPROT:ERMA_STAAN P0A0H2 rRNA adenine N-6-methyltransferase (2.1.1.184) (Erythromycin resistance protein) (Macrolide-lincosamide-streptogramin B resistance protein)
		
		#Sequence will be buried somewhere in the GFFs and we'll have to check for the name, then the seq is under translation=...
		
		'''
		
		primary_target_name = ">UNIPROT:"
		gff_search_name = ""
		
		
		#We'll have to get the proteins listed under this annotation; this stores those prots
		my_translations = {}
		#Look into these
		additionals = []
		#Parse the file
		for line in annot_data:
			if line.startswith("ID"):
				segs = line.strip().split()
				gff_search_name = segs[1]
				primary_target_name += segs[1]
				primary_target_name += " "
				primary_target_name += self.prot
			if line.startswith("DE"):
				segs = line.strip().split()
				if segs[1] == "RecName:":
					usable = ' '.join(segs[2:])
					usable = usable[5:(len(usable)-1)]
					usable = "(" + usable + ")"
					primary_target_name += " "
					primary_target_name += usable
				if segs[1] == "AltName:":
					usable = ' '.join(segs[2:])
					usable = usable[5:(len(usable)-1)]
					usable = "(" + usable + ")"
					primary_target_name += " "
					primary_target_name += usable
					
				if segs[1].startswith("EC="):
					usable = segs[1][3:(len(segs[1]) - 1)]
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
		
		gffs_name = os.path.normpath(self.output + self.prot + "/gffs") + "/"
		if not os.path.exists(gffs_name):
			os.mkdir(gffs_name)
			
		coords_name = os.path.normpath(self.output + self.prot + "/coords") + "/"
		if not os.path.exists(coords_name):
			os.mkdir(coords_name)
			
		genomes_name = os.path.normpath(self.output + self.prot + "/genomes") + "/"
		if not os.path.exists(genomes_name):
			os.mkdir(genomes_name)
		
		for p in additionals:
			
			one_file = self.download_gffs(p, gff_search_name)
			formatted_results = []
			successful_tags = []
			
			fh = open(coords_name + p +"_coords.txt", "w")
			for item_set in one_file:
				if item_set[0] in additionals or item_set[1] in my_translations:
					print(*item_set, file = fh)
					successful_tags.append(p)
					
			fh.close()
			
			#Get the FNA from the sequences
			print("    Collecting reference genomes for protein " + str(self.index), "of", str(self.total)+" (step 3 of 3)")
			for tag in successful_tags:
				try:
					sequence = requests.get("https://www.ebi.ac.uk/Tools/dbfetch/dbfetch/embl/"+tag+"/fasta")
					sequence = sequence.content.decode()
					fh = open(genomes_name + tag +"_fasta.txt", "w")
					fh.write(sequence)
					fh.close()
				except:
					pass
			print("    Protein", str(self.index), "complete!")
		
		return None
		
	def download_gffs(self, prot, target):
		#self.output + self.prot
		url_base = "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch/embl/"
		url_base += prot
		url_base += "/gff3"
		
		print("    Collecting additional gene annotations for protein " + str(self.index), "of", str(self.total)+" (step 2 of 3)")
		
		gff_data = requests.get(url_base).content.decode()
		
		output_gff = self.output + self.prot + "/gffs/" + prot + ".gff3"
		
		fh = open(output_gff, "w")
		fh.write(gff_data)
		fh.close()
		
		gff_data = gff_data.split("\n")
		
		#here's where we search for the target
		
		
		coordinate_info = self.gff_to_coords(gff_data)
		
		return coordinate_info
		
	def check_and_extract_primary_target_seq(self):
		target_seq = None
		
		#If seq not present, then none and we handle that up above; else the sequence was found and we want it.
		return target_seq
		
	#Assumes a downloaded GFF3 file that has been processed by download_gffs so that it is a list of strings with each item in the list corresponding to a line in the original GFF file.	
	def gff_to_coords(self, gff_data):
		results = []
		
		for line in gff_data:
			#Header lines. Skip em.
			if line.startswith("#"):
				continue
			segs = line.split("\t")
			#Skip line if the record is incomplete.
			if len(segs) < 9:
				continue
				
			#Set these to none; they'll be filled or the none value will be used as a check.
			prot_id, tran_id = None, None
			#The annotation chunk is a semicolon separated line of unknown length.	
			annotation = segs[8].split(";")
			for item in annotation:
				#The labels seem very backwards to me, but they match the original code.
				if item.startswith("protein_id="):
					#Tran ID is always the first thing and always starts with 'protein_id='. This extracts the segment after 'ID='
					tran_id = item[11:]
				if item.startswith("db_xref=UniProtKB"):
					prot_id = item.split(":")[1]

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
			
			results.append((prot_id, tran_id, from_coord, to_coord, strand))
			
			#my_prot_ids.append(prot_id)
			#my_tran_ids.append(tran_id)
			#my_start_coordinates.append(from_coord)
			#my_end_coordinates.append(to_coord)
			#my_strands.append(strand)
			#return [formatted, prot_id, tran_id, from_coord, to_coord, strand]
				
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
			dl = download_manager(prot, self.posdir, my_index = count, max_index = total)
			self.worklist.append(dl)
			count += 1
		for prot in self.negative_set:
			dl = download_manager(prot, self.negdir, my_index = count, max_index = total)
			self.worklist.append(dl)
			count += 1
		
	def execute_downloads(self):
		pool = multiprocessing.Pool(self.threads)
		
		pool.map(do_download, self.worklist)
		
		pool.close()
		pool.join()
		

#This is effectively the main function here. There's probably a better way to package all of this, but it works.
pos_list = sys.argv[1]
neg_list = sys.argv[2]
threads = sys.argv[3]

dl = uniprot_downloader(pos_list, neg_list, threads, "rdl")
dl.prepare_directories()
dl.prepare_downloaders()
dl.execute_downloads()


