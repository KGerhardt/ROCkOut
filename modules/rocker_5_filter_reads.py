import sys
import os
import multiprocessing
import subprocess


def import_model(file):
	multiple_alignment = None
	
	ma_prots = []
	ma_seqs = []
	model = {}
	
	in_model = True
	
	fh = open(file)
	
	for line in fh:
		if line.startswith("position_in_protein"):	
			continue
		if line.startswith("#####"):
			in_model = False
			continue
			
		if in_model:
			segs = line.strip().split()
			#Here's the components:
			midpoint = int(segs[0])
			if midpoint not in model:
				model[midpoint] = {}
				
			readlen = int(segs[1])				
			discriminant = float(segs[2])
			
			model[midpoint][readlen] = discriminant
			
		if not in_model:
			if line.startswith(">"):
				ma_prots.append(line.strip()[1:])
			else:
				ma_seqs.append(line.strip())
	
	fh.close()
	
	if len(ma_seqs) > 0:
		multiple_alignment = dict(zip(ma_prots, ma_seqs))
	
	return model, multiple_alignment
	
def convert_legacy_model(old, new):
	vers_info = []
	multiple_alignment = {}
	model = {}
	
	in_model = False
	in_MA = False
	
	fh = open(file)
	
	for line in fh:
		#Comments at top - skip
		if not in_model:
			if line.startswith("#"):
				vers_info.append(line.strip())
				continue
			else:
				in_model = True
		
		#This isn't part of the description or multiple alignment, e.g. it is a rocker model's direct info.
		if line.startswith("#:>"):
			in_MA = True
			
		if in_MA:
			#Is MA sequence name
			if line.startswith("#:>"):
				#Remove the identifier
				name = line.strip()[2:]
				multiple_alignment[name] = ""
				continue
			#Is MA seqline
			if line.startswith("#:"):
				#Remove the identifier
				seq = line.strip()[2:]
				#Name will always be defined, here
				multiple_alignment[name] += seq
				continue
		else:
			segs = line.strip().split()
			#Here's the components:
			start = int(segs[0])
			end = int(segs[1])
			#Maybe this is count of bases aligned?
			count = int(segs[2])
			#Smaller integers. Ask Miguel.
			unknown = int(segs[3])
			
			discriminant = float(segs[4])

	fh.close()
	
	return model, multiple_alignment, vers_info
	
#Takes a blast input
def filter_blast(blast_in, filt_blast_out, model):
	
	fh = open(blast_in)
	out = open(filt_blast_out, "w")
	
	for line in fh:
		segs =  line.strip().split("\t")
		query = segs[0]
		target = segs[1]
		
		#Skip if the target has no hit.
		if target not in model:
			continue
		
		#Locations of aln in ref.
		from_bp = int(segs[8])
		to_bp = int(segs[9])
		
		actual_readlen = from_bp - to_bp + 1
		
		bits = float(segs[11])
		
		on_target = query.split(";")[5]
		
		midpt = int((from_bp + to_bp) / 2.0)
		
		if midpt in model:
			if actual_readlen in model[midpt]:
				if bits >= model[midpt][actual_readlen]:
					print(query, target, from_bp, to_bp, bits, on_target, midpt)
	
	out.close()
	fh.close()
	
f = sys.argv[1]
rocker_model = import_model(f)

'''
	#This blurb was for emulating the original rocker's lencorr function, but it's no longer needed.
	#Default case: read is long enough, and so we aren't penalizing it.
	if actual_readlen >= exp_readlen:
		corrected_bit_score = bits
	else:
		bits_per_aa = bits/actual_readlen
		#This happens when the read is shorter than "expected" readlen AND we have a lencorr value AND the penalty is high.
		if len_corr_penalty > 1.0:
		#Figure out how many extra bases are needed
		extra = max([0, (actual_readlen * (max_correction + 1.0) - exp_readlen)])
		max_tri = max_correction * actual_readlen * bits_per_aa / 2
		if max_correction > 0:
			tanTheta = bits_per_aa/(max_corr*actual_readlen)
		else:
			tanTheta = 0
		extra_tri = tanTheta * extra * extra / 2
		corrected_bit_score = bits + (max_tri - extra_tri)
		#Shorter than expected readlen AND we have a length correction AND the penalty is normal.
		else:
		miss = min([exp_readlen - actual_readlen, max_correction * actual_readlen])
		corrected_bit_score = bits + (bits_per_aa * miss * (1.0-penalty))
	
	#The actual filtering decision: print if the corrected bs is greater than the model's bit score threshold
	#at the midpoint of the subject read.
	if corrected_bit_score > model[midpt]:
		print(line)
	
	Original ruby bitscore code
	
	Expected readlen is defined by the model, e.g. the 125 at the header of the AmoA model.
	exp_readlen = data.signatures[:l].to_i
	
	o[:lencorr_max] ||= 0.4
	
	def correct_bs(bh, readlen, exp_readlen, max_corr, penalty)
	bs = bh.bits
	return bs if @o[:lencorr].nil? or readlen.nil? or readlen >= exp_readlen
	bits_per_aa = bs.to_f / readlen
	if penalty.nil? or penalty > 1.0
	  extra = [0.0, readlen * (max_corr + 1.0) - exp_readlen].max
	  max_tri = max_corr * readlen * bits_per_aa / 2
	  tanTheta = max_corr > 0.0 ? bits_per_aa / (max_corr * readlen) : 0.0
	  extra_tri = extra * extra * tanTheta / 2
	  return bs + (max_tri - extra_tri)
	else
	  miss = [exp_readlen - readlen, max_corr * readlen].min
	  return bs + (bits_per_aa * miss * (1.0 - penalty))
	end
	  end
	
'''
	
