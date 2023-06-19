import urllib.request as urlreq
import plotly
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import sys
import os

#file = "multiple_alignment/complete_multiple_alignment_aa.fasta"

def read_fasta(file):
	cur_seq = ""
	cur_id = ""
	cur_def = ""
	max_seqlen = 0
	
	seq_dict = {}
	defline_dict = {}
	
	fh = open(file)
	
	for line in fh:
		line = line.strip()
		if line.startswith(">"):
			if len(cur_seq) > 0:
				max_seqlen = max([len(cur_seq), max_seqlen])
				seq_dict[cur_id] = cur_seq
				defline_dict[cur_id] = cur_def
				cur_seq = ""
			
			cur_def = line.split()
			#Get the seqid and remove '>'
			cur_id = cur_def[0][1:]
			
			#If there's additional stuff, record it.
			if len(cur_def) > 0:
				cur_def = " ".join(cur_def[1:])
		else:
			cur_seq += line	
	fh.close()
	
	#Final iteration doesn't happen otw
	if len(cur_seq) > 0:
		max_seqlen = max([len(cur_seq), max_seqlen])
		seq_dict[cur_id] = cur_seq
		defline_dict[cur_id] = cur_def
		cur_seq = ""
	
	return seq_dict, defline_dict, max_seqlen
	

def amino_acid_color_dict(list = 1):
	color_dict = {"wht":"#FFFFFF",
				"teal": "#50BBBB",
				"red":  "#F45044",
				"lav":  "#A0B8F4",
				"purp": "#D076D0",
				"grn":  "#40FF40",
				"ylw":  "#FFFF40",
				"sal":  "#F4AC76",
				"pink": "#F4A0A0"}
	
	aa_nums = [0., 0.05, 0.1 , 0.15, 0.2 , 0.25, 0.3 , 0.35, 0.4 , 0.45, 0.5 , 0.55, 0.6 , 0.65, 0.7 , 0.75, 0.8 , 0.85, 0.9 , 0.95, 1.  ]
	aa_list = ['-', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
	
	color_list_1 = ["wht", "lav", "pink", "purp", "purp", "lav", "sal", "teal", "lav", "red", "lav", "lav", 
					"grn", "ylw", "grn", "red", "grn", "grn", "lav", "lav", "teal"]
	
	final_cols = []
	
	'''
	plotly doc on colorscale
	colorscale – Sets the colorscale. The colorscale must be an array containing arrays mapping a normalized 
	value to an rgb, rgba, hex, hsl, hsv, or named color string. At minimum, a mapping for the lowest (0) 
	and highest (1) values are required. For example, [[0, 'rgb(0,0,255)'], [1, 'rgb(255,0,0)']]. To control 
	the bounds of the colorscale in color space, use zmin and zmax. Alternatively, colorscale may be a palette 
	name string of the following list: Blackbody,Bluered,Blues,C ividis,Earth,Electric,Greens,Greys,Hot,Jet,
	Picnic,Portl and,Rainbow,RdBu,Reds,Viridis,YlGnBu,YlOrRd.
	
	We're gonna make a color value for each 1-21 aa incl. '-' as a whitespace character. The z values will eventually use this
	'''
	
	if list == 1:
		for aa, num, col in zip(aa_list, aa_nums, color_list_1):
			hex = color_dict[col]
			#final_cols[aa] = [num, hex]
			final_cols.append([num, hex])
			
	return final_cols

	
def convert_seqs_to_xyz(seqs, labels, total_size):
	x = []
	for i in range(0, total_size):
		x.append(i)
	
	y = []
	z = []
	
	text = []
	row = 0

	aa_list = ['-', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
	aa_nums = [0., 0.05, 0.1 , 0.15, 0.2 , 0.25, 0.3 , 0.35, 0.4 , 0.45, 0.5 , 0.55, 0.6 , 0.65, 0.7 , 0.75, 0.8 , 0.85, 0.9 , 0.95, 1.  ]

	aa_to_num = dict(zip(aa_list, aa_nums))
	
	#for the z value, we need a colorscale match...
	
	for seqid in seqs:
		spacer_to_add = 0
		#Find the max seqlen.
		full_name = seqid + " " + labels[seqid]
		#y.append(full_name)
		y.append(seqid)
		cur_seq = seqs[seqid]
		cur_sl = len(cur_seq)
		
		split_seq = [aa for aa in cur_seq]
		
		if cur_sl < total_size:
			split_seq.extend(["-"]*(total_size - cur_sl))
			
		text.append(split_seq)
		
		next_z = [aa_to_num[aa] for aa in split_seq]
		
		z.append(next_z)
		
	
	return x, y, z, text
	

file = sys.argv[1]
output = sys.argv[2]
	
s, d, max_sl = read_fasta(file)

cols = amino_acid_color_dict()

x, y, z, t = convert_seqs_to_xyz(s, d, max_sl)

fig = go.Figure(data=go.Heatmap(
		x = x,
		y = y,
		z = z,
		text = t,
		#texttemplate="%{text}",
		colorscale = cols,
		showscale=False),
		)
		
if 75 > max_sl:
	start_zoom = max_sl
else:
	start_zoom = 75
	
fig.update_layout(	
	xaxis=dict(
			#range=[0, start_zoom]
		)
)
#fig.update_layout(dragmode='pan')
#fig.update_layout(selectdirection='h')

config = {
	'scrollZoom': False,
		'toImageButtonOptions': {
		'format': 'svg', # one of png, svg, jpeg, webp
		'filename': 'custom_image',
		'height': 1080,
		'width': 1920,
		'scale': 1 # Multiply title/legend/axis/canvas sizes by this factor
  }
}

#Takes config options. Write SVG here
fig.write_html(output, config = config)
