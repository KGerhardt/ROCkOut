class best_hit_filter:
	def __init__(self, file = None):
		self.input_reads = file
		
	def collect_bh(self):
		with open(self.input_reads) as fh:
			for line in fh:
				segs = line.strip().split("\t")
				read_id = segs[0]
				tgt = segs[1]
				bitscore = float(segs[11])