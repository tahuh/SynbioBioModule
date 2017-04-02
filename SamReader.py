#!/usr/bin/python

class sub_recorder(object):
	def __init__(self , data):
		self.query = data[0]
		self.flag = data[1]
		self.rname = data[2]
		self.pos = int(data[3])
		self.mapq = int(data[4])
		self.cigar = data[5]
		self.rnext = data[6]
		self.pnext = data[7]
		self.slen = int(data[8])
		self.seq = data[9]
		try:
			self.qual  = data[10]
		except IndexError:
			self.qual = None
		try:
			self.misc = data[11]
		except IndexError:
			self.misc = None
	def build_sam_block(self):
		
		forms = []
		forms.append(self.query)
		forms.append(self.flag)
		forms.append(self.rname)
		forms.append(self.pos)
		forms.append(self.mapq)
		forms.append(self.cigar)
		forms.append(self.rnext)
		forms.append(self.pnext)
		forms.append(self.slen)
		forms.append(self.seq)
		if not self.qual:
			forms.append(self.qual)
		if not self.misc:
			forms.append(self.misc)
		return "\t".join(map(lambda x : str(x) ,forms))
class SamReader(object):

	def __init__(self, samfile) :
		self.sam = samfile
		self.headers = []
		self.header_read = False
	def open(self):
		self.sam_file = open(self.sam)

	def read_header(self):
		for line in self.sam_file:
			if line.startswith("@"):
				self.headers.append(line.rstrip())
				continue
			else:
				self.sam_file.seek(0)
				break
		self.header_read = True

	def read_record(self):
		for line in self.sam_file:
			
			if line.startswith("@") : 
				if self.header_read == False:
					self.headers.append(line.rstrip())
					continue
				else:
					continue

			data_field = {}
			data = line.rstrip().split("\t")

			rec = sub_recorder(data)

			yield rec

	def show_header(Self):
		if self.headers == [] :
			raise AttributeError("No header read. Please use read_header method")
		else:
			return self.headers

	def rewind(self):
		self.sam_reader.seek(0)

	def close(self):
		self.sam_file.close()

	def cigar_alpha(self,cigar):
		alphas = re.split("[MIDNSHPX=]+" , cigar)[:-1]
		return alphas

	def cigar_chars(self, cigar):
		chars = re.split("[0-9]+" , cigar)[1:]
		return chars
