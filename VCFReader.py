#!/usr/bin/python

"""
VCFReader.py

A module that operates for VCF files

Author : Thomas Sunghoon Heo

Usage below

from VCFReader import VCFReader

vcffile_name = "some_vcf_file.vcf"
vcf_reader = VCFReader(vcffile_name)
vcf_reader.open()
for record in vcf_reader.parse():
	# Do something with VCF record using methods of VCFRecord class
vcf_reader.close()
"""

import gzip
def format_checker(fqname) :
	F = open(fqname,"rb")
	bytes = F.read(3)
	if bytes == "\x1f\x8b\x08" :
		F.seek(0)
		return gzip.GzipFile(fqname , fileobj=F)
	else:
		F.seek(0)
		return F

class VCFRecord:
	def __init__(self, line, lin_num=None):
		data = line.rstrip().split("\t")
		self.chrom = data[0]
		self.pos = int(data[1])
		self.id = data[2]
		self.ref = data[3]
		self.alt = data[4]
		self.qual = data[5]
		self.filter = data[6]
		self.info = {}
		num_info = len(data[7:])
		self.num_samples = num_info
		self.line_no = lin_num
		for i in range(num_info):
			self.info[i+1] = {}
			info = data[7+i]
			semi_colons = info.split(";")
			for item in semi_colons:
				x = item.split("=") ## VCF file info are separated by equality sign(=), and those are unique
				if len(x) == 2:
					k = x[0]
					v = x[1]
					self.info[i+1][k] = v
				elif len(x) == 1:
					self.info[i+1][x[0]] = None
				else:
					# Scan for the first = letter
					xloc = item.find('=')
					k = ''
					v = ''
					for j in range(xloc):
						k += item[j]
					for j in range(xloc+1 , len(item)):
						v += item[j]
					self.info[i+1][k] = v
	def __str__(self):
		array = [self.chrom, str(self.pos), self.id, self.ref, self.alt, self.qual, self.filter]
		for i in range(self.num_samples):
			mades = []
			for k, v in self.info[i+1].items():
				if v == None:
					make = k
				else:
					made = k + '=' + v
				mades.append(made)
			mades = ";".join(mades)
			array.append(mades)
		return "\t".join(array)

	def get_chrom(self):
		return self.chrom
	def get_pos(self):
		return self.pos
	def get_id(self):
		return self.id
	def get_ref(self):
		return self.ref
	def get_alt(self):
		return self.alt
	def get_qual(self):
		return self.qual
	def get_filter(self):
		return self.filter
	def get_variant_error_prob(self):
		if self.qual == '.':
			return None
		return 10 ** ((-1)*float(self.qual))
	def get_whole_info(self):
		return self.info
	def get_sample_info(self, index):
		return self.info[index]
	def get_num_samples(self):
		return self.num_samples
	def get_sample_tag_info(self, index, tag):
		try:
			return self.info[index][tag]
		except KeyError:
			if self.line_no == None:
				raise KeyError("No such tag in this VCF record\n"%(tag))
			else:
				raise KeyError("No such tag %s in this VCF record at VCF line %d\n"%(tag, self.line_no))
class VCFReader:
	def __init__(self, fname=None):
		self.fname = fname
		self.fptr = None
		self.headers = []
		self.chrom_line = None
	def open(self, fname=None):
		if self.fname == None:
			if fname == None:
				raise AttributeError("File name is not specified")
			else:
				self.fname = fname
		self.fptr = format_checker(self.fname)
		
	def close(self):
		self.fptr.close()
	def rewind(self):
		self.fptr.seek(0)

	def readline(self):
		self.read_headers()
		line = self.fptr.readline()
		return VCFRecord(line)

	def read_headers(self):
		for line in self.fptr:
			if line[0] == '#':
				if line.statswith('#CHROM'):
					if self.chrom_line != None:
						continue
					else:
						self.chrom_line = line.rstrip()
				else:
					if len(self.headers) != 0 :
						continue
					else:
						self.headers.append(line.rstrip())
	def parse(self):
		line_no = 0
		for line in self.fptr:
			if line[0] == '#':
				if line.startswith('#CHROM'):
					self.chrom_line = line.rstrip()
					continue
				else:
					self.headers.append(line.rstrip())
			else:
				line_no += 1
				yield VCFRecord(line, lin_num=line_no)
