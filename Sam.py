#!/usr/bin/python

"""
Sam.py

Stores Sam record and offers some manipulation functions

This only works for FASTQ read mapped SAM file

Author : Sunghoon Heo
"""

import re
class SamHeader:
	def __init__(self):
		self.headers = []
	def add(self, header):
		self.headers.append(header)
class SamRecord:
	def __init__(self, line):
		if "\n" in line : line = line[:-1]
		data = line.split("\t")
		self.qname = data[0]
		self.flag = int(data[1])
		self.rname = data[2]
		self.pos = int(data[3])
		self.mapq = int(data[4])
		self.cigar = data[5]
		self.rnext = data[6]
		self.pnext = int(data[7])
		self.tlen = int(data[8])
		self.seq = data[9]
		self.qual = data[10]
		self.aux = {}
		for d in data[11:]:
			k = d.split(":")[0]
			v = d.split(":")[-1]
			self.aux[k] = v

	def cigar_alphas(self):
		return re.split("[0-9]+" , self.cigar)[1:]
	def cigar_nums(self):
		X = re.split("[MINDNSHPX=]+" , self.cigar)[:-1]
		M = map(lambda x : int(x) , X)
		return M
		
	def get_tag(self, tag):
		try:
			return self.aux[k]
		except KeyError:
			return None
		
class Sam:
	def __init__(self , fname):
		self.fname = fname
		self.file_ptr = None
		self.headers = SamHeader()

	def open(self):
		self.file_ptr = open(self.fname)

	def close(self):
		self.file_ptr.close()

	def get_headers_of_sam(self):
		return self.headers.headers

	def parse_record(self):
		for line in self.file_ptr:
			if line.startswith("@") : 
				self.headers.add(line[1:-1])
			else:
				yield SamRecord(line)

	
		
