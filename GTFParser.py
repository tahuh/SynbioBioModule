#!/usr/bin/python

"""
GTFParser.py

A parser for GTF file

Author : Sunghoon Heo
"""

import gzip

class GTFContinaer:
	def __init__(self, gtf_line):
		### GTF file fields
		self.chrom = None ## self.seqname = None
		self.source = None
		self.feature = None
		self.start = None
		self.end = None
		self.score = None
		self.strand = None
		self.frame = None
		self.attribute = None
		self.process_gtf_line(gtf_file)

	def process_gtf_line(self , l):
		if "\n" in line:
			data = line[:-1].split("\t")
		else:
			data = line.tstrip().split("\t")
		self.chrom = data[0]
		self.source = data[1]
		self.feature = data[2]
		self.start = int(data[3])
		self.end = int(data[4])
		self.score = float(data[5])
		self.strand = data[6]
		self.frame = int(data[7])
		self.attribute = {}
		tmp = data[8].split(";")
		for t in tmp:
			if t == '' : continue
			tt = t.split(" ")
			key = tt[0]
			value = tt[1].replace('"' , '')
			self.attribute[key] = value

	def gene_id(self):
		try:
			return self.attribute["gene_id"]
		except:
			return None
	def gene_name(self):
		try:
			return self.attribute["gene_name"]
		except:
			return None

	def transcript_name(self):
		try:
			return self.attribute["transcript_name"]
		except:
			return None

	def transcript_id(self):
		try:
			return self.attribute["transcript_id"]
		except:
			return None

	def gene_biotype(self):
		try:
			return self.attribute["gene_biotype"]
		except :
			return None

	def gene_version(self):
		try:
			return self.attribute["gene_version"]
		except:
			return None
class GTFParser:
	def __init__(self, fname):
		self.fname = fname
		self.fh = None
		self.lines_read = 0
		self.bytes_read = []
	def open(self):
		if self.fname.endswith(".gz"):
			self.fh = gzip.open(self.fname)
		else:
			self.fh = open(self.fname)

	def readline(self):
		line = self.fh.readline()
		self.lines_read += 1
		self.bytes_read.append(len(line))
		return GTFContiner(line)

	def jump_nlines_forward(self, n):
		if n >= self.lines_read:
			self.fh.seek(0)
		else:
			tmp = self.bytes_read[::-1]
			byte = 0
			for i in range(n):
				byte += tmp[i]
			self.fh.seek(byte , 1)
			self.lines_read -= n
	def jump_nlines_backward(self, n):
		for i in range(n):
			line = self.fh.readline()
			if line == '' : break ## EOF
			self.lines_read += 1
			self.bytes_read.append(len(line))
	def close(self):
		self.fh.close()
