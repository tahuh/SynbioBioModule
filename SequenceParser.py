#!/usr/bin/python

"""
SequenceParser.py

This code parse FASTA/Q files for some benefit

Code only works for Illumina sequencing platform result

Author : Sunghoon Heo

How To Use

from SequenceParser import FQParser ## FASTQ
from SequenceParser import FAParser ## FASTA

fqname = "test.fq" ## fastq 
gzip_fname = "test.fq.gz" ## gzip fq
fasta_name = "test.fa"


### FASTQ synopsis

### Method 1
fqparser = FQParser()
fqparser.open(fqname)

for id , seq, qual in fqparser.parse() :
	### Process for your own

fqparser.close()

### Method 2

fqparser = FQParser(fqname)
fqparser.open()

for id , seq , qual in fqparser.parse() :
	### Process for your own

fqparser.close()

Note : One can use gzip_fname instead of fqname for gzip compressed file.
       we have same usage
### FASTQ synopsis

### Method 1
faparser = FAParser()
faparser.open(fasta_name)
for id , description, seq in faparser.parse() :
	### Process on your own

faparser.close()

### Method 2 

faparser = FAParser(fasta_name)
faparser.open()

for id , description, seq in faparser.parse():
	### Process on your own

faparser.close()

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

class FQParser(object):
	def __init__(self, fqname=None, format="r"):
		if not fqname : ## None
			self.fqname = None
		else:
			self.fqname = fqname

	def open(self, fqname=None):
		if self.fqname == None:
			if fqname == None:
				raise AttributeError("File name must be specified")
			else:
				self.fqname = fqname

		self.fqfile = format_checker(self.fqname)
	def rewind(self):
		self.fqfile.seek(0)
	def close(self) :
		self.fqfile.close()
	def next(self):
		id = self.fqfile.readline()
		seq = self.fqfile.readline()
		self.fqfile.readline()
		qual = self.fqfile.readline()
		return id, seq, qual
	def parse(self) :
		id = ''
		seq = ''
		qual = ''
		for i, line in enumerate(self.fqfile) :
			if i &3 == 0 : # Case when header
				id = line.rstrip()
			elif i & 3 == 1 :
				seq = line.rstrip()
			elif i & 3 == 2:
				continue
			elif i & 3 == 3 :
				qual = line.rstrip()
				yield id, seq, qual


class FAParser(object):
	def __init__(self, faname=None):
		if not faname :
			self.faname = None
		else:
			self.faname = faname

	def open(self, faname=None) :
		if self.faname == None:
			if faname == None:
				raise AttributeError("Fasta File name must be specified")

			else:
				self.faname = faname

		self.fafile = open(self.faname)
	def rewind(self):
		self.fafile.seek(0)
	def close(self) :
		self.fafile.close()

	def parse(self) :
		id = ''
		desc = ''
		seq = ''
		seq_trail = []
		for line in self.fafile:
			if line.startswith(">") :
				if seq_trail :
					yield id , desc , "".join(seq_trail)
				id = line.rstrip().split()[0][1:]
				desc = line.rstrip()[1:]
				seq_trail = []

			else:
				seq_trail.append(line.rstrip())
		if seq_trail:
			yield id , desc , "".join(seq_trail)
