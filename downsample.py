#!/usr/bin/python

"""
downsample.py

This code downsamples NGS fastq files

usage:

downsample.py \
--file1 read1.fq \
--file2 read2.fq \
--out1 read1_down.fq \
--out2 read2_out.fq \
--rate 0.5 

Author : Sunghoon Heo

Dependency : No dependency
"""
import argparse
import gzip
import random

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

	def close(self) :
		self.fqfile.close()

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

parser = argparse.ArgumentParser()

parser.add_argument("--file1" , help="fastq file 1" , type=str, required=True)
parser.add_argument("--file2" , help="fastq file 2" , type=str, required=True)
parser.add_argument("--out1" , help="Outfile 1" , type=str, required=True)
parser.add_argument("--out2" , help="Outfile 2" , type=str, required=True)
parser.add_argument("--rate" , type=float, help="Sampling rate. this value much will be sample. default=0.4" , default=0.4)

args = parser.parse_args()


file1 = args.file1
file2 = args.file2
out1 = args.out1
out2 = args.out2
rate = args.rate

print "Settings"
print "FILE1 : %s"%(file1)
print "FILE2 : %s"%(file2)
print "OUT1 : %s"%(out1)
print "OUT2 : %s"%(out2)
print "RATE : %f"%(rate)

parser1 = FQParser(file1)
parser1.open()
parser2 = FQParser(file2)
parser2.open()
out1 = open(args.out1 , "w")
out2 = open(args.out2 , "w")
for (id , seq ,qual) , (id2,seq2,qual2) in zip(parser1.parse() , parser2.parse()):
	if random.random() < rate:
		out1.write(id + "\n" + seq + "\n+\n" + qual + "\n")
		out2.write(id2 + "\n" + seq2 + "\n+\n" + qual2 + "\n")
parser1.close()
parser2.close()
out1.close()
out2.close()
