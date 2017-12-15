#!/usr/bin/python

"""
Seq.py

This source code offers some utilities for
DNA sequence manipulation

Author : Sunghoon Heo
"""

BASE_CONV = { x : y , for x , y in zip("ATGCN" , "TACGN") }

class Seq:
	def __init__(self, seq=None):
		self.seq = seq
		if self.seq == None:
			self.lseq = 0
		else:
			self.lseq = len(seq)

	def setup_sequence(self,seq):
		self.seq = seq
		self.lseq = len(seq)

	def length(self):
		return self.lseq

	def seq(self):
		return seq

	def complement(self, new = None):
		if self.seq == None: return None
		else:
			if new == None:
				s = []
				for b in self.seq:
					s.append(BASE_CONV[b])
					return "".join(s)
			else:
				s = []
				for b in new :
					s.append(BASe_CONV[b])
					return "".join(s)

	def revserse_complement(self, new =None):
		comp = self.complement(new=new)
		if comp == None: return None
		else: return comp[::-1]

	def dna_to_rna(self,seq=None):
		if seq == None:
			if self.seq == None:
				return None
			else:
				s = []
				for b in self.seq:
					if b =='T': s.append('U')
					else: s.append(b)
				return "".join(s)
		else:
			s = []
			for b in seq:
				if b == 'T': s.append('U')
				else: s.append(b)
			return "".join(s)
