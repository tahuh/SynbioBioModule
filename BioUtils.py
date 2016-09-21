#!/usr/bin/python

"""
BioUtils.py

Biopython like utils

Custom code

Author : Sunghoon Heo
"""

import CodonConstants

codon_table = CodonConstants.codon_table

base_dict = {"A" : "T" , "C'" : "G" , "G" : "C" , "T" : "A" , "R" : "Y" , "Y" : "R" , "N" : "N"}


class Seq(object):
	def __init__(self,seq):
		self.seq = seq
		
	def complement(self) :
		comp = []
		for b in self.seq :
			comp.append(base_dict[b])
			
		return ''.join(comp)
		
	def reverse_complement(self) :
		rev_comp = self.complement()
		return rev_comp[::-1]
		
	def translate(self) :
		if len(self.seq) % 3 != 0 :
			remainder = 3 - len(self.seq) % 3
			poly_n = "N" * remainder
			self.seq += poly_n
		aa = []
		for i in xrange(len(self.seq)):
			codon = self.seq[i:i+1]
			if codon not in codon_table :
				aa.append("*")
			else:
				aa.append(codon_table[codon])
				
		return ''.join(aa)