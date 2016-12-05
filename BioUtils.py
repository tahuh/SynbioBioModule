#!/usr/bin/python

"""
BioUtils.py

Biopython like utils

Custom code

Author : Sunghoon Heo
"""

import random
import CodonConstants

codon_table = CodonConstants.codon_table # acid : codons
rc_codon_table = CodonConstants.REVERSED_CODON # Codon : acid

base_dict = {"A" : "T" , "C" : "G" , "G" : "C" , "T" : "A" , "R" : "Y" , "Y" : "R" , "N" : "N"}


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
			codon = self.seq[i:i+3]
			if codon not in rc_codon_table :
				aa.append("*")
			else:
				aa.append(rc_codon_table[codon])
				
		return ''.join(aa)
	def to_dna(self) :
		"""operates on amino acid sequence"""
		seq = self.seq
		ret = []
		for acid in seq :
			codons = codon_table[acid]
			selected_codon = random.choice(codons)
			ret.append(selected_codon)
		return "".join(ret)
		
