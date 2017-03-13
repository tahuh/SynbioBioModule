#!/usr/bin/python

"""
ChipOligoDesigner.py

This code is advanced version of chip oligo design for members of
Sythetic Biology Lab.

This program translates protein sequence in to DN sequence using 
human codon usage


Author : Sunghoon Heo
"""
import random
import sys
import time

import CodonConstants as Constants
from BioUtils import Seq
from SequenceParser import FAParser

### GLOBAL CONSTANTS

BSA_FWD = "GGTCTC"
BSA_REV = "GAGACC"


class Designer(object) :
	def __init__(self,opts):
		self.codon_table = Constants.codon_table
		self.reverse_codon = Constants.REVERSED_CODON
		self.human_codon_usage = Constants.HUMAN_CODON_USAGE
		
		### Data field parsed from options
		
		self.left_flank = opts.left_flank ### left flanking to match multiple of tiling coverage
		self.right_flank = opts.right_flank ### right flanking to match multiple of tiling coverage
		self.left_insert = opts.left_insert ### Left insert that added to assembled DNA sequence
		self.right_insert = opts.right_insert ### Right insert that added to assembled DNA sequence
		self.faname = opts.faname
		self.oname = opts.oname
		self.frq = opts.frq
		self.coverage = opts.tile_cov
		self.slide_size = opts.slide_size
		self.final_add_left = opts.chip_left
		self.final_add_right = opts.chip_right
		self.dna_file = open("./log.log" , "w")
		
		
		if self.final_add_left == None:
			self.final_add_left = ''
		if self.final_add_right == None:
			self.final_add_right = ''
		if self.left_flank == None :
			self.left_flank = ''
		if self.right_flank == None:
			self.right_flank = ''
		
		if self.frq == None : self.frq = 0.10
		
		self.records = [] ### Will be stored in tuple format
		
		self._parse_fasta()
		
		self.usage_library = {}
		
		for aa in self.codon_table:
			vals = self.codon_table[aa]
			self.usage_library[aa] = []
			for codon in vals:
				frq = self.human_codon_usage[codon]
				FRQ = int(100 * frq)
				if frq > self.frq :
					i = 1
					while i <= FRQ :
						self.usage_library[aa].append(codon)
						i += 10
	
	def design(self) :
		
		back_translated = self._reverse_translation()
		cons_removed = self._remove_consecutive_codons(back_translated)
		homo_rmv = self._remove_homopolymer(cons_removed)
		enzyme_removed = self._remove_basI_site(homo_rmv)
		
		for id , dna in enzyme_removed :
			self.dna_file.write(">" + str(id) + "\n" + dna + "\n")
			
		inserted = self._insert_dna(enzyme_removed)
		slice_and_add = self._slice_and_add(inserted)
		self.result = slice_and_add
		
		
	def _parse_fasta(self):
		fa_parser = FAParser()
		
		fa_parser.open(self.faname)
		
		for id , _ , seq in fa_parser.parse() :
			self.records.append((id,seq.upper()))
			print id
		fa_parser.close()
		
	def _insert_dna(self , enz_removed_set):
		inserted = []
		for id , dna in enz_removed_set :
			ins = self.__insert_dna(dna)
			inserted.append((id, ins))
			
		return inserted
		
	def __insert_dna(self, dna) :
		new = self.left_insert + dna + self.right_insert
		if len(dna) % (self.slide_size/self.coverage) == 0 :
			return new
			
		i = 0
		
		while True:
			left_add = self.left_flank[i]
			new = left_add + new
			if len(new) % (self.slide_size/self.coverage) == 0 :
				break
			else:
				right_add =self.right_flank[i]
				new = new + right_add
				if len(new) % (self.slide_size/self.coverage) == 0 :
					break
				else:
					i += 1
					continue
					
		return new
		
		
	def _slice_and_add(self , inserted_set):
		rets = []
		for id , inserted_dna in inserted_set :
			fragments = self.__fragmentation(inserted_dna)
			ret = []
			for slice in fragments :
				new = self.final_add_left + BSA_FWD + slice + BSA_REV + self.final_add_right
				ret.append(new)
				
			rets.append((id,ret))
			
		return rets
		
	def _is_homopolymer(self, window):
		n = 1
		t= window[0]
		for b in window[1:] :
			if b == t : n += 1
			
		if n == len(window) :
			return True
		else:
			return False
		
	def _remove_homopolymer(self,con_removed) :
		""" Scan over DNA sequence and remove holompolyers """
		R = []
		for id, dna in con_removed :
			new_dna = self.__remove_homopolymer(dna)
			R.append((id,new_dna))
			
		return R
	def __remove_homopolymer(self,_dna):
		codon_table = self.codon_table
		reverse_codon = self.reverse_codon
		dna = _dna
		for x in xrange(0,len(dna)-5,3):
			window = dna[x:x+5]
			if self._is_homopolymer(window):
				front_codon = dna[x:x+3]
				a = reverse_codon[front_codon]
				codon_list = codon_table[a]
				for c in codon_list :
					if c!= front_codon:
						dna1 = dna[:x]
						dna2 = dna[x+3:]
						dna_temp = dna1 + c + dna2
						temp_window = dna_temp[x:x+5]
						if not self._is_homopolymer(temp_window):
							dna = dna_temp
			else:
				continue
				
		return dna
	def __fragmentation(self , dna) :
		fragments = []
		for i in xrange(0,len(dna) - (self.slide_size - 1), self.slide_size/self.coverage):
			fragments.append(dna[i:i+self.slide_size])
		return fragments
	def _reverse_translation(self) :
		reverse_translated = []
		for id , seq in self.records :
			iid , dna = self.__reverse_translation(id,seq)
			
			reverse_translated.append((iid, dna))
			
		return reverse_translated
		
	def __reverse_translation(self,id, seq):
		""" reverse translation using codon optimized property """
		id = id
		seq = seq
		dna = ''
		for acid in seq :
			codons = self.usage_library[acid]
			random.shuffle(codons)
			codon = random.choice(codons)
			dna += codon
		return id , dna
	
	def _remove_consecutive_codons(self , translated_list):
		rmv_set = []
		for id , dna in translated_list:
			
			dna_origin = Seq(dna).translate()
			rmv_seq = self.__remove_consecutive_codons(dna)
			rmv_translated = Seq(rmv_seq).translate()
			C = 0
			Z = zip(dna_origin , rmv_translated)
			for z in Z :
				if z[0] == z[1] :
					C += 1
					
			assert len(dna_origin) == C
			rmv_set.append((id , rmv_seq))
			
		return rmv_set
		
	def __remove_consecutive_codons(self, dna) :
		removed = dna
		for i in xrange(len(dna) - 18) :
			end = removed[i:i+18]
			seq_obj = Seq(end)
			trans_obj = seq_obj.translate()
			consecutive = list(trans_obj)
			tmp = consecutive[0]
			mark = 0
			for ele in consecutive :
				new = ele
				if new == tmp:
					mark += 1
			if mark == 6 :
				last = end[-3:]
				orogin = self.codon_table[last]
				for ele in origin :
					if ele != last :
						dna1 = dna[:i+15]
						dna2 = dna[i+18:]
						removed = dna1 + ele + dna2
						break
		return removed
	def _remove_basI_site(self , cons_removed_set):
		S = []
		for id , dna in cons_removed_set:
			iid = id
			basi_removed = self.__remove_basI_site(dna)
			S.append((iid,basi_removed))
			
		return S
	def __remove_basI_site(self,_dna) :
		dna = _dna
		for x in xrange(0,len(dna)-5,1) :
			small = dna[x:x+6]
			if small == BSA_FWD:
				if x % 3 == 0 :
					front = dna[x-3:x]
					codon = dna[x:x+3]
					rear = dna[x+3:x+6]
					acid = self.reverse_codon[codon]
					codons = self.codon_table[acid]
					
					for c in codons :
						N = front + c
						M = c + rear
						if c != codon and N != BSA_FWD and N != BSA_REV and M != BSA_FWD and M != BSA_REV:
							
							dna = dna[:x] + c + dna[x+3:]
							
							break
				elif x % 3 == 1 :
					front = dna[x-1:x+2]
					codon = dna[x+2:x+5]
					rear = dna[x+5:x+8]
					acid = self.reverse_codon[codon]
					codons = self.codon_table[acid]
					
					for c in codons :
						N = front + c
						M = c + rear
						if c != codon and N != BSA_FWD and N != BSA_REV and M != BSA_FWD and M != BSA_REV:
							dna = dna[:x+2] + c + dna[x+5:]
							break
							
				elif x % 3 == 2 :
					front = dna[x-2:x+1]
					codon = dna[x+1:x+4]
					rear = dna[x+4:x+7]
					acid = self.reverse_codon[codon]
					codons = self.codon_table[acid]
					
					for c in codons :
						N = front + c
						M = c + rear
						if c != codon and N != BSA_FWD and N != BSA_REV and M != BSA_FWD and M != BSA_REV:
							dna = dna[:x+1] + c + dna[x+4:]
							break
			elif small == BSA_REV:
				if x % 3 == 0 :
					front = dna[x-3:x]
					codon = dna[x:x+3]
					rear = dna[x+3:x+6]
					acid = self.reverse_codon[codon]
					codons = self.codon_table[acid]
					
					for c in codons :
						N = front + c
						M = c + rear
						if c != codon and N != BSA_FWD and N != BSA_REV and M != BSA_FWD and M != BSA_REV:
							dna = dna[:x] + c + dna[x+3:]
							break
				elif x % 3 == 1 :
					front = dna[x-1:x+2]
					codon = dna[x+2:x+5]
					rear = dna[x+5:x+8]
					acid = self.reverse_codon[codon]
					codons = self.codon_table[acid]
					
					for c in codons :
						N = front + c
						M = c + rear
						if c != codon and N != BSA_FWD and N != BSA_REV and M != BSA_FWD and M != BSA_REV:
							dna = dna[:x+2] + c + dna[x+5:]
							break
							
				elif x % 3 == 2 :
					front = dna[x-2:x+1]
					codon = dna[x+1:x+4]
					rear = dna[x+4:x+7]
					acid = self.reverse_codon[codon]
					codons = self.codon_table[acid]
					
					for c in codons :
						N = front + c
						M = c + rear
						if c != codon and N != BSA_FWD and N != BSA_REV and M != BSA_FWD and M != BSA_REV:
							dna = dna[:x+1] + c + dna[x+4:]
							break
			else:
				continue
		return dna
		
		
