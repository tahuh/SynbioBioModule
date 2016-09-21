#!/usr/bin/python

# Constant field
# codon_usage = """TTT F 0.45 16.9(336562)  TCT S 0.18 14.6(291040)  TAT Y 0.43 12.0(239268)  TGT C 0.45  9.9(197293)
# TTC F 0.55 20.4(406571)  TCC S 0.22 17.4(346943)  TAC Y 0.57 15.6(310695)  TGC C 0.55 12.2(243685)
# TTA L 0.07  7.2(143715)  TCA S 0.15 11.7(233110)  TAA * 0.28  0.7(14322)  TGA * 0.52  1.3(25383)
# TTG L 0.13 12.6(249879)  TCG S 0.06  4.5(89429)  TAG * 0.20  0.5(10915)  TGG W 1.00 12.8(255512)
# CTT L 0.13 12.8(253795)  CCT P 0.28 17.3(343793)  CAT H 0.41 10.4(207826)  CGT R 0.08  4.7(93458)
# CTC L 0.20 19.4(386182)  CCC P 0.33 20.0(397790)  CAC H 0.59 14.9(297048)  CGC R 0.19 10.9(217130)
# CTA L 0.07  6.9(138154)  CCA P 0.27 16.7(331944)  CAA Q 0.25 11.8(234785)  CGA R 0.11  6.3(126113)
# CTG L 0.41 40.3(800774)  CCG P 0.11  7.0(139414)  CAG Q 0.75 34.6(688316)  CGG R 0.21 11.9(235938)
# ATT I 0.36 15.7(313225)  ACT T 0.24 12.8(255582)  AAT N 0.46 16.7(331714)  AGT S 0.15 11.9(237404)
# ATC I 0.48 21.4(426570)  ACC T 0.36 19.2(382050)  AAC N 0.54 19.5(387148)  AGC S 0.24 19.4(385113)
# ATA I 0.16  7.1(140652)  ACA T 0.28 14.8(294223)  AAA K 0.42 24.0(476554)  AGA R 0.20 11.5(228151)
# ATG M 1.00 22.3(443795)  ACG T 0.12  6.2(123533)  AAG K 0.58 32.9(654280)  AGG R 0.20 11.4(227281)
# GTT V 0.18 10.9(216818)  GCT A 0.26 18.6(370873)  GAT D 0.46 22.3(443369)  GGT G 0.16 10.8(215544)
# GTC V 0.24 14.6(290874)  GCC A 0.40 28.5(567930)  GAC D 0.54 26.0(517579)  GGC G 0.34 22.8(453917)
# GTA V 0.11  7.0(139156)  GCA A 0.23 16.0(317338)  GAA E 0.42 29.0(577846)  GGA G 0.25 16.3(325243)
# GTG V 0.47 28.9(575438)  GCG A 0.11  7.6(150708)  GAG E 0.58 40.8(810842)  GGG G 0.25 16.4(326879)"""

codon_table = {
    'A': ('GCT', 'GCC', 'GCA', 'GCG'),
    'C': ('TGT', 'TGC'),
    'D': ('GAT', 'GAC'),
    'E': ('GAA', 'GAG'),
    'F': ('TTT', 'TTC'),
    'G': ('GGT', 'GGC', 'GGA', 'GGG'),
    'I': ('ATT', 'ATC', 'ATA'),
    'H': ('CAT', 'CAC'),
    'K': ('AAA', 'AAG'),
    'L': ('TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'),
    'M': ('ATG',),
    'N': ('AAT', 'AAC'),
    'P': ('CCT', 'CCC', 'CCA', 'CCG'),
    'Q': ('CAA', 'CAG'),
    'R': ('CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'),
    'S': ('TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'),
    'T': ('ACT', 'ACC', 'ACA', 'ACG'),
    'V': ('GTT', 'GTC', 'GTA', 'GTG'),
    'W': ('TGG',),
    'Y': ('TAT', 'TAC'),
    '*': ('TAA', 'TAG', 'TGA'),
}

REVERSED_CODON = {
	'CTT': 'L', 'ATG': 'M', 'AAG': 'K', 'AAA': 'K', 
	'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R', 
	'CCT': 'P', 'ACT': 'T', 'AGC': 'S', 'ACA': 'T', 
	'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I', 
	'CTG': 'L', 'CTA': 'L', 'CTC': 'L', 'CAC': 'H', 
	'ACG': 'T', 'CAA': 'Q', 'AGT': 'S', 'CAG': 'Q', 
	'CCG': 'P', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G', 
	'TGT': 'C', 'CGA': 'R', 'CCA': 'P', 'TCT': 'S', 
	'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C', 
	'GGG': 'G', 'TAG': '*', 'GGA': 'G', 'TAA': '*', 
	'GGC': 'G', 'TAC': 'Y', 'GAG': 'E', 'TCG': 'S', 
	'TTA': 'L', 'GAC': 'D', 'TCC': 'S', 'GAA': 'E', 
	'TCA': 'S', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A', 
	'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'TTC': 'F', 
	'GTT': 'V', 'GCT': 'A', 'ACC': 'T', 'TGA': '*', 
	'TTG': 'L', 'CGT': 'R', 'TGG': 'W', 'CGC': 'R'
}
HUMAN_CODON_USAGE = {
	'CTT': 0.13, 'ACC': 0.36, 'ACA': 0.28, 'AAA': 0.42, 
	'ATC': 0.48, 'AAC': 0.54, 'ATA': 0.16, 'AGG': 0.2, 
	'CCT': 0.28, 'ACT': 0.24, 'AGC': 0.24, 'AAG': 0.58, 
	'AGA': 0.2, 'CAT': 0.41, 'AAT': 0.46, 'ATT': 0.36, 
	'CTG': 0.41, 'CTA': 0.07, 'CTC': 0.2, 'CAC': 0.59, 
	'ACG': 0.12, 'CAA': 0.25, 'AGT': 0.15, 'CCA': 0.27, 
	'CCG': 0.11, 'CCC': 0.33, 'TAT': 0.43, 'GGT': 0.16, 
	'TGT': 0.45, 'CGA': 0.11, 'CAG': 0.75, 'TCT': 0.18, 
	'GAT': 0.46, 'CGG': 0.21, 'TTT': 0.45, 'TGC': 0.55, 
	'GGG': 0.25, 'TAG': 0.2, 'GGA': 0.25, 'TGG': 1.0, 
	'GGC': 0.34, 'TAC': 0.57, 'TTC': 0.55, 'TCG': 0.06, 
	'TTA': 0.07, 'TTG': 0.13, 'CGT': 0.08, 'GAA': 0.42, 
	'TAA': 0.28, 'GCA': 0.23, 'GTA': 0.11, 'GCC': 0.4, 
	'GTC': 0.24, 'GCG': 0.11, 'GTG': 0.47, 'GAG': 0.58, 
	'GTT': 0.18, 'GCT': 0.26, 'TGA': 0.52, 'GAC': 0.54, 
	'TCC': 0.22, 'TCA': 0.15, 'ATG': 1.0, 'CGC': 0.19 
}