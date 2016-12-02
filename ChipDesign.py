#!/usr/bin/python

"""
ChipDesign.py

code availity : Prof. Bang's members ( Synthetic biology lab. Yonsei Univ. Korea )

Author : Sunghoon Heo
"""

from argparse import ArgumentParser
from Designer import Designer

parser = ArgumentParser()

parser.add_argument("--left_flank" , type=str, required=True, help= "Left flanking sequence")
parser.add_argument("--right_flank" , type=str, required=True, help= "Right flanking sequence")
parser.add_argument("--left_insert" , type=str, required=True, help= "Left insert sequence")
parser.add_argument("--right_insert" , type=str, required=True, help= "Right insert sequence")
parser.add_argument("--faname" , type=str, required=True, help = "FASTA name with full path")
parser.add_argument("--oname" , type=str, required=True, help="Output name with full path")
parser.add_argument("--freq" , type=float, required=False, help="freq cut for human codon usage" , default=0.10)

parser.add_argument("--tile_cov",type=int , required = False, help= "Tile coverage such as 10 ,.." , default=10)

parser.add_argument("--slide_size",type=int,required = False , help= "Tile window(insert size). Must be multiple of tile_cov" , default = 70)

parser.add_argument("--chip_left" , type=str, required=False, help="Left addition for Chip oligo")
parser.add_argument("--chip_right" , type=str, required=False, help="Right addition for Chip oligo")

opts = parser.parse_args()

designer = Design(opts)

optdict = vars(opts)

print "Designer starts to design"

designer.design()

ofile = open(opts.oname, "w")
for id , ret_set in designer.records:
	length = len(ret_set)
	for k in xrange(len(length)) :
		ofile.write(id + "_" + str(k+1) + "\t" + ret_set[k] + "\n")
		
ofile.close()

print "Done design"
