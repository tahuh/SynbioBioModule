#!/usr/bin/python

"""
BarcodeGenerator.py

Generates barcode sequences with length n

Usage : python BarcodeGenerator.py -n 8 -o output.fasta

Author : Sunghoon Heo
"""
import itertools
import argparse
impot sys
parser = argparse.ArgumentParser()
parser.add_argument("-n" ,"--mers" ,  type=int, required=True, help="Length of barcode")
parser.add_argument("-o", "--outfile", type=str ,required=True, help="Output fasta")

args =parser.parse_args()

o = open(args.outfile , "w")
sys.stdout.write("Generate\n")
combs = itertools.combinations_with_replacement("ATGC", args.mers)
i = 0
for c in combs:
    o.write(">barcode-" + str(i) + "\n")
    seq = ''.join(c)
    o.write(seq + "\n")
print "Done"
o.close()
