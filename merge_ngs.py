#!/usr/bin/python

"""
merge_ngs.py

Merge NGS files in a single file

usage

python merge_ngs.py --list <listfile> --outfile out.fq

"""
import gzip
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--list", type=str, required=True, help="A file contains list of fastq files to merge. Each file in a row with full path")
parser.add_argument("--outfile", type=str, required=True, help="Output target file")

args = parser.parse_args()

def format_checker(fname):
    F = open(fname , "rb")
    bytes = F.read(3)
    if bytes == '\x1f\x8b\x08':
        F.seek(0)
        return gzip.GzipFile(fname , fileobj=F)
    else:
        F.seek(0)
        return F

O = open(args.outfile , "w")

for line in open(args.list):
    fname = line[:-1]
    handle = format_checker(fname)
    for lin in handle :
        O.write(lin + "\n")

O.close()
