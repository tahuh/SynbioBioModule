#!/usr/bin/python
"""
count_num_seqs.py

Count the number of sequences in fasta/fastq file

Author : Thomas Sunghoon Heo
"""
import argparse
import sys
parser = argparse.ArgumentParser()
parser.add_argument("--infile", type=str, required= True, help="infile")
parser.add_argument("--type", type=str,default="fasta", help="infile type(fasta/fastq) default=fasta)

if args.type == 'fasta':
    head = '>'
elif args.type == 'fastq':
    head = '@'
else:
    raise AttributeError("--type argument must be fasta for fastq. others are not accepted\n")
cnt = 0
with open(args.infile) as F:
    for line in F:
        if line[0] == head:
            cnt += 1
sys.stdout.write("File %s contains %d sequences\n"%(args.infile , cnt))
