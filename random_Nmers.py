"""
random_Nmers.py

Generates random N mers

Author : Thomas Sunghoon Heo
"""

import sys
import itertools
import time

argv = sys.argv
argc = len(argv)
usage = """random_Nmers.py [N] [outfilename]
Commandline argument description
N            INT     The number of 'N' mers for user to choose
outfilename  STR     Outfile name to write down result of random N mers"""

# Main function
def main():
	if argc < 3:
		sys.stdout.write(usage + "\n")
		exit(-1)
	sys.stdout.write("Generating.....\n")
	t1 = time.time()
	nmers = [ x for x in list(itertools.product("ATGC", repeat=int(argv[1]))) ]
	o = open(argv[2] , "w")
	for mer in nmers:
		o.write(mer + "\n")
	o.close()
	t2 = time.time() - t1
	sys.stdout.write("Done generating random %d mers in %.3f sec\n"%(int(argv[1]), t2))
if __name__ == "__main__":
	main()
