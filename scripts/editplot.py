#! /usr/bin/python
import re
import os
import sys

def plot(dat, outname):
	infile = open(dat, 'r')
	outfile = open(outname, 'w')

	exon = 0
	for line in infile:
		line = line.rstrip()
		locate = line.split()[1]
		num = int(line.split()[0])
		if not re.search('downstream|upstream|ncRNA', line):
			if re.search('exonic|UTR3|UTR5', line):
				exon += num
			else:
				outfile.write('%s\t%s\n' %(locate, num))
	outfile.write('exonic\t%s\n' %exon)
	infile.close()
	outfile.close()

if len(sys.argv) > 1:
	plot(sys.argv[1], sys.argv[2])
else:
	print 'Usage: plot <dat.txt> <outfile>'
	sys.exit(0)
