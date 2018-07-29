#! /usr/bin/python

import re
import os
import sys
import gzip

def rmN(Args):
	if not Args:
		print 'Usage: rmN <in.fq> <max_N_ratio><outDir>'
		print 'Max_N_ratio [float] , e.g. 0.05'
		sys.exit(0)
	if re.search('gz$|gzip$', Args[0]):
		infile = gzip.open(Args[0], 'rb')
	else:
		infile = open(Args[0], 'r')
	maxN = float(Args[1])
	out_dir = Args[-1]
	fqname = os.path.basename(Args[0]).replace('.gz', '')
	fqname = fqname.replace(fqname.split('.')[-1], 'fastq')
	remain = open(Args[-1] + '/remained_' + fqname, 'w')
	#remove = open(Args[-1] + '/removed_N_' + fqname, 'w')
	remove_log  = open(Args[-1] + '/removed_N_' + fqname + '.log', 'w')
	remove_reads = 0
	total_reads = 0
	while True:
		readid = infile.readline().rstrip()
		if len(readid) == 0:
			break
		seq = infile.readline().rstrip()
		mark = infile.readline().rstrip()
		qual = infile.readline().rstrip()
		seqlen = len(seq)
		total_reads += 1
		N_count = seq.count('N')
		N_ratio = float(N_count) / seqlen
		if N_ratio <= maxN:
			remain.write('%s\n%s\n%s\n%s\n' %(readid, seq, mark, qual))
		else:
	#		remove.write('%s\n%s\n%s\n%s\n' %(readid, seq, mark, qual))
			remove_reads += 1 
	remove_rate = float(remove_reads) / total_reads * 100
	remove_log.write('Total Reads\tRemoved N Reads\t%Removed N Rate\n')
	remove_log.write('%s\t%s\t%.2f\n' %(total_reads, total_reads - remove_reads, remove_rate))

if __name__ == '__main__':
	if len(sys.argv) > 1:
		rmN(sys.argv[1:])
	else:
		print 'Usage: rmN <in.fq> <max_N_ratio><outDir>'
		print 'Max_N_ratio [float] , e.g. 0.05'
		sys.exit(0)



