#! /usr/bin/python

import re
import os
import sys
import gzip

def rmUpaired(Args):
	if not Args:
		print 'Usage: rmUpaired <fq1.in> <fq2.in> <fq1.out> <fq2.out>' 
	infile1 = open(Args[0], 'r')
	infile2 = open(Args[1], 'r')
	outfile1=open(Args[2], 'w')
	outfile2=open(Args[3], 'w')
	#fq1_tag=[]
	fq1_reads=[]
	fq2_reads=[]
	total_reads1 = 0
	total_reads2 = 0
	while True:
		readid = infile1.readline().rstrip()
		if len(readid) == 0:
			break
		seq = infile1.readline().rstrip()
		mark = infile1.readline().rstrip()
		qual = infile1.readline().rstrip()
		fq1_reads.append(readid)
		#print readid
		#outfile.write('%s\n' %(fq1_tag[total_reads]))
		total_reads1 += 1
	#print fq1_reads
	print 'total reads from fq_1 is: ',total_reads1
	
	while True:
		readid = infile2.readline().rstrip()
		if len(readid) == 0:
			break
		seq = infile2.readline().rstrip()
		mark = infile2.readline().rstrip()
		qual = infile2.readline().rstrip()
		fq2_reads.append(readid)
		#print readid
		#outfile.write('%s\n' %(fq1_tag[total_reads]))
		total_reads2 += 1
	#print fq2_reads
	print 'total reads from fq_2 is: ',total_reads2
	infile1.close()
	infile2.close()
	fq1_reads=set(fq1_reads)
	fq2_reads=set(fq2_reads)
	
	fq_reads=fq1_reads.intersection(fq2_reads)
	#fq_reads=list(fq_reads)
	#fq_reads=[val for val in fq1_reads if val in fq2_reads]
	#print fq_reads
	print 'after rmUnpaired total reads is: ',len(fq_reads)
	
	infile1 = open(Args[0], 'r')
	infile2 = open(Args[1], 'r')
	
	n=0
	while True:
		readid = infile1.readline().rstrip()
		if len(readid) == 0:
			break
		seq = infile1.readline().rstrip()
		mark = infile1.readline().rstrip()
		qual = infile1.readline().rstrip()
		
		if (readid in fq_reads) :
			n=n+1
		#	readid+='/1'
		#	print str(n)+readid+'\n'
			outfile1.write('%s\n%s\n%s\n%s\n' %(readid,seq,mark,qual))
	outfile1.close()		
	n=0
	while True:
		readid = infile2.readline().rstrip()
		if len(readid) == 0:
			break
		seq = infile2.readline().rstrip()
		mark = infile2.readline().rstrip()
		qual = infile2.readline().rstrip()
		
		if (readid in fq_reads) :
			n=n+1
		#	readid+='/2'
		#	print str(n)+readid+'\n'
			outfile2.write('%s\n%s\n%s\n%s\n' %(readid,seq,mark,qual))
	outfile2.close()
	
	
	
if __name__ == '__main__':
	if len(sys.argv) > 1:
		rmUpaired(sys.argv[1:])
	else:
		
		sys.exit(0)



