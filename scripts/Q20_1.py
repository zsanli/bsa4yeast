import re
import os
import sys
import gzip

def Qstat(Args):
	if not Args:
		print 'Usage: Qstat <in.fq> <outDir> <phred>'
		sys.exit(0)
		
	infile = open(Args[0], 'r')
	out_dir=Args[1]
	outname = 'Q20_Q30_%s' %os.path.basename(Args[0])
	outfile = open(out_dir + '/' + outname + '.txt', 'w')
	min_qual = int(Args[2])

	Q20 = 0
	Q30 = 0
	Total_base = 0
	read = 0
	GC = 0

	while True:
		seqid = infile.readline().rstrip()
		if len(seqid) == 0:
			break
		seq = infile.readline().rstrip()
		GC_num = seq.count('G') + seq.count('C')
		GC += GC_num
		mark = infile.readline().rstrip()
		qual = infile.readline().rstrip()
		Total_base += len(qual)
		read += 1
		for qu in qual:
			qual_score = ord(qu) - min_qual
			if qual_score >= 20:
				Q20 += 1
			if qual_score >= 30:
				Q30 += 1
	#print "%s\t%s" %(GC,Total_base)
	#sys.exit(0)
	pGC = GC/float(Total_base) *100

	PQ20 = float(Q20) / Total_base * 100
	PQ30 = float(Q30) / Total_base * 100
	outfile.write('Minimum Quality: %s\n' % min_qual)
	outfile.write('Total reads\t%s\n' %read)
	outfile.write('Total bases\t%s\n' %format(Total_base, ','))
	outfile.write('GC%%\t%.2f\n' %pGC)
	outfile.write('Q20 bases\t%s\n' %format(Q20,','))
	outfile.write('Q20 bases%%\t%.2f\n' %PQ20)
	outfile.write('Q30 bases\t%s\n' %format(Q30,','))
	outfile.write('Q30 bases%%\t%.2f\n' %PQ30)

	infile.close()
	outfile.close()

if __name__ == '__main__':
	if len(sys.argv) > 1:
		Qstat(sys.argv[1:])
	else:
		print 'Usage: Qstat <in.fq> <outDir>'
		sys.exit(0)


	
