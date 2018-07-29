from operator import itemgetter, attrgetter
import sys

def sortCover(args):
	infile=open(args[0],'r')
	outfile=open(args[1],'w')
	content_infile=[]
	while True:
		line = infile.readline().rstrip()
		if not line:
			break
		
		line_list=line.split('\t')
		line_list[1]=int(line_list[1])
		content_infile.append(line_list)
	
	content_infile=sorted (content_infile, key=itemgetter(0,1))
	
	for content in content_infile:
		content[1]=str(content[1])
		line='\t'.join(content)		
		line=line+'\n'
		outfile.write(line)
	infile.close()
	outfile.close()

if __name__ == '__main__':
	if len(sys.argv) > 1:
		sortCover(sys.argv[1:])
	else:
		print 'Usage: sortCover <in> <out> '
		sys.exit(0)

	
