import sys

inFile = open(sys.argv[1],'r')
line=inFile.readline().strip().split('\t')
line.append('marker')
#line.append('HetMarker')
print '\t'.join([str(x) for x in line])
for line in inFile:
    data = line.strip().split('\t')
    marker=''
    if (data[-4]=='Pass' and data[-2]=='Pass' and data[-3]!=data[-1]) :
        marker='MPass'
    else:
        marker='No_Pass'


    data.append(marker)
    print '\t'.join([str(x) for x in data])

