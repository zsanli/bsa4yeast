import sys

inFile = open(sys.argv[1],'r')
line=inFile.readline().strip().split('\t')
print '\t'.join([str(x) for x in line])
chr_data={'I':[],'II':[],'III':[],'IV':[],'V':[],'VI':[],'VII':[],'VIII':[],'IX':[],'X':[],'XI':[],'XII':[],'XIII':[],'XIV':[],'XV':[],'XVI':[]}
for line in inFile:
    data = line.strip().split('\t')
    if data[0] in chr_data.keys():
        chr_data[data[0]].append(line.strip())

for i in ['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI']:
    if len(chr_data[i])>0 :
        print '\n'.join([str(x) for x in chr_data[i]])

