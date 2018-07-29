import sys

def main(sortedFile,n_file):
    inFile = open(sortedFile,'r')
    line=inFile.readline().strip().split('\t')
    line.append('H_A0')
    line.append('H_A1')
    line.append('L_A0')
    line.append('L_A1')
    #line.append('HetMarker')
    print '\t'.join([str(x) for x in line])
    for line in inFile:
        data = line.strip().split('\t')
        h_nucl=data[18:22]
        l_nucl=data[26:30]
        ha0=h_nucl[int(data[35])-1]
        ha1=h_nucl[int(data[37])-1]
        la0=l_nucl[int(data[35])-1]
        la1=l_nucl[int(data[37])-1]
    	
        data.append(ha0)
        data.append(ha1)
        data.append(la0)
        data.append(la1)
        print '\t'.join([str(x) for x in data])
          
if __name__=='__main__':
    if len(sys.argv)>1:
        main(sys.argv[1],int(sys.argv[2]))
    else:
        exit(0)