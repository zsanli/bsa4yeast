import sys

def main(sortedFile,n_file,bulk_design):
    inFile = open(sortedFile,'r')
    line=inFile.readline().strip().split('\t')
    if bulk_design==2:
        for i in range((n_file-2)/2):
            line.append('H'+str(i)+'_A0')
            line.append('H'+str(i)+'_A1')
        for i in range((n_file-2)/2):
            line.append('L'+str(i)+'_A0')
            line.append('L'+str(i)+'_A1')
        #line.append('HetMarker')
        print '\t'.join([str(x) for x in line])
        for line in inFile:
            data = line.strip().split('\t')
            for j in range((n_file-2)/2):
                h_nucl=data[18+j*8:18+j*8+4]
                data.append(h_nucl[int(data[n_file*8+3])-1])
                data.append(h_nucl[int(data[n_file*8+5])-1])
            for k in range((n_file-2)/2):
                l_nucl=data[18+(n_file-2)*4+k*8:18+(n_file-2)*4+k*8+4]
                data.append(l_nucl[int(data[n_file*8+3])-1])
                data.append(l_nucl[int(data[n_file*8+5])-1])
            print '\t'.join([str(x) for x in data])
    else:
        for i in range(n_file-2):
            line.append('B'+str(i)+'_A0')
            line.append('B'+str(i)+'_A1')
        print '\t'.join([str(x) for x in line])
        for line in inFile:
            data = line.strip().split('\t')
            for j in range(n_file-2):
                b_nucl=data[18+j*8:18+j*8+4]
                data.append(b_nucl[int(data[n_file*8+3])-1])
                data.append(b_nucl[int(data[n_file*8+5])-1])
            print '\t'.join([str(x) for x in data])
          
if __name__=='__main__':
    if len(sys.argv)>1:
        main(sys.argv[1],int(sys.argv[2]),int(sys.argv[3]))
    else:
        exit(0)