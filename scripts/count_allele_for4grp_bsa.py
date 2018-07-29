'''
@author: zhizhang
'''
import sys

def main(mpile,n_file):
    inFile = open(mpile,'r')
    print "chr\tbp","\tA\tG\tC\tT\tdel\tins\tinserted\tambiguous"*n_file
    for line in inFile:
        out=[]
        data = line.strip().split('\t')
        chr=data[0]
        bp = data[1]
        ref = data[2].upper()
        out=[chr,bp]
        for j in range(n_file):
            bases1 = data[3*j+4].upper()
            types1 = {'A':0,'G':0,'C':0,'T':0,'-':0,'+':[],'X':[]}
        
            i = 0
            while i < len(bases1):
                base= bases1[i]
                if base == '*' and len(bases1)==1:
                    types1['A']=0
                    types1['G']=0
                    types1['C']=0
                    types1['T']=0
                    types1['-']=0
                    types1['+']=[]
                    types1['X']=[]
                    i+=1
                    break
                elif base == '^' or base == '$':
                    i += 1
                elif base == '-':
                    i += 1
                elif base == '*':
                    types1['-'] += 1
                elif base == '+':
                    i += 1
                    try:
                        addNum = int(bases1[i])
                        addSeq = ''
                        for a in range(addNum):
                            i += 1
                            addSeq += bases1[i]
                            types1['+'].append(addSeq)
                    except ValueError:
                        pass
                elif base == '.' or base == ',':
                        types1[ref] += 1
                else:
                    if types1.has_key(base) and base != 'X':
 #                       print base
 #                       print types1
 #                       print types1[base]
 #                       print line
                        types1[base] += 1
                    else:
                        types1['X'].append(base)
                i += 1
            adds1 = '.'
            if len(types1['+']) > 0:
                adds1 = ','.join(types1['+'])
            amb1 = '.'
            if len(types1['X']) > 0:
                amb1 = ','.join(types1['X'])
            out1= [types1['A'],types1['G'],types1['C'],types1['T'],types1['-'],len(types1['+']),adds1,amb1]
            out=out+out1 
        print '\t'.join([str(x) for x in out])       
            
if __name__ == '__main__':
    if len(sys.argv)>1:
        main(sys.argv[1],int(sys.argv[2]))
