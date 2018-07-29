#!/usr/bin/env python
#---------------import--------------------#
import re
import sys
import os
import commands
import argparse
import collections
import operator

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('-i','--inputfile1',help = 'avi input file1')
    parser.add_argument('-j','--inputfile2',help = 'avi input file2')
    parser.add_argument('-b','--bedfile',help = 'bed file')
    parser.add_argument('-o','--outputfile',help = 'output file')
    args=parser.parse_args()
    inputfile1=args.inputfile1
    inputfile2=args.inputfile2
    bedfile=args.bedfile
    outputfile=args.outputfile
    
    chr_len=collections.OrderedDict()
    f_chr=open(bedfile,'r')
    chr=[]
    for line in f_chr:
	    i=line.rstrip().split('\t')
	    chr_len[i[0]]=int(i[2])
	    chr.append(i[0])
    print chr_len
    print chr
    inb1=open(inputfile1,'r')
    inb2=open(inputfile2,'r')
    outb=open(outputfile,'a')
    target=[]
    
    for i in inb1:
    	i_item=i.rstrip().split('\t')
    	#print i_item
    	i_item.append('P1')
    	target.append(i_item)
    	
    for i in inb2:
    	i_item=i.rstrip().split('\t')
    	i_item.append('P2')
    	target.append(i_item)
    #map(lambda x: int(x[4]),target)
    for i in target:
    	i[4]=int(i[4])
    
    
    
    target.sort(key = operator.itemgetter(3, 4))
#    print target
    print len(target)
    i=0
    while True:
    #for i in range(len(target)):
        print i
        if i==len(target)-1:
            if target[i][-1]=='P1':
                outb.write(target[i][3]+'\t'+str(target[i][4])+'\t'+target[i][6]+'\t'+target[i][7]+'\t'+target[i][6]+'\t'+target[i][1]+'\t'+target[i][2]+'\n')
            elif target[i][-1]=='P2':
                outb.write(target[i][3]+'\t'+str(target[i][4])+'\t'+target[i][6]+'\t'+target[i][6]+'\t'+target[i][7]+'\t'+target[i][1]+'\t'+target[i][2]+'\n')
            break
        elif target[i][4]==target[i+1][4]  and target[i][-1]=='P1':
            outb.write(target[i][3]+'\t'+str(target[i][4])+'\t'+target[i][6]+'\t'+target[i][7]+'\t'+target[i+1][7]+'\t'+target[i][1]+'\t'+target[i][2]+'\n')
            i=i+2
        elif target[i][4]==target[i+1][4]  and target[i][-1]=='P2':
            outb.write(target[i][3]+'\t'+str(target[i][4])+'\t'+target[i][6]+'\t'+target[i+1][7]+'\t'+target[i][7]+'\t'+target[i][1]+'\t'+target[i][2]+'\n')
            i=i+2
    	elif  target[i][4]<>target[i+1][4] and  target[i][-1]=='P1':
    	    outb.write(target[i][3]+'\t'+str(target[i][4])+'\t'+target[i][6]+'\t'+target[i][7]+'\t'+target[i][6]+'\t'+target[i][1]+'\t'+target[i][2]+'\n')
            i=i+1
    	elif  target[i][4]<>target[i+1][4] and  target[i][-1]=='P2':
            outb.write(target[i][3]+'\t'+str(target[i][4])+'\t'+target[i][6]+'\t'+target[i][6]+'\t'+target[i][7]+'\t'+target[i][1]+'\t'+target[i][2]+'\n')
            i=i+1
        else:
    		pass
    		
if __name__=='__main__':
    if len(sys.argv)>1:
        main()
    else:
        print 'Type "-h" or "--help" for more help'
        sys.exit(0)