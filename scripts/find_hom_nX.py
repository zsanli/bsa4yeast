import sys
n=int(sys.argv[2])
inFile = open(sys.argv[1],'r')
line=inFile.readline().strip().split('\t')
line.append('hom1')
line.append('pos1')
line.append('hom2')
line.append('pos2')
#line.append('HetMarker')
print '\t'.join([str(x) for x in line])
for line in inFile:
    data = line.strip().split('\t')
#    het_marker=''
    data1=[int(x) for x in data[2:6]]
    data2=[int(x) for x in data[10:14]]
    #print data1
    #print data2
    hom1=''
    pos1=''
    hom2=''
    pos2=''
#    print data1
  #  print data1
   # print data2
    if (int(data1[0])>n and int(data1[1])==0 and int(data1[2])==0 and int(data1[3])==0) :
        hom1='Pass'
        pos1=1
    elif (int(data1[0])==0 and int(data1[1])>n and int(data1[2])==0 and int(data1[3])==0):
        hom1='Pass'
        pos1=2
    elif (int(data1[0])==0 and int(data1[1])==0 and int(data1[2])>n and int(data1[3])==0):
        hom1='Pass'
        pos1=3
    elif (int(data1[0])==0 and int(data1[1])==0 and int(data1[2])==0 and int(data1[3])>n) :
        hom1='Pass'
        pos1=4
    else:
        hom1='No_Pass'
        pos1=0



    if (int(data2[0])>n and int(data2[1])==0 and int(data2[2])==0 and int(data2[3])==0):
        hom2='Pass'
        pos2=1
    elif (int(data2[0])==0 and int(data2[1])>n and int(data2[2])==0 and int(data2[3])==0):
        hom2='Pass'
        pos2=2
    elif  (int(data2[0])==0 and int(data2[1])==0 and int(data2[2])>n and int(data2[3])==0):
        hom2='Pass'
        pos2=3
    elif (int(data2[0])==0 and int(data2[1])==0 and int(data2[2])==0 and int(data2[3])>n):
        hom2='Pass'
        pos2=4
    else:
        hom2='No_Pass'
        pos2=0


    data.append(hom1)
    data.append(pos1)
    data.append(hom2)
    data.append(pos2)
    print '\t'.join([str(x) for x in data])

