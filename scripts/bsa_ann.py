#!/usr/bin/env python
#---------------import--------------------#
import argparse
import pandas as pd
import sys

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('-i','--inputfile1',help = 'avi input file1')
    parser.add_argument('-j','--inputfile2',help = 'avi input file2')
    parser.add_argument('-o','--outputfile',help = 'output file')
    args=parser.parse_args()
    inputfile1=args.inputfile1
    inputfile2=args.inputfile2
    outputfile=args.outputfile
    t1=pd.read_table(inputfile1, header=None, names=['P1_line', 'P1_variant_type', 'P1_anno', 'P1_chr',  'P1_start', 'P1_stop', 'P1_ref', 'P1_alt', 'P1_blank1', 'P1_blank2','P1_DP'])
    t2=pd.read_table(inputfile2, header=None, names=['P2_line', 'P2_variant_type', 'P2_anno', 'P2_chr',  'P2_start', 'P2_stop', 'P2_ref', 'P2_alt', 'P2_blank1', 'P2_blank2','P2_DP'])
    gene_id_P1=[];exon_P1=[];nt_P1=[];protein_P1=[];gene_id_link_P1=[]
    for attr1 in t1['P1_anno'].values:
        results = attr1.split(',')
        if results:
            gene_id=[];exon=[];nt=[];protein=[];gene_id_link=[]
            for i in results[:-1]:
                k = i.split(':')
                if len(k)!=5:
                    gene_id.append(k[0])
                    exon.append(k[2])
                    gene_id_link.append('<a href="http://www.yeastgenome.org/locus/' + k[0] + '/overview">' + k[0] + '</a>')
                else:
                    gene_id.append(k[0])
                    exon.append(k[2])
                    nt.append(k[3])
                    protein.append(k[4])
                    gene_id_link.append('<a href="http://www.yeastgenome.org/locus/' + k[0] + '/overview">' + k[0] + '</a>')
            gene_id_P1.append(';'.join(gene_id));exon_P1.append(';'.join(exon));nt_P1.append(';'.join(nt));protein_P1.append(';'.join(protein));gene_id_link_P1.append(';'.join(gene_id_link))
        else:
            gene_id_P1.append('NaN');exon_P1.append('NaN');nt_P1.append('NaN');protein_P1.append('NaN');gene_id_link_P1.append('NaN')
    t1['gene_id_P1']=gene_id_P1;t1['exon_P1']=exon_P1;t1['nt_P1']=nt_P1;t1['protein_P1']=protein_P1;t1['gene_id_link_P1']=gene_id_link_P1
    gene_id_P2=[];exon_P2=[];nt_P2=[];protein_P2=[];gene_id_link_P2=[]
    for attr1 in t2['P2_anno'].values:
        results = attr1.split(',')
        if results:
            gene_id=[];exon=[];nt=[];protein=[];gene_id_link=[]
            for i in results[:-1]:
                k = i.split(':')
                if len(k)!=5:
                    gene_id.append(k[0])
                    exon.append(k[2])
                    gene_id_link.append('<a href="http://www.yeastgenome.org/locus/' + k[0] + '/overview">' + k[0] + '</a>')
                else:
                    gene_id.append(k[0])
                    exon.append(k[2])
                    nt.append(k[3])
                    protein.append(k[4])
                    gene_id_link.append('<a href="http://www.yeastgenome.org/locus/' + k[0] + '/overview">' + k[0] + '</a>')
            gene_id_P2.append(';'.join(gene_id));exon_P2.append(';'.join(exon));nt_P2.append(';'.join(nt));protein_P2.append(';'.join(protein));gene_id_link_P2.append(';'.join(gene_id_link))
        else:
            gene_id_P2.append('NaN');exon_P2.append('NaN');nt_P2.append('NaN');protein_P2.append('NaN');gene_id_link_P2.append('NaN')
    t2['gene_id_P2']=gene_id_P2;t2['exon_P2']=exon_P2;t2['nt_P2']=nt_P2;t2['protein_P2']=protein_P2;t2['gene_id_link_P2']=gene_id_link_P2

    t1['P1_index'] =t1[['P1_chr','P1_blank1','P1_start']].astype(str).sum(axis=1)
    t2['P2_index'] =t2[['P2_chr','P2_blank1','P2_start']].astype(str).sum(axis=1)
    t1.index=t1['P1_index'].astype(str)
    t2.index=t2['P2_index'].astype(str)
    t1.index.is_unique
    t2.index.is_unique
    t1.index.duplicated
    t2.index.duplicated
    t11=t1[~t1.index.duplicated()]
    t22=t2[~t2.index.duplicated()]
    t33 = pd.concat([t11, t22], axis=1)
    t33['chr'] =[i.split('.')[0] for i in t33.index.astype(str).tolist()]
    t33['coord'] =[i.split('.')[1] for i in t33.index.astype(str).tolist()]
    t33['chr_number']  = t33['chr']     
    t33['chr_number']=t33['chr_number'].replace(["I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI"], range(1,17))
    aa=t33['coord'].astype(int).tolist()
    t33['coord_number']=aa
    t4=t33.sort(['chr_number', 'coord_number'], ascending=[True, True])
    t4['P1_chr']=t4['P1_chr'].fillna(t4['P2_chr'])
    t4['P1_start']=t4['P1_start'].fillna(t4['P2_start'])
    t4['P1_stop']=t4['P1_stop'].fillna(t4['P2_stop'])
    t4['P1_ref']=t4['P1_ref'].fillna(t4['P2_ref'])
    t4['P1_alt']=t4['P1_alt'].fillna(t4['P2_ref'])
    t4['P2_chr']=t4['P2_chr'].fillna(t4['P1_chr'])
    t4['P2_start']=t4['P2_start'].fillna(t4['P1_start'])
    t4['P2_stop']=t4['P2_stop'].fillna(t4['P1_stop'])
    t4['P2_ref']=t4['P2_ref'].fillna(t4['P1_ref'])
    t4['P2_alt']=t4['P2_alt'].fillna(t4['P1_ref'])
    t5=t4[~(t4['P1_alt']==t4['P2_alt'])]
    t5['P1_start']=t5['P1_start'].astype(int)
    t5['P1_stop']=t5['P1_stop'].astype(int)
    t5['P2_start']=t5['P2_start'].astype(int)
    t5['P2_stop']=t5['P2_stop'].astype(int)
    t5['P1_anno']=t5['P1_anno'].fillna('-')
    t5['P2_anno']=t5['P2_anno'].fillna('-')
    t5['P1_variant_type']=t5['P1_variant_type'].fillna('-')
    t5['P2_variant_type']=t5['P2_variant_type'].fillna('-')
    t5['gene_id_P1']=t5['gene_id_P1'].fillna('-')
    t5['exon_P1']=t5['exon_P1'].fillna('-')
    t5['nt_P1']=t5['nt_P1'].fillna('-')
    t5['protein_P1']=t5['protein_P1'].fillna('-')
    t5['gene_id_link_P1']=t5['gene_id_link_P1'].fillna('-')
    t5['gene_id_P2']=t5['gene_id_P2'].fillna('-')
    t5['exon_P2']=t5['exon_P2'].fillna('-')
    t5['nt_P2']=t5['nt_P2'].fillna('-')
    t5['protein_P2']=t5['protein_P2'].fillna('-')
    t5['gene_id_link_P2']=t5['gene_id_link_P2'].fillna('-')
    t6=t5[['chr',    'coord','P1_start','P1_stop','P1_ref','P1_alt','P2_alt','P1_variant_type',    'P1_anno','gene_id_P1',    'exon_P1','nt_P1','protein_P1','gene_id_link_P1','P2_variant_type',    'P2_anno','gene_id_P2',    'exon_P2','nt_P2','protein_P2','gene_id_link_P2'    ]]
    t6['P1_SNAP']=t6['P1_SNAP'].fillna('-')
    t6['P2_SNAP']=t6['P2_SNAP'].fillna('-')
    t6.to_csv(outputfile,sep='\t',index=False)
    
if __name__=='__main__':
    if len(sys.argv)>1:
        main()
    else:
        print 'Type "-h" or "--help" for more help'
        sys.exit(0)
