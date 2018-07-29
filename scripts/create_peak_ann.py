#!/usr/bin/env python
#---------------import--------------------#
import argparse
import pandas as pd
import sys

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('-i','--anno_file',help = 'annotation file')
    parser.add_argument('-j','--peak_file',help = 'peak file')
    parser.add_argument('-o','--outputfile',help = 'output file')
    args=parser.parse_args()
    anno_file=args.anno_file
    peak_file=args.peak_file
    outputfile=args.outputfile
    
    tt = pd.read_table(anno_file, header=0, index_col=False)
    peak_table=pd.read_table(peak_file, header=None, index_col=False)
    peak_table.columns=['chr','start','end']
    peak_table_anno=pd.DataFrame()
    for i in range(len(peak_table)):
        tt_filter=tt.loc[(tt['chr'] == peak_table.loc[i]['chr']) & tt['coord'].isin(range(peak_table.loc[i]['start']-1000,peak_table.loc[i]['end']+1000))]
        if not tt_filter.empty:
            tt_filter.loc[:,('chr_peak')]=peak_table.loc[i]['chr']
            tt_filter.loc[:,('start_peak')]=peak_table.loc[i]['start']
            tt_filter.loc[:,('end_peak')]=peak_table.loc[i]['end']
            peak_table_anno=peak_table_anno.append(tt_filter,ignore_index=True)
    peak_table_anno.to_csv(outputfile,sep='\t',index=False)

if __name__=='__main__':
    if len(sys.argv)>1:
        main()
    else:
        print 'Type "-h" or "--help" for more help'
        sys.exit(0)
        
        
