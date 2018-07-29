#!/usr/bin/env python
#---------------import--------------------#
import argparse
import pandas as pd
import sys

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('-i','--anno_file',help = 'annotation file')
    parser.add_argument('-j','--region_file',help = 'region file')
    parser.add_argument('-o','--outputfile',help = 'output file')
    args=parser.parse_args()
    anno_file=args.anno_file
    region_file=args.region_file
    outputfile=args.outputfile
    
    tt = pd.read_table(anno_file, header=0, index_col=False)
    region_table=pd.read_table(region_file, header=None, index_col=False)
    region_table.columns=['chr','start','end']
    region_table_anno=pd.DataFrame()
    for i in range(len(region_table)):
        tt_filter=tt.loc[(tt['chr'] == region_table.loc[i]['chr']) & tt['coord'].isin(range(region_table.loc[i]['start'],region_table.loc[i]['end']))]
        if not tt_filter.empty:
            tt_filter.loc[:,('chr_region')]=region_table.loc[i]['chr']
            tt_filter.loc[:,('start_region')]=region_table.loc[i]['start']
            tt_filter.loc[:,('end_region')]=region_table.loc[i]['end']
            region_table_anno=region_table_anno.append(tt_filter,ignore_index=True)
    region_table_anno.to_csv(outputfile,sep='\t',index=False)

if __name__=='__main__':
    if len(sys.argv)>1:
        main()
    else:
        print 'Type "-h" or "--help" for more help'
        sys.exit(0)
        
        
