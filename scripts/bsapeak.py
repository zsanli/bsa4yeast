#!/usr/bin/python
import bsautil
import numpy as np
import bsadraw
import sys
import os
#print bsautil.compute_half_sample([1,2,3,4,5])
#print bsautil.half_sample_mode([1,2,3,4,5])
#print bsautil.identify_outliers(np.array([1,2,3,4,5]))
#RESULTSDIR=os.path.join(os.path.dirname(os.path.dirname(__file__)),'results')
def main(G_file,output_folder,prefix):
	gstat=[]
	infile=open(G_file,'r')
	while True:
		line =infile.readline().rstrip()
		if not line:
			break
		slast=float(line.split('\t')[-1])
		gstat.append(slast)
#	print '-'*20
#	print gstat
#	print '-'*20
	threshold_file=open(os.path.join(output_folder,prefix+'_threshold.txt'),'w')
	mu,delta= bsautil.robust_lognormal_mean_var(np.array(gstat))
	pcutoff, gcutoff= bsautil.estimate_smoothG_cutoff_lognorm_theory(np.array(gstat),mu,delta)
#	pcutoff, gcutoff=bsautil.estimate_smoothG_cutoff_lognorm_robust(np.array(gstat))
#	print infile
	gstat_file=bsadraw.load_Gstats(G_file)
	peak=bsautil.peaks(gstat_file,gcutoff)
#	print peak
#	print dir(peak)
	threshold_file.write('pcutoff\t'+'gcutoff\n'+str(pcutoff)+'\t'+str(gcutoff)+'\n')
	threshold_file.close()
	#print peak.__getitem__
	peak_file=open(os.path.join(output_folder,prefix+'_peak.txt'),'w')
	region_file=open(os.path.join(output_folder,prefix+'_region.txt'),'w')
	sym_table={0:'I',1:'II',2:'III',3:'IV',4:'V',5:'VI',6:'VII',7:'VIII',8:'IX',9:'X',10:'XI',11:'XII',12:'XIII',13:'XIV',14:'XV',15:'XVI'}
	line1=''
	line2=''
	for (i,j) in peak.items():
	#	print j
	#	print j.__class__
		if len(j)<>0:	
			for kk in j:
#				print kk.apex
#				print kk.left
#				print kk.right
				line1=sym_table[i]+'\t'+str(kk.apex)+'\t'+str(kk.apex)+'\n'
				line2=sym_table[i]+'\t'+str(kk.left)+'\t'+str(kk.right)+'\n'
				peak_file.write(line1)
				region_file.write(line2)
	peak_file.close()
	region_file.close()
	return pcutoff, gcutoff,os.path.join(output_folder,prefix+'_peak.txt'),os.path.join(output_folder,prefix+'_region.txt')
if __name__=='__main__':
	if len(sys.argv)>1:
		pcutoff, gcutoff,peak_file,region_file=main(sys.argv[1],sys.argv[2],sys.argv[3])
		print pcutoff, gcutoff,peak_file,region_file

		
		