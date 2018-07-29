#!/usr/bin/python
#---------------import--------------------#
import re
import sys
import os,glob
import commands
import argparse
import os.path

#---------------BIN path setting--------------------#
BIN=os.path.dirname(sys.argv[0])
if not re.search('/',BIN):
    BIN=os.getcwd()

fastqc=os.path.join(os.path.dirname(__file__),'FastQC.app/Contents/MacOS/fastqc')
fastx_trimmer=os.path.join(os.path.dirname(__file__),'fastx_trimmer')
fastq_quality_filter=os.path.join(os.path.dirname(__file__),'fastq_quality_filter')
Q20=os.path.join(os.path.dirname(__file__),'Q20_1.py')
rmN=os.path.join(os.path.dirname(__file__),'rmN.py')
ANNOVAR=os.path.join(os.path.dirname(__file__),'annotate_variation.pl')
scdb=os.path.join(os.path.dirname(os.path.dirname(__file__)),'genome/sacCer3')
convert2annovar=os.path.join(os.path.dirname(__file__),'convert2annovar.pl')
tableannovar=os.path.join(os.path.dirname(__file__),'table_annovar.pl')
ReadMap = os.path.join(os.path.dirname(__file__),'ReadMap.R')
sortCover=os.path.join(os.path.dirname(__file__),'sortCover2.py')
matePair=os.path.join(os.path.dirname(__file__),'rmUnpaired1.py')
cutadapt=os.path.join(os.path.dirname(__file__),'cutadapt')
coverDepth=''
GATK=os.path.join(os.path.dirname(__file__),'GenomeAnalysisTK.jar')
editPlot=os.path.join(os.path.dirname(__file__),'editplot.py')

n_pattern=re.compile(r'\s+(\d+)')
#---------------general function--------------------#
def getFile(path, regex):
    pattern='%s/*%s'%(path,regex)
    files = glob.glob(pattern)
    return files

def mkDir(path):
    if not os.path.exists(path):
        os.mkdir(path)

#---------------NGS analysis pipeline--------------------#
#---------------data cleaning--------------------#
#---------------trim short reads--------------------#
def Trim(fqlist,phred,lastkeep,out_dir):
    out_clean=out_dir+'/01.DataCleaning'
    print '-'*60+'\n'
    print 'Cleaning fastq is in:',out_clean
    mkDir(out_clean)
    trim_shell=''
    phred=int(phred)
    if len(fqlist)==2:
        fq_1=fqlist[0]
        fq_2=fqlist[1]
        trim_1=out_clean+'/trimmed_'+os.path.basename(fq_1)
        trim_2=out_clean+'/trimmed_'+os.path.basename(fq_2)
        trim_shell='%s -v -f 1 -l %s -i %s -o %s -Q %s \n'%(fastx_trimmer,lastkeep,fq_1,trim_1,phred)
        trim_shell+='%s -v -f 1 -l %s -i %s -o %s -Q %s \n'%(fastx_trimmer,lastkeep,fq_2,trim_2,phred)
    elif len(fqlist)==1:
        fq_1=fqlist[0]
        trim_1=out_clean+'/trimmed_'+os.path.basename(fq_1)
        trim_shell='%s -v -f 1 -l %s -i %s -o %s -Q %s \n'%(fastx_trimmer,lastkeep,fq_1,trim_1,phred)
    else:
        print >>stderr, 'There must be one or two fastq files!!!'
        sys.exit(0)
    return trim_shell

#---------------QC check--------------------#
def QC(fqlist,phred,out_dir):
    out_QC=out_dir+'/00.QC/'
    mkDir(out_QC)
    qc_shell=''
    if len(fqlist)==2:
        fq_1=fqlist[0]
        fq_2=fqlist[1]
        qc_shell='%s %s %s -o %s\n'%(fastqc,fq_1,fq_2,out_QC)
        qc_shell+='python2.7 %s %s %s %s \n'%(Q20,fq_1,out_QC,phred)
        qc_shell+='python2.7 %s %s %s %s \n'%(Q20,fq_2,out_QC,phred)
    elif len(fqlist)==1:
        fq_1=fqlist[0]
        qc_shell='%s %s -o %s\n'%(fastqc,fq_1,out_QC)
        qc_shell+='python2.7 %s %s %s %s \n'%(Q20,fq_1,out_QC,phred)
    else:
        print >>stderr,'There must be one or two fq files!!!'
        sys.exit(0)
    return qc_shell

#--------------cut adapter-------------------#
def rmAdapt(fq,adapter,overlap,min_len,out_dir):
    outfq=out_dir+'/rmadapter_'+os.path.basename(fq)
    rmadapt_log=out_dir+'/run.adapter.'+os.path.basename(fq)+'.log'
    rmadapt_shell='%s -a %s -O %s -m %s %s 1>%s 2>>%s '%(cutadapt,adapter,overlap,min_len,fq,outfq,rmadapt_log)
    return rmadapt_shell

#--------------cut N-------------------#
def removeN(fqlist,N_cutoff,out_dir):
    out_clean=out_dir+'/01.DataCleaning'
    mkDir(out_clean)
    if len(fqlist)==2:
        fq_1=fqlist[0]
        fq_2=fqlist[1]
        rmN_shell=''
        rmN_shell='python2.7 %s %s %s %s \n'%(rmN,fq_1,N_cutoff,out_clean) 
        rmN_shell+='python2.7 %s %s %s %s \n'%(rmN,fq_2,N_cutoff,out_clean) 
    elif len(fqlist)==1:
        fq_1=fqlist[0]
        rmN_shell=''
        rmN_shell='python2.7 %s %s %s %s \n'%(rmN,fq_1,N_cutoff,out_clean)
    else:
        print >>stderr,'There must be one or two fq files!!!'
        sys.exit(0)
    return rmN_shell

#--------------filterQ-------------------#
def filterQ(fqlist,phred, minQ,pminQ,out_dir):
    out_clean=out_dir+'/01.DataCleaning'
    mkDir(out_clean)
    phred=int(phred)
    if len(fqlist)==2:
        fq_1=fqlist[0]
        fq_2=fqlist[1]
        fq_1_basename=os.path.basename(fq_1)
        fq_2_basename=os.path.basename(fq_2)
        filtQ_1 = out_clean + '/filtQ_' + fq_1_basename
        filtQ_2 = out_clean + '/filtQ_' + fq_2_basename
        filtQ_1_log=filtQ_1+'.log'
        filtQ_2_log=filtQ_2+'.log'
        filterQ_shell=''
        filterQ_shell = '%s -v -q %s -p %s -i %s -o %s -Q %s >%s \n' %(fastq_quality_filter,minQ, pminQ, fq_1, filtQ_1,phred, filtQ_1_log)
        filterQ_shell += '%s -v -q %s -p %s -i %s -o %s -Q %s  >%s\n' %(fastq_quality_filter,minQ, pminQ, fq_2, filtQ_2,phred, filtQ_2_log)
    elif len(fqlist)==1:
        fq_1=fqlist[0]
        fq_1_basename=os.path.basename(fq_1)
        filtQ_1 = out_clean + '/filtQ_' + fq_1_basename
        filtQ_1_log=filtQ_1+'.log'
        filterQ_shell=''
        filterQ_shell = '%s -v -q %s -p %s -i %s -o %s -Q %s >%s \n' %(fastq_quality_filter,minQ, pminQ, fq_1, filtQ_1,phred, filtQ_1_log)
    else:
        print >>stderr,'There must be one or two fq files!!!'
        sys.exit(0)
    return filterQ_shell


#--------------remove unpaired reads---------------#
def rmUnpaired(fqlist,out_dir):
    out_clean=out_dir+'/01.DataCleaning'
    mkDir(out_clean)
    fq_1=fqlist[0]
    fq_2=fqlist[1]
    fq_1_out=out_clean+'/clean_'+os.path.basename(fq_1)
    fq_2_out=out_clean+'/clean_'+os.path.basename(fq_2)
    rmUnpaired_shell=''
    rmUnpaired_shell='python2.7 %s\t%s\t%s\t%s\t%s \n'%(matePair,fq_1,fq_2,fq_1_out,fq_2_out) 
    return rmUnpaired_shell

#----------------QC statistics------------------#
def qcStat(fqlist,out_dir):
    out_clean=out_dir+'/01.DataCleaning'
    Rstrip=open(out_clean+'/qcStat.R','w')
    fqlist=fqlist
    if len(fqlist)==2:
        if getFile(os.path.dirname(out_clean),'01.DataCleaning/rmadapter*.fq') or getFile(os.path.dirname(out_clean),'01.DataCleaning/rmadapter*.fastq'):
            fq_rmadapter=getFile(os.path.dirname(out_clean),'01.DataCleaning/rmadapter*.fq') or getFile(os.path.dirname(out_clean),'01.DataCleaning/rmadapter*.fastq')
            n_fq_rmadapter_left=re.findall(n_pattern,commands.getoutput('wc -l %s'%(fq_rmadapter[0])))[0]     
            n_fq_rmadapter_right=re.findall(n_pattern,commands.getoutput('wc -l %s'%(fq_rmadapter[1])))[0]
        else:
            n_fq_rmadapter_left=n_fq_rmadapter_right=0
        if getFile(os.path.dirname(out_clean),'01.DataCleaning/remained*.fq') or getFile(os.path.dirname(out_clean),'01.DataCleaning/remained*.fastq'):
            fq_remained=getFile(os.path.dirname(out_clean),'01.DataCleaning/remained*.fq') or getFile(os.path.dirname(out_clean),'01.DataCleaning/remained*.fastq')
            n_fq_remained_left=re.findall(n_pattern,commands.getoutput('wc -l %s'%(fq_remained[0])))[0]      
            n_fq_remained_right=re.findall(n_pattern,commands.getoutput('wc -l %s'%(fq_remained[1])))[0]     
        else:
            n_fq_remained_left=n_fq_remained_right=0
        if getFile(os.path.dirname(out_clean),'01.DataCleaning/filtQ*.fq') or getFile(os.path.dirname(out_clean),'01.DataCleaning/filtQ*.fastq'):
            fq_filtQ=getFile(os.path.dirname(out_clean),'01.DataCleaning/filtQ*.fq') or getFile(os.path.dirname(out_clean),'01.DataCleaning/filtQ*.fastq') 
            n_fq_filtQ_left=re.findall(n_pattern,commands.getoutput('wc -l %s'%(fq_filtQ[0])))[0]     
            n_fq_filtQ_right=re.findall(n_pattern,commands.getoutput('wc -l %s'%(fq_filtQ[1])))[0]    
        else:
            n_fq_filtQ_left=n_fq_filtQ_right=0
        if getFile(os.path.dirname(out_clean),'01.DataCleaning/clean*.fq') or getFile(os.path.dirname(out_clean),'01.DataCleaning/clean*.fastq'):
            fq_clean=getFile(os.path.dirname(out_clean),'01.DataCleaning/clean*.fq') or getFile(os.path.dirname(out_clean),'01.DataCleaning/clean*.fastq')
            n_fq_clean_left=re.findall(n_pattern,commands.getoutput('wc -l %s'%(fq_clean[0])))[0]      
            n_fq_clean_right=re.findall(n_pattern,commands.getoutput('wc -l %s'%(fq_clean[1])))[0]      
        else:
            n_fq_clean_left=n_fq_clean_right=0
        script='pdf(file="%s/qcStatistics.pdf")\n'% out_clean
        script+='barplot(c(%s,%s,%s,%s,%s,%s,%s,%s),names.arg=c("rmadapter_left","remained_left","filtQ_left","clean_left","rmadapter_right","remained_right","filtQ_right","clean_right"),las=2)\n'%(int(n_fq_rmadapter_left),int(n_fq_remained_left),int(n_fq_filtQ_left),int(n_fq_clean_left),int(n_fq_rmadapter_right),int(n_fq_remained_right),int(n_fq_filtQ_right),int(n_fq_clean_right))
        script+='dev.off()\n'
        Rstrip.write('%s\n'%script)
        Rstrip.close()
        Rplot=out_clean+'/qcStat.R'
        qcStat_shell='Rscript %s\n'% Rplot
    elif len(fqlist)==1:
        if getFile(os.path.dirname(out_clean),'01.DataCleaning/rmadapter*.fq') or getFile(os.path.dirname(out_clean),'01.DataCleaning/rmadapter*.fastq'):
            fq_rmadapter=getFile(os.path.dirname(out_clean),'01.DataCleaning/rmadapter*.fq') or getFile(os.path.dirname(out_clean),'01.DataCleaning/rmadapter*.fastq')
            n_fq_rmadapter_left=re.findall(n_pattern,commands.getoutput('wc -l %s'%(fq_rmadapter[0])))[0]       
        else:
            n_fq_rmadapter_left=0
        if getFile(os.path.dirname(out_clean),'01.DataCleaning/remained*.fq') or getFile(os.path.dirname(out_clean),'01.DataCleaning/remained*.fastq'):
            fq_remained=getFile(os.path.dirname(out_clean),'01.DataCleaning/remained*.fq') or getFile(os.path.dirname(out_clean),'01.DataCleaning/remained*.fastq')
            n_fq_remained_left=re.findall(n_pattern,commands.getoutput('wc -l %s'%(fq_remained[0])))[0]       
        else:
            n_fq_remained_left=0
        if getFile(os.path.dirname(out_clean),'01.DataCleaning/filtQ*.fq') or getFile(os.path.dirname(out_clean),'01.DataCleaning/filtQ*.fastq'):
            fq_filtQ=getFile(os.path.dirname(out_clean),'01.DataCleaning/filtQ*.fq') or getFile(os.path.dirname(out_clean),'01.DataCleaning/filtQ*.fastq') 
            n_fq_filtQ_left=re.findall(n_pattern,commands.getoutput('wc -l %s'%(fq_filtQ[0])))[0]    
        else:
            n_fq_filtQ_left=0
        if getFile(os.path.dirname(out_clean),'01.DataCleaning/clean*.fq') or getFile(os.path.dirname(out_clean),'01.DataCleaning/clean*.fastq'):
            fq_clean=getFile(os.path.dirname(out_clean),'01.DataCleaning/clean*.fq') or getFile(os.path.dirname(out_clean),'01.DataCleaning/clean*.fastq')
            n_fq_clean_left=re.findall(n_pattern,commands.getoutput('wc -l %s'%(fq_clean[0])))[0]    
        else:
            n_fq_clean_left=0
        print "****************"
        print "n_fq_rmadapter_left",n_fq_rmadapter_left,"n_fq_remained_left",n_fq_remained_left,"n_fq_filtQ_left",n_fq_filtQ_left,"n_fq_clean_left",n_fq_clean_left
        print "****************"
        script='pdf(file="%s/qcStatistics.pdf")\n'% out_clean
        script+='barplot(c(%s,%s,%s,%s),names.arg=c("rmadapter_left","remained_left","filtQ_left","clean_left"),las=2)\n'%(int(n_fq_rmadapter_left),int(n_fq_remained_left),int(n_fq_filtQ_left),int(n_fq_clean_left))
        script+='dev.off()\n'
        Rstrip.write('%s\n'%script)
        Rstrip.close()
        Rplot=out_clean+'/qcStat.R'
        qcStat_shell='Rscript %s\n'% Rplot
    else:
        print >>stderr,'There must be one or two fq files!!!'
        sys.exit(0)
    return qcStat_shell
    
#----------------mapping------------------#
def mapping(fqlist,rgid,species,ref, mismatch, gap_open, gap_extention, insize, out_dir):
    out_aln=out_dir+'/02.Alignment'
    mkDir(out_aln)
    mapping_shell=''
    if len(fqlist)==2:
        fq_1=fqlist[0]
        fq_2=fqlist[1]
        fq_1_basename=os.path.basename(fq_1)
        fq_2_basename=os.path.basename(fq_2)
        species=species
#       fq_1_basename=fq_1_basename.replace('clean_filtQ_remained_rmadapter_trimmed_','')
#       fq_2_basename=fq_2_basename.replace('clean_filtQ_remained_rmadapter_trimmed_','')
        sai1 =out_aln+'/'+fq_1_basename+'.sai'
        sai2 =out_aln+'/'+fq_2_basename+'.sai'
        sam=out_aln+'/'+fq_1_basename+'.sam'
        shell_aln_1 = 'bwa aln  -n %s -o %s -e %s -m 100000 -t 6 -q 20  %s %s>%s \n' %(mismatch, gap_open, gap_extention, ref, fq_1, sai1)
        shell_aln_2 = 'bwa aln  -n %s -o %s -e %s -m 100000 -t 6 -q 20  %s %s>%s\n' %(mismatch, gap_open, gap_extention, ref, fq_2, sai2)
        shell_sam = 'bwa sampe -r "@RG\\tID:%s\\tLB:%s\\tPL:illumina\\tSM:%s\\tPI:%s" -a %s %s %s %s %s %s>%s\n' %(rgid,rgid,rgid, insize,insize,  ref,sai1, sai2,fq_1, fq_2, sam)
        mapping_shell = shell_aln_1 + shell_aln_2 + shell_sam
        bam = sam.replace('.sam', '.bam')
        sortbam_prefix = sam.replace('.sam', '.sort')
        sortbam = sam.replace('.sam', '.sort.bam')
        rmdupbam = sortbam.replace('.sort.bam', '.rmdup.sort.bam')
        uniqbam = sortbam.replace('.sort.bam', '.uniq.rmdup.sort.bam')
        rmduplog = rmdupbam + '.log'
        shell_bam = 'samtools view -Sb %s > %s\n' %(sam, bam)
        shell_sort = 'samtools sort %s >%s\n' % (bam, sortbam)
        shell_index_sort = 'samtools index %s\n' % sortbam
        shell_stat_raw = 'samtools flagstat %s >>%s\n' % (sortbam, rmduplog)
        shell_rmdup = 'samtools rmdup -sS %s %s\n' % (sortbam, rmdupbam)
        shell_stat_rmdup = 'samtools flagstat %s >>%s\n' %(rmdupbam, rmduplog)
        shell_uniq = "samtools view -h %s |awk '$0~/^@/||$0~/XT:A:U/||$0~/XT:A:M/{print $0}'|samtools view -Sb - >%s\n" %(rmdupbam, uniqbam)
        shell_stat_uniq = 'samtools flagstat %s >>%s\n' %(uniqbam, rmduplog)
        shell_index_uniq = 'samtools index %s\n' % uniqbam
        shell_samtools = shell_bam + shell_sort + shell_index_sort + shell_stat_raw + shell_rmdup + shell_stat_rmdup + shell_uniq + shell_stat_uniq + shell_index_uniq
    elif len(fqlist)==1:
        fq_1=fqlist[0]
        fq_1_basename=os.path.basename(fq_1)
        species=species
#       fq_1_basename=fq_1_basename.replace('clean_filtQ_remained_rmadapter_trimmed_','')
        sai1 =out_aln+'/'+fq_1_basename+'.sai'
        sam=out_aln+'/'+fq_1_basename+'.sam'
        shell_aln_1 = 'bwa aln  -n %s -o %s -e %s -m 100000 -t 6 -q 20  %s %s>%s \n' %(mismatch, gap_open, gap_extention, ref, fq_1, sai1)
        shell_sam = 'bwa samse -r "@RG\\tID:%s\\tLB:%s\\tPL:illumina\\tSM:%s\\tPI:%s"  %s %s %s >%s\n' %(rgid,rgid,rgid, insize, ref,sai1,fq_1, sam)
        mapping_shell = shell_aln_1 + shell_sam
        bam = sam.replace('.sam', '.bam')
        sortbam_prefix = sam.replace('.sam', '.sort')
        sortbam = sam.replace('.sam', '.sort.bam')
        rmdupbam = sortbam.replace('.sort.bam', '.rmdup.sort.bam')
        uniqbam = sortbam.replace('.sort.bam', '.uniq.rmdup.sort.bam')
        rmduplog = rmdupbam + '.log'
        shell_bam = 'samtools view -Sb %s > %s\n' %(sam, bam)
        shell_sort = 'samtools sort %s >%s\n' % (bam, sortbam)
        shell_index_sort = 'samtools index %s\n' % sortbam
        shell_stat_raw = 'samtools flagstat %s >>%s\n' % (sortbam, rmduplog)
        shell_rmdup = 'samtools rmdup -sS %s %s\n' % (sortbam, rmdupbam)
        shell_stat_rmdup = 'samtools flagstat %s >>%s\n' %(rmdupbam, rmduplog)
        shell_uniq = "samtools view -h %s |awk '$0~/^@/||$0~/XT:A:U/||$0~/XT:A:M/{print $0}'|samtools view -Sb - >%s\n" %(rmdupbam, uniqbam)
        shell_stat_uniq = 'samtools flagstat %s >>%s\n' %(uniqbam, rmduplog)
        shell_index_uniq = 'samtools index %s\n' % uniqbam
        shell_samtools = shell_bam + shell_sort + shell_index_sort + shell_stat_raw + shell_rmdup + shell_stat_rmdup + shell_uniq + shell_stat_uniq + shell_index_uniq
    else:
        print >>stderr,'There must be one or two fq files!!!'
        sys.exit(0)    
    return mapping_shell+shell_samtools
    
#def callVcf(ref,out_dir):  
#    out_var = out_dir  + '/03.Variants'
#    out_aln = out_dir + '/02.Alignment'
#    mkDir(out_var)
#    sam=getFile(out_aln,'.sam')
#    sam=sam[0]
#    bam = sam.replace('.sam', '.bam')
#    sortbam_prefix = sam.replace('.sam', '.sort')
#    sortbam = sam.replace('.sam', '.sort.bam')
#    rmdupbam = sortbam.replace('.sort.bam', '.rmdup.sort.bam')
#    uniqbam = sortbam.replace('.sort.bam', '.uniq.rmdup.sort.bam')
#    print '-'*60
#    print 'clean bam file is:',uniqbam
#    print '-'*60
#    bamname = os.path.basename(uniqbam)
#    mpileup = out_var + '/' + bamname.replace('.bam', '.mpileup')
#    vcfFilt = out_var + '/' + bamname.replace('.bam', '.vcf')
#    realign_interval=out_var + '/' + bamname.replace('.bam', '.realign.intervals')
#    bam_mpileup = 'samtools mpileup -Os -f %s %s > %s \n ' %(ref, uniqbam, mpileup)
##  vcf_mpileup = 'samtools mpileup -m 5 -u -f %s %s |bcftools view -cgv - | %s varFilter -D 100 - >%s\n' %(ref, rmdupbam, vcfutils, vcfFilt)
#    #vcf_mpileup = 'java -jar %s -T UnifiedGenotyper --genotype_likelihoods_model BOTH -rf BadCigar -R %s -I %s -o %s\n' %(GATK, ref, uniqbam, vcfFilt)
#    vcf_mpileup='java -jar %s -T RealignerTargetCreator -rf BadCigar -R %s -I %s -o %s\n'%(GATK,ref,uniqbam,realign_interval)
#    realign_interval_bam=out_var + '/' + bamname.replace('.bam', '.realign.intervals.bam')
#    vcf_mpileup+='java -jar %s -T IndelRealigner -R %s -targetIntervals %s -I %s -o %s \n'%(GATK,ref,realign_interval,uniqbam,realign_interval_bam)
#    known_VCF=out_var + '/'+'known.vcf'
#    known_INDEL=out_var+'/'+'known_INDEL.vcf'
#    known_SNP=out_var+'/'+'known_SNP.vcf'
#    knownsite_shell='java -jar %s -T UnifiedGenotyper --genotype_likelihoods_model BOTH -rf BadCigar -R %s -I %s -o %s \n'%(GATK,ref,uniqbam,known_VCF)
#    knownsite_shell+='java -jar %s -T UnifiedGenotyper --genotype_likelihoods_model INDEL -rf BadCigar -R %s -I %s -o %s \n'%(GATK,ref,uniqbam,known_INDEL)
#    knownsite_shell+='java -jar %s -T UnifiedGenotyper --genotype_likelihoods_model SNP -rf BadCigar -R %s -I %s -o %s \n'%(GATK,ref,uniqbam,known_SNP)
#    recal_table=out_var + '/' + bamname.replace('.bam', '.recal.table')
#    recal_post_table=out_var + '/' + bamname.replace('.bam', '.recal.post.table')
#    BQSR_shell='java -Xmx2g -jar %s -T BaseRecalibrator -R %s -I %s -knownSites %s -o %s \n'%(GATK,ref,realign_interval_bam,known_VCF,recal_table)
#    BQSR_shell+='java -Xmx2g -jar %s -T BaseRecalibrator -R %s -I %s  -knownSites %s -BQSR %s  -o %s \n'%(GATK,ref,realign_interval_bam,known_VCF,recal_table,recal_post_table)
#    recal_pdf=out_var + '/' + bamname.replace('.bam', '.pdf')
#    BQSR_shell+='java -Xmx2g -jar %s -T AnalyzeCovariates -l DEBUG -R %s -before %s -after %s -plots %s \n'%(GATK,ref,recal_table,recal_post_table,recal_pdf)
#    recal_bam=out_var + '/' + bamname.replace('.bam', '.recal.bam')
#    BQSR_shell+='java -Xmx2g -jar %s -T PrintReads -R %s -I %s -BQSR %s -o %s \n'%(GATK,ref,realign_interval_bam,recal_table,recal_bam)
#    GVCF=out_var + '/' + bamname.replace('.bam', '.hc.pairHMM.gvcf')
#    BQSR_shell+='java -Xmx8g -jar %s -T HaplotypeCaller -rf BadCigar -R %s -I %s -pairHMM VECTOR_LOGLESS_CACHING -nct 2 --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o %s \n'%(GATK,ref,recal_bam,GVCF)
#    VCF=out_var + '/' + bamname.replace('.bam', '.hc.pairHMM.vcf')
#    gvcf2vcf_shell='java -Xmx2g -jar %s -T GenotypeGVCFs -R %s -o %s -V %s \n'%(GATK,ref,VCF,GVCF)
#    callVcf_shell = bam_mpileup + vcf_mpileup+knownsite_shell+BQSR_shell+gvcf2vcf_shell    
#    return callVcf_shell
    
#----------------Coverage------------------#
def coverGenome(bam, seqlen, genomebed,windowsize,out_dir):
    out_cover=out_dir+'/04.Coverage'
    cover_shell = ''
    mkDir(out_cover)
    cover_genome =out_cover + '/' + os.path.basename(bam).replace('bam', 'cover.genome')
    cover_window=out_cover+'/'+os.path.basename(bam).replace('bam','cover.window')
    window_bed=out_cover+'/'+os.path.basename(bam).replace('bam','window.bed')
    make_window_shell='bedtools makewindows -b %s  -w %s  >%s\n'%(genomebed,windowsize,window_bed)
    cover_shell= 'bedtools coverage -abam %s -b %s >%s\n' %(bam, genomebed, cover_genome)
    cover_shell+='bedtools coverage -abam %s -b %s > %s\n' %(bam,window_bed,cover_window)   
    sort_cover=out_cover+'/'+os.path.basename(bam).replace('bam','sorted.cover')
    pdf=out_cover+'/'+os.path.basename(bam).replace('bam','cover.pdf')
    sort_shell ='python2.7 %s %s %s\n'% (sortCover,cover_window,sort_cover)
    plot_cover_shell = 'Rscript %s --args %s %s\n' %(ReadMap, sort_cover, pdf)
    return make_window_shell+cover_shell+sort_shell+plot_cover_shell
    
#def annotate(vcf,annodb,species,out_dir):
#    out_anno=out_dir+'/05.Annotation'
#    mkDir(out_anno)
#    anno_shell=''
#    #vcf_file=getFile(out_dir +'/03.Variants','vcf')[0]    
#    vcf_file_basename=os.path.basename(vcf)
#    avinput_file=out_dir +'/05.Annotation/'+vcf_file_basename.replace('vcf','avinput')
#    #convert_annovar_shell='perl %s --includeinfo --allallele -format vcf4 %s > %s \n' %(convert2annovar,vcf,avinput_file)
#    convert_annovar_shell='perl %s -format vcf4 %s > %s \n'%(convert2annovar,vcf,avinput_file)
#    anno_shell= 'perl %s %s %s -buildver %s\n' %(ANNOVAR, avinput_file, annodb, species)
#    varfuction=avinput_file+'.variant_function'
#    dat = out_anno + '/dat.txt'
#    staread = "awk '{print $1}' %s |sort |uniq -c > %s\n" %(varfuction, dat)
#    editfile = out_anno + '/plot.dat'
#    editplot_shell= 'python2.7 %s %s %s\n' %(editPlot, dat, editfile)
#    anno_shell =  convert_annovar_shell + anno_shell + staread + editplot_shell
#    Rstrip = open(out_anno + '/plotpie.R', 'w')
#    script = 'library("plotrix")\n'
#    script += "pdf('%s/readsmap.pdf')\n"  %(out_anno)
#    script += 'dat=read.table("%s")\n' %editfile
#    script += 'ratio=sprintf("%.2f",100*dat[,2]/sum(dat[,2]))\n'
#    script += 'ratio=paste(ratio,"%",sep="")\n'
#    script += 'label=paste(dat[,1],ratio,sep="\\n")\n'
#    script += 'pie3D(dat[,2],col=rainbow(10),main="", border="black", labels=label,font=2,labelcex=1,,explode=0.1,radius=0.95)\n'
#    Rstrip.write('%s\n' %script)
#    Rstrip.close()
#    Rstrip = out_anno + '/plotpie.R'
#    Rplot = 'Rscript %s\n' %Rstrip
#    anno_shell += Rplot
#    return anno_shell

#---------------main function-------------#
def main():  
    parser=argparse.ArgumentParser()
    parser.add_argument('-i','--indir',help = 'fastq files path')
    parser.add_argument('-o','--outdir',help = 'output path')
    parser.add_argument('-T','--trim',help = 'trim given fastq short reads as --lastkeep defined')
    parser.add_argument('-l','--lastkeep',help = 'with "trim" option, last bases to keep',default=100,type=int)
    parser.add_argument('-P','--phred',help = 'phred score used in platform [33]',default=33,type=int)    
    parser.add_argument('-Q','--qccheck',help='do quality check [true]',default='true')
    parser.add_argument('-r','--rmadapt',help='remove adapter [true]',default='true')
    parser.add_argument('-L','--ladapter',help='left adapter [AGATCGGAAGAGC]',default='AGATCGGAAGAGC')
    parser.add_argument('-R','--radapter',help='right adapter [AGATCGGAAGAGC]',default='AGATCGGAAGAGC')
    parser.add_argument('-O','--overlap',help='If the overlap between the read and adapter is shorter than the overlap length, the read will NOT be modified. [6]',default=6, type=int)
    parser.add_argument('-m','--minlen', help = 'Discard trimmed reads that are shorter than "minlen" [75]', default = 75,type=int)
    parser.add_argument('-N', '--removeN', help='remove "N" bases [true]', default = 'true')
    parser.add_argument('-c','--Ncutoff', help='with "removeN" option, N cutoff [0.1]', default = "0.1",type=float)
    parser.add_argument('-F', '--filtQ', help = 'Filters sequences based on quality [true]', default = 'true')
    parser.add_argument( '--minQ', help = 'Minimum quality score to keep [20]', default = 20,type=int)
    parser.add_argument( '--pminQ', help = 'Minimum percent of bases [80]', default = 80,type=int)
    parser.add_argument( '-U','--rmUnpaired', help = 'remove unpaied reads [true]', default = 'true')   
    parser.add_argument('-q','--qcStat',help='generate QC statistic plot for seq',default='true')
    parser.add_argument('-M', '--map', help='read mapping [true]', default = 'true')
    parser.add_argument('--ref', help='fasta to use as reference')
    parser.add_argument('--insertsize', help = 'insert size [400]', default = 400,type=int)
    parser.add_argument( '--mismatch', help = 'max #diff (int) or missing prob under 0.02 err rate (float) [0.04]', default = 0.04,type=float)
    parser.add_argument( '--gapopen', help = 'maximum number or fraction of gap opens [1]', default = 1,type=int)
    parser.add_argument( '--gapextention', help = 'maximum number of gap extensions, -1 for disabling long gaps [-1]', default = -1,type=int)
    parser.add_argument('-G', '--genome', help = 'coverage of genome region [false]', default ='true')
    parser.add_argument('-E', '--exome', help = 'coverage of exome and near 200bp region [false]', default ='false')
    parser.add_argument( '--readlen', help = 'read length [100]', default = 100,type=int)
    parser.add_argument( '--callVcf', help = 'generate mileup and Vcf file', default = 'true')
    parser.add_argument( '--genomebed', help = 'coverage genome bed file')      
    parser.add_argument( '--windowsize', help = 'coverage window size [500000]', default = 500000,type=int)
    parser.add_argument( '--species', help = 'species', default = 'hg19')
    args=parser.parse_args()
    in_dir=args.indir
    out_dir=args.outdir
    mkDir(out_dir)
    RGID=os.path.basename(in_dir)
    run_log=open(out_dir+'/reseq_login.sh','a')
    fqlist=getFile(in_dir,'.fq') or getFile(in_dir,'.fastq')
    fqlist=sorted(fqlist)
    print '-'*60
    print 'Original fastq list is :',fqlist
    
    if re.search('true',args.trim,re.I):
        shell_trim=Trim(fqlist,args.phred,args.lastkeep,out_dir)
        run_log.write('-'*60+'\n')
        run_log.write('Trim starting at %s\n' % commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"'))
        run_log.write('-'*60+'\n')
        run_log.close()
        print '-'*60+'\n'
        print 'Trim starting at %s\n' % commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"')
        print '-'*60+'\n'
        status=os.system(shell_trim)
        if status==0:
            run_log=open(out_dir+'/reseq_login.sh','a')
            run_log.write(shell_trim)
            run_log.write('-'*60+'\n')
            run_log.write('Trim finished at %s\n'%commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"'))
            run_log.write('-'*60+'\n')
            run_log.close()
            print '-'*60+'\n'
            print 'Trim finished at %s\n'%commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"')  
            print '-'*60+'\n'
        else:
            run_log=open(out_dir+'/reseq_login.sh','a')
            run_log.write ('-'*60+'\n')
            run_log.write ('Trim Error!!! at %s\n'%commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"'))
            run_log.write  ('-'*60+'\n')
            run_log.close()
            print '-'*60+'\n'
            print 'Trim Error!!! at %s\n'%commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"')
            print '-'*60+'\n'
        out_clean=out_dir+'/01.DataCleaning'
        if getFile(out_clean,'trimmed*fq') or getFile(out_clean,'trimmed*fastq'):
            fqlist=getFile(out_clean,'trimmed*fq') or getFile(out_clean,'trimmed*fastq')
        print '-'*60+'\n'
        print 'Trimmed fastq file is:',fqlist

    if re.search('true',args.qccheck,re.I):
        qc_shell=QC(fqlist,args.phred,out_dir)
        run_log=open(out_dir+'/reseq_login.sh','a')
        run_log.write('-'*60+'\n')
        run_log.write('QC starting at %s\n'%commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"'))
        run_log.write('-'*60+'\n')
        run_log.close()
        print '-'*60+'\n'
        print 'QC starting at %s\n'%commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"')
        print '-'*60+'\n'
        status=os.system(qc_shell)
        if status==0:
            run_log=open(out_dir+'/reseq_login.sh','a')
            run_log.write('-'*60+'\n')
            run_log.write(qc_shell)
            run_log.write('-'*60+'\n')
            run_log.write('QC finished at %s \n'%commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"'))
            run_log.write('-'*60+'\n')
            run_log.close()
            print '-'*60+'\n'
            print 'QC finished at %s \n'%commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"')
            print '-'*60+'\n'
        else:
            run_log=open(out_dir+'/reseq_login.sh','a')
            run_log.write('-'*60+'\n')
            run_log.write('%s\n' % qc_shell)
            run_log.write('-'*60+'\n')
            run_log.write('QC Error!!!\n%s\n' %commands.getoutput('date'))
            run_log.write('-'*60+'\n')
            run_log.close()
            print '-'*60+'\n'
            print 'QC Error!!!\n%s\n' %commands.getoutput('date')
            print '-'*60+'\n'
#   sys.exit(0) 

    if re.search('true', args.rmadapt, re.I):
        leftadapt=args.ladapter
        rightadapt=args.radapter
        min_len=args.minlen
        overlen=int(args.overlap)
        out_qc=out_dir+'/01.DataCleaning'
        mkDir(out_qc)
        rmadapt_shell=''
        if len(fqlist)==2:
            rmadapt_shell=rmAdapt(fqlist[0],leftadapt,overlen,min_len,out_qc)+' \n'
            rmadapt_shell+=rmAdapt(fqlist[1],rightadapt,overlen,min_len,out_qc)+'\n'
            run_log=open(out_dir+'/reseq_login.sh','a')
            run_log.write('-'*60+'\n')
            run_log.write('rmadpt starting at %s\n' %commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"'))
            run_log.write('-'*60+'\n')
            run_log.close()
            print '-'*60+'\n'
            print 'rmadpt starting at %s\n' %commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"')
            print '-'*60+'\n'
        elif len(fqlist)==1:
            rmadapt_shell=rmAdapt(fqlist[0],leftadapt,overlen,min_len,out_qc)+'\n'
            run_log=open(out_dir+'/reseq_login.sh','a')
            run_log.write('-'*60+'\n')
            run_log.write('rmadpt starting at %s\n' %commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"'))
            run_log.write('-'*60+'\n')
            run_log.close()
            print '-'*60+'\n'
            print 'rmadpt starting at %s\n' %commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"')
            print '-'*60+'\n'
        else:
            print >>stderr,'There must be one or two fq files!!!'
            sys.exit(0)     
        status = os.system(rmadapt_shell)
        if status==0:
            run_log=open(out_dir+'/reseq_login.sh','a')
            run_log.write('-'*60+'\n')
            run_log.write(rmadapt_shell)
            run_log.write('-'*60+'\n')
            run_log.write('rmadapt finished at %s \n'%commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"'))
            run_log.write('-'*60+'\n')
            run_log.close()
            print '-'*60+'\n'
            print 'rmadapt finished at %s \n'%commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"')
            print '-'*60+'\n'
        else:
            run_log=open(out_dir+'/reseq_login.sh','a')
            run_log.write('-'*60+'\n')
            run_log.write('%s\n' % rmadapt_shell)
            run_log.write('-'*60+'\n')
            run_log.write('rmadpt Error!!!\n%s\n' %commands.getoutput('date'))
            run_log.write('-'*60+'\n')
            run_log.close()
            print '-'*60+'\n'
            print 'rmadapt Error!!!\n%s\n' %commands.getoutput('date')
            print '-'*60+'\n'
        out_clean=out_dir+'/01.DataCleaning'
        if getFile(out_clean,'rmadapter*fq') or getFile(out_clean,'rmadapter*fastq'):
            fqlist=getFile(out_clean,'rmadapter*fq') or getFile(out_clean,'rmadapter*fastq')
        print '-'*60+'\n'
        print 'rmadapter fastq file is:',fqlist
    
    if re.search('true', args.removeN, re.I):
        #print 'removeN starting at %s\n'
        N_cutoff=args.Ncutoff
        #print  N_cutoff
        removeN_shell=removeN(fqlist,N_cutoff,out_dir)
        #print removeN_shell
        run_log=open(out_dir+'/reseq_login.sh','a')
        run_log.write('-'*60+'\n')
        run_log.write('removeN starting at %s\n' %commands.getoutput('date'))
        run_log.write('-'*60+'\n')
        run_log.close()
        print '-'*60+'\n'
        print 'removeN starting at %s\n' %commands.getoutput('date')
        print '-'*60+'\n'
        status = os.system(removeN_shell)
        if status == 0:
            run_log=open(out_dir+'/reseq_login.sh','a')
            run_log.write('-'*60+'\n')
            run_log.write(removeN_shell)
            run_log.write('-'*60+'\n')
            run_log.write('rmoveN finished at %s \n'%commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"'))
            run_log.write('-'*60+'\n')
            run_log.close()
            print '-'*60+'\n'
            print 'rmoveN finished at %s \n'%commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"')
            print '-'*60+'\n'
        else:
            run_log=open(out_dir+'/reseq_login.sh','a')
            run_log.write('-'*60+'\n')
            run_log.write('%s\n' % removeN_shell)
            run_log.write('-'*60+'\n')
            run_log.write('rmoveN Error!!!\n%s\n' %commands.getoutput('date'))
            run_log.write('-'*60+'\n')
            run_log.close()
            print '-'*60+'\n'
            print 'rmoveN Error!!!\n%s\n' %commands.getoutput('date')
            print '-'*60+'\n'
        out_clean = out_dir + '/01.DataCleaning'
        if getFile(out_clean, 'remained*fq') or getFile(out_clean, 'remained*fastq'):
            fqlist = getFile(out_clean, 'remained*fq') or getFile(out_clean, 'remained*fastq')
        print '-'*60
        print 'removeN remained fq file is: ',fqlist
 
    if re.search('true', args.filtQ, re.I):
        phred=args.phred
        minQ=args.minQ
        pminQ=args.pminQ
        filterQ_shell=filterQ(fqlist,phred,minQ,pminQ,out_dir)
       # print filterQ_shell
        run_log=open(out_dir+'/reseq_login.sh','a')
        run_log.write('-'*60+'\n')
        run_log.write('filterQ starting at %s\n' %commands.getoutput('date'))
        run_log.write('-'*60+'\n')
        run_log.close()
        print '-'*60+'\n'
        print 'filterQ starting at %s\n' %commands.getoutput('date')
        print '-'*60+'\n'
        status = os.system(filterQ_shell)
        if status == 0:
            run_log=open(out_dir+'/reseq_login.sh','a')
            run_log.write('-'*60+'\n')
            run_log.write(filterQ_shell)
            run_log.write('-'*60+'\n')
            run_log.write('filterQ finished at %s \n'%commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"'))
            run_log.write('-'*60+'\n')
            run_log.close()
            print '-'*60+'\n'
            print 'filterQ finished at %s \n'%commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"')
            print '-'*60+'\n'
        else:
            run_log=open(out_dir+'/reseq_login.sh','a')
            run_log.write('-'*60+'\n')
            run_log.write('%s\n' % filterQ_shell)
            run_log.write('-'*60+'\n')
            run_log.write('filterQ Error!!!\n%s\n' %commands.getoutput('date'))
            run_log.write('-'*60+'\n')
            run_log.close()
            print '-'*60+'\n'
            print 'filterQ Error!!!\n%s\n' %commands.getoutput('date')
            print '-'*60+'\n'
        out_clean = out_dir + '/01.DataCleaning'
        if getFile(out_clean, 'filtQ*fq') or getFile(out_clean, 'filtQ*fastq'):
            fqlist = getFile(out_clean, 'filtQ*fq') or getFile(out_clean, 'filtQ*fastq')
        print '-'*60
        print 'filterQ fqlist is: ',fqlist

    if re.search('true', args.rmUnpaired, re.I):
        if len(fqlist)==2:
            rmUnpaired_shell=rmUnpaired(fqlist,out_dir)
            run_log=open(out_dir+'/reseq_login.sh','a')
            run_log.write('-'*60+'\n')
            run_log.write('rmUnpaired starting at %s\n' %commands.getoutput('date'))
            run_log.write('-'*60+'\n')
            run_log.close()
            print '-'*60+'\n'
            print 'rmUnpaired starting at %s\n' %commands.getoutput('date')
            print '-'*60+'\n'
            status = os.system(rmUnpaired_shell)
            if status == 0:
                run_log=open(out_dir+'/reseq_login.sh','a')
                run_log.write('-'*60+'\n')
                run_log.write(rmUnpaired_shell)
                run_log.write('-'*60+'\n')
                run_log.write('rmUnpaired finished at %s \n'%commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"'))
                run_log.write('-'*60+'\n')
                run_log.close()
                print '-'*60+'\n'
                print 'rmUnpaired finished at %s \n'%commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"')
                print '-'*60+'\n'
            else:
                run_log=open(out_dir+'/reseq_login.sh','a')
                run_log.write('-'*60+'\n')
                run_log.write('%s\n' % rmUnpaired_shell)
                run_log.write('-'*60+'\n')
                run_log.write('rmUnpaired Error!!!\n%s\n' %commands.getoutput('date'))
                run_log.write('-'*60+'\n')
                run_log.close()
                print '-'*60+'\n'
                print 'rmUnpaired Error!!!\n%s\n' %commands.getoutput('date')
                print '-'*60+'\n'
            out_clean = out_dir + '/01.DataCleaning'
            if getFile(out_clean, 'clean*fq') or getFile(out_clean, 'clean*fastq'):
                fqlist = getFile(out_clean, 'clean*fq') or getFile(out_clean, 'clean*fastq')
            print '-'*60
            print 'remove Unpaired fq file is: ',fqlist
        else:
            print '-'*60
            print 'remove Unpaired fq file is: ',fqlist
            pass
 
    if re.search('true', args.qcStat, re.I):
        qcStat_shell=qcStat(fqlist,out_dir)
        run_log=open(out_dir+'/reseq_login.sh','a')
        run_log.write('-'*60+'\n')
        run_log.write('QC statistics starting at %s\n' %commands.getoutput('date'))
        run_log.write('-'*60+'\n')
        run_log.close()
        print '-'*60+'\n'
        print 'QC statistics starting at %s\n' %commands.getoutput('date')
        print '-'*60+'\n'
        status = os.system(qcStat_shell)
        if status == 0:
            run_log=open(out_dir+'/reseq_login.sh','a')
            run_log.write('-'*60+'\n')
            run_log.write(qcStat_shell)
            run_log.write('-'*60+'\n')
            run_log.write('qcStat finished at %s \n'%commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"'))
            run_log.write('-'*60+'\n')
            run_log.close()
            print '-'*60+'\n'
            print 'qcStat finished at %s \n'%commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"')
            print '-'*60+'\n'
        else:
            run_log=open(out_dir+'/reseq_login.sh','a')
            run_log.write('-'*60+'\n')
            run_log.write( qcStat_shell)
            run_log.write('-'*60+'\n')
            run_log.write('qcStat Error!!!\n%s\n' %commands.getoutput('date'))
            run_log.write('-'*60+'\n')
            run_log.close()
            print '-'*60+'\n'
            print 'qcStat Error!!!\n%s\n' %commands.getoutput('date')
            print '-'*60+'\n'
 
    if re.search('true', args.map, re.I):
        ref = args.ref
        mis = args.mismatch
        gapo = args.gapopen
        gape = args.gapextention
        insize = args.insertsize
        species=args.species
        mapping_shell = mapping(fqlist,RGID,species,ref, mis, gapo, gape, insize, out_dir)
        run_log=open(out_dir+'/reseq_login.sh','a')
        run_log.write('-'*60+'\n')
        run_log.write('mapping starting at %s\n' %commands.getoutput('date'))
        run_log.write('-'*60+'\n')
        run_log.close()
        print '-'*60+'\n'
        print 'mapping starting at %s\n' %commands.getoutput('date')
        print '-'*60+'\n'
        status = os.system(mapping_shell)
        if status == 0:
            run_log=open(out_dir+'/reseq_login.sh','a')
            run_log.write('-'*60+'\n')
            run_log.write(mapping_shell)
            run_log.write('-'*60+'\n')
            run_log.write('mapping finished at %s \n'%commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"'))
            run_log.write('-'*60+'\n')
            run_log.close()
            print '-'*60+'\n'
            print 'mapping finished at %s \n'%commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"')
            print '-'*60+'\n'
        else:
            run_log=open(out_dir+'/reseq_login.sh','a')
            run_log.write('-'*60+'\n')
            run_log.write( mapping_shell)
            run_log.write('-'*60+'\n')
            run_log.write('mapping Error!!!\n%s\n' %commands.getoutput('date'))
            run_log.write('-'*60+'\n')
            run_log.close()
            print '-'*60+'\n'
            print 'mapping Error!!!\n%s\n' %commands.getoutput('date')
            print '-'*60+'\n'
        out_aln = out_dir + '/02.Alignment'
        if getFile(out_aln, 'uniq.rmdup.sort.bam'):
            bam_file = getFile(out_aln, 'uniq.rmdup.sort.bam')
        print '-'*60
        print 'clean bam file is: ',bam_file

#   if re.search('true', args.map, re.I):
#       if not args.ref:
#           print '-f option is needed!!!\n'
#           sys.exit(0)
#       ref = args.ref
#       callVcf_shell = callVcf(ref,out_dir)
#       run_log=open(out_dir+'/reseq_login.sh','a')
#       run_log.write('='*60+'\n')
#       run_log.write('calling Vcf starting at %s\n' %commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"'))
#       run_log.write('='*60+'\n')
#       run_log.close()
#       print '='*60+'\n'
#       print 'calling Vcf starting at %s\n' %commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"')
#       print '='*60+'\n'
#       status = os.system(callVcf_shell)
#       if status == 0:
#           run_log=open(out_dir+'/reseq_login.sh','a')
#           run_log.write('#%s\n' %  '\n#'.join(callVcf_shell.split('\n')[:-1]))
#           run_log.write('='*60+'\n')
#           run_log.write('calling Vcf finished at %s\n' %commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"'))
#           run_log.write('='*60+'\n')
#           run_log.close()
#           print '='*60+'\n'
#           print 'calling Vcf finished at %s\n' %commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"')
#           print '='*60+'\n'
#       else:
#           run_log=open(out_dir+'/reseq_login.sh','a')
#           run_log.write('%s\n' % callVcf_shell)
#           run_log.write('='*60+'\n')
#           run_log.write('calling Vcf Error!!!\n%s\n' %commands.getoutput('date'))
#           run_log.write('='*60+'\n')
#           run_log.close()
#           out_aln = out_dir + '/02.Alignment'
#           out_var = out_dir  + '/03.Variants'
#           vcf_file = getFile(out_var, '.vcf')
#           print '-'*60+'\n'
#           print 'Vcf file is: ',vcf_file
#           print '-'*60+'\n'

    if re.search('true', args.genome,re.I):
        genome_bed = args.genomebed
        window_size=args.windowsize
        bam=getFile(out_dir+'/02.Alignment','.uniq.rmdup.sort.bam')[0]
        seq_len=args.readlen
        if not genome_bed:
            print 'The genomebed file is needed for genome option!!!'
            sys.exit(0)     
        coverGenome_shell = coverGenome(bam, seq_len, genome_bed,window_size,out_dir)
        run_log=open(out_dir+'/reseq_login.sh','a')
        run_log.write('='*60+'\n')
        run_log.write('calculating coverage  starting at %s\n' %commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"'))
        run_log.write('='*60+'\n')
        run_log.close()
        print '='*60+'\n'
        print 'calculating coverage starting at %s\n' %commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"')
        print '='*60+'\n'
        status = os.system(coverGenome_shell)
        if status == 0:
            run_log=open(out_dir+'/reseq_login.sh','a')
            run_log.write('-'*60+'\n')
            run_log.write(coverGenome_shell)
            run_log.write('-'*60+'\n')
            run_log.write('calculating coverage finished at %s \n'%commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"'))
            run_log.write('-'*60+'\n')
            run_log.close()
            print '-'*60+'\n'
            print 'calculating coverage finished at %s \n'%commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"')
            print '-'*60+'\n'
        else:
            run_log=open(out_dir+'/reseq_login.sh','a')
            run_log.write('-'*60+'\n')
            run_log.write( coverGenome_shell)
            run_log.write('-'*60+'\n')
            run_log.write('calculating coverage Error!!!\n%s\n' %commands.getoutput('date'))
            run_log.write('-'*60+'\n')
            run_log.close()
            print '-'*60+'\n'
            print 'calculating coverage Error!!!\n%s\n' %commands.getoutput('date')
            print '-'*60+'\n'
        out_coverage = out_dir + '/04.Coverage'
        if getFile(out_coverage, '.cover.genome'):
            coverage_file = getFile(out_coverage, '.cover.genome')
        print '-'*60
        print 'coverage file is: ',coverage_file

#   if re.search('true',args.genome,re.I):
#       species= args.species
#       vcf=getFile(out_dir +'/03.Variants','.hc.pairHMM.vcf')[0]
#       anno_shell = annotate(vcf,annodb, species,out_dir)
#       run_log=open(out_dir+'/reseq_login.sh','a')
#       run_log.write('='*60+'\n')
#       run_log.write('annotation  starting at %s\n' %commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"'))
#       run_log.write('='*60+'\n')
#       run_log.close()
#       print '='*60+'\n'
#       print 'annotation starting at %s\n' %commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"')
#       print '='*60+'\n'
#       status = os.system(anno_shell)
#       if status == 0:
#           run_log=open(out_dir+'/reseq_login.sh','a')
#           run_log.write('-'*60+'\n')
#           run_log.write(anno_shell)
#           run_log.write('-'*60+'\n')
#           run_log.write('annotation finished at %s \n'%commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"'))
#           run_log.write('-'*60+'\n')
#           run_log.close()
#           print '-'*60+'\n'
#           print 'annotation finished at %s \n'%commands.getoutput('date "+%Y-%m-%d%:%H:%M:%S"')
#           print '-'*60+'\n'
#       else:
#           run_log=open(out_dir+'/reseq_login.sh','a')
#           run_log.write('-'*60+'\n')
#           run_log.write( anno_shell)
#           run_log.write('-'*60+'\n')
#           run_log.write('annotation Error!!!\n%s\n' %commands.getoutput('date'))
#           run_log.write('-'*60+'\n')
#           run_log.close()
#           print '-'*60+'\n'
#           print 'annotation Error!!!\n%s\n' %commands.getoutput('date')
#           print '-'*60+'\n'
#           out_anno = out_dir + '/05.Annotation'
#           var_function_file = getFile(out_anno, 'variant_function')
#           print '-'*60
#           print 'function annotation file is: ',var_function_file

if __name__=='__main__':
    if len(sys.argv)>1:
        main()
    else:
        print 'Type "-h" or "--help" for more help'
        sys.exit(0)
