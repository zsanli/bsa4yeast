#!/usr/bin/env python
"""
#
#------------------------------------------------------------------------------
#                         University of Minnesota
#         Copyright 2012, Regents of the University of Minnesota
#------------------------------------------------------------------------------
# Author:
#
#  James E Johnson
#  Jesse Erdmann
#
#------------------------------------------------------------------------------
"""

"""
Generate a VCF file from a samtools pileup
filtering on read coverage, base call quality, and allele frequency

Pileup Format
http://samtools.sourceforge.net/pileup.shtml
Columns:
chromosome, 1-based coordinate, reference base, the number of reads covering the site, read bases, base qualities

VCF Format
http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
The header line names the 8 fixed, mandatory columns. These columns are as follows:
CHROM POS ID REF ALT QUAL FILTER INFO
"""

import sys,re,os.path
import optparse
from optparse import OptionParser

vcf_header =  """\
##fileformat=VCFv4.0
##source=pileup_to_vcf.pyV1.2
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">
##INFO=<ID=SAF,Number=.,Type=Float,Description=\"Specific Allele Frequency\">
##FILTER=<ID=DP,Description=\"Minimum depth of %s\">
##FILTER=<ID=SAF,Description=\"Allele frequency of at least %s with base quality minimum %d\">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\
"""

def __main__():
  #Parse Command Line
  parser = optparse.OptionParser()
  # files
  parser.add_option( '-i', '--input', dest='input', help='The input pileup file (else read from stdin)' )
  parser.add_option( '-o', '--output', dest='output', help='The output vcf file (else write to stdout)' )
  # filters
  parser.add_option( '-c', '--min_coverage', type='int', default='1', dest='min_coverage', help='The minimum read coverage depth to report a variant [default 1]' )
  parser.add_option( '-r', '--report_depth', type='choice', choices=['source','ref','qual','all'], default='ref', dest='report_depth', help='Coverage depth as: source - as reported in the ipleup source, ref - reads rcovering base, qual - reads with base call above min_base_quality, all - reads spanning base (includes indels)' )
  parser.add_option( '-b', '--min_base_quality', type='int', default='0', dest='min_base_quality', help='The minimum base quality for a base call to be counted' )
  parser.add_option( '-f', '--min_allele_freq', type='float', default='.5', dest='min_allele_freq', help='The minimum frequency of an allele for it to be reported (default .5)' )
  parser.add_option( '-m', '--allow_multiples', action="store_true", dest='allow_multiples', default=False, help='Allow multiple alleles to be reported' )
  parser.add_option( '-s', '--snps_only', action="store_true", dest='snps_only', default=False, help='Only report SNPs, not indels' )
  # ID to use 
  parser.add_option( '-I', '--id', dest='id', default=None, help='The value for the VCF ID field' )
  # select columns
  parser.add_option( '-C', '--chrom_col', type='int', default='1', dest='chrom_col', help='The ordinal position (starting with 1) of the chromosome column' )
  parser.add_option( '-P', '--pos_col', type='int', default='2', dest='pos_col', help='The ordinal position (starting with 1) of the position column' )
  parser.add_option( '-R', '--ref_col', type='int', default='3', dest='ref_col', help='The ordinal position (starting with 1) of the reference base  column' )
  parser.add_option( '-D', '--coverage_col', type='int', default='4', dest='coverage_col', help='The ordinal position (starting with 1) of the read coverage depth column' )
  parser.add_option( '-B', '--base_call_col', type='int', default='5', dest='base_call_col', help='The ordinal position (starting with 1) of the base call column' )
  parser.add_option( '-Q', '--base_qual_col', type='int', default='6', dest='base_qual_col', help='The ordinal position (starting with 1) of the base quality column' )
  # debug
  parser.add_option( '-d', '--debug', action="store_true", dest='debug', default=False, help='Print Debug messages' )
  # 
  (options, args) = parser.parse_args()

  debug = options.debug if options.debug else False

  #Coverage must have at least a total depth of min_coverage to be considered
  min_coverage = options.min_coverage
  #Count depth as: ref - reads covercovering base, all - reads spanning base, qual - reads with base call above min_base_quality
  dp_as = options.report_depth
  #Base qualities must be at least the min_base_quality value to be counted
  min_base_quality = options.min_base_quality
  #Coverage must be more than this pct in support of the variant, e.g. variant support / total coverage > min_allele_freq
  min_allele_freq = options.min_allele_freq
  #Allow multiple variants to be reported for a position, min_allele_freq must be less than 0.5 for this to have effect.
  allow_multiples = options.allow_multiples
  # Whether to report indels 
  report_indels = not options.snps_only

  # Check for 1-based column options, and subtract 1 for python array index
  #Which column (count starting with 1) contains the chromosome (Default assumes 6-col pileup or 12-col filtered pileup)
  chrom_col = options.chrom_col - 1
  #Which column (count starting with 1) contains the position (Default assumes 6-col pileup or 12-col filtered pileup)
  pos_col = options.pos_col - 1
  #Which column (count starting with 1) contains the reference base (Default assumes 6-col pileup or 12-col filtered pileup)
  ref_col = options.ref_col - 1
  #Which column (count starting with 1) contains the total coverage for the position (Default assumes 6-col pileup or 12-col filtered pileup)
  coverage_col = options.coverage_col - 1
  #Which column (count starting with 1) contains the base calls for the position (Default assumes 6-col pileup or 12-col filtered pileup)
  base_call_col = options.base_call_col - 1
  #Which column (count starting with 1) contains the base call qualities for the position (Default assumes 6-col pileup or 12-col filtered pileup)
  base_qual_col = options.base_qual_col - 1

  # pileup input 
  if options.input != None:
    try:
      inputPath = os.path.abspath(options.input)
      inputFile = open(inputPath, 'r')
    except Exception, e:
      print >> sys.stderr, "failed: %s" % e
      exit(2)
  else:
    inputFile = sys.stdin
  # vcf output 
  if options.output != None:
    try:
      outputPath = os.path.abspath(options.output)
      outputFile = open(outputPath, 'w')
    except Exception, e:
      print >> sys.stderr, "failed: %s" % e
      exit(3)
  else:
    outputFile = sys.stdout

  vcf_id = options.id if options.id else "."

  indel_len_pattern = '([1-9][0-9]*)'
  ref_skip_pattern = '[<>]'

  SAF_header = "SAF";

  print >> outputFile, vcf_header % (min_coverage,min_allele_freq,min_base_quality)

  try:  
    for linenum,line in enumerate(inputFile):
      if debug: print >> sys.stderr, "%5d:\t%s" % (linenum,line)
      fields = line.split('\t')
      sdp = int(fields[coverage_col]) # count of all reads spanning ref base from source coverage column
      # Quick check for coverage
      if dp_as != 'all' and sdp < min_coverage :  
        continue
      bases = fields[base_call_col]
      m = re.findall(ref_skip_pattern,bases) 
      # Quick check for coverage, pileup coverage - reference skips
      rdp = sdp - (len(m) if m else 0)
      if dp_as == 'ref'  or dp_as == 'qual':
        if rdp < min_coverage:
          continue
      quals = fields[base_qual_col]
      chrom = fields[chrom_col]
      pos = int(fields[pos_col])
      ref = fields[ref_col]
      deletions = dict()
      max_deletion = 0
      insertions = dict()
      snps = {'A':0,'C':0,'G':0,'T':0,'N':0}
      bi = 0	# bases index
      qi = 0	# quals index
      adp = 0	# read coverage depth
      qdp = 0	# read coverage depth with base call quality passing minimum
      variants = dict() # variant -> count
      mc = 0 # match count
      while bi < len(bases):
        base = bases[bi] 
        # if debug: print >> sys.stderr, "%3d\t%s\t%3d\t%s %3d" % (bi,base,qi,quals[qi],ord(quals[qi]) - 33)
        bi += 1
        if base in '<>' :	#reference skip (between paired reads)
          if debug: print >> sys.stderr, "%3d\t%s\t%3d\t%s %3d" % (bi-1,base,qi,quals[qi],ord(quals[qi]) - 33)
          qi += 1
        elif base in '.,' : # match reference on forward/reverse strand
          if debug: print >> sys.stderr, "%3d\t%s\t%3d\t%s %3d" % (bi-1,base,qi,quals[qi],ord(quals[qi]) - 33)
          mc += 1
          adp += 1
          qual = ord(quals[qi]) -33
          qi += 1
          # count match 
          if qual >= min_base_quality:
            qdp += 1
        elif base in 'ACGTNacgtn' : # SNP on forward/reverse strand
          if debug: print >> sys.stderr, "%3d\t%s\t%3d\t%s %3d" % (bi-1,base,qi,quals[qi],ord(quals[qi]) - 33)
          adp += 1
          qual = ord(quals[qi]) -33
          qi += 1
          # record snp variant
          if qual >= min_base_quality:
            qdp += 1
            snps[base.upper()] += 1
        elif base in '+-' : # indel 
          adp += 1
          m = re.match(indel_len_pattern,bases[bi:]) 
          indel_len_str = m.group(0) # get the indel len string
          bi += len(indel_len_str) # increment the index by the len of the len string
          indel_len = int(indel_len_str)
          indel = bases[bi:bi+indel_len].upper() # get the indel bases
          if debug: print >> sys.stderr, "%3d\t%s%s" % (bi-2,base,indel)
          bi += indel_len # increment the index by the len of the indel 
          if base == '+':
            # record insertion variant
            if indel in insertions: 
              insertions[indel] += 1
            else:
              insertions[indel] = 1
          elif base == '-': 
            # record deletion variant
            max_deletion = max(indel_len,max_deletion)
            if indel in deletions: 
              deletions[indel] += 1
            else:
              deletions[indel] = 1
        elif base == '*' : # a deletion from a prior base
          if debug: print >> sys.stderr, "%3d\t%s\t%3d\t%s %3d" % (bi-1,base,qi,quals[qi],ord(quals[qi]) - 33)
          qi += 1
        elif base == '^' : #start of read (followed by read quality char)
          if debug: print >> sys.stderr, "%3d\t%s%s" % (bi-1,base,bases[bi])
          read_mapping_qual = ord(bases[bi]) - 33 
          bi += 1
        elif base == '$' : #end of read
          if debug: print >> sys.stderr, "%3d\t%s" % (bi-1,base)
          pass # no associated quality in either bases or quals
      dp = rdp if dp_as == 'ref' else qdp if dp_as == 'qual' else sdp if dp_as == 'source' else adp
      if debug: print >> sys.stderr, "match count: %d" % mc
      if debug: print >> sys.stderr, "snps: %s\tins: %s\tdel: %s" % (snps,insertions,deletions)
      if dp < max(1,min_coverage) :
        continue
      # output any variants that meet the depth coverage requirement
      vcf_pos = pos
      vcf_ref = ref
      suffix = ''
      alts = []
      safs = []
      # If we have any deletions to report, need to modify the ref string
      if debug: print >> sys.stderr, "max deletions: %s" % deletions
      if report_indels:
        for k,v in deletions.items():
          saf = v/float(dp)
          if saf >= min_allele_freq:
            if len(k) > len(suffix):
              suffix = k
        vcf_ref = (ref + suffix).upper()
      if debug: print >> sys.stderr, "snps: %s" % snps
      for k,v in snps.items():
        if debug: print >> sys.stderr, "snp: %s %d" % (k,v)
        saf = v/float(dp)
        if saf >= min_allele_freq:
          alts.append(k + suffix)
          safs.append(saf)
      if debug: print >> sys.stderr, "insertions: %s" % insertions
      if report_indels:
        for k,v in insertions.items():
          saf = v/float(dp)
          if saf >= min_allele_freq:
            alts.append(ref + k + suffix)
            safs.append(saf)
        if debug: print >> sys.stderr, "deletions: %s" % deletions
        for k,v in deletions.items():
          saf = v/float(dp)
          if saf >= min_allele_freq:
            alts.append(vcf_ref[:len(vcf_ref) - len(k)])   # TODO alt will be a substring of vcf_ref,  test this
            safs.append(saf)
      if len(alts) > 0:
        vcf_qual = "." 
        vcf_filter = "PASS"
        # if not allow_multiples, report only the most freq alt
        if not allow_multiples and len(alts) > 1:
          max_saf = max(safs)
          for i,v in enumerate(safs):
            if v == max_saf:
              vcf_alt = alts[i]
          info = "SAF=%f;DP=%d" % (max_saf,dp)
        # if allow_multiples
        else:
          vcf_alt = ','.join(alts)
          isafs = []
          for saf in safs:
            isafs.append("SAF=%f" % saf)
          info = "%s;DP=%d" % (','.join(isafs),dp)
        print >> outputFile,"%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s" % (chrom,vcf_pos,vcf_id,vcf_ref,vcf_alt,vcf_qual,vcf_filter,info)
  except Exception, e:
    print >> sys.stderr, "failed: %s" % e
    exit(1)

if __name__ == "__main__": __main__()

