bsa4yeast
===========================
This is the pipeline code in bsa4yeast.lcsb.uni.lu.  
****
|Author|Zhi Zhang|
|---|---
|E-mail|zhi.zhang@uni.lu

****
## Table of Contents
* [run with map file](#map)
* [run with bam file](#bam)
* [run with fastq file](#fastq)
* [scripts](./scripts)

### map
-----------
run with map file
```Bash
python2.7 bsacalc.py -L myresult.LMap -H myresult.HMap -w 33750 -k tricube -o myresult.G
python2.7 bsapeak.py myresult.G . myresult
python2.7 bsadraw.py -g myresult.G -c sc_chr.length --ylim 0 30  -o myresult.pdf --noraw --threshold 4 -n 0
paste myresult.HMap | awk '{ s=$3;p=$4; for(i=2; i<=NF/4; i+=1) {s+=$i*4-1;p+=$i*4}; print $1,$2, s,p }'|awk '{if ($3+$4==0)print $1,$2,$3/($3+$4+1);else print $1,$2,$3/($3+$4);}' > myresult_Hmap_freq.txt
paste myresult.LMap | awk '{ s=$3;p=$4; for(i=2; i<=NF/4; i+=1) {s+=$i*4-1;p+=$i*4}; print $1,$2, s,p }'|awk '{if ($3+$4==0)print $1,$2,$3/($3+$4+1);else print $1,$2,$3/($3+$4);}' > myresult_Lmap_freq.txt
paste myresult_Hmap_freq.txt  myresult_Lmap_freq.txt  |awk '{print $1,$2,$3-$6}' > myresult_freq.txt
Rscript --vanilla freq_diff.r myresult_freq.txt myresult_freq.pdf
```

### bam
-----------
run with bam file
```Bash
samtools mpileup -P ILLUMINA  -f  sc64.fa sake_strain_P1.bam white_tecc_strain_P2.bam Bulk_H.bam Bulk_L.bam > myresult.mpileup
python2.7 count_allele_for4grp_bsa.py myresult.mpileup 4 > myresult.count
python2.7 find_hom_nX.py myresult.count 10 > myresult.10xhom
python2.7 find_marker.py  myresult.10xhom > myresult.marker
python2.7 sort_chr.py myresult.marker > myresult.sorted
python2.7 find_4cell.py myresult.sorted 4 2 > myresult.4cell
grep "MPass" myresult.4cell > myresult.MPass
cut -f1,2,40,41 myresult.MPass > myresult.HMap
cut -f1,2,42,43 myresult.MPass > myresult.LMap
python2.7 bsacalc.py -L myresult.LMap -H myresult.HMap -w 33750 -k tricube -o myresult.G
python2.7 bsapeak.py myresult.G . myresult
python2.7 bsadraw.py -g myresult.G -c sc_chr.length --ylim 0 30  -o myresult.pdf --noraw --threshold 4 -n 0
paste myresult.HMap | awk '{ s=$3;p=$4; for(i=2; i<=NF/4; i+=1) {s+=$i*4-1;p+=$i*4}; print $1,$2, s,p }'|awk '{if ($3+$4==0)print $1,$2,$3/($3+$4+1);else print $1,$2,$3/($3+$4);}' > myresult_Hmap_freq.txt
paste myresult.LMap | awk '{ s=$3;p=$4; for(i=2; i<=NF/4; i+=1) {s+=$i*4-1;p+=$i*4}; print $1,$2, s,p }'|awk '{if ($3+$4==0)print $1,$2,$3/($3+$4+1);else print $1,$2,$3/($3+$4);}' > myresult_Lmap_freq.txt
paste myresult_Hmap_freq.txt  myresult_Lmap_freq.txt  |awk '{print $1,$2,$3-$6}' > myresult_freq.txt
Rscript --vanilla freq_diff.r myresult_freq.txt myresult_freq.pdf
samtools mpileup -E -P ILLUMINA  -f sc64.fa sake_strain_P1.bam > sake_strain_P1.mpileup
samtools mpileup -E -P ILLUMINA  -f sc64.fa white_tecc_strain_P2.bam > white_tecc_strain_P2.mpileup
python2.7 pileup_to_vcf.py -i sake_strain_P1.mpileup -o sake_strain_P1.vcf
python2.7 pileup_to_vcf.py -i white_tecc_strain_P2.mpileup -o white_tecc_strain_P2.vcf
perl convert2annovar.pl -format vcf4 sake_strain_P1.vcf > sake_strain_P1.avi
perl convert2annovar.pl -format vcf4 white_tecc_strain_P2.vcf > white_tecc_strain_P2.avi
perl annotate_variation.pl -buildver sacCer3 sake_strain_P1.avi sacCer3
perl annotate_variation.pl -buildver sacCer3 white_tecc_strain_P2.avi sacCer3
python2.7 bsa_ann3.py -i <P1_exonic_variant_function_file> -j <P2_exonic_variant_function_file> -se scoreNonSynSNPS.pl -si myresult.SNAP -so myresult.SNAP.out -o myresult_P1_P2_combined_exonic_variant_function.txt
python2.7 create_peak_ann.py -i myresult_P1_P2_combined_exonic_variant_function.txt -j myresult_peak.txt  -o myresult_peak_ann.txt
python2.7 create_region_ann.py -i myresult_P1_P2_combined_exonic_variant_function.txt -j myresult_region.txt  -o myresult_region_ann.txt
```

### fastq
-----------
run with fastq file
```Bash
bwa mem -R "@RG\tID:USER\tLB:USER\tPL:illumina\tSM:USER\tPI:USER" sc64.fa sake_strain_P1_1.fastq sake_strain_P1_2.fastq > sake_strain_P1.sam
samtools view -Sb sake_strain_P1.sam > sake_strain_P1.bam
samtools sort sake_strain_P1.bam sake_strain_P1.bam.sort 

bwa mem -R "@RG\tID:USER\tLB:USER\tPL:illumina\tSM:USER\tPI:USER" sc64.fa white_tecc_strain_P2_1.fastq white_tecc_strain_P2_2.fastq > white_tecc_strain_P2.sam
samtools view -Sb white_tecc_strain_P2.sam > white_tecc_strain_P2.bam
samtools sort white_tecc_strain_P2.bam white_tecc_strain_P2.bam.sort 

bwa mem -R "@RG\tID:USER\tLB:USER\tPL:illumina\tSM:USER\tPI:USER" sc64.fa Bulk_H.fastq > Bulk_H.sam
samtools view -Sb Bulk_H.sam > Bulk_H.bam
samtools sort Bulk_H.bam Bulk_H.bam.sort 

bwa mem -R "@RG\tID:USER\tLB:USER\tPL:illumina\tSM:USER\tPI:USER" sc64.fa Bulk_L.fastq > Bulk_L.sam
samtools view -Sb Bulk_L.sam > Bulk_L.bam
samtools sort Bulk_L.bam Bulk_L.bam.sort 

samtools mpileup -P ILLUMINA  -f  sc64.fa sake_strain_P1.bam.sort white_tecc_strain_P2.bam.sort Bulk_H.bam.sort Bulk_L.bam.sort > myresult.mpileup
python2.7 count_allele_for4grp_bsa.py myresult.mpileup 4 > myresult.count
python2.7 find_hom_nX.py myresult.count 10 > myresult.10xhom
python2.7 find_marker.py  myresult.10xhom > myresult.marker
python2.7 sort_chr.py myresult.marker > myresult.sorted
python2.7 find_4cell.py myresult.sorted 4 2 > myresult.4cell
grep "MPass" myresult.4cell > myresult.MPass
cut -f1,2,40,41 myresult.MPass > myresult.HMap
cut -f1,2,42,43 myresult.MPass > myresult.LMap
python2.7 bsacalc.py -L myresult.LMap -H myresult.HMap -w 33750 -k tricube -o myresult.G
python2.7 bsapeak.py myresult.G . myresult
python2.7 bsadraw.py -g myresult.G -c sc_chr.length --ylim 0 30  -o myresult.pdf --noraw --threshold 4 -n 0
paste myresult.HMap | awk '{ s=$3;p=$4; for(i=2; i<=NF/4; i+=1) {s+=$i*4-1;p+=$i*4}; print $1,$2, s,p }'|awk '{if ($3+$4==0)print $1,$2,$3/($3+$4+1);else print $1,$2,$3/($3+$4);}' > myresult_Hmap_freq.txt
paste myresult.LMap | awk '{ s=$3;p=$4; for(i=2; i<=NF/4; i+=1) {s+=$i*4-1;p+=$i*4}; print $1,$2, s,p }'|awk '{if ($3+$4==0)print $1,$2,$3/($3+$4+1);else print $1,$2,$3/($3+$4);}' > myresult_Lmap_freq.txt
paste myresult_Hmap_freq.txt  myresult_Lmap_freq.txt  |awk '{print $1,$2,$3-$6}' > myresult_freq.txt
Rscript --vanilla freq_diff.r myresult_freq.txt myresult_freq.pdf
samtools mpileup -E -P ILLUMINA  -f sc64.fa sake_strain_P1.bam > sake_strain_P1.mpileup
samtools mpileup -E -P ILLUMINA  -f sc64.fa white_tecc_strain_P2.bam > white_tecc_strain_P2.mpileup
python2.7 pileup_to_vcf.py -i sake_strain_P1.mpileup -o sake_strain_P1.vcf
python2.7 pileup_to_vcf.py -i white_tecc_strain_P2.mpileup -o white_tecc_strain_P2.vcf
perl convert2annovar.pl -format vcf4 sake_strain_P1.vcf > sake_strain_P1.avi
perl convert2annovar.pl -format vcf4 white_tecc_strain_P2.vcf > white_tecc_strain_P2.avi
perl annotate_variation.pl -buildver sacCer3 sake_strain_P1.avi sacCer3
perl annotate_variation.pl -buildver sacCer3 white_tecc_strain_P2.avi sacCer3
python2.7 bsa_ann3.py -i <P1_exonic_variant_function_file> -j <P2_exonic_variant_function_file> -se scoreNonSynSNPS.pl -si myresult.SNAP -so myresult.SNAP.out -o myresult_P1_P2_combined_exonic_variant_function.txt
python2.7 create_peak_ann.py -i myresult_P1_P2_combined_exonic_variant_function.txt -j myresult_peak.txt  -o myresult_peak_ann.txt
python2.7 create_region_ann.py -i myresult_P1_P2_combined_exonic_variant_function.txt -j myresult_region.txt  -o myresult_region_ann.txt
```

