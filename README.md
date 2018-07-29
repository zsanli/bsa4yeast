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

```
### bam
-----------
run with bam file
```Bash
samtools mpileup -P ILLUMINA  -f  sc64.fa sake_strain_P1.bam white_tecc_strain_P2.bam Bulk_H.bam Bulk_L.bam > myresult.mpileup
```

### fastq
-----------
run with fastq file
```Bash

```

