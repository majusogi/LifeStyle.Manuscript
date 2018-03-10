#!/bin/bash

00.HiSeq.data


1. To merge lanes of read 1 and 2: 
zcat MG67_S*_L00*_R2_001.fastq.gz > MG67.2.fastq

2. To count the No. of raw reads: 

for i in *.1.fastq; do l=$( cat $i | paste - - - - | wc -l | sed -e 's/ *//') ; echo $(basename $i .1.fastq) $l; done | tr ' ' '\t' > No.fastq.raw.Quito.txt



3. To remove human reads: 

Flag -q of 0 and 1 specify fasta and fastq input

dir="/nv/hp10/mjsg3/shared3/projects/EcoZUR/00.raw.MGs/Quito.controls/00.Raw.data"

for i in $(cat Quito.list.txt); do

./bmtagger.sh -b hs_ref_GRCh38.bitmask -x hs_ref_GRCh38.srprism -T tmp -q 1 -1 $dir/$i.1.fastq -2 $dir/$i.2.fastq -o $dir/$i.bmtagger -X ;

done



