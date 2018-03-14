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

4. Change format of the name to MIGA style: 

For read 1:
for f in *_1.fastq; do mv "$f" "$(awk -F '.' '{print $1".1.fastq" }' <<<"$f")"; done

for read 2:
for f in *_2.fastq; do mv "$f" "$(awk -F '.' '{print $1".2.fastq" }' <<<"$f")"; done


***Count Number of trimmed reads: 
for i in *.1.fasta; do l=$(grep -c '>' $i); echo $i $l ; done | tr ' ' '\t' > trimmed.Quito.control.txt


4. Nonpareil: 

for i in *.1.fasta; do ~/shared3/bin/nonpareil -s $i -b $(basename $i .1.fasta) -d 0.7 -t 30 -R 40000 -T alignment; done



5. MASH DISTANCES

Fast genome and metagenome distance estimation using MinHash
module unload gcc/4.6.2
module load gcc/4.9.0
module load capnproto-c++/0.5.3
module load autoconf/2.69
module load zlib/1.2.8
module load boost/1.57.0
module load mash/1.0.2

sketching first

For sequences to be compared with mash, they must first be sketched, which creates vastly reduced representations of them. This will happen automatically if mash dist is given raw sequences. However, if multiple comparisons will be performed, it is more efficient to create sketches with mash sketch first and provide them to mash dist in place of the raw sequences. mash sketch creates a sketch file, which is a reduced representation of a sequence or set of sequences (based on min-hashes) that can be used for fast distance estimations. Input can be fasta or fastq files.

mash sketch -l list.txt -o ~/data3/from-pmicro1/norovirus/metagenomes/12.MASH -k 25 -s 10000

options



8. MetaPhlAn: Metagenomic Phylogenetic Analysis (MetaPhlAn version 1.7.7 )

module load python/2.7
module load bowtie2/2.1.0











