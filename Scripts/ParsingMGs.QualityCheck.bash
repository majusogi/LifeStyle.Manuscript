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

###sketching first
Sketch size corresponds to the number of (non-redundant) min-hashes that are kept.
For sequences to be compared with mash, they must first be sketched, which creates vastly reduced representations of them. 
This will happen automatically if mash dist is given raw sequences. However, if multiple comparisons will be performed, 
it is more efficient to create sketches with mash sketch first and provide them to mash dist in place of the raw sequences. 
mash sketch creates a sketch file, which is a reduced representation of a sequence or set of sequences (based on min-hashes) that can be used for fast distance estimations. 
Input can be fasta or fastq files: Coupled.Reads.fna

Folder: 03.MASH.dist
for i in /nv/hp10/mjsg3/shared3/projects/EcoZUR/MGs.lifestyle/Rios.controls/data/04.trimmed_fasta/*.CoupledReads.fa; 
do mash sketch -k 25 -s 10000 $i -o $(basename $i .CoupledReads.fa); done

options
-l List input. Each file contains a list of sequence files, one per line.
-k Kmer size. Hashes will be based on strings of this many nucleotides 
-s Sketch size. Each sketch will have at most this many non-redundant min-hashes (Default: 1000)


###Calculating distances:

mash dist <query1> <query2>

To automate this analysis, I used a bash script ‘dist.mash’. 

#!/bin/bash

CTG=$(ls *.msh)

for i in $CTG ; do
   for j in $CTG ; do
      [[ "$i" == "$j" ]] && break
      o="dist/$(basename $i .msh)-$(basename $j .msh)"
      [[ -s $o.dist ]] || mash dist  $i  $j   > $o.dist
   done
done


Given k-mer sets A and B, the MinHash algorithm provides an estimation of the Jaccard index:
The Jaccard index is a useful measure of global sequence similarity because it correlates with ANI.






6. MetaPhlAn: Metagenomic Phylogenetic Analysis (MetaPhlAn version 1.7.7 )

module load python/2.7
module load bowtie2/2.1.0

./metaphlan2.py 4-4M.CoupledReads.fa --input_type fasta  --bowtie2out metagenome.bowtie2.bz2  --nproc 16 > profiled_metagenome.txt

./metaphlan2.py --input_type fasta environmental.CoupledReads.fa --tax_lev f --bowtie2out envri.bowtie2.bz2 --nproc 16 > profiled_envi.MG.txt



to classify at different levels using the envri.bowtie2.bz2 already computed:

./metaphlan2.py --input_type bowtie2out Cas27.bt2out  --tax_lev f > pfile.family.cas27.txt



for Multiple samples: 

for f in */*.bowtie2.bz2; do /nv/hp10/mjsg3/data/tools/metaphlan2/metaphlan2.py $f --input_type bowtie2out --tax_lev g --nproc 16 > ${f%.bowtie2.bz2}_prof.genus.txt; done





2. merge metaphaln tables. The resulting table contains relative abundances with microbial clades as rows and samples as columns.

utils/merge_metaphlan_tables.py /09.MetaPhlAn/profiled_samples/*.txt > 09.MetaPhlAn/output/merged_abundance_table.txt

/nv/hp10/mjsg3/data/tools/metaphlan2/utils/merge_metaphlan_tables.py casP7/profiled_casP7.txt casP4/profiled_casP4.txt casP5/profiled_casP5.txt casP8/profiled_casP8.txt > merged_abundance_table.txt



Visualization with Krona: 

./metaphlan2krona.py -p ../profiled_MG.27.txt -k krona.27.txt

 /Users/mjsg3/Software/KronaTools-2.4/scripts/ImportText.pl  LC1.krona.txt -o krona.27.html












