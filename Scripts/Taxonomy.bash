## Taxonomic annotation of MGs: 

## Metaphlan 

module load python/2.7
module load bowtie2/2.1.0


./metaphlan.py 4-4M.CoupledReads.fa --bowtie2db /nv/hp10/avpg3/nsegata-metaphlan/bowtie2db/mpa -t rel_ab --input_type multifasta --nreads 10101132 --tax_lev s --bowtie2out 4-4M.bt2out > 09.MetaPhlAn/profiled_samples/4-4M.txt

for i in $(ls *.bmtagger.out.fa); do metaphlan2.py $i --input_type fasta --nproc 16 > ~/data3/from-pmicro1/norovirus/metagenomes/09.MetaPhlAn/$(basename $i bmtagger.out.fa)_profile.txt; done

--tax_lev TAXONOMIC_LEVEL

--nproc 16 threads
./metaphlan2.py 4-4M.CoupledReads.fa --input_type fasta  --bowtie2out metagenome.bowtie2.bz2  --nproc 16 > profiled_metagenome.txt


###to classify at different levels using the envri.bowtie2.bz2 already computed:

./metaphlan2.py --input_type bowtie2out Cas27.bt2out --tax_lev f > pfile.family.cas27.txt

for Multiple samples: 

for f in */*.bowtie2.bz2; do /nv/hp10/mjsg3/data/tools/metaphlan2/metaphlan2.py $f --input_type bowtie2out --tax_lev g 
--nproc 16 > ${f%.bowtie2.bz2}_prof.genus.txt; done



2. merge metaphaln tables. The resulting table contains relative abundances with microbial clades as rows and samples as columns.

utils/merge_metaphlan_tables.py /09.MetaPhlAn/profiled_samples/*.txt > 09.MetaPhlAn/output/merged_abundance_table.txt

/nv/hp10/mjsg3/data/tools/metaphlan2/utils/merge_metaphlan_tables.py Rios.cases/profi*.txt Quito.cases/profi*.txt > cases_abundance_table.txt



3.Plot results as heatmap
module unload python/2.7
module load anaconda2


/nv/hp10/mjsg3/data/tools/metaphlan2/utils/metaphlan_hclust_heatmap.py --minv 0.1 --tax_lev 'o' -s lin --in 09.MetaPhlAn/merged_table.txt --out ~/09.MetaPhlAn/plot

./metaphlan_hclust_heatmap.py -x 0.01 -y 0.05 -c bbcry -d braycurtis --minv 0.1 --top 20 --tax_lev 'g' -s lin --in /09.MetaPhlAn/merged_table.txt --out ~/09.MetaPhlAn/plot_top20g

 /nv/hp10/mjsg3/data/tools/metaphlan2/utils/metaphlan_hclust_heatmap.py --minv 0.1 --top 30 --tax_lev 's' -s lin --in cases_abundance_table.txt --out cases_top30spe
 
 Options: 
 -d {euclidean,minkowski,cityblock,seuclidean)... default: "braycurtis" 
 -x 0.01:  Width of heatmap cells. 
 -y 0.05: Height of heatmap cells
 --minv 0.1 : Minimum value to display.
 

*** ConStrains:
Metagenomic strain-level population genomics

Before run it, we have to run first Metaphlan and get the output:  profiled_metagenome.txt

Load Modules:
module load gcc/4.9.0
module load mkl/11.2
module load bowtie2/2.1.0
module load samtools/1.0
module load openmpi/1.8
module load python/2.7

Run Constrains, but paste the metaphlan output in the same location with the fasta files:

python /nv/hp10/mjsg3/data/tools/constrains/ConStrains.py -m /nv/hp10/mjsg3/data/tools/metaphlan2/metaphlan2.py -c test.conf -o test.results -t 10  

--min_cov=FLOAT     Minimum coverage of a species in a sample to be
                    considered [default: 10, range: 5+].
 
## OJO: modify format of the fastq name!!!

for f in *.2.clipped.fastq; do mv "$f" "$(awk -F '.' '{print $1".2.fastq" }' <<<"$f")"; done


test.conf format: 
//
sample: sample_1
fq1: ./fq/sample_1.1.fq
fq2: ./fq/sample_1.2.fq
metaphlan: ./fq/sample_1.metaphlan.txt
//
sample: sample_2
fq1: ./fq/sample_2.1.fq
fq2: ./fq/sample_2.2.fq
metaphlan: ./fq/sample_2.metaphlan.txt


## Output DIR: Results

Intra_sp_rel_ab.profiles     # this is a tabular file with strain relative abundance within species. See header for details.
Overall_rel_ab.profiles      # this is a tabular file with strain relative abundance in overall samples. See header for details.
uniGcode                     # this is the directory for all ".uniGcode" files, which are genotypes for strains.


# Species   strain_ID   masked_samples  sample_1   sample_2
Escherichia_coli    str-1   NA  53.252835   37.212245
Escherichia_coli    str-2   NA  46.747165   62.787755
Salmonella_typhi    str-1   1   15.492194   41.212442
Salmonella_typhi    str-2   1   38.313441   21.483291

This means there are two species that passed the minimum coverage requirement for strain inference, and the relative abundance of 
each strain (E.coli has 2 strains and S.typhi has 3) are listed in sample_1 and sample_2 columns.



*****In the "uniGcode/" directory, you will find a few "*.uniGcode" files with names indicating the species. 
The format looks like the example shown below:
# *: not covered base; -: uncertain base
#pid    position    ref str-1   str-2
p0387   1   A   *   *   # <- insufficient mapped reads for inference
p0387   2   T   *   *
......




*** Functional Profiling with HUMAnN2:

1.1 Requirements    
    module load bowtie2/2.3.2 

First I downloaded two databases: UniRef50 and chocophlan
humann2_databases --download chocophlan full $INSTALL_LOCATION
humann2_databases --download uniref uniref50_diamond $INSTALL_LOCATION

humann2 --input examples/demo.fastq --output ensayo --diamon /nv/hp10/mjsg3/data/tools --metaphlan /nv/hp10/mjsg3/data/tools/metaphlan2



### MaxBin:

Bowtie2 = module load bowtie2/2.1.0 
FragGeneScan = 
Hmmer3 =  module load hmmer/3.1b1
IDBA-UD  =  module load idba/1.1.1 

perl /nv/hp10/mjsg3/data/tools/MaxBin-2.1.1/run_MaxBin.pl -contig all.Contigs.cassava.fna -abund myout.contig.tmp.reads.abund3 -out Maxbin.results -thread 16

-contig all.Contigs.cassava.fna -reads_list reads_list_file.txt -out myout -prob_threshold 0.5 -markerset 100 -thread 16






















