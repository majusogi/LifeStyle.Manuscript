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
 













