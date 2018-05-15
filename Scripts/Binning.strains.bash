Binning using MaxBin: 


Bowtie2 = module load bowtie2/2.1.0 
FragGeneScan = 
Hmmer3 =  module load hmmer/3.1b1
IDBA-UD  =  module load idba/1.1.1 


OJO:
module unload perl5.10.1
module load perl/5.14.4

perl run_MaxBin.pl -contig <assembly file name>  -out <output file name> -reads <libraries fastA> -thread 8 -plotmarker <file.pdf>

perl /nv/hp10/mjsg3/data/tools/MaxBin-2.2.4/run_MaxBin.pl -contig /nv/hp10/mjsg3/shared3/projects/EcoZUR/MGs.lifestyle/Rios.controls/data/05.assembly/MG59.LargeContigs.fna 
-reads /nv/hp10/mjsg3/shared3/projects/EcoZUR/MGs.lifestyle/Rios.controls/data/04.trimmed_fasta/MG59.CoupledReads.fa -out MG59 -thread 20


CheckM in cluster:

module purge
module load python/2.7
module load prodigal/2.6.1                          
module load ocaml/4.01.0    
module load pplacer/1.1alpha15
#module load pgi/14.10
module load hmmer/3.1b1

-x : extension: fasta/fna
checkm tree -x fna $dir/05.taxonomy $dir/05.taxonomy/checkM -t 16;

checkm tree_qa $dir/05.taxonomy/checkM > $dir/05.taxonomy/checkM.taxonomy.txt;

Adding new scripts

