### Stats in 16S rRNA amplicons: 

Adonis:

map=read.table("maping.Adonis.mod.txt", header=TRUE, sep="\t")

adonis(braycurtis.d~Location + New.Age + race + parto + education	+ sanitation_casa_cat	+ water_casa_cat + 
treat_casa + Gender, data=map, permutations=999)

Multiplication of Factors:

 adonis(braycurtis.d~Location* New.Age * race * parto * education * sanitation_casa_cat	* water_casa_cat *
 treat_casa * Gender, data=map, permutations=999)
 
