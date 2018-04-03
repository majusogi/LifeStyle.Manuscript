## Non-parametric or  T-test

1. Crear subset of the data: 

Quito_case=subset(table, location=="Quito" & status=="case")
Rios_case=subset(table, location=="Rios" & status=="case")
Quito_control=subset(table, location=="Quito" & status=="control")
Rios_control=subset(table, location=="Rios" & status=="control")

# Generate BoxPlot
gb<- ggplot(table, aes(x=loc.stat,y=Shannon, fill=status))+geom_boxplot() + theme_classic() + xlab("") +
ylab("Shannon Index") + theme(text = element_text(size=13)) 
col = c("darkred","darkblue")
gb + scale_fill_manual(values = col) 

# To sort order based on medians: 
library(forcats)

ggplot(iris, aes(x = fct_reorder(Species, Sepal.Width, fun = median, .desc =TRUE), y = Sepal.Width)) + geom_boxplot()

ggplot(tabla, aes(x=fct_reorder(zone, shannon, fun = median, .desc =TRUE) ,y=shannon))+ 
geom_boxplot(fill="#4271AE") + theme_classic() + theme(text = element_text(size=9))



# Generate Geo.point with different variables: 
> ggplot(tabla, aes(x=remoteness ,y=shannon))+ geom_point(aes(colour = factor(zone))) + theme_classic() + 
theme(text = element_text(size=9))




## Primero test for equality of variances:
var.test(Shannon ~ location, table, alternative = "two.sided")

t.test(Quito_control$phylo.div, Rios_control$phylo.div, var.equal = T)



## Create Violin Plots: 

# Run this function to calculate the mean and dev.std:

> data_summary <- function(x) {
+   m <- mean(x)
+   ymin <- m-sd(x)
+   ymax <- m+sd(x)
+   return(c(y=m,ymin=ymin,ymax=ymax))
+ }

diver<-read.csv("diversity.nonpareil.csv", sep = ",", header = T)

ggplot(diver, aes(x=diver$location, y=diver$diversity, fill=diver$location)) +  geom_violin() +
scale_fill_brewer(palette="Blues") + theme_classic() + stat_summary(fun.data=data_summary) +  xlab("") +
ylab("") + theme(text = element_text(size=13))


2. ADoNis: PERMANOVA

#FROM a distance matrix: 
Necisto dos sets de datos: 
1: matrix
2: metadata

#option 1: is to add metadata to the matrix in a single file

location	query	MG59	MG60	MG61
Rios	MG59	0	0.108981	0.0647395
Rios	MG60	0.108981	0	0.0912695
Rios	MG61	0.0647395	0.0912695	0
Rios	MG62	0.0986788	0.0740981	0.0670813

curtis<-read.table('BrayCurtis.matrix.tsv', sep="\t", header = T, row.names=1)

#option2 upload two separate files

PermANOVA:
adonis_location = adonis(matrix_data ~ location, inputFile)
adonis_location = adonis(otus_dist ~ location, metadata)

adonis_location # take a look at results; there are no summary() or plot() methods included

library(vegan)

























