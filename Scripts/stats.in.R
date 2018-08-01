## Non-parametric or  T-test

1. Crear subset of the data: 

Quito_case=subset(table, location=="Quito" & status=="case")
Rios_case=subset(table, location=="Rios" & status=="case")
Quito_control=subset(table, location=="Quito" & status=="control")
Rios_control=subset(table, location=="Rios" & status=="control")

# Generate BoxPlot
col = c("darkred","darkblue")

## To plot only Controls
col = c("darkorange","forestgreen")
ggplot(B, aes(x=location,y=Shannon, fill=location))+geom_boxplot()+ theme_classic() + xlab("") +ylab("Shannon Index") + 
theme(text = element_text(size=13)) + scale_fill_manual(values = col) + xlim(-5000, 5000)

# To sort order based on medians: 
library(forcats)

ggplot(tabla, aes(x=fct_reorder(zone, shannon, .fun = median, .desc =TRUE) ,y=shannon))+ 
geom_boxplot(fill="#4271AE") + theme_classic() + theme(text = element_text(size=9))


# Generate Geo.point with different variables: 
> ggplot(tabla, aes(x=remoteness ,y=shannon))+ geom_point(aes(colour = factor(zone))) + theme_classic() + 
theme(text = element_text(size=9))


# To Order according with the Remoteness value: 
Bdata$factor<-factor(Bdata$zone, levels=c("Colon Eloy", "Maldonado", "Timbire", "Atahualpa", "playa de oro", "Santo Domingo",
                                           "Selva Alegre", "Zapallo Grande", "San Francisco"))

boxplot(Bdata$shannon ~ Bdata$factor)

## Set up colors/ Colored By Zone
library(RColorBrewer)
NumberOfLevels<-length(levels(Bdata$factor))
mycolors<-brewer.pal(n=NumberOfLevels, name="Set1")
boxplot(Bdata$shannon ~ Bdata$factor, col=mycolors)

To contorl font size:
par(cex.lab=1.5) # is for y-axis

par(cex.axis=1.5) # is for x-axis


## Statistics: 
kruskal.test(shannon ~ zone, data = Bdata)



## Primero test for equality of variances:
var.test(Shannon ~ location, table, alternative = "two.sided")

#####If variances are not equal: 
A modification of the t-test known as Welchs test is said to correct for this problem by 
estimating the variances, and adjusting the degrees of freedom to use in the tes


t.test(Quito_control$phylo.div, Rios_control$phylo.div)

t.test(x, y, alternative = "two.sided", var.equal = FALSE)


###If variances are equal: 
Agregar var.equal = T
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

























