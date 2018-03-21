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




