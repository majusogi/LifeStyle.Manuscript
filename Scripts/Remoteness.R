

> ggplot(algo3, aes(x=algo3$average, y=algo3$Chao1, color=algo3$zone)) + geom_point() + theme_classic() + 
theme(text = element_text(size=13)) + ylab("Chao1") + xlab("") 
