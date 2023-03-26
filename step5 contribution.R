rm(list = ls()) 
options=(stringsAsFactors=F)
load(file = "NES.Rdata")
load(file = "immune NES.Rdata")
load(file = "cor.Rdata")

NES <- NES[match(colnames(cor),rownames(NES)),]


# #Calculating the correlation coefficient for each cancer
 cancer_cor <- vector("numeric",30)
 for (i in 1:30) {
   a <- cor[,i]
   b <- NES[i,]
   cancer_cor[[i]] <- cor(a,b,method = "pearson")
 }
 names(cancer_cor) <- colnames(cor)
# cancer_cor
 mean(cancer_cor)
 
#Calculating the correlation coefficient for each herbal ingredient
 herb_cor <- vector("numeric",496)
 for (i in 1:496) {
   a <- cor[i,]
   b <- NES[,i]
   herb_cor[[i]] <- cor(a,b,method = "pearson")
 }
 names(herb_cor) <- rownames(cor)
# herb_cor
 mean(herb_cor)
 

 


#herb bar chart
herb_cor <- as.data.frame(herb_cor)
rownames(herb_cor) <- rownames(cor)
library(ggplot2)
ggplot(herb_cor, aes(herb_cor)) +
  geom_histogram(bins = 12,color = "black",fill = "#00cd00")+
  labs(x = 'correlation coeffient',y = 'count')+
  geom_vline(xintercept = 0.3420564,linetype = "dashed",color="black",size = 1)+
  geom_vline(xintercept = 0.3527726,linetype = "dashed",color="black",size = 1)+
  annotate(geom = "text", fontface = "bold", color="black",
           x = 0.21,y =110,
           label = 'drug average', size=6)+
  scale_x_continuous(breaks=seq(-0.4, 0.9, 0.1))+
  theme_classic()+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1))+
  theme(axis.ticks.x=element_line(color="black",size=1,lineend = 1))+
  theme(axis.line.y=element_line(linetype=1,color="black",size=1))+
  theme(axis.ticks.y=element_line(color="black",size=1,lineend = 1))+
  theme(axis.title.x = element_text(size = 20, color = "black", face = "bold", vjust = 0.1, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 20, color = "black", face = "bold"))+#, vjust = 0.5, hjust = 0.5))
  theme(axis.text.x = element_text(size = 15, color = "black"))+
  theme(axis.text.y = element_text(size = 15, color = "black"))





#cancer bar chart
cancer_cor <- as.data.frame(cancer_cor)
rownames(cancer_cor) <- colnames(cor)
ggplot(cancer_cor, aes(cancer_cor)) +
geom_histogram(bins = 10,color = "black",fill = "#cecd00")+
  labs(x = 'correlation coeffient',y = 'count')+
  geom_vline(xintercept = 0.5299466,linetype = "dashed",color="black",size = 1)+
  #geom_vline(xintercept = 0.3527726,linetype = "dashed",color="black",size = 1)+
  annotate(geom = "text", fontface = "bold", color="black",
           x = 0.6, y=8,
           label = 'cancer average', size=6)+
  scale_x_continuous(breaks=seq(0, 0.9, 0.1))+
  theme_classic()+
  theme(axis.line.x=element_line(linetype=1,color="black",size=1))+
  theme(axis.ticks.x=element_line(color="black",size=1,lineend = 1))+
  theme(axis.line.y=element_line(linetype=1,color="black",size=1))+
  theme(axis.ticks.y=element_line(color="black",size=1,lineend = 1))+
  theme(axis.title.x = element_text(size = 20, color = "black", face = "bold", vjust = 0.1, hjust = 0.5))+
  theme(axis.title.y = element_text(size = 20, color = "black", face = "bold"))+#, vjust = 0.5, hjust = 0.5))
  theme(axis.text.x = element_text(size = 15, color = "black"))+
  theme(axis.text.y = element_text(size = 15, color = "black"))


#save(cancer_cor,herb_cor,file = "NES&corr.Rdata")
#save(cancer_cor,herb_cor,file = "NES&corr1.Rdata")








