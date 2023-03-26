rm(list = ls()) 
options=(stringsAsFactors=F)
#Loading cancer expression profiles
load(file = "expression group_list.Rdata")
load(file = "difference.Rdata")

#Loading herbal ingredients Expression profiles
load(file = "logfc20210824.Rdata")

#Reading immune signature
install.packages("openxlsx")
library(openxlsx)                                 
signature <- read.xlsx("1-s2.0-S1043661818316074-mmc1.xlsx")
colnames(signature) <- signature[1,]
signature <- signature[-1,]
#Collation format
signature <- as.data.frame(t(as.matrix(signature)))
colnames(signature) <- signature[1,]
signature <- signature[-1,]
signature$a <- c(1:250)

library(data.table)
signature <- melt(signature,id="a",na.rm = FALSE,variable.name="Geneset_name",value.name="Gene_symbol")
signature <- signature[,-1]

#Constructing a deg list of the results of the differential gene expression for each cancer
deg <- list(deg[[1]],deg[[2]],deg[[3]],deg[[4]],deg[[5]],deg[[6]],deg[[7]],deg[[8]],deg[[9]],
             deg[[10]],deg[[11]],deg[[12]],deg[[13]],deg[[14]],deg[[15]],deg[[16]],deg[[17]],
             deg[[18]],deg[[19]],deg[[20]],deg[[21]],deg[[22]],deg[[23]],
            deg[[25]],deg[[26]],deg[[27]],deg[[28]],deg_CHOL,deg_LAML,deg_READ)
y <- c("ACC","BLCA","LGG","BRCA","CESC","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC",
       "LUAD","LUSC","OV","PAAD","PCPG","PRAD","SARC","SKCM","STAD","THYM","THCA","UCS","UCEC","CHOL","LAML","READ")
names(deg) <- y


#Calculating GSEA results for each cancer
BiocManager::install("clusterProfiler")
library(clusterProfiler)
cancer_NES <- vector("list",30)
cancer_GSEA <- vector("list",30)
cancer_gsea <- vector("list",30)
for (i in 1:30) {
  print(i)
  logFC <- deg[[i]][,1]
  names(logFC) <- rownames(deg[[i]])
  logFC <-  sort(logFC,decreasing = T)
  cancer_gsea[[i]] <- GSEA(logFC,TERM2GENE = signature,pvalueCutoff=1,nPermSimple = 100000)
  cancer_GSEA[[i]] <- cancer_gsea[[i]]@result
  names(cancer_GSEA)[[i]]=x[i]
}



#Matching the row names of each cancer GSEA result so that the row names are in the same order for each cancer result  
for (i in 1:29) { 
  cancer_GSEA[[i+1]] <- cancer_GSEA[[i+1]][match(rownames(cancer_GSEA[[1]]),rownames(cancer_GSEA[[i+1]])),]
}
#Calculating the NES value for each cancer 
for (i in 1:30) {   
  cancer_NES[[i]] <- cancer_GSEA[[i]]$NES
  names(cancer_NES)[[i]]=y[i]
}

#Combining the NES values for each cancer
pancancer_NES <- do.call(cbind,cancer_NES)
rownames(pancancer_NES) <- rownames(cancer_GSEA[[1]])


#Calculating GSEA results for each herbal ingredinet
herb_nes <- vector("list", 496)
herb_gsea <- vector("list",496)
herb_GSEA <- vector("list",496)
for (i in 1:496) {
  print(i)
  herb_logfc <- logfc[[i+2]]
  names(herb_logfc) <- logfc[[1]]
  herb_logfc <-  sort(herb_logfc,decreasing = T)
  herb_gsea[[i]] <- GSEA(herb_logfc,TERM2GENE = signature,pvalueCutoff=1,nPermSimple = 1000000)
  herb_GSEA[[i]] <- herb_gsea[[i]]@result
  names(herb_GSEA)[[i]]=colnames(logfc)[i+2]
} 


for (i in 1:495) {
  #Matching the row names of each herbal ingredient GSEA result so that the row names are in the same order for each herbal ingredient result  
  herb_GSEA[[i+1]] <- herb_GSEA[[i+1]][match(rownames(herb_GSEA[[1]]),rownames(herb_GSEA[[i+1]])),]
}
  
#Calculating the NES value for each herbal ingredient 
for (i in 1:496) { 
  herb_nes[[i]] <- herb_GSEA[[i]]$NES
  names(herb_nes)[[i]]=colnames(logfc)[i+2]
}

#Combining the NES values of each herbal ingredinet
herb_NES <- do.call(cbind,herb_nes)
rownames(herb_NES) <- rownames(herb_GSEA[[1]])

#Extraction of signatures with FDR <= 0.05 in the GSEA result for cancer 
#Count the number of signatures significantly associated with cancer
b <- vector("list",30)
up <- vector("list",30)
for (i in 1:30) {
  b[[i]] <- cancer_GSEA[[i]][which(cancer_GSEA[[i]]$qvalues<=0.05),]
  up[[i]] <- sum(b[[i]]$NES>0)
  names(up)[[i]]=y[i]
}

down <- vector("list",30)
for (i in 1:30) {
  down[[i]] <- sum(b[[i]]$NES<0)
  names(down)[[i]]=y[i]
}

cancer_up <- do.call(rbind,up)
colnames(cancer_up) <- "up"
cancer_down <- do.call(rbind,down)
colnames(cancer_down) <- "down"
cancer_FDR <- merge(cancer_down,cancer_up,by="row.names")
rownames(cancer_FDR) <- cancer_FDR[,1]
cancer_FDR <- cancer_FDR[,-1]
write.csv(cancer_FDR,file = "cancer_FDR.csv")

#Extraction of signatures with FDR <= 0.05 in the GSEA result for herbal ingredients 
#Counting the number of signatures significantly associated with herbal ingredinets
a <- vector("list",496)
up <- vector("list",496)
for (i in 1:496) {#Counting the number of NES > 0
  a[[i]] <- herb_GSEA[[i]][which(herb_GSEA[[i]]$qvalues<=0.05),]
  up[[i]] <- sum(a[[i]]$NES>0)
  names(up)[[i]]=colnames(logfc)[i+2]
}

down <- vector("list",496)
for (i in 1:496) {#Counting the number of NES < 0
  down[[i]] <- sum(a[[i]]$NES<0)
  names(down)[[i]]=colnames(logfc)[i+2]
}

herb_up <- do.call(rbind,up)
colnames(herb_up) <- "up"
herb_down <- do.call(rbind,down)
colnames(herb_down) <- "down"
herb_FDR <- merge(herb_down,herb_up,by="row.names")
rownames(herb_FDR) <- herb_FDR[,1]
herb_FDR <- herb_FDR[,-1]
write.csv(herb_FDR,file = "herb_FDR.csv")


#save(cancer_GSEA,cancer_NES,herb_GSEA,herb_nes,signature,logfc,file = "GSEA.Rdata")
#save(pancancer_NES,herb_NES,file = "immune NES.Rdata")
#save(deg,file = "diff.Rdata")


#Extraction of signatures(FDR <= 0.05) from the GSEA results for cancers
immune_pathway_cancer <- vector("list",30)
for (i in 1:30) {
  immune_pathway_cancer[[i]] <- cancer_GSEA[[i]][which(cancer_GSEA[[i]]$qvalues<=0.05),]$ID
  names(immune_pathway_cancer)[[i]]=y[i]
}

#Counting the number of cancers significantly associated with each signature
number_immune_pathway_cancer <- vector("list",608)
for (i in 1:608) {
  a <- cancer_GSEA$ACC$ID[[i]]
  b <- vector("numeric",30)
  for (x in 1:30) {
    b[[x]] <- a%in%immune_pathway_cancer[[x]]
  }
  number_immune_pathway_cancer[[i]] <- sum(b)
  names(number_immune_pathway_cancer)[[i]] <- cancer_GSEA$ACC$ID[[i]]
}

number_immune_pathway_cancer <- as.data.frame(do.call(rbind,number_immune_pathway_cancer)) 
colnames(number_immune_pathway_cancer) <- "total"

library(openxlsx)                                 
signature <- read.xlsx("1-s2.0-S1043661818316074-mmc1.xlsx")
colnames(signature) <- signature[1,]
signature <- signature[-1,]
rownames(signature) <- c(1:608)

number_immune_pathway_cancer <- number_immune_pathway_cancer[match(signature$Geneset_name,rownames(number_immune_pathway_cancer)),,drop=FALSE]
write.csv(number_immune_pathway_cancer,file = "cancer_number_immune_pathway.csv")



#Extraction of signature(FDR <= 0.05) from the GSEA results for herbal ingredients
immune_pathway_herb <- vector("list",30)

for (i in 1:496) {
  immune_pathway_herb[[i]] <- herb_GSEA[[i]][which(herb_GSEA[[i]]$qvalues<=0.05),]$ID
  names(immune_pathway_herb)[[i]]=names(herb_GSEA)[[i]]
}

#Counting the number of herbal ingredients significantly associated with each signature
number_immune_pathway_herb <- vector("list",608)
for (i in 1:608) {
  a <- herb_GSEA$Myricetrin$ID[[i]]
  b <- vector("numeric",30)
  for (x in 1:30) {
    b[[x]] <- a%in%immune_pathway_herb[[x]]
  }
  number_immune_pathway_herb[[i]] <- sum(b)
  names(number_immune_pathway_herb)[[i]] <- herb_GSEA$Myricetrin$ID[[i]]
}

number_immune_pathway_herb <- as.data.frame(do.call(rbind,number_immune_pathway_herb)) 
colnames(number_immune_pathway_herb) <- "total"

number_immune_pathway_herb <- number_immune_pathway_herb[match(signature$Geneset_name,rownames(number_immune_pathway_herb)),,drop=FALSE]
write.csv(number_immune_pathway_herb,file = "herb_number_immune_pathway.csv")




