#Extraction of the first 250 up- and down-regulated genes in cancer differential genes

rm(list = ls()) 
options=(stringsAsFactors=F)
load(file = "diff.Rdata")
#BiocManager::install("reshape2")
library(reshape2)

deg$TGCT <- NULL


#The top 250 up-regulated and down-regulated gene were selected as the cancer signatures
y <- c("ACC","BLCA","LGG","BRCA","CESC","COAD","DLBC","ESCA","GBM","HNSC","KICH","KIRC","KIRP","LIHC",
       "LUAD","LUSC","OV","PAAD","PCPG","PRAD","SARC","SKCM","STAD","THYM","THCA","UCS","UCEC","CHOL","LAML","READ")
up_signature <- vector("list",30)
down_signature <- vector("list",30)
for (i in 1:30) {
  #Ranking of differential genes from smallest to largest p-value
  deg[[i]] <- deg[[i]][order(deg[[i]]$P.Value),]
  #Among the highly expressed drugs, the 250 with the smallest p-value were selected
  up <- head(deg[[i]][which(deg[[i]]$change=="up"),],250)
  up_signature[[i]] <- rownames(up)
  names(up_signature)[[i]] <- y[i]
  #Of the low-expressing drugs, the 250 with the smallest p-value were selected
  down <- head(deg[[i]][which(deg[[i]]$change=="down"),],250)
  down_signature[[i]] <- rownames(down)
  names(down_signature)[[i]] <- y[i]
}

#Merge up and down cancer signature lists into data frame
up_signature <- do.call(cbind,up_signature)
up_signature <- as.data.frame(up_signature)
down_signature <- do.call(cbind,down_signature)
down_signature <- as.data.frame(down_signature)
#Constructing up-regulated cancer signatures
up_signature$a <- c(1:250)
up_signature <- melt(up_signature,id="a",na.rm = FALSE,variable.name="Geneset_name",value.name="Gene_symbol")
up_signature <- up_signature[,-1]

#Constructing down-regulated cancer signatures
down_signature$a <- c(1:250)
down_signature <- melt(down_signature,id="a",na.rm = FALSE,variable.name="Geneset_name",value.name="Gene_symbol")
down_signature <- down_signature[,-1]
