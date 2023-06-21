rm(list = ls()) 
options=(stringsAsFactors=F)
load(file = "logfc20210824.Rdata")
load(file = "diff.Rdata")


# BiocManager::install("clusterProfiler",force = TRUE)
library(clusterProfiler)

deg$TGCT <- NULL


#Use all genes of herbal ingredients as background genes
bg_gene <- logfc[[1]]


#Random generation of signatures from background genes
ran_signature <- vector("list",1000)
for (j in 1:1000) {
  gene <- vector("list",608)
  for (i in 1:608) {# randomly generated 608 signatures and data fusion
    gene[[i]] <- sample(bg_gene,size=250,replace=F)#Randomly generated 250 genes in the background genes to form a signature
  }
  message("608 list created!")
  gene <- do.call(cbind,gene)
    #colnames(ran_signature[[i]]) <- c(1:608)
  ran_signature[[j]] <- as.data.frame(gene)
  ran_signature[[j]]$a <- c(1:250)
  ran_signature[[j]] <- melt(ran_signature[[j]],id="a",na.rm = FALSE,variable.name="Geneset_name",value.name="Gene_symbol")
  ran_signature[[j]] <- ran_signature[[j]][,-1]
  message(j,"created!")
}
  
save(ran_signature,file = "ran_signature.Rdata")

