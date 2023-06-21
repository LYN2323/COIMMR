rm(list = ls()) 
options=(stringsAsFactors=F)
#load(file = "logfc20210824.Rdata")
load(file = "diff.Rdata")
#loading 1000 randomly generated signatures
load(file = "ran_signature.Rdata")

library(clusterProfiler)
library(furrr)
plan(multisession, workers = availableCores()-48)

options(future.globals.maxSize= 3221222500)


#GSEA of 1000signature in cancer expression profiles
system.time(
  result1 <- furrr::future_map(1:30, function(i) {
    logFC <- deg[[i]][,1]
    names(logFC) <- rownames(deg[[i]])
    logFC <-  sort(logFC,decreasing = T)
    
    rv <- purrr::map(1:1000, function(x) {
      cancer_GSEA <- suppressMessages(GSEA(logFC,TERM2GENE = ran_signature[[x]],pvalueCutoff=1,nPermSimple = 100000)@result)
      z <- cancer_GSEA[,c(1,5)]
      z
    })
    rv
  }, .progress = TRUE)
)

#save(result1,file = "cancer_GSEA1000.Rdata")

#Harmonizing the row names of GSEA results for each cancer
aa <- result1
for (i in 1:30) {
  for (x in 1:1000) {
    aa[[i]][[x]] <- aa[[i]][[x]][match(rownames(aa[[1]][[1]]),rownames(aa[[i]][[x]])),]
  }
}



#Present the NES values for each cancer to form a data frame
cancer_NES <- vector("list",1000)
for (i in 1:1000) {
  bb <- data.frame(matrix(nrow = 608,ncol = 30))
  for (x in 1:30) {
   bb[[x]] <- aa[[x]][[i]]$NES 
  }
  rownames(bb) <- rownames(aa[[1]][[1]])
  cancer_NES[[i]] <- bb
  rm(bb)
}

for (i in 1:1000) {
  colnames(cancer_NES[[i]]) <- y
}

#save(cancer_NES,result1,file = "cancer_NES_1000.Rdata")











