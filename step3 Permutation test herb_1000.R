rm(list = ls())
options(stringsAsFactors=F)
load(file = "logfc20210824.Rdata")
load(file = "ran_signature.Rdata")

library(clusterProfiler)
library(furrr)
plan(multisession, workers = availableCores()-48)

options(future.globals.maxSize= 3221222500)

#GSEA of 1000signature in expression profiles of herbal ingredients 
system.time(
  result <- furrr::future_map(1:496, function(i) {
    logFC <- logfc[[i+2]]
    names(logFC) <- logfc[[1]]
    logFC <-  sort(logFC,decreasing = T)

    rv <- purrr::map(1:1000, function(x) {
      herb_GSEA <- suppressMessages(GSEA(logFC,TERM2GENE = ran_signature[[x]],pvalueCutoff=1,nPermSimple = 10000)@result)
      z <- herb_GSEA[,c(1,5)]
      z
    })
    rv
  }, .progress = TRUE)
)

#save(result,file = "herb_GSEA1000.Rdata")


#Unify the row names of all elements of the herb_NES list
for (i in 1:1000) {
    herb_NES[[i]] <- herb_NES[[i]][match(rownames(herb_NES[[1]]),rownames(herb_NES[[i]])),]
}

for (i in 1:1000) {
  colnames(herb_NES[[i]]) <- colnames(logfc)[c(3:498)]
}


#Harmonisation of row names for each drug GSEA result
bb <- result
for (i in 1:496) {
  for (x in 1:1000) {
    bb[[i]][[x]] <- bb[[i]][[x]][match(rownames(bb[[1]][[1]]),rownames(bb[[i]][[x]])),]
  }
}



#Present the NES values for each herbal ingredient to form a data frame
herb_NES <- vector("list",1000)
for (i in 1:1000) {
  dd <- data.frame(matrix(nrow = 608,ncol = 496))
  for (x in 1:496) {
    dd[[x]] <- bb[[x]][[i]]$NES 
  }
  rownames(dd) <- rownames(bb[[1]][[1]])
  herb_NES[[i]] <- dd
  rm(dd)
}

for (i in 1:1000) {
  colnames(herb_NES[[i]]) <- colnames(logfc)[3:498]
}



#save(result,file = "herb_GSEA1000.Rdata")
#save(herb_NES,file = "herb_NES1000.Rdata")




















