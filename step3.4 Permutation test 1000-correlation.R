rm(list = ls()) 
options=(stringsAsFactors=F)
load(file = "cor.Rdata")
load(file = "herb_NES1000.Rdata")
load(file = "cancer_NES_1000.Rdata")
#Immune

#Harmonisation of line names for cancer and herbal small molecule NES values
for (i in 1:1000) {
  herb_NES[[i]] <- herb_NES[[i]][match(rownames(cancer_NES[[1]]),rownames(herb_NES[[i]])),]
}


#计算相关性
corr <- vector("list",1000)
for (i in 1:1000) {
  corr[[i]] <- cor(herb_NES[[i]],cancer_NES[[i]],method = "pearson",use = "complete.obs")
}

save(corr,file = "immune 1000corr.Rdata")












