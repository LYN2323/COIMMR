rm(list = ls()) 
options=(stringsAsFactors=F)
load(file = "cor.Rdata")
load(file = "herb_NES1000.Rdata")
load(file = "cancer_NES_1000.Rdata")
#免疫

#统一癌症和中药小分子NES值的行名
for (i in 1:1000) {
  herb_NES[[i]] <- herb_NES[[i]][match(rownames(cancer_NES[[1]]),rownames(herb_NES[[i]])),]
}


#计算相关性
corr <- vector("list",1000)
for (i in 1:1000) {
  corr[[i]] <- cor(herb_NES[[i]],cancer_NES[[i]],method = "pearson",use = "complete.obs")
}

save(corr,file = "immune 1000corr.Rdata")












