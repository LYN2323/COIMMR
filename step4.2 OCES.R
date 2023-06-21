rm(list = ls()) 
options=(stringsAsFactors=F)
load(file = "cancer_signature.Rdata")
load(file = "logfc20210824.Rdata")

#Cancer up signature enriched to herbal ingredients induced expression profile
up_nes <- vector("list", 496)
up_gsea <- vector("list",496)#GSEA的总的结果
up_GSEA <- vector("list",496)#GSEA提取result
for (i in 1:496) {
  print(i)
  herb_logfc <- logfc[[i+2]]
  names(herb_logfc) <- logfc[[1]]
  herb_logfc <-  sort(herb_logfc,decreasing = T)
  up_gsea[[i]] <- GSEA(herb_logfc,TERM2GENE = up_signature,pvalueCutoff=1,nPermSimple = 100000)
  up_GSEA[[i]] <- up_gsea[[i]]@result
  names(up_GSEA)[[i]]=colnames(logfc)[i+2]
} 

#Separate enrichment analysis for Trigonelline
herb_logfc <- logfc$Trigonelline
names(herb_logfc) <- logfc[[1]]
herb_logfc <-  sort(herb_logfc,decreasing = T)
Trigonelline <- GSEA(herb_logfc,TERM2GENE = up_signature,pvalueCutoff=1,nPermSimple = 1000000)
Trigonelline <- Trigonelline@result
#Separate enrichment analysis for Roburic acid
herb_logfc <- logfc$`Roburic acid`
names(herb_logfc) <- logfc[[1]]
herb_logfc <-  sort(herb_logfc,decreasing = T)
Roburic_acid <- GSEA(herb_logfc,TERM2GENE = up_signature,pvalueCutoff=1,nPermSimple = 1000000)
Roburic_acid <- Roburic_acid@result

up_GSEA$Trigonelline <- Trigonelline
up_GSEA$`Roburic acid` <- Roburic_acid

for (i in 1:495) {
  #Match the row names of each herbal ingredient up_GSEA result so that the row names of each GSEA knotted herbal ingredient are in the same order
  up_GSEA[[i+1]] <- up_GSEA[[i+1]][match(rownames(up_GSEA[[1]]),rownames(up_GSEA[[i+1]])),]
}

#Calculate the up_NES value for each herbal ingredient  
for (i in 1:496) { 
  up_nes[[i]] <- up_GSEA[[i]]$NES
  names(up_nes)[[i]]=colnames(logfc)[i+2]
}

#Combine the up_NES values for each herbal ingredient
up_NES <- do.call(cbind,up_nes)
rownames(up_NES) <- rownames(up_GSEA[[1]])


#Cancer down signature enrichment to herbal ingredisents induced expression profiles
down_nes <- vector("list", 496)
down_gsea <- vector("list",496)
down_GSEA <- vector("list",496)
for (i in 1:496) {
  print(i)
  herb_logfc <- logfc[[i+2]]
  names(herb_logfc) <- logfc[[1]]
  herb_logfc <-  sort(herb_logfc,decreasing = T)
  down_gsea[[i]] <- GSEA(herb_logfc,TERM2GENE = down_signature,pvalueCutoff=1,nPermSimple = 100000)
  down_GSEA[[i]] <- down_gsea[[i]]@result
  names(down_GSEA)[[i]]=colnames(logfc)[i+2]
} 

for (i in 1:495) {
  #Match the row names of each herbal ingredients down_GSEA result so that the row names of each herbal ingredients GSEA result are in the same order
  down_GSEA[[i+1]] <- down_GSEA[[i+1]][match(rownames(down_GSEA[[1]]),rownames(down_GSEA[[i+1]])),]
}

#Calculate the down_NES value for each herbal ingredient  
for (i in 1:496) { 
  down_nes[[i]] <- down_GSEA[[i]]$NES
  names(down_nes)[[i]]=colnames(logfc)[i+2]
}

#Combine the down_NES values for each herbal ingredient
down_NES <- do.call(cbind,down_nes)
rownames(down_NES) <- rownames(down_GSEA[[1]])



##Merge data frames and subtract values 
down_NES <- down_NES[match(rownames(up_NES),rownames(down_NES)),]
a <- cbind(up_NES,down_NES)
#a <- as.data.frame(NES)

NES <- vector("list",496)
for(i in 1:496){
  NES[[i]] <- a[,i]-a[,i+496]
  names(NES)[[i]] <- colnames(a)[i]
}

NES <- do.call(cbind,NES)
rownames(NES) <- rownames(a)

save(NES,file = "NES.Rdata")


NES_ <- t(as.matrix(NES))
#NES1 <- NES_[c(1:248),]
#NES2 <- NES_[c(249:496),]

pheatmap::pheatmap(NES_,
                   scale = "row",
                   color = colorRampPalette(c("#00008B", "white", "#ff9900"))(256),
                   cellheight = 2.5,
                   cellwidth = 10,
                   border_color = NA,
                   cluster_rows = T,
                   cluster_cols = T,
                   fontsize_row = 3,
                   fontsize_col = 6)
                   #show_rownames = F)
