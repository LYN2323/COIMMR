rm(list = ls()) 
options=(stringsAsFactors=F)
load(file = "cor.Rdata")
load(file = "immune 1000corr.Rdata")



new_new <- furrr::future_map(1:30, function(x) {
  cv <- purrr::map(1:496, function(n) {
    rv <- purrr::map_dbl(1:1000, function(i) {
      cc <- corr[[i]][n,x]
    })
    rv
  })
})

#Position of cor rows with correlation coefficient < 0
cor_number2 <- vector("list",30)
for (i in 1:30) {
  cor_number2[[i]] <-which(cor[,i]<0)
  names(cor_number2)[[i]] <- colnames(cor)[i]
  
}

#Correlation coefficient of <0
cor3 <- vector("list",30)
for (i in 1:30) {
  cor3[[i]] <- cor[cor_number2[[i]],i]
  names(cor3)[[i]] <- colnames(cor)[i]
}

#Corresponding correlation coefficients to <0 in the permutation test
new3 <- vector("list",30)
for (i in 1:30) {
  names(new3)[[i]] <- colnames(cor)[i]
  if (length(cor_number2[[i]])>0){
    for (x in 1:length(cor_number2[[i]])) {
      new3[[i]][[x]] <- new_new[[i]][[cor_number2[[i]][[x]]]]  
    }
  }
}

#p-value for correlation coefficient < 0  
pp3 <- furrr::future_map(1:30, function(x) {
  cv <- purrr::map_chr(1:length(cor3[[x]]), function(n) {
    sum(new3[[x]][[n]]<=cor3[[x]][[n]])/1000
  })
  cv
})

#Calculating fdr values
#Converted to vectors
c_pp3 <- unlist(pp3)
#Calculating fdr values
fdr3 <- p.adjust(c_pp3,method = "BH")
#Converted to list
lengthh3 <- lengths(pp3)

FDR3 <- vector("list",30)
for (i in 2:30) {
  FDR3[[1]] <- fdr3[1:sum(lengthh3[seq(1,1,1)])]
  FDR3[[i]] <- fdr3[(sum(lengthh3[seq(1,i-1,1)])+1):sum(lengthh3[seq(1,i,1)])]
  names(FDR3)[[i]] <- colnames(cor)[i]
}

#save(a,new,pp,file = "immune1000_result.Rdata")
a3 <- vector("numeric",30)
for (i in 1:30) {
  a3[[i]] <- sum(FDR3[[i]]<0.05)
}
names(a3) <- colnames(cor)
sort(a3)

#position of FDRs significantly associated with the first 3 cancers
site4 <- vector("list",30)
for (i in 1:30) {
  site4[[i]] <- which(FDR3[[i]]<0.05)
  names(site4)[[i]] <- colnames(cor)[i]
}

site5 <- vector("list",30)
for (i in 1:30) {
  names(site5)[[i]] <- colnames(cor)[i]
  if (length(site4[[i]])>0){
    for (x in 1:length(site4[[i]])) {
      site5[[i]][[x]] <- cor_number2[[i]][[site4[[i]][[x]]]]  
    } 
  }
}


#Convert the elements of a list into a vector
for (i in 1:30) {
  site5[[i]] <- unlist(site5[[i]])
}

#Presenting the first 3 rows of cancer in position
name_top3 <- names(sort(a3,decreasing = T))[1:3]
site_top3 <- vector("numeric",3)
for (i in 1:3) {
  site_top3[i] <- which(colnames(cor)==name_top3[[i]])
}
#Position of the cor rows selected for significant correlation in the top 3 cancers
site6 <- vector("list",3)
for (i in 1:3) {
  site6[[i]] <- site5[[site_top3[[i]]]]
  names(site6)[[i]] <- colnames(cor)[site_top3[[i]]]
}

#Which drugs are significantly associated with the top 3 cancers
top3_herbname <- vector("list",3)
for (i in 1:3) {
  top3_herbname[[i]] <- rownames(cor)[site6[[i]]]
  names(top3_herbname)[[i]] <- names(site6)[[i]]
}



#Correlation coefficient data frame for the top 3 cancers
newcor2 <- cor[,c(3,9,23)]
STAD_herb <- vector("numeric",239)
for (i in 1:239) {
  STAD_herb[i] <- newcor2[c(site6[[3]][[i]]),3]
  names(STAD_herb)[i] <- rownames(newcor2)[site6[[3]][[i]]]
}
STAD_top48 <- sort(STAD_herb)[1:48]#前10%的药物
#save(STAD_top24,file = "STAD_top24.Rdata")



load(file = "NES.Rdata")
NES <- as.data.frame(t(as.matrix(NES)))
#NES values for drugs proposed for STADtop48
STAD_NES <- NES$STAD
names(STAD_NES) <- rownames(NES)

#save(STAD_herb1,file = "LGG_STAD_top24nes<0.Rdata")



#The  contribution for STAD top48 drugs
load(file = "NES&corr.Rdata")
STAD_herb_NEScor <- herb_cor[[1]]
names(STAD_herb_NEScor) <- rownames(herb_cor)
STAD_herb_NEScor <- sort(STAD_herb_NEScor,decreasing = T)[1:48]

#save(STAD_herb_NEScor,file = "LGG_STAD_TOP24_nescor.Rdata")


STAD_herb <-Reduce(intersect,list(names(STAD_top48),names(STAD_herb1),names(STAD_herb_NEScor)))
 
STAD_herbrank <- data.frame(a=STAD_top48[STAD_herb],B=STAD_herb1[STAD_herb],C=STAD_herb_NEScor[STAD_herb])
#save(STAD_herbrank,file = "immune herb rank1111.Rdata")


