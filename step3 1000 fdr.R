rm(list = ls()) 
options=(stringsAsFactors=F)
load(file = "cor.Rdata")
load(file = "immune 1000corr.Rdata")

#Organize the correlation coefficients in the permutation test, each element being 10,000 correlation coefficients for the corresponding position of cor
new <- furrr::future_map(1:496, function(x) {
  cv <- purrr::map(1:30, function(n) {
    rv <- purrr::map_dbl(1:1000, function(i) {
      cc <- corr[[i]][x,n]
    })
    rv
  })
})


#Position of columns with correlation coefficient < 0
cor_number <- vector("list",496)
for (i in 1:496) {
  cor_number[[i]] <-which(cor[i,]<0)
  names(cor_number)[[i]] <- rownames(cor)[i]
}

#Correlation coefficient of <0
cor1 <- vector("list",496)
for (i in 1:496) {
  cor1[[i]] <- cor[i,cor_number[[i]]]
  names(cor1)[[i]] <- rownames(cor)[i]
}

#Corresponding correlation coefficients  <0 in the permutation test
new1 <- vector("list",496)
for (i in 1:496) {
  
  if (length(cor_number[[i]])>0){
    for (x in 1:length(cor_number[[i]])) {
      new1[[i]][[x]] <- new[[i]][[cor_number[[i]][[x]]]] 
    }
  }
}

#p-value for correlation coefficient < 0 
new1[[433]] <- NA  
cor1[[433]] <- NA
pp1 <- furrr::future_map(1:496, function(x) {
  cv <- purrr::map_chr(1:length(cor1[[x]]), function(n) {
    sum(new1[[x]][[n]]<=cor1[[x]][[n]])/1000
  })
  cv
})

#Converted to vectors
c_pp <- unlist(pp1)
#Calculating FDR values
fdr1 <- p.adjust(c_pp,method = "BH")
#Converted to list
lengthh1 <- lengths(pp1)

FDR1 <- vector("list",496)
for (i in 2:496) {
  FDR1[[1]] <- fdr1[1:sum(lengthh1[seq(1,1,1)])]
  FDR1[[i]] <- fdr1[(sum(lengthh1[seq(1,i-1,1)])+1):sum(lengthh1[seq(1,i,1)])]
  names(FDR1)[[i]] <- rownames(cor)[i]
}



#save(a,new,pp,file = "immune1000_result.Rdata")
a1 <- vector("numeric",496)
for (i in 1:496) {
  a1[[i]] <- sum(FDR1[[i]]<0.05)
}
sum(a1[1:432],a1[434:496])
names(a1) <- rownames(cor)
sort(a1,decreasing = T)
names(sort(a1,decreasing = T))[1:50]



#Position of correlation coefficient > 0
cor_number1 <- vector("list",496)
for (i in 1:496) {
  cor_number1[[i]] <-which(cor[i,]>0)
  names(cor_number1)[[i]] <- rownames(cor)[i]
}
#Correlation coefficient for >0
cor2 <- vector("list",496)
for (i in 1:496) {
  cor2[[i]] <- cor[i,cor_number1[[i]]]
  names(cor2)[[i]] <- rownames(cor)[i]
}
#Corresponding correlation coefficient  >0 in the permutation test
new2 <- vector("list",496)
for (i in 1:496) {
  
  if (length(cor_number1[[i]])>0){
    for (x in 1:length(cor_number1[[i]])) {
      new2[[i]][[x]] <- new[[i]][[cor_number1[[i]][[x]]]]  
    }
  }
}
#p-value for correlation coefficient > 0 
new2[[494]] <- NA  
cor2[[494]] <- NA
new2[[346]] <- NA  
cor2[[346]] <- NA
#p-value for correlation coefficient > 0 
pp2 <- furrr::future_map(1:496, function(x) {
  cv <- purrr::map_chr(1:length(cor2[[x]]), function(n) {
    sum(new2[[x]][[n]]>=cor2[[x]][[n]])/1000 #计算p值
  })
  cv
})

#Converted to vectors
c_pp2 <- unlist(pp2)
#Calculating FDR
fdr2 <- p.adjust(c_pp2,method = "BH")
#Converted to list
lengthh2 <- lengths(pp2)

FDR2 <- vector("list",496)
for (i in 2:496) {
  FDR2[[1]] <- fdr2[1:sum(lengthh2[seq(1,1,1)])]
  FDR2[[i]] <- fdr2[(sum(lengthh2[seq(1,i-1,1)])+1):sum(lengthh2[seq(1,i,1)])]
  names(FDR2)[[i]] <- rownames(cor)[i]
}


a2 <- vector("numeric",496)
for (i in 1:496) {
  a2[[i]] <- sum(FDR2[[i]]<0.05)
}
sum(a2[1:345],a2[347:493],a2[495:496])
names(a2) <- rownames(cor)



##cancer
#Organize the correlation coefficients in the permutation test, each element being 10,000 correlation coefficients for the corresponding position of cor
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

#Corresponding correlation coefficients<0 in the substitution test
new3 <- vector("list",30)
for (i in 1:30) {
  names(new3)[[i]] <- colnames(cor)[i]
  if (length(cor_number2[[i]])>0){
    for (x in 1:length(cor_number2[[i]])) {
      new3[[i]][[x]] <- new_new[[i]][[cor_number2[[i]][[x]]]]  # new 里面的相关系数小于0的，计算P值
    }
  }
}

#p-value for correlation coefficient < 0
pp3 <- furrr::future_map(1:30, function(x) {
  cv <- purrr::map_chr(1:length(cor3[[x]]), function(n) {
    sum(new3[[x]][[n]]<=cor3[[x]][[n]])/1000#(cor[x,n]<0))/1000  #计算p值
  })
  cv
})

#Calculating FDR
#Converted to vectors
c_pp3 <- unlist(pp3)
#Calculating FDR
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
#sum(a1[1:432],a1[434:496])
names(a3) <- colnames(cor)
sort(a3)
#names(a3)[which(a3>240)]


#save(a1,a2,a3,file = "immune1000_result.Rdata")
