rm(list = ls()) 
options=(stringsAsFactors=F)
load(file = "immune NES.Rdata")


#Calculating correlation coefficients for cancer and small molecule enrichment fractions
pancancer_NES <- pancancer_NES[match(rownames(herb_NES),rownames(pancancer_NES)),]

#Calculating the correlation coefficient
cor <- cor(herb_NES,pancancer_NES,method = "pearson")

library(dendextend)
library(ggsci)
library(rio)
library(circlize)
library(ComplexHeatmap)


standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}


cor <- standarize.fun(cor)

### Heatmap
col_range = c(-1, 0, 1)
col_fun <- colorRamp2(
  col_range, 
  c("#bbc622", "white", "#3b2d8e")
)  
  
  ht <- Heatmap(
    matrix = cor,
    col = col_fun,
    name = "ht1",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_column_names = T,
    show_row_names = F,
    show_heatmap_legend = FALSE,
    row_names_rot = 45,
    row_names_gp = gpar(
      fontsize = 5
    ),
    column_names_rot = 45,
    column_names_gp = gpar(
      fontsize = 10
    )
  )
  
  lgd <- Legend(
    col_fun = col_fun, 
    title = "", 
    at = col_range, 
    labels = c("-1","","1"), 
    direction = "horizontal",
    legend_width = unit(1, "in"),
    border = FALSE
  )
  
draw(ht, heatmap_legend_list = list(lgd), heatmap_legend_side = c("bottom"))
  






