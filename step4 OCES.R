rm(list = ls()) 
options=(stringsAsFactors=F)
load(file = "NES.Rdata")
NES_ <- t(as.matrix(NES))


pheatmap::pheatmap(NES_,
                   scale = "row",
                   color = colorRampPalette(c("#00008B", "white", "#ff9900"))(256),
                   cellheight = 2.5,
                   cellwidth = 2.5,
                   border_color = NA,
                   cluster_rows = T,
                   cluster_cols = T,
                   fontsize_row = 3,
                   fontsize_col = 6)



standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}


NES_ <- standarize.fun(NES_,halfwidth = 2) 



col_range = c(-2, 0, 2)
col_fun <- colorRamp2(
  col_range, 
  c("#00008B", "white", "#ff9900")
)


col_fun <- colorRamp2(
  c(-3, 0, 3), 
  c("#344747", "#ffffff", "#e3181f")
)

col_fun <- colorRamp2(
  c(-3, 0, 3), 
  c("#344747", "#ffffff", "#e3181f")
)


ht <- Heatmap(
  matrix = NES_,
  col = col_fun,
  name = "ht1",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_column_names = T,
  show_row_names = F,
  show_heatmap_legend = FALSE,
  row_names_rot = 45,
  #column_split = 3,
  #row_split = NA,
  row_gap = unit(3, "mm"),
  column_gap = unit(3, "mm"),
  border = TRUE,
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
  labels = c("low","","high"), 
  direction = "horizontal",
  legend_width = unit(1, "in"),
  border = FALSE
)

draw(ht, heatmap_legend_list = list(lgd), heatmap_legend_side = c("bottom"))





