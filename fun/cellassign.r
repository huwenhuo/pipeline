library(cellassign)
library(Seurat)
library(data.table)
setwd('/ifs/work/solitlab/huw/solit/study/hiseq/10xgenomic/')
load('cranger/sq/sq_oo_fil.RData')
source('~/pipeline/fun/param.r')
cellassign.marker.mtx = marker_list_to_mat(cellassign.marker)
mtx = GetAssayData(sq.oo.fil, slot='counts')
cm = intersect(row.names(cellassign.marker.mtx), row.names(mtx))
mtx.sel = mtx[cm, ]
mtx.sel = as.matrix(mtx.sel)
cellassign.marker.mtx = cellassign.marker.mtx[cm, ]
mtx.s = apply(mtx, 2, sum)
sq.cass = cellassign(t(mtx.sel), marker_gene_info=cellassign.marker.mtx, s=mtx.s) 

