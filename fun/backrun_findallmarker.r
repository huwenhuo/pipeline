.libPaths('/home/huw/anaconda3/envs/r4/lib/R/library/')
#library(reticulate)
library(argparse)
library(Seurat)
library(future)
library(data.table)

plan('multisession', workers = 5)

parser <- ArgumentParser(description='findallmarker')
parser$add_argument('-i', '--fname.in')
parser$add_argument('-o', '--fname.out')
args <- parser$parse_args()
cat(args$fname.oo, '\n')

oo = SeuratDisk::LoadH5Seurat(args$fname.in)
oo.marker = FindAllMarkers(object = oo, only.pos = F, assay = 'RNA')
oo.marker = as.data.table(oo.marker, keep.rownames=T)
save(oo.marker, file = args$fname.out)
