.libPaths('/home/huw/anaconda3/envs/r4/lib/R/library/')
#library(reticulate)
library(argparse)
library(Seurat)
library(future)
library(data.table)

plan('multisession', workers = 5)

parser <- ArgumentParser(description='merge.cca')
parser$add_argument('-i', '--fname.in')
parser$add_argument('-o', '--fname.out')
args <- parser$parse_args()
cat(args$fname.oo, '\n')

source('/home/huw/pipeline/work/g10_all_.r')

oo = SeuratDisk::LoadH5Seurat(args$fname.in)
oo.singler = singler(oo, nn = 300)
oo.singler.dt = as.data.table(oo.singler, keep.rownames=T)
save(oo.singler.dt, file=args$fname.out)
