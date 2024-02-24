.libPaths('/home/huw/anaconda3/envs/r4/lib/R/library/')
#library(reticulate)
library(argparse)
library(Seurat)
library(future)
library(data.table)

options(future.globals.maxSize= 20*1024*1024^2) # 50G
plan('multisession', workers = 5)

parser <- ArgumentParser(description='backrun')
parser$add_argument('-x', '--fname.args', default = NULL)
parser$add_argument('-f', '--fun.call', default = NULL)
parser$add_argument('-i', '--fname.in', default = NULL)
parser$add_argument('-o', '--fname.out', default = NULL)
args <- parser$parse_args()
cat(args$fname.args, '\n')

source('/home/huw/pipeline/work/g10_all_.r')

#arg.lst = list(fun.call = 'findallmarker.fun', fname.in = fname.in, args = list(fname.out = fname.out))

#arg.lst = list(fun.call = 'merge.fun', args = list(sel = sel, meta = meta, merge.way = 'cca', cranger.dir = 'cranger2', if.singlet = F, ncol.cut = 120, 
#						   n.features = 2000, fname.out = fname.out))

if(!is.null(args$fname.args)){
	load(args$fname.args)
	do.call(arg$fun.call, args$args)
	cat('+++++++++++++++++++++\n')
	cat(arg$fun.call, ' finished.\n')
	cat('+++++++++++++++++++++\n')
}

if(is.null(args$fname.args)){
	if(tools::file_ext(args$fname.in) == 'h5seurat'){
		cat('loading Seurat object', args$fname.in, '\n')
		oo  = SeuratDisk::LoadH5Seurat(args$fname.in)
		do.call(args$fun.call, oo, fname.out = args$fname.out)
		cat('+++++++++++++++++++++\n')
		cat(args$fun.call, ' finished.\n')
		cat('+++++++++++++++++++++\n')
	}else{
		stop('nothing done, check the arguments')
	}
}
