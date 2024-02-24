
options(future.globals.maxSize= 50*1024*1024^2) # 50G

library(Seurat)
library(dplyr)
library(future)
library(data.table)
library(ggplot2)
library(celldex)
library(SingleR)
library(SeuratDisk)
library(ComplexHeatmap)
library(log4r)

plan('sequential')
plan('multisession', workers = 10)

setwd('/juno/work/solitlab/huw/solit/study/hiseq/10xgenomic')
source('~/pipeline/work/g10_all_.r')

logfile = paste0(getwd(), '/log'); logfile
logger =  create.logger()
logfile(logger) = logfile
level(logger) = 'INFO' 

source('~/pipeline/work/g10_all_.r')
sels = c(
	'DS008',
	'DS010',
	'DS015_1',
	'DS015_2',
	'DS018',
	'DS019',
	'DS049',
	'DS052',
	'DS054L',
	'DS054P',
	'DS063',
	'DS043',
	'DS056',
	'DS066',
	'DS067',
	'DS068',
	'DS069'
)

sel = 'DS054P'
cwd = getwd()
future.apply::future_lapply(sel, function(sel){setwd(cwd); monocle3.fun(sel = sel); 'done'})  
