
get.geneset.genes = function(geneset=''){
	#geneset = 'HOLLERN_SQUAMOUS_BREAST_TUMOR'
	gs = system(paste0('grep ', geneset, ' /home/huw/program/gmt/v6.2/msigdb.v6.2.symbols.gmt'), intern=T)
	gs = unlist(strsplit(gs, '\t'))
	gs = gs[3:length(gs)]
	gs
}
