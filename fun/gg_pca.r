gg.pca = function(mtx, rn = NULL){
	tmp = apply(mtx, 1, sd)
	tmp = mtx[order(tmp, decreasing = T), ]
	tmp = tmp[1:1000,]
	col = c(rep(2, ncol(mtx)), rep(3, ncol(rsem_shen$counts)))
	names(col)=c(colnames(blca.cc), colnames(rsem_shen$counts))
	col[grep('p', names(col))] = 4
	length(col)
	tmp = prcomp(t(tmp))
	tmp = tmp$x
	tmp = as.data.frame(tmp)

	tmp$samp_class= "TCGA BLCA"
	tmp$samp_class[grep("Org", row.names(tmp))] = "organoids"
	tmp$samp_class[grep("tumor", row.names(tmp))] = "tumor"
	tmp$col = 2
	tmp$col[grep("Org", row.names(tmp))] = 3
	tmp$col[grep("tumor", row.names(tmp))] = 4

	pdf("results/pcav2.pdf", width=8, height=6)
	ggplot(tmp, aes(PC1, PC2, label = samp_class,  color=samp_class)) + geom_point() + theme(legend.position="right") 
	dev.off()
}

gg.pca = function(mtx, ann, topn = NULL, basemean = NULL, rn = NULL){
	mtx=log2(res.cc+1)
	#ann = anno.old
	if(!is.null(topn)){
		row.sd = apply(mtx, 1, sd)
		row.sd = row.sd[order(row.sd, decreasing=T)]
		row.mn = rowSums(mtx)
		names(row.sd[1:topn]) -> sel.id ## top varied genes
	}
	## basemean 
	if(!is.null(basemean)){
		names(row.mn[row.mn > basemean & names(row.mn) %in% sel.id]) -> sel.id
	}
	sel.id
	mtx = mtx[row.names(mtx) %in% sel.id, ]
	dim(mtx)
	pca = prcomp(t(mtx), scale. = TRUE)
	pca.sum = summary(pca) -> pca.sum
	pca = pca$x[, 1:3]
	pca = as.data.table(as.data.frame(pca), keep.rownames = T)
	pca

	pca.sum
	pca.sum.key = pca.sum$importance
	xlab=paste0('PC1 ', round(100*pca.sum.key[2, 1], 0), '%') 
	ylab=paste0('PC2 ', round(100*pca.sum.key[2, 2], 0), '%') 

	source('~/pipeline/g10_all_.r')
	pca = merge.df.fun(pca, ann, if.returndt = T)
	pca

	##gg = ggplot(pca, aes(PC1, PC2, label = grp,  color=Histology)) + geom_point() + theme(legend.position="right")  +
	##	labs(color = 'Histology')
	##ggsave(gg, file='res/gg_pca.pdf', width=5, height=4)

	gg = ggpubr::ggscatter(pca, x='PC1', y='PC2', color=grp) 
	if(!is.null(fname)){
		ggsave(gg, file = fname, width = ww, height = hh )
	}
	gg
}

cal.pca = function(mtx, rsum = 100, topn = 1000, grp = 'grp'){
	if(!exists('loadfun')){ source('~/pipeline/fun/loadfun.r') }
	if(inherits(mtx, 'data.table')){
		    source('~/pipeline/fun/dt_to_df.r')
		    mtx = dt.to.df(mtx)
	}
	r.s = rowSums(mtx)
	mtx = mtx[r.s > rsum, ]
	if(!is.null(topn) & is.numeric(topn) & topn > 0){
		rsd = apply(mtx, 1, sd)
		mtx = mtx[order(rsd, decreasing=T), ]
		mtx = mtx[1:topn, ]
	}
	pca1 = prcomp(t(mtx), scale. = TRUE)
	pca1$x[, 1:3] -> pca.x
	pca.x = as.data.frame(pca.x)
	summary(pca1) -> pca.sum
	pca.sum
	pca.sum.key = pca.sum$importance
	pca.sum$importance
	xlab=paste0('PC1 ', round(100*pca.sum.key[2, 1], 0), '%') ; xlab
	ylab=paste0('PC2 ', round(100*pca.sum.key[2, 2], 0), '%') 
	if(!exists('df.to.dt')){ 
		loadfun('df.to.dt')
	}
	pca.x = df.to.dt(pca.x)
	list(mtx = pca.x, xlab=xlab, ylab=ylab)
}
