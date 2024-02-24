consensus_redraw=function(oo, optk=3, fname, ww = 10, hh = 10, anno.col=NULL){
	fm = oo[[optk]][[4]]
	hc = oo[[optk]][['consensusTree']] #=hclust( as.dist( 1 - fm ), method= 'average');
	ct = oo[[optk]][['consensusClass']]; ct

	if(is.null(anno.col)){
		tmp = oo[[optk]]
		anno.col = as.data.frame(tmp[[3]])
		colnames(anno.col) = 'cluster'
		anno.col$cluster = paste0('c', anno.col$cluster)
	}

	pc = fm 
	row.names(pc) = names(ct)
	colnames(pc) = names(ct)
	pc = pc[hc$order,hc$order] 

	breaks = c( seq(from = 0, to = 1, length.out = 9));breaks = unique(breaks);breaks
	library(RColorBrewer)
	clr = brewer.pal(length(breaks + 1), 'Blues'); clr
	#library(ColorRampPalette)
	#clr = colorRampPalette(clr)

	onco.col = 'forestgreen'; lonco.col= adjustcolor('dodgerblue', alpha.f = .7)
	anno.color = list(
			  cluster =c(c1 = 'slateblue1', c2 = 'orange3', c3 = 'turquoise'), 
			  FGFR3 =  c('Oncogenic' = onco.col, 'Likely Oncogenic' = lonco.col), 
			  FGFR2 =  c('Oncogenic' = onco.col, 'Likely Oncogenic' = lonco.col), 
			  TP53 =   c('Oncogenic' = onco.col, 'Likely Oncogenic' = lonco.col), 
			  BRAF =   c('Oncogenic' = onco.col, 'Likely Oncogenic' = lonco.col), 
			  KRAS =   c('Oncogenic' = onco.col, 'Likely Oncogenic' = lonco.col), 
			  HRAS =   c('Oncogenic' = onco.col, 'Likely Oncogenic' = lonco.col), 
			  ERBB2 =  c('Oncogenic' = onco.col, 'Likely Oncogenic' = lonco.col), 
			  PIK3CA = c('Oncogenic' = onco.col, 'Likely Oncogenic' = lonco.col) )

	onco.col = 'forestgreen'; lonco.col= adjustcolor('dodgerblue', alpha.f = .7)
	anno.color = list(
			  cluster =c(c1 = 'slateblue1', c2 = 'orange3', c3 = 'turquoise'), 
			  pT.Stag2 =c(pT4 = onco.col, T4 = onco.col, pT3 = adjustcolor(onco.col, alpha.f=.7), 
				      pT3b = adjustcolor(onco.col, alpha.f=.7), 
				      pT2 = adjustcolor(onco.col, alpha.f=.4), 
				      pT1 = adjustcolor(onco.col, alpha.f=.4), 
				      pTa = adjustcolor(onco.col, alpha.f=.4)), 
			  pN.Stag2 =c('pN1/2' = onco.col), 
			  Smoking.Status.at.Diagnosis=c(Active = onco.col, Former = adjustcolor(onco.col, alpha.f=.7), Never = adjustcolor(onco.col, alpha.f=.4)), 
			  FGFR3 =  c('Oncogenic' = onco.col, 'Likely Oncogenic' = lonco.col), 
			  FGFR2 =  c('Oncogenic' = onco.col, 'Likely Oncogenic' = lonco.col), 
			  TP53 =   c('Oncogenic' = onco.col, 'Likely Oncogenic' = lonco.col), 
			  BRAF =   c('Oncogenic' = onco.col, 'Likely Oncogenic' = lonco.col), 
			  KRAS =   c('Oncogenic' = onco.col, 'Likely Oncogenic' = lonco.col), 
			  HRAS =   c('Oncogenic' = onco.col, 'Likely Oncogenic' = lonco.col), 
			  ERBB2 =  c('Oncogenic' = onco.col, 'Likely Oncogenic' = lonco.col), 
			  PIK3CA = c('Oncogenic' = onco.col, 'Likely Oncogenic' = lonco.col) )

	#anno.color = list( c1 = 'slateblue1', c2 = 'orange3', c3 = 'turquoise', 'Oncogenic' = onco.col, 'Likely Oncogenic' = lonco.col)
	head(anno.col)
	tmp = anno.col[, c(1,2,3, 7,6,8,9,10,5)]
	hp = pheatmap(pc, color=clr, scale='none', show_rownames = F, breaks = breaks, border_color = NA, filename = fname, width=ww, height=hh,
		      main='',  annotation_col = tmp, annotation_colors = anno.color, cluster_rows = F, cluster_cols = F, show_colnames = F, fontsize_row = 6)
	scp(fname)

	pc

}


consensus_cluster = function(xx, pdffile, fdr='consensuscluster', optK, maxK, reps, isScaled=F){
	source('/home/huw/program/ConsensusClusterPlus/R/ConsensusClusterPlus.R')
	if(!isScaled){
		cn = colnames(xx)
		xx= t(apply(xx, 1, scale))
		colnames(xx) = cn
	}

	tf = apply(!is.na(as.matrix(xx)), 1, all)
	xx = xx[tf, ]

	if(missing(maxK)){maxK = 5}
	if(missing(reps)){reps = 500000}
	if(missing(pdffile)){cat("missing pdf file name\n"); return ;}

	pdf(pdffile)
	res = ConsensusClusterPlus(xx, maxK=5,reps=1000,pItem=0.8,pFeature=1, 
				   title=fdr,clusterAlg="hc",distance="pearson",seed=1262118388.71279)
	dev.off()

	##https://www.biostars.org/p/198789/
	if(!missing(optK)){
		Kvec = 2:Kmax
		x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
		PAC = rep(NA,length(Kvec)) 
		names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK
		for(i in Kvec){
			M = results[[i]]$consensusMatrix
			Fn = ecdf(M[lower.tri(M)])
			PAC[i-1] = Fn(x2) - Fn(x1)
		}
		# The optimal K
		optK = Kvec[which.min(PAC)]
		cat("optK is ", optK, "\n");
	}

	optClass = res[[optK]]$consensusClass
	res

}
