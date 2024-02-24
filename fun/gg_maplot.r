 wi = 6
 hi = 4
 fdr= 0.05
 log2fc = 1
 filename=''
 topn = 10
 genes = ''
 labx = 'log2 Fold Change'
 lab=F
 laby="-log10(adjusted p value)"
 max.overlaps = Inf
gg_maplot = function(dt, wi = 6, hi = 4, fdr= 0.05, log2fc = 1, filename='', topn = 10, genes = '', labx = 'log2 Fold Change', lab=F, laby="-log10(adjusted p value)", max.overlaps = Inf){
	require(ggplot2)
	require(ggrepel)
	require(dplyr)
	dt = dt[!is.na(padj),] 
	dt[, yy := -log10(padj)]
	dt = dt[order(padj),]
	#dt[, yy := ifelse(abs(yy) <= 1, 0, sign(yy)*log10(abs(yy)))]
	#dt[, yy := sqrt(yy1)] 
	dt[, Significant := 'NS']
	dt[log2FoldChange > log2fc  & padj < fdr, Significant := 'UP']
	dt[log2FoldChange < -log2fc & padj < fdr, Significant := 'DN']
	tt = dt[F,]
	if(topn > 0){
		tt1 = dt %>% filter(padj < fdr) %>% filter(log2FoldChange > 0) %>% head(n=topn) ; tt1
		tt2 = dt %>% filter(padj < fdr) %>% filter(log2FoldChange < 0) %>% head(n=topn) ; tt2
		tt = rbind(tt2, tt1)
	}
	if(length(genes) > 0){
		tt3 = dt[symbol %in% genes,]
		tt = rbind(tt, tt3)
	}
	dt[symbol %in% tt$symbol , Significant := 'Spot']
	tt[, Significant := 'Spot']
	myColors = c(NS = adjustcolor('black', .7), UP = adjustcolor('red', alpha.f = .7), DN = adjustcolor('blue', alpha.f=.7), Spot = adjustcolor('green', .9))
	myColors
	gg = ggplot() + 
		geom_point(data=dt, aes(x = log2FoldChange, y = yy, color=Significant), size=1)  + 
		geom_point(data=tt, aes(x = log2FoldChange, y = yy, color=Significant), size=1)  + 
		theme(legend.position="none")  + 
		scale_colour_manual(name = "Significance", values = myColors) + 
		scale_y_continuous(trans = "log1p") + 
		theme_bw(base_size = 8) + theme(legend.position='none') + ylab(laby) + xlab (labx) 
	if(lab == T){
		gg = gg + geom_text_repel(data=tt, aes(x = log2FoldChange, y = yy, label = symbol),  max.overlaps = Inf,
					   size=2, segment.color = 'grey50') 
	}
	if(filename != ''){
		ggsave(gg, file=filename, width=wi, height=hi)
	}
	gg
}
