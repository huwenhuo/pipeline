scan.dir = function(pat = NULL, dd = NULL, if.fastq = T){
	if(is.null(dd)){ stop('no folder provided') }
	pat = '07697'
	dd = '/igo/delivery/share/alahmadh/'
	fd = dir(dd, pattern = pat)
	fd
	if(if.fastq){
		fd.list = sapply(fd, function(x){system(paste0('find -L ', dd, '/', x, ' -name *.fastq.gz'), intern=T)})
		fd.list = unlist(fd.list)
		fd.list
		fd.dt = data.table(fq = fd.list)
		fd.dt[, fastq.dir := paste0(dirname(fq), '/')]
		fd.dt[, sample.name := sub('.*Sample_(.*)_IGO.*', '\\1', basename(fastq.dir))]
		fd.dt = fd.dt[!duplicated(fastq.dir)]
		fd.dt
		fd.dt[grep('X1CF', sample.name), ]
	}


}

heatmap1.fun = function(mtx, ann = NULL, ann.row = NULL, ha = NULL, ra = NULL, ww = 10, hh = 8, fname = NULL, ...){
	par.list = rlang::dots_list(...)

	mtx = rowname.fun(mtx)
	mtx = mtx.scale.fun(mtx)
	if(!is.null(ann)){ 
		ann = dt.to.df(ann) 
		ov = intersect(colnames(mtx), rownames(ann))
		mtx = mtx[, ov]
		ann = ann[ov, ]
		ha = ComplexHeatmap::HeatmapAnnotation(df = ann, annotation_name_side = "left" )
		par.chp$top_annotation  = ha
	}
	if(!is.null(ann.row)){ 
		ann.row = dt.to.df(ann.row) 
		ov = intersect(rownames(mtx), rownames(ann.row))
		mtx = mtx[ov, ]
		ann.row = ann.row[ov, ]
		ra = ComplexHeatmap::rowAnnotation(df = ann.row)
		par.chp$left_annotation = ra
	}
	par.list[['matrix']] = mtx
	par.list[['top_annotation']] = ha
	par.list[['left_annotation']] = ra
	gg = do.call('Heatmap', args = par.list)

	message("plotting ", fname)
	pdf(file=fname, width = ww, height = hh)
	ComplexHeatmap::draw(gg)
	dev.off()
}

mtx.scale.fun = function(mtx, if.log2 = T, if.scale.row = T, if.percent.col = F, mm = 2, if.sig = T, res = NULL, genes = NULL){
	if(if.sig == T & !is.null(res)){ mtx = mtx[row.names(mtx) %in% res[abs(log2FoldChange) > 1 & padj < 0.05, rn], ] }
	mtx = rowname.fun(mtx)
	if(!is.null(genes)){
		mtx = mtx[rownames(mtx) %in% genes, ]
	}
	if(if.log2 == T){
		if(max(mtx) > 100){ mtx = log2(mtx + 1) }
	}
	if(if.percent.col){ 
		mtx = RcmdrMisc::colPercents(mtx) 
		mtx = mtx[1:(nrow(mtx)-2), ]
	}
	if(if.scale.row & max(mtx) > 5 ){ # in case the rows has been scaled
		cn = colnames(mtx)
		mtx = t(scale(t(mtx)))
		colnames(mtx) = cn
		if(!is.null(mm)){
			mtx[mtx > mm] = mm
			mtx[mtx < -mm] = -mm
		}
	}
	mtx = mtx[!is.na(mtx[, 1]), ]
	mtx
}

rowname.fun = function(mtx, rname = NULL, rn.out = 'symbol'){
	if(grepl('ENSG', row.names(mtx)[1])){ row.names(mtx) = substr(row.names(mtx), 17, length(row.names(mtx))) }
	mtx = mtx[!duplicated(rownames(mtx)), ]
	mtx
}

dotplot.group.fun = function(mtx, genes = NULL, group.by = NULL, ann = NULL, ii.ww = 3, ii.hh = 3, fname = NULL){
	mtx = rowname.fun(mtx)
	mtx = mtx.scale.fun(mtx, if.scale.row = F)

	dd = melt(mtx[rownames(mtx) %in% genes, ])
	dd = merge(dd, ann[, c('rn', group.by), with = F], by.x = 'Var2', by.y = 'rn')
	dd = as.data.table(dd)
	dd = dd[!is.na(group.by), ]
	dd[, val := value]
	
	gg.list = lapply(genes, function(xx){ggpubr::ggboxplot(dd[Var1 == xx, ], x = group.by, y = 'val', title = xx, ylab = 'log2 #reads', 
							       add = 'jitter', fill = group.by) +
			 ggpubr::stat_compare_means() + ggpubr::rotate_x_text() + ggpubr::rremove('legend')})
	gg = gridExtra::marrangeGrob(gg.list, ncol=3, nrow = 2, as.table = T)
	if(!is.null(fname)){
		ggplot2::ggsave(gg, file = fname, width = ii.ww * 3, height = ii.hh * 2)
		message(fname, ' saved')
		dd
	}else{
		gg
	}
}

geneset.heatmap.fun = function(fl.dt = NULL, res.cc = NULL, res = NULL, geneset.name.list = NULL, topn = 5, geneset.list = c('h.all', 'c7.all', 'c3.tft'), bbase = NULL){
	if(is.null(fl.dt)){ fl.dt = gsea.read.fun() }
	tmp = copy(fl.dt)

	if(!is.null(base)){
		tmp = tmp %>% filter(base %in% bbase)
	}

	if(is.null(geneset.name.list)){
		if(!is.null(topn)){
			tmp = tmp %>% group_by(geneset) %>% top_n(topn, wt = qval)
		}
		if(!is.null(geneset.list)){
			tmp = tmp %>% filter(geneset %in% geneset.list)
		}
	}else{
		tmp = tmp %>% filter(name %in% geneset.name.list)
	}

	geneset.name.list = tmp$name
	geneset.name.list
	message('selected geneset:\n',  paste(geneset.name.list, collapse = '\n'))
	lapply(geneset.name.list, function(geneset.name){
		       message(' dealting with geneset ', geneset.name)
		       genes = get.geneset.genes(geneset.name)
		       if(!is.null(res)){
			       genes = unique(res[symbol %in% genes,][padj < 0.05,  symbol])
			       genes
		       }
		       if(length(genes) < 5){ message(geneset.name, ' has ', length(genes), ', skipped')}else{
			       genes
			       tmp = rowname.fun(res.cc)
			       mtx = log2(tmp[rownames(tmp) %in% genes, ] + 1)
			       mtx = mtx.scale.fun(mtx)
			       mtx = mtx[, dsn$sample.name.2]
			       head(mtx)
			       fname = paste0('res/heatmap_geneset_', geneset.name, '.pdf'); fname
			       tmp = data.frame(grp = dsn$grp, row.names = dsn$sample.name.2)
			       ww = ncol(mtx) * 0.5 + 2
			       hh = nrow(mtx) * 0.13 + 1 
			       message(ncol(mtx), ' x ', nrow(mtx), ' in mtx')
			       heatmap.fun(mtx, ann = tmp, fname = fname, show_column_dend = F, show_row_dend = F, cluster_columns = F, column_split = tmp$grp, hh =hh, ww = ww) 
			       fname
		       }
		      })
}

read.gsea.fun = function(idir = 'res/gsea', gsea.list = NULL){
	
	if(is.null(gsea.list)){
   		gsea.list = system(paste0('find ', idir, ' -name "*gsea_report_for*.xls"'), intern=T)
		gsea.list
	}

	lapply(gsea.list, fread) -> fl.list
	names(fl.list) = gsea.list
	as.data.table(rbindlist(fl.list, idcol = TRUE)) -> fl.dt
	fl.dt[,  base := sub('.*gsea\\.(.*).rnk.*', '\\1', .id)]
	fl.dt[, geneset := sub(".*rnk.(.*).v6.2.symbols.gmt.G.*", "\\1", .id)]
	fl.dt[, tag := paste0( base, ':', base)]
	fl.dt[NES < 0, updn := 'DN']
	fl.dt[NES > 0, updn := 'UP']
	fl.dt[, qval := as.numeric(`FDR q-val`)]
	fl.dt[, name := NAME]
	fl.dt[, NES := as.numeric(NES)]
	fl.dt = fl.dt[order(qval), ]
	fl.dt[, col.fill := ifelse(sign(NES)>0, 'up', 'dn')]
	fl.dt
}

sum.gsea = function(idir = 'res/gsea', plot=T, gsea.list = NULL){
	
	fl.dt = read.gsea.fun(idir = idir, gsea.list = gsea.list)

	require(dplyr)
	require(gridExtra)
	require(ggforce)

	tmp = fl.dt  %>% filter(qval < 0.05) %>% group_by(base, geneset, sign(NES)) %>% 
		top_n(10, wt = abs(NES)) %>% as.data.table()
	tmp

	for(ii in unique(tmp$base)){
		for(jj in unique(tmp$geneset)){
			tt = tmp[base == ii & geneset == jj, ]
			gg = ggplot(tt, aes(x=reorder(name, NES), y=NES, fill = col.fill)) + 
				geom_bar(width=0.6, height=.05, stat = 'identity') 
			gg = gg+ theme(axis.text.x = element_text(angle=-90, vjust = 0.5, hjust=0)) + 
				ylab("NES") + xlab("Gene set names") + coord_flip() 
			fname = paste0('res/gsea/gsea_sum_', ii, '_', jj, '.pdf'); fname
			hh = 1 + .2*nrow(tt)
			ggsave(gg, file=fname, width=15, height=hh, limitsize=F) # 
			message(fname, ' plotted')
		}
	}
	fl.dt
}

run.gsea.fun = function(rnkfile, cmdfile=NULL, gmtfile=NULL, cwd=NULL, nperm=1000, mem=6, cpu=5, if.run=T, if.svg=F, if.cluster = T){
	require(tools)
	if(missing(rnkfile)){
		print("rnkfile must provided\n")
		return
	}else{
		if(!all(file.exists(rnkfile))){
			return("rnkfile is not all available")
		}
	}
	rnkfile = normalizePath(rnkfile)
	rnkfile

	if(missing(cwd)){ cwd = getwd(); }
	if(missing(gmtfile)){
		target = c("c1.all.v6.2.symbols.gmt" ,
			 "c2.cgp.v6.2.symbols.gmt" ,
			 "c2.cp.biocarta.v6.2.symbols.gmt" ,
			 "c2.cp.kegg.v6.2.symbols.gmt" ,
			 "c2.cp.reactome.v6.2.symbols.gmt" ,
			 "c2.cp.v6.2.symbols.gmt" ,
			 "c3.mir.v6.2.symbols.gmt" ,
			 "c3.tft.v6.2.symbols.gmt" ,
			 "c4.cgn.v6.2.symbols.gmt" ,
			 "c4.cm.v6.2.symbols.gmt" ,
			 "c5.bp.v6.2.symbols.gmt" ,
			 "c5.cc.v6.2.symbols.gmt" ,
			 "c5.mf.v6.2.symbols.gmt" ,
			 "c6.all.v6.2.symbols.gmt" ,
			 "c7.all.v6.2.symbols.gmt" ,
			 "h.all.v6.2.symbols.gmt")
		target = paste0('/home/huw/program/gmt/v6.2/', target)
		target
		file.exists(target)
		file.size(target)
	}else if(file_ext(gmtfile) == 'gmt'){
		target = gmtfile
	}else{## choose a single gene set
		#gmtfile = 'JAEGER_METASTASIS_UP'
		#tmpfile = tempfile(pattern = "tmp", tmpdir = tempdir(), fileext = ".gmt")
		tmpfile = normalizePath('tmp_file.gmt')
		all.gmt = '/home/huw/program/gmt/msigdb.v5.0.symbols.gmt'
		system(paste0('grep ', gmtfile, ' ', all.gmt, ' > ', tmpfile))
		if(file.size(tmpfile) < 10){
			return('gmt file has problem')
		}
		target = tmpfile
		#system(paste0('rm ', tmpfile))
	}

	if(!all(file.exists(target)) ){ ## if gmt file is not at current working dir
		xx = !file.exists(target)
		target[xx] = paste0("/home/huw/program/gmt/v6.2/", target[xx]) # look into this dir
	}
	if(!all(file.exists(target))){ # if still cannot find gmt file
		cat("gmt file can not find\n")
		cat(target[!file.exists(target)], sep="\t")
		cat("\n")
		return
	}
	targetname = basename(target)

	java="/home/huw/program/jdk7/bin/java"
	java = '/opt/common/CentOS_6/java/jdk1.8.0_31/bin/java'
	xmx = paste0(" -Xmx", mem, "g")
	cmd0 = paste0(java,  xmx, " -cp /home/huw/program/gsea-3.0.jar xtools.gsea.GseaPreranked")
	cmd0 = paste0(cmd0, " -collapse false -mode Max_probe -norm meandiv -nperm ", nperm, " -scoring_scheme classic -create_svgs ", if.svg)
	cmd0 = paste0(cmd0, " -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp ")
	cmd0 = paste0(cmd0, " -set_max 500 -set_min 15 -zip_report false -gui false")


	for(j in 1:length(rnkfile)){
		rnkfilej = rnkfile[j]
		cmdfile = sub(".rnk$", "_gsea.sh", rnkfilej);
		if( cmdfile == rnkfilej ) {
			cmdfile = paste0(rnkfilej, "_gsea.sh")
		}
		cat(NULL, file=cmdfile)
		outdir = file.path(normalizePath(dirname(rnkfilej)), "gsea"); outdir
		cmd1 = paste0(cmd0, " -out ",  outdir, " ")
		cmd1
		#source('~/program/configure.R')
		for(i in 1:length(target)){
			targeti = target[i]
			cmd2 = paste0(cmd1, " -gmx ", targeti)

			targetnamei = targetname[i]
			cmd2 = paste0(cmd2, " -rnk ", rnkfilej)
			label = basename(rnkfilej)
			cmd2 = paste0(cmd2, " -rpt_label gsea.", label, ".", targetnamei)
			jobname = paste0(basename(targeti), '.', basename(rnkfilej))
			if(if.cluster == T){
				cmd2 = paste0(bsub.head(jobname = jobname, hmem = mem, mem=mem, cpu=cpu), '"', cmd2, ' "')
			}
			cat("## now dealing with ", rnkfilej, "\t\t", targeti, "\n", file=cmdfile, append=T)
			cat(cmd2, "\n")
			cat(cmd2, "\n", file=cmdfile, append=T)
			if(if.run == T){
				system(cmd2)
			}
		}
		cat("command saved to ",  cmdfile, "\n")
	}
}
run_gsea = run.gsea.fun

heatmap.fun = function(mtx, ann, fname = NULL, hh = 6, ww = 5, ...) {
	ha = HeatmapAnnotation(df = ann)
	ch = Heatmap(mtx, top_annotation = ha, ...)
	#show_row_names=F, show_column_dend=T, show_row_dend=T, name='z-score', top_annotation=ha, column_split = ann$grp, km=km, cluster_columns = F, ...)
	if(!is.null(fname)){
		pdf(file=fname, width=ww, height=hh)
		draw(ch)
		dev.off()
		message(fname, ' saved!')
	}
	ch
}

get.geneset.genes = function(geneset=''){
	#geneset = 'HOLLERN_SQUAMOUS_BREAST_TUMOR'
	gs = system(paste0('grep ', geneset, ' /home/huw/program/gmt/v6.2/msigdb.v6.2.symbols.gmt'), intern=T)
	gs = unlist(strsplit(gs, '\t'))
	gs = gs[3:length(gs)]
	gs
}

gg.maplot.fun = function(mtx, ww = 6, hh = 4, fdr= 0.05, log2fc = 1, fname='', topn = 10, genes = '', 
			 labx = 'log2 Fold Change', if.lab=F, laby="-log10(adjusted p value)", max.overlaps = Inf, gtitle = NULL){
	require(ggplot2)
	require(ggrepel)
	require(dplyr)
	mtx = mtx[!is.na(padj),] 
	mtx[, yy := -log10(padj)]
	#mtx[, yy := ifelse(abs(yy) <= 1, 0, sign(yy)*log10(abs(yy)))]
	#mtx[, yy := sqrt(yy1)] 
	mtx[, Significant := 'NS']
	mtx[log2FoldChange > log2fc  & padj < fdr, Significant := 'UP']
	mtx[log2FoldChange < -log2fc & padj < fdr, Significant := 'DN']
	tt = mtx[F,]
	if(topn > 0){
		tt = mtx %>% filter(Significant %in% c('UP', 'DN')) %>% group_by(Significant) %>% dplyr::top_n(10, wt = -padj)
		tt
	}
	if(length(genes) > 0){
		tt1 = mtx[symbol %in% genes,]
		tt = rbind(tt, tt1)
		tt = as.data.table(tt)
		tt = tt[!duplicated(rn), ]
		tt
	}
	tt = as.data.table(tt)
	mtx[rn %in% tt$rn, Significant := 'Spot']
	tt[, Significant := 'Spot']
	myColors = c(NS = adjustcolor('black', .7), UP = adjustcolor('red', alpha.f = .7), DN = adjustcolor('blue', alpha.f=.7), Spot = adjustcolor('green', .9))
	myColors
	gg = ggplot() + 
		geom_point(data=mtx, aes(x = log2FoldChange, y = yy, color=Significant), size=1)  + 
		#geom_point(data=tt, aes(x = log2FoldChange, y = yy, color=Significant), size=1)  + 
		theme(legend.position="none")  + 
		scale_colour_manual(name = "Significance", values = myColors) + 
		scale_y_continuous(trans = "log1p") + 
		theme_bw(base_size = 8) + theme(legend.position='none') + ylab(laby) + xlab (labx) 
	if(if.lab == T){
		gg = gg + geom_text_repel(data=tt, aes(x = log2FoldChange, y = yy, label = symbol),  max.overlaps = Inf,
					   size=2, segment.color = 'grey50') 
	}
	if(!is.null(gtitle)){
		gg = gg + ggtitle(gtitle)
	}
	if(fname!= ''){
		ggsave(gg, file=fname, width=ww, height=hh)
		message(fname, ' saved!')
	}
	gg
}

gg.pca.fun = function(mtx, ann, topn = NULL, basemean = NULL, rn = NULL, fname = NULL, hh = 3, ww = 4){
	if(max(mtx) > 30){ mtx=log2(mtx+1) }
	#ann = anno.old
	if(!is.null(topn)){
		row.sd = apply(mtx, 1, sd)
		row.sd = row.sd[order(row.sd, decreasing=T)]
		sel.id = names(row.sd[1:topn])
		if(!is.null(basemean)){
			row.mn = rowSums(mtx)
			sel.id = names(row.mn[row.mn > basemean & names(row.mn) %in% sel.id])
		}
		mtx = mtx[row.names(mtx) %in% sel.id, ]
	}
	pca = prcomp(t(mtx), scale. = TRUE)
	pca.sum = summary(pca)
	pca = pca$x[, 1:3]
	pca = as.data.table(as.data.frame(pca), keep.rownames = T)

	pca.sum.key = pca.sum$importance
	xlab=paste0('PC1 ', round(100*pca.sum.key[2, 1], 0), '%') 
	ylab=paste0('PC2 ', round(100*pca.sum.key[2, 2], 0), '%') 

	source('~/pipeline/work/g10_all_.r')
	pca = merge.df.fun(pca, ann, if.returndt = T)

	gg = ggpubr::ggscatter(pca, x='PC1', y='PC2', color=grp, xlab = xlab, ylab = ylab)

	if(!is.null(fname)){
		ggsave(gg, file = fname, width = ww, height = hh )
		message(fname, ' saved')
	}
	gg
}

oncoplot.topn.gene.fun = function(maf, n.topn.gene = 20){
	if(class(maf) == 'MAF'){maf = maf@data}
	tbl.dt = as.data.table(as.data.frame(table(maf$Hugo_Symbol, maf$simple.onco)))[Freq > 0, ]
	setnames(tbl.dt, c('Hugo_Symbol', 'simple.onco', 'Freq'))
	topn.gene = as.character(tbl.dt[order(simple.onco, -Freq), ][1:ifelse(nrow(tbl.dt) < n.topn.gene, nrow(tbl.dt), n.topn.gene), Hugo_Symbol])
	topn.gene
}

named.vector.fun = function(dt.dt, name.c = 'simple.vc.onco', val.c = 'color'){
	if(length(colnames(dt.dt))< 2){
		stop('not enough columns in the data frame provided')
	}
	if(val.c %in% colnames(dt.dt)){
		tmp = unlist(dt.dt[, val.c, with=F])
	}else{
		tmp = unlist(dt.dt[, 2])
	}
	if(name.c %in% colnames(dt.dt)){
		nn = unlist(dt.dt[, name.c, with=F])
	}else{
		nn = unlist(dt.dt[, 1])
	}
	names(tmp) = nn
	tmp
}

list.len.max = function(l1, l2){
        oo = c();
        for(i in 1:length(l1)){
                oo = c(oo, max(length(l1[[i]]), length(l2[[i]])));
        }
        return(oo)
}
list.intersect = function(l1, l2){
        oo = list();
        for(i in 1:length(l1)){
                oo[[i]] = intersect(l1[[i]], l2[[i]]);
        }
        return(oo)
}
list.len.min = function(l1, l2){
        oo = c();
        for(i in 1:length(l1)){
                oo = c(oo, min(length(l1[[i]]), length(l2[[i]])));
        }
        return(oo)
}
list.intersect.len = function(l1, l2){
        oo = c();
        for(i in 1:length(l1)){
                oo = c(oo, length(intersect(l1[[i]], l2[[i]])));
        }
        return(oo)
}

list.setdiff.len = function(l1, l2){
        oo = c();
        for(i in 1:length(l1)){
                oo = c(oo, length(intersect(l1[[i]], l2[[i]])));
        }
        return(oo)
}

list.setdiff = function(l1, l2){
        oo = list();
        for(i in 1:length(l1)){
		nn = names(l1)[i];
		paste0(nn, "_sub") -> nn
                oo[[nn]] = intersect(l1[[i]], l2[[i]]);
        }
        return(oo)
}

## write to rnk file for gsea
## xx is dataframe(genesymbol, log2FoldChange, pvalue )
write_rnk = function(data = xx, fname = NULL, pval.cutoff=1, neg = F){
	tmp = data;
	if(neg){tmp[, log2FoldChange := -log2FoldChange]}
	if(!('symbol' %in% colnames(tmp))){
		if('gene' %in% colnames(tmp)){
			tmp$symbol = tmp$gene
		}else{
			tmp$symbol = sub(".*_", "", row.names(tmp))
		}
	}

	if(!('log2FoldChange' %in% colnames(tmp)) & 'avg_log2FC' %in% colnames(tmp)){
		tmp[, log2FoldChange := avg_log2FC]
	}

	if(!('padj' %in% colnames(tmp)) & 'p_val_adj' %in% colnames(tmp)){
		tmp[,  padj := p_val_adj]
	}

	if(!('pvalue' %in% colnames(tmp)) & 'p_val' %in% colnames(tmp)){
		tmp[,   pvalue := p_val]
	}

	tmp = tmp[tmp$symbol != '' & !is.na(tmp$log2FoldChange) & !is.na(tmp$pvalue) & tmp$pvalue < pval.cutoff, ] 
	tmp$gsea = sign(tmp$log2FoldChange) * (- log10(tmp$pvalue))
	tmp = tmp[order(tmp$pvalue),]
	#strsplit(as.character(tmp[,1]), "///") -> tt;
	#lapply(tt, length) -> tt.len;
	#unlist(tt.len) -> tt.len;
	#Map(function(x,y){rep(x,each=y)}, tmp[,2], tt.len) -> xx;
	#unlist(xx) -> xx;
	#data.frame(tt = unlist(tt), xx) -> tmp2;
	duplicated(tmp[,"symbol"]) -> dup;
	tmp2 = tmp[!dup,];
	if(is.infinite(tmp2$gsea[1])){
		top = tmp2$gsea[is.infinite(tmp2$gsea)]
		mm = max(abs(tmp2$gsea[!is.infinite(tmp2$gsea)])) +1
		top_seq = seq(from = mm, by = 1, length.out = length(top))
		top_seq = top_seq * sign(tmp$log2FoldChange[length(top):1])
		tmp2$gsea[is.infinite(tmp2$gsea)] = top_seq
	}
	tmp2
	write.table(tmp2[,c("symbol", "gsea")], file= fname, quote=F, sep="\t", row.names=F, col.names=F);
}

bed_to_granges <- function(file){
   df <- read.table(file,
                    header=F,
                    stringsAsFactors=F)
 
   if(length(df) > 6){
      df <- df[,-c(7:length(df))]
   }
 
   if(length(df)<3){
      stop("File has less than 3 columns")
   }
 
   header <- c('chr','start','end','id','score','strand')
   names(df) <- header[1:length(names(df))]
 
   if('strand' %in% colnames(df)){
      df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
   }
 
   library("GenomicRanges")
 
   if(length(df)==3){
      gr <- with(df, GRanges(chr, IRanges(start, end)))
   } else if (length(df)==4){
      gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
   } else if (length(df)==5){
      gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
   } else if (length(df)==6){
      gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
   }
   return(gr)
}

### load gmt file as a list
load_as_list = function(x){

	x <- scan(x, what="", sep="\n");
	y <- strsplit(x, "[[:space:]]+");
	names(y) <- sapply(y, `[[`, 1);
	y <- lapply(y, `[`, -1);
	y <- lapply(y, `[`, -1);
	return(y);
}

generate_breaks = function(x){
        rnge = range(x)
        t1 = seq(from = -1, to = 1, length.out = 8)
        t2 = seq(from = -2, to = -1, length.out = 8)
        t3 = seq(from = 1, to = 2, length.out = 8)
        t4 = seq(from = -10, to = -2, length.out = 4)
        t5 = seq(from = 2, to = 10, length.out = 4)
        tt = unique(c(t1, t2, t3, t4, t5, rnge))
        tt = tt[order(tt)] 
        tt[tt >= floor(min(rnge)) & tt <= ceiling(max(rnge))]     
}

generate_z_score = function(x) {
	sweep(x, 1, apply(x, 1, median, na.rm=T))
}

cbioMaf = function(hmaf, if.all = F){
	require(data.table)
	#source('~/pipeline/fun/param.r')

	maf.ori = copy(hmaf)
	maf.ori[, simple.vc := Variant_Classification ]
	maf.ori[Variant_Classification == "5'Flank", simple.vc := 'Promoter Mutation']
	maf.ori[Variant_Classification %in% c("Frame_Shift_Del", "Frame_Shift_Ins"), simple.vc := 'Truncating Mutation']
	maf.ori[Variant_Classification %in% c("In_Frame_Del", "In_Frame_Ins"), simple.vc := 'Inframe']
	maf.ori[Variant_Classification %in% c("Missense_Mutation"), simple.vc := 'Missense Mutation']
	maf.ori[Variant_Classification %in% c("Nonsense_Mutation"), simple.vc := 'Truncating Mutation']
	maf.ori[Variant_Classification %in% c("Splice_Region", "Splice_Site"), simple.vc := 'Splice Mutation']
	maf.ori[Variant_Classification %in% c("3'Flank", "RNA", "Intron", "IGR", "5'UTR", "3'UTR", "Silent"), simple.vc := 'ex']
	maf.ori = maf.ori[simple.vc != 'ex', ]
	maf.ori

	unique(maf.ori$Tumor_Sample_Barcode)

	maf.ori[oncogenic %in% onco.sel, simple.onco := 'putative driver']
	maf.ori[!(oncogenic %in% onco.sel), simple.onco := 'unknown significance']
	maf.ori[, simple.vc.onco := paste0(simple.vc, '(', simple.onco, ')')]
	#maf.ori[!duplicated(paste(Variant_Classification, simple.onco)), .(Tumor_Sample_Barcode, Hugo_Symbol, HGVSp_Short, simple.vc, Variant_Classification, Variant_Type, oncogenic, simple.onco, simple.vc.onco)]

	maf.ori[, Variant_Classification_old := Variant_Classification]
	maf.ori[, Variant_Classification := simple.vc.onco]

	ll = simple.vc[order(simple.vc$level), simple.vc.onco]
	maf.ori[, simple.vc.onco := factor(simple.vc.onco, levels = ll)]
	maf.ori

}
