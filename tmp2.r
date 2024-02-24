
## tempo wes
{

	cmoid = fread('res/CMOID_v2.xls.txt', sep="\t")
	cmoid[histology != 'v', HIS := histology]
	cmoid

	sdt = fread('CCS_YMJ4BKK8/sample_data.txt')
	sdt = sdt[, .(sample, purity, ploidy, MSIscore, Number_of_Mutations, TMB)]
	sdt[, Tumor_Sample_Barcode := sub('__.*', '', sample)]
	sdt[, sname := Tumor_Sample_Barcode ]
	sdt

	sdt = merge(sdt, cmoid, by.x='sname', by.y = 'CMOID')
	sdt[HIS == 'NOS', Type := 'NOS']
	sdt[HIS != 'NOS', Type := 'NXS']
	sdt[, patient := substr(sname, 5, 10)]
	sdt[, patient.type := paste0(patient, '_', Type)]
	sdt[, sample := substr(sname, 12, 15)]
	sdt[, HIS := relevel(factor(HIS), ref='NOS')]
	sdt[, Type := relevel(factor(Type), ref='NOS')]
	wes.info = as.data.table(sdt)

	wes.dir = 'CCS_YMJ4BKK8'
	maf.files = system(paste0('find -L ', wes.dir, ' -name "*somatic.final.maf"'), intern = T)
	maf.files
	wes.maf = lapply(maf.files, fread)
	wes.maf = rbindlist(wes.maf)
	wes.maf[grep('SMAR', Hugo_Symbol), cc.sel, with=F]
	wes.maf[is.na(oncogenic), oncogenic := 'Unknown']
	wes.maf[oncogenic == '', oncogenic := 'Unknown']
	table(wes.maf$oncogenic)
	head(wes.maf$Tumor_Sample_Barcode)
	wes.maf[Hugo_Symbol == 'ERBB2', .(Hugo_Symbol, HGVSp_Short, Tumor_Sample_Barcode, oncogenic)]
	fname = 'wes_maf.txt'
	fwrite(wes.maf, file=fname,sep="\t")


	## cnTable
	{
		#cnTable Custom copynumber data if gistic results are not available. Input file or a data.frame should contain three columns in aforementioned order with gene name, Sample name and copy number status (either 'Amp' or 'Del'). Default NULL.

		cnv.dt = fread(paste0(wes.dir, '/cna_genelevel.txt'))
		gg =  oncoplot.topn.gene.fun(mm.o)
		cnv.dt = cnv.dt[gene %in% gg, ]
		cnv.dt = cnv.dt[, sample.name := sub('__.*', '', sample) ]
		cnv.dt[tcn >3, cn := 'AMP']
		cnv.dt[tcn <2, cn := 'DEL']
		cnv.dt = cnv.dt[!is.na(cn), .(gene, sample.name, cn)]
		setnames(cnv.dt, c('gene name' ,'sample name', 'cn status'))
		cnv.dt

	}

	fwrite(wes.maf, file='wes_maf.txt', sep="\t")
	wes.maf = fread('wes_maf.txt')

	source('~/pipeline/fun/cbioMaf.R')

	wes.maf.x = cbioMaf(wes.maf)
	wes.maf.o = read.maf(wes.maf.x, vc_nonSyn = simple.vc$simple.vc.onco, cnTable = cnv.dt)
	fname = 'oncoplot_mcp_cnv.pdf'
	pdf(file = fname, width=8, height=7)
	gg =  oncoplot.topn.gene.fun(mm.o)
	oncoplot(wes.maf.o, genes = gg, annotationFontSize = 1, anno_height = 1, 
		 clinicalFeatures= c('Type', 'HIS'), annotationDat = wes.info, sortByAnnotation=T,  writeMatrix = T,
		 colors = named.vector.fun(simple.vc), removeNonMutated=F, showTumorSampleBarcodes = F, drawColBar = F)
	dev.off()

	head(res.cc)
	wes.maf$

 	plotname = colnames(read.csv('onco_matrix.txt', sep="\t"))
 	plotname
 	setkey(wes.info, 'Tumor_Sample_Barcode')
 	tmp = wes.info[plotname, .(Tumor_Sample_Barcode, TMB, Type)]
	tmp[, Type := factor(Type, levels = c('NOS', 'NXS'))]
 	tmp[TMB > 50, TMB := 50]
 	fname = 'oncoplot_tmb.pdf'
 	gg = ggpubr::ggbarplot(tmp, x = 'Tumor_Sample_Barcode', y = 'TMB', lab.size=0, fill = 'Type') + ggpubr::rotate_x_text() + ggplot2::ylim(0, 50)
 	ggplot2::ggsave(gg, file=fname, width=12, height=3.5)

	gg = ggpubr::ggdotplot(wes.info, x='Type', y='TMB', add='mean_sd',  error.plot = "crossbar", add.params = list(width = 0.5), color = 'Type') + ggpubr::stat_compare_means(method = "t.test")
	fname = 'tmb_bar_v2.pdf'
	ggplot2::ggsave(gg, file = fname, width=2.5, height=4)


	tmp = wes.maf[Hugo_Symbol %in% impct$symbol, ]
	tmp.c = cbioMaf(tmp)
	tmp.o = read.maf(tmp.c$maf, vc_nonSyn = tmp.c$vc.key) 
	fname = 'res/oncoplot_tempo_wes_all.pdf'
	pdf(file = fname, width=9, height=8)
	oncoplot(tmp.o, colors = tmp.c$color.table, removeNonMutated = T, showTumorSampleBarcode = T, barcode_mar = 10, annotationDat = wes.info, clinicalFeatures = c('Type', 'HIS'), sortByAnnotation=T)
	dev.off()

	tmp = wes.maf[Hugo_Symbol %in% impct$symbol, ]
	tmp = merge(tmp, sdt[, .(sname, Type)], by.x = 'Tumor_Sample_Barcode', by.y = 'sname', all.x=T)
	tmp[, patient.type := paste0(patient,'_', Type)]
	tmp[, Tumor_Sample_Barcode := sub('C_', '', patient.type)]
	tmp[, tag := paste0(patient, '_', Type, '_', Start_Position)]
	tmp = tmp[!duplicated(tag), ]
	tmp.c = cbioMaf(tmp)
	tmp.o = read.maf(tmp.c$maf, vc_nonSyn = tmp.c$vc.key) 
	fname = 'res/oncoplot_tempo_wes_patient.pdf'
	sdt[, patient.type := paste0(patient, '_', Type)]
	tmp.info = sdt[!duplicated(patient.type), ]
	tmp.info[, Tumor_Sample_Barcode := patient.type]
	tmp.info
	pdf(file = fname, width=6, height=8)
	oncoplot(tmp.o, colors = tmp.c$color.table, removeNonMutated = T, showTumorSampleBarcode = T, barcode_mar = 10, annotationDat = tmp.info, clinicalFeatures = c('Type', 'HIS'), sortByAnnotation=T)
	dev.off()

	tmp = wes.maf[oncogenic != '', ][!is.na(oncogenic), ]
	table(tmp$oncogenic)
	(tmp$oncogenic)
	tmp.c = cbioMaf(tmp)
	tmp.o = read.maf(tmp.c$maf, vc_nonSyn = tmp.c$vc.key) 
	fname = 'oncoplot_oncogenic.pdf'
	pdf(file = fname, width=6, height=80)
	oncoplot(tmp.o, colors = tmp.c$color.table, removeNonMutated = T, showTumorSampleBarcode = T, barcode_mar = 10)
	dev.off()

	## IMPACT505
	impct = fread('/work/ci/resources/roslin_resources/targets/IMPACT505/b37/IMPACT505_b37_targets.bed')
	impct = impct[grep('Tiling', V5, invert=T), ]
	impct[, symbol := sub('_target.*', '', V5)]
	impct = impct[!duplicated(symbol), ]
	impct = impct[grep('snp', symbol, invert=T), ]
	impct[, symbol := sub('_.*', '', symbol)]
	impct = impct[!duplicated(symbol), ]
	impct

	tmp = wes.maf[Hugo_Symbol %in% c(impct$symbol), ]#[oncogenic != 'Unknown', ]
	tmp.c = cbioMaf(tmp)
	tmp.o = read.maf(tmp.c$maf, vc_nonSyn = tmp.c$vc.key) 
	fname = 'mcp_oncoplot_tempo_wes_impact505_top20.pdf'
	pdf(file = fname, width=6, height=8)
	gg = oncoplot(tmp.o, colors = tmp.c$color.table, top = 20, removeNonMutated = F, showTumorSampleBarcode = T, barcode_mar = 10, writeMatrix = T)
	dev.off()

	## use this one
	cca.sel = c('ARID1A', 'CDH1', 'EP300', 'KMT2C', 'SMARCA4', 'TP53', 'KMT2B', 'KRAS', 'MGA', 'NF1', 'NF2', 'PBRM1', 'PPP2R1A', 'SMAD4', 'BRCA1', 'FBXW7', 'LZTR1', 'NTHL1', 'PIK3CA', 'PIK3R1', 'SETD2', 'STK11')
	tmp = wes.maf[Tumor_Sample_Barcode %in% info[cca == T, sname], ]
	tmp.c = cbioMaf(tmp)
	tmp.o = read.maf(tmp.c$maf, vc_nonSyn = tmp.c$vc.key) 
	fname = 'oncoplot_cca_sel_gene.pdf'
	pdf(file = fname, width=5.5, height=8)
	gg = oncoplot(tmp.o, genes = cca.sel, colors = tmp.c$color.table, removeNonMutated = F, showTumorSampleBarcode = T, barcode_mar = 10, writeMatrix = T)
	dev.off()

	wes.maf[Tumor_Sample_Barcode == 's_C_JJ9RNX_P003_d', cc.sel, with=F][oncogenic != 'Unknown', ]
	wes.maf[Tumor_Sample_Barcode == 's_C_JJ9RNX_P004_d', cc.sel, with=F][oncogenic != 'Unknown', ]
	wes.maf[Tumor_Sample_Barcode == 's_C_JJ9RNX_P005_d', cc.sel, with=F][oncogenic != 'Unknown', ]


	wes.maf[grep('SMAR', Hugo_Symbol), cc.sel, with=F]
	wes.maf[grep('7MM', Tumor_Sample_Barcode), cc.sel, with=F][HGVSp_Short != '', ][Variant_Classification == 'Nonsense_Mutation', ]
	(wes.maf[grep('7MM', Tumor_Sample_Barcode), cc.sel, with=F][HGVSp_Short != '', ][Hugo_Symbol %in% cca.gg, ])
	wes.maf[grep('7MM', Tumor_Sample_Barcode), cc.sel, with=F][Hugo_Symbol %in% cca.gg]
	wes.maf[grep('JJ', Tumor_Sample_Barcode), cc.sel, with=F][Hugo_Symbol %in% cca.gg]
	wes.maf[oncogenic != 'Unknown', cc.sel, with=F][grep('7MM', Tumor_Sample_Barcode), ]

	dd = read.table('onco_matrix.txt', sep="\t")
	sorder = colnames(dd)
	sorder

	cca.g0 = rownames(dd)
	## Whole-genome sequencing revealed novel prognostic biomarkers and promising targets for therapy of ovarian clear cell carcinoma
	cca.g1 = c('NBPF20', 'NBPF10', 'NBPF14', 'NEFH', 'ARID1A', 'PIK3CA', 'ARID1B', 'PRR23A', 'HYDIN', 'USP17L5', 'HLA-DRB1', 'KRTAP4-5', 'WRNIP1', 
		   'HLA-A', 'DSPP', 'ATXN1', 'TTN', 'GPRIN2', 'FOXE1', 'HGC6.3', 'CREG2', 'MUC2', 'FKRP', 'COL18A1', 'PUS7', 'TBC1D3', 'MESP2', 'ABCA1', 
		   'TBP', 'ATN1', 'SRRD', 'MUC12', 'MUC4', 'ZNF717', 'NBPF11', 'C2orf81', 'GIGYF2')
	## genomics of endometrial clear cell carcinoma
	cca.g2 = c('TP53', 'PPP2R1A', 'PIK3CA', 'FBXW7', 'ARID1A', 'SPOP', 'PIK3R1', 'KRAS', 'HNF1A', 'KMT2D', 'HLA-A', 'KMT2C', 'ATM', 'SMARCA4', 'KEAP1', 
		   'MAP3K1', 'BCOR1', 'AMER1', 'PMS2', 'APC', 'NF1', 'EPHA5', 'ERBB2', 'MPL', 'JAK2', 'FAT4', 'PTEN', 'DICER1', 'ERBB3', 'CCNE1', 'CEBPA', 
		   'DAXX', 'ERBB2', 'MEF2B', 'WT1', 'MCL1', 'NOTCH2', 'MET', 'AKT2')
	## Analysis of gene expression signatures identifies prognostic and functionally distinct ovarian clear cell carcinoma subtypes
	cca.g3 = c('ARID1A', 'ARID1B', 'SMARCA4', 'KMT2A', 'KMT2B', 'KMT2C', 'KMT2D', 'KRAS', 'PPP2R1A', 'MET', 'PIK3CA', 'PIK3R1', 'PTEN', 'AKT1', 'AKT2', 
		  'ERBB2', 'MYC', 'ZNF217', 'CCNE1', 'CDKN2A', 'CDKN2B', 'IL6', 'TP53', 'NF1', 'CDK12', 'RB1', 'BRCA1', 'BRCA2')
	intersect(cca.g0, cca.g1)
	intersect(cca.g0, cca.g2)
	intersect(cca.g0, cca.g3)
	intersect(cca.g1, cca.g2)
	intersect(cca.g1, cca.g3)
	intersect(cca.g1, cca.g3)
	wes.maf[Hugo_Symbol %in% c(cca.g0, cca.g1, cca.g2, cca.g3), .(Tumor_Sample_Barcode, Hugo_Symbol, oncogenic, HGVSp_Short)][order(Tumor_Sample_Barcode), ][grep('FP3', Tumor_Sample_Barcode), ]
	wes.maf[grep('SMAR', Hugo_Symbol), cc.sel, with=F]


	## mutation burden
	dim(dd)
	intersect(colnames(dd), sdt$Tumor_Sample_Barcode)
	sdt
	tmp = sdt[cca == T, ]
	tmp[, Tumor_Sample_Barcode := factor(Tumor_Sample_Barcode, levels=colnames(dd))]
	gg.tmb = ggplot(tmp, aes(x = Tumor_Sample_Barcode,  TMB)) + geom_bar(stat='identity') + theme(axis.text.x=element_text(angle=90, vjust = 0.5, hjust = 1))
	fname = 'tmb.pdf'
	ggsave(gg.tmb, file = fname, width=6, height=3)

	## copy number var
	{

		cnv.dt = fread(paste0(wes.dir, '/cna_genelevel.txt'))
		cnv.dt = cnv.dt[spans_segs == T, ][gene %in% impct$symbol, ][gene %in% c(rownames(dd), mut.2017), ]
		cnv.dt = cnv.dt[cn_state != 'INDETERMINATE', ]
		cnv.dt = cnv.dt[cn_state != 'DIPLOID or CNLOH', ]
		cnv.dt[grep('^CNLOH', cn_state), ]
		cnv.dt[, sname := sub('.*__', '', sample)]
		cnv.dt[, gene := factor(gene)]
		cnv.dt[, sname := factor(sname)]
		cnv.dt[, cnv := cn_state]
		cnv.dt[grep('LOSS', cn_state), cnv := 'LOSS']
		cnv.dt[grep('GAIN', cn_state), cnv := 'GAIN']
		cnv.dt[grep('AMP', cn_state), cnv := 'GAIN']
		cnv.dt[grep('TETRA', cn_state), cnv := 'TETRAPLOID']
		cnv.dt[grep('HOMDEL', cn_state), cnv := 'LOSS']
		cnv.dt[, gene.nn := .N, by = 'gene']
		cnv.dt[, cnv.nn := .N, by = 'sname']
		cnv.dt[, gene := factor(gene, levels = cnv.dt[order(gene.nn, decreasing = F), ][!duplicated(gene), gene])]
		#cnv.dt[, sname := factor(sname, levels = colnames(dd))]
		cnv.dt[, gene.i :=  as.numeric(gene)]
		cnv.dt[, sname.i :=  as.numeric(sname)]
		cnv.dt[, ll := 0]
		cnv.dt[, bb := 0]
		cnv.dt[, x1 := sname.i]
		cnv.dt[, x2 := x1 + 0.8]
		cnv.dt[, y1 := gene.i]
		cnv.dt[, y2 := y1 + 0.8]
		cnv.dt

		gg = ggplot(cnv.dt, aes(sname, gene, fill = cnv)) + geom_tile() + theme(axis.text.x = element_text(angle=90))
		fname = 'cnv.pdf'
		ggsave(gg, file=fname, width=8, height=6)

		gg = ggplot(cnv.dt, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, fill = cnv)) + geom_rect() + theme_void() +
			geom_text(aes(ll, gene.i, label=gene)) +
			geom_text(aes(sname.i, bb, label=sname), angle=90) 
		fname = 'cnv.pdf'
		ggsave(gg, file=fname, width=8, height=6)



		cnv.dt = fread(paste0(wes.dir, '/cna_facets_run_info.txt'))
		cnv.dt[grep('hisens', Sample), .(Sample, Ploidy, Purity)]


	}

	## mutation signature
	{

		s.dt = fread(paste0(wes.dir, '/sample_data.txt'))
		s.dt$sample
		cn = c('sample', grep('SBS', names(s.dt), value=T)); cn
		tmp = s.dt[, ..cn]
		tmp = melt(tmp)
		tmp[grep('observed', variable), param := 'sval']
		tmp[grep('pvalue', variable), param := 'pval']
		tmp[, sig.name := sub('\\..*', '', variable)]
		tmp = dcast(tmp, sample + sig.name ~ param)
		tmp = tmp[pval < 0.05 & sval > 0.10, ]
		tmp[, sname := sub('__.*', '', sample)]
		tmp = tmp[sname %in% sorder, ]
		tmp
		hh = tmp[!duplicated(sname),  ]
		hh[, sig.name := 'Others']
		hh[, pval := 0]
		hh[, sval := 0]
		sigs = rbind(hh, tmp)
		sigs = sigs[!(sig.name %in% c('SBS43', 'SBS49', 'SBS57')), ] ## possible sequencing artefact
		for(i in 1:nrow(hh)){
			i.val = 1 - sum(sigs[sig.name != 'Others' & sname == hh$sname[i], sval])
			sigs[sig.name == 'Others' & sname == hh$sname[i], sval := i.val]
		}
		sigs[sig.name == 'Others', ]
		sigs = merge(sigs, sig.table, by = 'sig.name')
		hh = sigs[!duplicated(sig.name), ]
		cc = hh$sig.col; names(cc) = hh$tag
		#sigs[, sname := factor(sname, levels = colnames(dd))]
		sigs[sig.name != 'Others', ll := sig.name]
		sig.plot = ggplot(data=sigs[order(sig.name)], aes(x=sname, y=sval, fill=tag, order=sig.name)) 
		sig.plot = sig.plot + geom_bar(stat = 'identity', position='stack') 
		sig.plot = sig.plot + scale_fill_manual(values=cc, name='Mutational Signature') 
		#sig.plot = sig.plot + scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), labels=c('0', '', '0.4', '', '0.8', '')) 
		sig.plot = sig.plot + labs(x='', y='Fraction of Mutations') 
		sig.plot = sig.plot + geom_text(aes(label = sig.name), position = position_stack(), vjust = 0, size = 2)
		#sig.plot = sig.plot + facet_grid(.id ~ patient.name, space="free_x", scales="free_x", drop = TRUE) 
		sig.plot = sig.plot + theme_classic(20) 
		sig.plot = sig.plot + theme(axis.title=element_text(size=8, face = "bold"), 
					    axis.text.x=element_text(size=9, angle=90, face = "bold", vjust=.5, hjust=1),
					    axis.text.y=element_text(size=9),
					    legend.text=element_text(size=9),
					    legend.title=element_text(size=9, face = "bold"),
					    axis.ticks.x=element_blank(),
					    axis.ticks.y=element_blank())
		sig.plot = sig.plot + theme(strip.text.x = element_text(angle = 90, size=9))
		#sig.plot = sig.plot + geom_text(data = nn, y = -0.1, aes(x=sample.name, label=N), inherit.aes=F, size=3)

		fname = 'res/signaturev2.pdf'
		ggsave(sig.plot, file=fname, width=15, height=6)


	}

}

## autospy
{

	mcp.maf.fil = wes.maf
	mcp.maf.fil[, IMPACT_468 := Hugo_Symbol %in% impct$symbol]
	mcp.maf.fil[, Variant_Bioportal := paste0(Variant_Classification, oncogenic)]
	mcp.maf.fil

	unique(mcp.maf.fil$Tumor_Sample_Barcode)
	pid = 'C_[[:alnum:]]+'
	stringr::str_extract(mcp.maf.fil$Tumor_Sample_Barcode, pid)
	stringr::str_extract(mcp.maf.fil$Tumor_Sample_Barcode, sid)

	sid = '_[P|M][[:alnum:]]+'
	mcp.maf.fil[, patient := stringr::str_extract(Tumor_Sample_Barcode, pid)]
	mcp.maf.fil[, sample := stringr::str_extract(Tumor_Sample_Barcode, sid)]
	mcp.maf.fil[, .(patient, sample)]
	dim(mcp.maf.fil)

	mcp.maf.fil[, t_var_freq := as.numeric(t_alt_count) / t_depth]
	mcp.maf.fil[, tm := paste0(Hugo_Symbol, " ", as.numeric(gsub("[^\\d]+", "", HGVSp_Short, perl=T)))] # 
	mcp.maf.fil[, TAG := paste0('chr', Chromosome, ':', Start_Position, '-', End_Position, ':', Reference_Allele, ':', Tumor_Seq_Allele2)]
	mcp.maf.fil[, variant := paste0(TAG, '::', Hugo_Symbol, ':', HGVSp_Short)]

	## tempo use new mcp.maf
	getwd()
	wes.dir ='res/autospy_tempo'
	dir.create(wes.dir)
	setwd(wes.dir)

	## 
	source('~/program/autospy/R/process_autopsy_maf.R')
	source('~/program/autospy/R/stratton_plot.R')

	### refilter
	filter_results <- filter_maf_report(mcp.maf.fil)
	tmp <- filter_results$maf
	nrow(tmp)

	write.table(tmp, 'filtered.maf', sep="\t", quote=F, row.names=F)

	### overlap plot
	make_mutation_overlap_plot(tmp, pid = pid, log = F, out = "pre_filter")

	## facets summary
	jc = '/ifs/work/solitlab/jessica/projects/squamous/analysis/tempo/Result_squamous_1/somatic/'
	arm = system(paste0('find ', jc, ' -name "*_hisens.cncf.txt"'), intern=T)
	arm.dt = lapply(arm, fread)
	arm.dt = rbindlist(arm.dt)
	arm.dt[, sname := sub('__.*', '', ID)]

	chr.hg19 = fread(cfg$genomeHg19chromsize)
	chr.hg19 = chr.hg19[!grep('_', V1), ]
	chr.hg19 = chr.hg19[!(V1 %in% c('chrY', 'chrM')), ]
	chr.hg19[V1 == 'chrX', V1 := 'chr23']
	chr.hg19[, chr := as.numeric(sub('chr', '', V1))]
	chr.hg19[, len := as.numeric(V2)]
	chr.hg19 = chr.hg19[order(chr), ]
	chr.hg19$cum = cumsum(c(0, chr.hg19$len[1:22]))
	chr.hg19

	arm.dt = merge(arm.dt, chr.hg19[, .(chr, len, cum)], by.x='chrom', by.y='chr', all.x=T)
	arm.dt

	ll = unique(arm.dt[order(sname), sname]); ll
	arm.dt[, pid := factor(sname, levels=ll)]
	arm.dt = arm.dt[tcn.em != 2, ]
	arm.dt

	arm.dt[, xmin := loc.start + cum]
	arm.dt[, xmax := loc.end   + cum]
	arm.dt[, ymin := as.numeric(pid)]
	arm.dt[, ymax := ymin + 0.8]

	bb = data.table(pid = factor(ll, levels=ll))
	bb[, xmin := 0]
	bb[, xmax := sum(chr.hg19$len)]
	bb[, ymin := as.numeric(pid)]
	bb[, ymax := ymin + 0.8]
	bb

	arm.dt[, tcn.col := tcn.em]
	arm.dt[tcn.em >= 5, tcn.col := 5]
	arm.dt[, tcn.col := as.character(tcn.col)]
	arm.dt$tcn.col

	# tempo ploidy
	pld = system(paste0('find ', jc, ' -name "*_OUT.txt" '), intern=T)
	pld.dt = lapply(pld, fread)
	pld.dt = rbindlist(pld.dt)
	pld.dt = pld.dt[grep('hisens', Sample), ]
	pld.dt[, pid := sub('__.*', '', Sample)]
	pld.dt[, pid := factor(pid, levels=ll)]
	pld.dt[, y := as.numeric(pid)]
	pld.dt[, x := sum(chr.hg19$len)]
	pld.dt

	limx=c(-0.2*max(bb$xmax), max(bb$xmax))
	bb.col = adjustcolor('black', alpha.f=.3)
	tcn.col.val = c('0'=adjustcolor('red', .4), '1'=adjustcolor('red', .2), '3' = adjustcolor('blue', .3),  '4' = adjustcolor('blue', .5), '5' = adjustcolor('blue', .8))
	gg = ggplot() + xlim(limx) +
		geom_rect(data=bb, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=bb.col) +
		geom_rect(data=arm.dt, aes(xmin = xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=tcn.col)) +
		scale_fill_manual(values=tcn.col.val) +
		geom_text(data=bb, aes(x=xmin, y=ymin+.4, label=pid), vjust=.5, hjust=1.05) +
		theme_void()  + geom_vline(xintercept = chr.hg19$cum, color=adjustcolor('#7d0229', alpha.f=.5)) +
		geom_text(data=chr.hg19, aes(y=0, x=0.5*len + cum, label=chr), hjust=.5) + 
		geom_text(data=pld.dt, aes(y=y, x=x, label=Ploidy), hjust=-0.2, vjust=-0.3, size=3)  
	fname = 'cnv_genome.pdf'
	ggsave(gg, file=fname, width=12, height=8)
	#ggsave(gg, file=fname, width=12, height=8,dpi=300)


	## DoFacets reanalysis
	cmd = paste0('/opt/local/singularity/3.3.0/bin/singularity exec /juno/work/taylorlab/cmopipeline/singularity_images/cmopipeline-temposig-0.2.2.img ')
	cmd = paste0(cmd, ' Rscript /usr/bin/facets-suite/run-facets-wrapper.R ')
	cmd = paste0(cmd, ' --cval 100 ')
	cmd = paste0(cmd, ' --snp-window-size 250 ')
	cmd = paste0(cmd, ' --normal-depth 35 ')
	cmd = paste0(cmd, ' --min-nhet 25 ')
	cmd = paste0(cmd, ' --purity-cval 500 ')
	cmd = paste0(cmd, ' --purity-min-nhet 25 ')
	cmd = paste0(cmd, ' --genome hg19 ')
	cmd = paste0(cmd, ' --counts-file ', )
	cmd = paste0(cmd, ' --sample-id ${tag} ')
	cmd = paste0(cmd, ' --directory ${outputDir} ')
	cmd = paste0(cmd, ' --facets-lib-path /usr/local/lib/R/site-library ')
	cmd = paste0(cmd, ' --seed 100 ')
	cmd = paste0(cmd, ' --everything ')
	cmd = paste0(cmd, ' --legacy-output T')


	sig.table

	cmd = 'python ~/program/autospy/inst/mutation-signatures/main.py --seed 100 ~/program/autospy/inst/mutation-signatures/Stratton_signatures29.txt filtered.maf signature.out'
	system(cmd)
	fread("signature.out") -> mut.sig
	fwrite(filter_results$report, file='report.txt', sep="\t", quote=F)
	ggsave(plot_mutation_signatures(mut.sig, pid, sid, fraction_threshold = 0.5), 
	       filename='mutation_singaure_decompoisition_combined.pdf', width=5, height=6)

	### patient wise
	patients = unique(mcp.maf.fil$patient)
	tumor_samples = unique(mcp.maf.fil$Tumor_Sample_Barcode)

	for(pat in patients){
		primary <- sampleID[patient == pat & `Sample Class` == 'Primary', `Sample ID`]
		primary
		pat.maf <- mcp.maf.fil[patient== pat,] 
		pat.maf
		print(paste("make_variant_classification_plot for patient ", pat))
		make_variant_classification_plot(pat.maf, out = pat)
		print("make_ccf_plots")
		make_ccf_plots(pat.maf, tumor_samples)
		print("make_stratton_plots")
		make_stratton_plots(pat.maf, tumor_samples, out = pat)
		print("make_binary_tree")
		print("make_mutation_signatures")
		write.table(pat.maf, 'pat.maf', sep="\t", quote=F, row.names=F)
		cmd = paste0('python ~/program/autospy/inst/mutation-signatures/main.py --seed 100 ')
		cmd = paste0(cmd, '~/program/autospy/inst/mutation-signatures/Stratton_signatures29.txt pat.maf signature.out')
		system(cmd)
		fread("signature.out") -> mut.sig
		ggsave(plot_mutation_signatures(mut.sig, pid, sid, fraction_threshold = 0.5), 
		       filename=paste0('mutation_singaure_decompoisition_', pat, '.pdf'), width=5, height=6)
		samples = unique(pat.maf$sample)
		sample_pairs <- combn(samples, 2, simplify = F)
		dir.create('ccf_2d_plots')
		for(sample_pairs in sample_pairs){
			make_ccf_2d(pat.maf,
				    sample_pairs,
				    out = pat,
				    directory = 'ccf_2d_plots')
		}
	}

	tmp = copy(mcp.maf.fil)
	tmp[, patient := sub('C_', '', patient)]
	pats = unique(tmp$patient)
	pats
	setwd('res/autospy_tempo')
	for(pat in pats){
		pat.maf = tmp[patient == pat, ]
		if(length(unique(pat.maf$sample))>1){
			pat.maf[, IMPACT_410 := IMPACT_468]
			pp = sdt[patient == pat & Primary == 'yes', sname]
			make_binary_tree(pat.maf, pp, hotspots = autospy::hotspots, vertical = TRUE, margins = TRUE, out = pat)
		}
	}
	primary = 's_C_JJ9RNX_P003_d'

	tmp = wes.maf
	tmp = maf
	t1 = tmp[Tumor_Sample_Barcode == 's_C_JJ9RNX_P003_d', Hugo_Symbol]
	t2 = tmp[Tumor_Sample_Barcode == 's_C_JJ9RNX_P004_d', Hugo_Symbol]
	t3 = tmp[Tumor_Sample_Barcode == 's_C_JJ9RNX_P005_d', Hugo_Symbol]
	intersect(t1, t2)
	intersect(t2, t3)
	wes.maf[Hugo_Symbol %in% intersect(t1, t3), cc.sel, with=F][order(Hugo_Symbol), ]
	wes.maf[Hugo_Symbol == 'EPHA2', cc.sel, with=F ]

}
