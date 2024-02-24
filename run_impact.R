library(data.table)
library(maftools)
library(ggplot2)
library(scales)
library(ggrepel)
library(future)
library(sqldf)
require(RSQLite)

lapply(paste('package:', names(sessionInfo()$otherPkgs), sep=""), detach, character.only=TRUE, unload=TRUE)


plan(multiprocess)
s.i %<-% { save() }
s.i

setwd('/ifs/work/solitlab/huw/solit/study/hiseq/impact')

## sqlite3
{
	con = dbConnect(RSQLite::SQLite(), '/ifs/work/solitlab/huw/sqlite/maf.sqlite')
	con.rna = dbConnect(RSQLite::SQLite(), '/ifs/work/solitlab/huw/sqlite/tcga_rna.sqlite')
	con.clin = dbConnect(RSQLite::SQLite(), '/ifs/work/solitlab/huw/sqlite/tcga_clin.sqlite')
	con.meth = dbConnect(RSQLite::SQLite(), '/ifs/work/solitlab/huw/sqlite/tcga_meth.sqlite')
	conn = dbConnect(RSQLite::SQLite(), '/ifs/work/solitlab/huw/sqlite/ccle.sqlite')
	dbDisconnect(con)
	dbWriteTable(con, 'impact', maf) # datable of impact maf file
	dbWriteTable(con, 'tcga', tcga.maf, overwite=T) # table of tcga maf file
	dbWriteTable(con, 'tcga_cnv', tcga.cnv, overwite=T) # table of tcga cnv
	dbWriteTable(con, 'tcga_tsb_disease', tmp) # table of tsb and disease name
	dbWriteTable(con, 'tcga_tsb_disease_mapt', tsb.mapt.df) # table of tsb and disease name only for patient has mapt mutations
	dbWriteTable(con, 'tcga_clin', tcga.clin) # tcga a simple clinical information
	dbWriteTable(con, 'laml_clin', laml.clin) # laml a simple clinical information
	dbWriteTable(con, 'laml_cnv', laml.cnv) # laml a simple cnvical information
	dbWriteTable(con, 'laml', laml.maf) # laml a simple mafical information

	dbSendStatement(con, 'create index tcga_disease on tcga(DISEASE)')

	skcm = dbGetQuery(con, 'select * from impact where CANCER_TYPE = "Melanoma"')
	dbGetQuery(con, 'select distinct CANCER_TYPE from impact')


}

## import impact data into sqlite
{
	maffile = 'msk-impact_20180625/data_mutations_merged_maf2maf_facets_clin_oncokb_prepped.txt'
	fread(maffile) -> maf
	maf[, vaf := t_alt_count / t_depth]
	maf = maf[vaf > 0.05 & t_alt_count > 3, ]
	dim(maf)
	colnames(maf)
	length(unique(imp.skcm$Tumor_Sample_Barcode))
	unique(maf$CANCER_TYPE)
	maf[CANCER_TYPE == 'Melanoma', ] -> imp.skcm
	t1 = imp.skcm[Hugo_Symbol == 'NF1', Tumor_Sample_Barcode]
	t2 = imp.skcm[Hugo_Symbol == 'BRAF' & HGVSp_Short == 'p.V600E', Tumor_Sample_Barcode]
	t22 = imp.skcm[Hugo_Symbol == 'BRAF', Tumor_Sample_Barcode]
	t3 = imp.skcm[Hugo_Symbol == 'TP53', Tumor_Sample_Barcode]
	t4 = imp.skcm[Hugo_Symbol == 'NRAS', Tumor_Sample_Barcode]
	## NF1 and BRAF V600E, no TP53 or Ras
	setdiff (intersect(t1, t2), c(t3, t4)) -> ss
	sel = imp.skcm[Tumor_Sample_Barcode %in% ss, .(Tumor_Sample_Barcode, Hugo_Symbol, HGVSp_Short, oncogenic)][Hugo_Symbol %in% c('TP53', 'NF1', 'BRAF', 'NRAS'), ]
	sel
	fwrite(sel, file='res/caseid_nP53_nRAS_yNf1_yBrafV600E.xls', sep="\t")
	sync('res')

	# no triple
	intersect(t3, intersect(t1, t2)) -> ss
	sel = imp.skcm[Tumor_Sample_Barcode %in% ss, .(Tumor_Sample_Barcode, Hugo_Symbol, HGVSp_Short, oncogenic)][Hugo_Symbol %in% c('TP53', 'NF1', 'BRAF'), ]
	imp.skcm[Tumor_Sample_Barcode %in% ss, .(Tumor_Sample_Barcode, Hugo_Symbol, HGVSp_Short, oncogenic)][Hugo_Symbol %in% c('TP53', 'NF1', 'BRAF', 'NRAS'), ]
	fwrite(sel, file='res/caseid_p53_nf1_braf.xls', sep="\t")


	# NF1 P53, no BRAF NRAS
	intersect(t1, t3) -> ss
	setdiff(ss, t22) -> ss
	setdiff(ss, t4) -> ss
	sel = imp.skcm[Tumor_Sample_Barcode %in% ss, .(Tumor_Sample_Barcode, Hugo_Symbol, HGVSp_Short, oncogenic)][Hugo_Symbol %in% c('TP53', 'NF1', 'BRAF'), ]
	unique(sel$Tumor_Sample_Barcode) 
	fwrite(sel, file='res/caseid_p53_nf1.xls', sep="\t")
	sync('res')

	# BRAF 600E, P53, no NF1, NRAS
	intersect(t2, t3) -> ss
	setdiff(ss, t1) -> ss
	setdiff(ss, t4) -> ss
	sel = imp.skcm[Tumor_Sample_Barcode %in% ss, .(Tumor_Sample_Barcode, Hugo_Symbol, HGVSp_Short, oncogenic)][Hugo_Symbol %in% c('TP53', 'NF1', 'BRAF'), ]
	unique(sel$Tumor_Sample_Barcode) 
	fwrite(sel, file='res/caseid_p53_BRAF.xls', sep="\t")
	sync('res')

	intersect(t1, t3) -> ss
	# no triple
	sel = imp.skcm[Tumor_Sample_Barcode %in% ss, .(Tumor_Sample_Barcode, Hugo_Symbol, HGVSp_Short, oncogenic)][Hugo_Symbol %in% c('TP53', 'NF1', 'BRAF'), ]
	fwrite(sel, file='res/caseid_p53_nf1.xls', sep="\t")

	dbWriteTable(con, 'impact', maf)
	dbSendStatement(con, 'create index impact_hs on impact(Hugo_Symbol)')
	dbSendStatement(con, 'create index impact_tsb on impact(Tumor_Sample_barcode)')
	dbSendStatement(con, 'alter table impact add flag int')

	maf.blca = maf[CANCER_TYPE == 'Bladder Cancer', ] 

	imp.clin.patient.file = 'msk-impact_20180625/data_clinical_patient.txt'
	fread(imp.clin.patient.file, skip=4) -> imp.clin.patient
	imp.clin.patient$pa:t
	cn = colnames(imp.clin.patient)
	cn
	cn[grep('STAGE', cn)]
	imp.clin.patient[, grep('STAGE', cn), with=F] # CDB_INIT_X_STAGE ?
	cn
	imp.clin.patient[, .(INITIAL_STAGE, STAGE)]
	dbWriteTable(con, 'imp_clin_patient', imp.clin.patient)

	kk.id = scan('kk.id', character()); kk.id
	kk.id = fread('kk.id', header=F)
	kk.id = merge(kk.id, imp.clin.patient, by.x = 'V1', by.y='PATIENT_ID')
	kk.id
	fwrite(kk.id, sep="\t", file='kk_id.xls')
	source('~/pipeline/fun/sync.r')
	scp('kk_id.xls')

	maf[substr(Tumor_Sample_Barcode, 1, 9) %in% kk.id, ] -> tmp
	fname = 'res/dmp_msi.xls'
	fwrite(tmp[!duplicated(Tumor_Sample_Barcode), .(Tumor_Sample_Barcode, MSI_SCORE, msi_high)], file=fname)
	scp(fname)

	imp.clin.sample.file = 'msk-impact_20180625/data_clinical_sample.txt'
	fread(imp.clin.sample.file, skip = 4) -> imp.clin.sample
	head(imp.clin.sample)
	dbWriteTable(con, 'imp_clin_sample', imp.clin.sample)

	cnvfile = '/ifs/res/taylorlab/ang46/data/allgenes_facets_calls/collapsed.txt' 
	impact.cnv = fread(cnvfile)
	save(impact.cnv, file='impact_cnv_allgene_facets_call.RData')

	dbWriteTable(con, 'impact_cnv', impact.cnv, overwite=T)
	dbSendStatement(con, 'create index tcga_cnv_hs on tcga_cnv(Hugo_Symbol)')
	dbSendStatement(con, 'create index tcga_cnv_tsb on tcga_cnv(Tumor_Sample_barcode)')
	dbSendStatement(con, 'alter table tcga_cnv add flag int')
	#rm(maf)

}

##TCGA data import into sqlite
{
	maffile = '/ifs/res/taylorlab/ang46/data/tcga/pancan_atlas/original_data/mc3.v0.2.8.PUBLIC.LAML_PATCH_prepped_facets_oncokb_hotspots_absolute.maf'
	tcga.maf = fread(maffile)
	tcga.maf[, vaf := t_alt_count / t_depth]
	setnames(tcga.maf, 'STRAND', 'STRAND2')

	dbWriteTable(con, 'tcga', tcga.maf, overwite=T)
	dbSendStatement(con, 'create index tcga_hs on tcga(Hugo_Symbol)')
	dbSendStatement(con, 'create index tcga_tsb on tcga(Tumor_Sample_barcode)')
	dbSendStatement(con, 'alter table tcga add flag int')

	tcga.clin.file = '/ifs/res/taylorlab/ang46/data/tcga/pancan_atlas/original_data/pancancer_sample_info.tsv'
	tcga.clin = fread(tcga.clin.file)
	save(tcga.clin, file='tcga_clin_info.RData')

	dbWriteTable(con, 'tcga_clin', tcga.clin)

	## copy number variation
	fd = '/ifs/res/taylorlab/tcga_wes_facets/merged_calls'
	tcga.cnv.files = list.files(fd)
	base = toupper(sub("_.*", "", tcga.cnv.files)) 
	tcga.cnv.files = paste0(fd, '/', tcga.cnv.files)
	file.exists(tcga.cnv.files)
	sapply(tcga.cnv.files, fread) -> tcga.cnv.lst
	head(tcga.cnv.lst[[1]])
	names(tcga.cnv.lst) = base
	rbindlist(tcga.cnv.lst, fill=T,  use.names = T, idcol=T) -> tcga.cnv.dt
	tcga.cnv.dt[, Tumor_Sample_Barcode := sub(".*(TCGA.*)_Dual.*", "\\1", Tumor_Sample_Barcode)]
	tcga.cnv.dt[grep('Dual', Tumor_Sample_Barcode, invert=T), Tumor_Sample_Barcode := sub(".*(TCGA.*)_his.*", "\\1", Tumor_Sample_Barcode)]
	tcga.cnv.dt[, patient := substr(Tumor_Sample_Barcode, 1, 12)]
	tail(tcga.cnv.dt)

	dbWriteTable(con, 'tcga_cnv', tcga.cnv.dt)

	dbSendStatement(con, 'drop table if exists tcga_cnv')
	dbSendStatement(con, 'drop index tcga_cnv_hs')
	dbSendStatement(con, 'drop index tcga_cnv_tsb')
	dbSendStatement(con, 'create index tcga_cnv_hs on tcga_cnv(Hugo_Symbol)')
	dbSendStatement(con, 'create index tcga_cnv_tsb on tcga_cnv(Tumor_Sample_barcode)')
	dbSendStatement(con, 'alter table tcga_cnv add flag int')

	dbGetQuery(con, 'select count(*) from tcga')
	dbGetQuery(con, 'select count(*) from tcga_cnv')
	dbGetQuery(con, 'select * from tcga_cnv limit 10')

	dbGetQuery(con, 'select name from sqlite_master where type="table" order by name')
	dbGetQuery(con, 'desc table sqlite_master order by name')


}

## TCGA basal vs uro
## TCGA squamous vs uro
{
	ba.vs.uro = fread('../tcga/tcga_basal_vs_uro.tsv')
	sq.vs.uro = fread('../tcga/tcga_sq_vs_uro.tsv')
	sq.vs.uro[symbol == 'NR3C1', ]
	ba.vs.uro[symbol == 'NR3C1', ]
	rm(sq.vs.uro)
	rm(ba.vs.uro)
	ls()

}


## mapt
{
	as.data.table(dbGetQuery(conn, 'select * from ccle_maf limit 10'))

	## tcga 
	{
		source('/home/huw/program/fun/cbioMaf.R')
		source('/home/huw/program/fun/sync.r')
		tcga.Gmapt.maf = as.data.table(dbGetQuery(con, 'select * from tcga where Hugo_Symbol == "MAPT"'))
		dim(tcga.Gmapt.maf)
		aa = table(tcga.Gmapt.maf$HGVSp_Short)
		aa[order(aa)]
		tcga.Gmapt.maf[, cc.1, with=F]
		tcga.Gmapt.maf$oncogenic
		tcga.Gmapt.maf.o = read.maf(tcga.Gmapt.maf)
		pdf('mt/tcga_lollipop.pdf', width=6, height=3)
		lollipopPlot(tcga.Gmapt.maf.o, gene='MAPT', cBioPortal=T)
		dev.off()
		sync('mt')
	}

	{# ccle mutations lollipop
		# others can see from ccle.r
		tmp = as.data.table(dbGetQuery(conn, 'select * from ccle_maf where Hugo_Symbol == "MAPT"'))
		tmp[, Tumor_Seq_Allele2 := Tumor_Seq_Allele1 ]
		tmp.o = read.maf(tmp)
		pdf('mt/ccle_lollipop.pdf', width=6, height=3)
		lollipopPlot(tmp.o, gene='MAPT', cBioPortal=T, labelPos = T)
		dev.off()
		sync('mt')
		tmp[grep('SKIN', Tumor_Sample_Barcode), Tumor_Sample_Barcode]


	}

	{# % of cases have mutations of mapt by disease in tcga
		tsb.df = dbGetQuery(con, 'select DISEASE, Tumor_Sample_Barcode, from tcga')
		tsb.df = as.data.table(tsb.df)
		tsb.df[!duplicated(Tumor_Sample_Barcode), ] -> tsb.df

		dbWriteTable(con, 'tcga_tsb_disease', tmp)

		tmp = dbGetQuery(con, 'select DISEASE, Tumor_Sample_Barcode from tcga where Hugo_Symbol = "MAPT"')
		tmp = as.data.table(tmp)
		tmp[!duplicated(Tumor_Sample_Barcode), ] -> tsb.mapt.df
		dbWriteTable(con, 'tcga_tsb_disease_mapt', tsb.mapt.df)

		tsb.df.n = tsb.df[, nn := .N, by = 'DISEASE'][!duplicated(DISEASE), ][, .(DISEASE, nn)]
		tsb.df.n
		tsb.mapt.df.n = tsb.mapt.df[, n := .N, by = 'DISEASE'][!duplicated(DISEASE), ][, .(DISEASE, n)]
		tsb.mapt.df.n
		tsb.n = merge(tsb.df.n, tsb.mapt.df.n, by = 'DISEASE') 
		tsb.n[, lbl := paste0(n, '/', nn)]
		tsb.n[, rr := 100*n/nn]
		tsb.n[, or := rr * n]
		tsb.n
		gg = ggplot(tsb.n, aes(x=reorder(DISEASE, or), y = nn)) +
			geom_bar(stat='identity') +
			theme(axis.text.x = element_text(angle=90, size=6)) +
			xlab('') + ylab('MAPT mutated cases (%)')
		ggsave(gg, file=(paste0('mt/tcga_cases_pertentage.pdf')))
		sync('mt')
	}

	cbio = cbioMaf(tcga.Smapt.maf)
	cbio.o = read.maf(cbio$maf, vc_nonSyn = cbio$vc.key)
	pdf('mt/tcga_mt_oncoplot.pdf', width=12, height=6)
	oncoplot(cbio.o, keepGeneOrder=T, fontSize=6, color = cbio$color.table, GeneOrderSort = F, drawRowBar=F, drawColBar=F, sortByMutation=T, writeMatrix=T, showTumorSampleBarcodes=F, SampleNamefontSize=6, top=30)
	dev.off()
	sync('mt')

	tmp = tcga.Smapt.maf[Hugo_Symbol %in% c(i468, 'MAPT'), ][oncogenic %in% onco.sel, ]
	tmp1 = tcga.Smapt.maf[Hugo_Symbol %in% c('MAPT'), ]
	rbind(tmp, tmp1) -> tmp
	rm(tmp1)

	cbio = cbioMaf(tmp)
	cbio.o = read.maf(cbio$maf, vc_nonSyn = cbio$vc.key)
	pdf('mt/tcga_Smapt_oncoplot_Gimpact.pdf', width=12, height=6)
	oncoplot(cbio.o, keepGeneOrder=T, fontSize=6, color = cbio$color.table, GeneOrderSort = F, drawRowBar=F, drawColBar=F, sortByMutation=T, writeMatrix=T, showTumorSampleBarcodes=F, SampleNamefontSize=6, top=30)
	dev.off()
	sync('mt')
	
	## look at bladder cancer mapt mutated patient
	{
		tcga.Smapt.maf[DISEASE == 'BLCA' & Hugo_Symbol %in% i468, ] -> tcga.blca.Smapt.maf
		tcga.blca.Smapt.maf
		cbio = cbioMaf(tcga.blca.Smapt.maf)
		cbio.o = read.maf(cbio$maf, vc_nonSyn = cbio$vc.key)
		pdf('mt/tcga_blca_Smapt_oncoplot_Gimpact.pdf', width=12, height=6)
		oncoplot(cbio.o, keepGeneOrder=T, fontSize=6, color = cbio$color.table, GeneOrderSort = F, drawRowBar=F, drawColBar=F, sortByMutation=T, writeMatrix=T, showTumorSampleBarcodes=F, SampleNamefontSize=6, top=30)
		dev.off()
		sync('mt')

		tcga.blca.Smapt.maf[Hugo_Symbol == 'NF1', cc.]

	}

}

## impact kdm6a mutation in blca analysis
## impact kdm6a mutation clonality
{
	unique(maf$CANCER_TYPE)
	unique(maf[CANCER_TYPE == 'Bladder Cancer', 'CANCER_TYPE_DETAILED'] )
	maf.blca = maf[CANCER_TYPE == 'Bladder Cancer', ] 
	ls()
	fread('msk_impact_blca_grade.csv') -> grade
	setkey(grade, 'id')
	grade
	setkey(maf.blca, 'Tumor_Sample_Barcode')
	maf.blca = merge(maf.blca, grade[, .(id, class)], by.x = 'Tumor_Sample_Barcode', by.y = 'id')
	#g1 = 'KDM6A'
	#g2 = 'STAG2'
	#t1 = maf.blca[Hugo_Symbol %in% g1 & oncogenic %in% onco.sel, Tumor_Sample_Barcode]
	#t2 = maf.blca[Hugo_Symbol %in% g2 & oncogenic %in% onco.sel, Tumor_Sample_Barcode]
	#t1;t2
	#tt = intersect(t1, t2)
	#tmp = maf.blca[Hugo_Symbol %in% c('CDKN1A', 'KDM6A', 'ARID1A', 'TP53') & oncogenic %in% onco.sel & Tumor_Sample_Barcode %in% tt, ]

	tmp = maf.blca[Hugo_Symbol %in% c('CDKN1A', 'KDM6A', 'ARID1A', 'TP53', 'STAG2') & oncogenic %in% onco.sel, ]
	tmp[, class := factor(class, levels=c('LG', 'HG', 'Met'))]
	gg = ggplot(tmp, aes(x=class, y = ccf_1copy, fill = class)) + 
		geom_bar(stat='summary', fun.y='mean', position='dodge') +
		stat_summary(fun.data = mean_se, geom = "errorbar", width=0.3, size=1) + 
		facet_wrap(~ Hugo_Symbol) + xlab('Clonality')
	ggsave(gg, file='res/impact_clonality.pdf', width=3.5, height=6)
	gg = ggplot(tmp, aes(x=class, y =vaf, fill = class)) + 
		geom_bar(stat='summary', fun.y='mean', position='dodge') +
		stat_summary(fun.data = mean_se, geom = "errorbar", width=0.3, size=1) + 
		facet_wrap(~ Hugo_Symbol) + xlab('Clonality')
	ggsave(gg, file='res/impact_clonality_vaf.pdf', width=3.5, height=6)
	sync('res')

	t.test(tmp[class == 'LG' & Hugo_Symbol == 'KDM6A', vaf], tmp[class == 'HG' & Hugo_Symbol == 'KDM6A', vaf]) 
	t.test(tmp[class == 'LG' & Hugo_Symbol == 'KDM6A', vaf], tmp[class == 'Met' & Hugo_Symbol == 'KDM6A', vaf]) 
	t.test(tmp[class == 'LG' & Hugo_Symbol == 'ARID1A', vaf], tmp[class == 'HG' & Hugo_Symbol == 'KDM6A', vaf]) 
	t.test(tmp[class == 'LG' & Hugo_Symbol == 'ARID1A', vaf], tmp[class == 'Met' & Hugo_Symbol == 'KDM6A', vaf]) 

	t.test(tmp[class == 'LG' & Hugo_Symbol == 'KDM6A', ccf_1copy], tmp[class == 'HG' & Hugo_Symbol == 'KDM6A', ccf_1copy]) 
	t.test(tmp[class == 'LG' & Hugo_Symbol == 'KDM6A', ccf_1copy], tmp[class == 'Met' & Hugo_Symbol == 'KDM6A', ccf_1copy]) 
	t.test(tmp[class == 'LG' & Hugo_Symbol == 'ARID1A', ccf_1copy], tmp[class == 'HG' & Hugo_Symbol == 'KDM6A', ccf_1copy]) 
	t.test(tmp[class == 'LG' & Hugo_Symbol == 'ARID1A', ccf_1copy], tmp[class == 'Met' & Hugo_Symbol == 'KDM6A', ccf_1copy]) 

	gg = ggplot(tmp, aes(x=class, y = ccf_1copy, fill = class, color = class)) + 
		geom_jitter(width=.2) +
		facet_wrap(~ Hugo_Symbol) + xlab('Clonality')
	ggsave(gg, file='res/impact_clonality_points.pdf', width=3.5, height=6)
	gg = ggplot(tmp, aes(x=class, y =vaf, fill = class, color = class)) + 
		geom_jitter(width=.2) +
		facet_wrap(~ Hugo_Symbol) + xlab('Clonality')
	ggsave(gg, file='res/impact_clonality_vaf_points.pdf', width=3.5, height=6)
	sync('res')
}


## impact pancancer melanoma odd ratio of two gene mutations
## and survival curve
## another try
## save(maf.skcm, file='impact_maf_skcm.RData')
{
	maffile = 'msk-impact_20180625/data_mutations_merged_maf2maf_facets_clin_oncokb_prepped.txt'
	fread(maffile) -> maf
	maf[, vaf := t_alt_count / t_depth]
	maf = maf[vaf > 0.05 & t_alt_count > 3, ]
	maf.skcm = maf[CANCER_TYPE == 'Melanoma', ] 
	maf.skcm[Hugo_Symbol == 'BRAF', HGVSp_Short]
	save(maf.skcm, file='impact_maf_skcm.RData')

	maf[Tumor_Sample_Barcode == 's_JC_ucc_032_TZ', ]
	utuc.clin = fread('~/utuc_cbe_colemanj_JC_ucc_clinical_data.tsv')
	maf[Tumor_Sample_Barcode %in% utuc.clin[, `Sample ID`], ]
	maf$Tumor_Sample_Barcode 

	## combine tcga and impact together
	tcga.maf.skcm
	ov = intersect(colnames(tcga.maf.skcm), colnames(maf.skcm))
	tmp.maf = rbind(tcga.maf.skcm[, ov, with=F], maf.skcm[, ov, with=F])
	tmp.maf

	#tmp.maf = copy(maf.skcm) ## impact only
	tmp.maf = tmp.maf[Hugo_Symbol == 'BRAF' & HGVSp_Short == 'p.V600E', Hugo_Symbol := 'BRAF:V600E']
	tmp.maf = tmp.maf[Hugo_Symbol == 'BRAF', Hugo_Symbol := 'BRAF:OTHERS']
	tmp.maf
	tmp.maf[Hugo_Symbol == 'BRAF', HGVSp_Short]

	onco.orig = onco.sel
	onco.sel = unique(tmp.maf$oncogenic)
	#genes = c('NF1', 'TP53', 'V600E', 'BRAF:NOS', 'NRAS', 'PTEN', 'CDKN2A'); #, 'KRAS', 'HRAS')
	genes = c('BRAF', 'NF1', 'BRAF:V600E', 'BRAF:OTHERS', 'TP53', 'NRAS', 'PTEN', 'CDKN2A', 'MAP2K1', 'KRAS', 'HRAS')
	dd = data.table(gene1 = character(),  gene2 = character(),  gene1.cnt = numeric(),  gene2.cnt = numeric(),  both.cnt = numeric(),  gene1.only = numeric(),  gene2.only = numeric(),  none.cnt = numeric(),  ttl = numeric())
	for( i in 1:(length(genes)-1)){
		for(j in (i+1):length(genes)){
			tmp1 = unique(tmp.maf[Hugo_Symbol == genes[i] & oncogenic %in% onco.sel, Tumor_Sample_Barcode])
			tmp2 = unique(tmp.maf[Hugo_Symbol == genes[j] & oncogenic %in% onco.sel, Tumor_Sample_Barcode])
			tmp1 = substr(tmp1, 1, 12)
			tmp2 = substr(tmp2, 1, 12)
			co = intersect(tmp1, tmp2)
			ttl = length(unique(tmp.maf$Tumor_Sample_Barcode))
			none = length(unique(tmp.maf[!(Hugo_Symbol %in% c(genes[i], genes[j])) & oncogenic %in% onco.sel, Tumor_Sample_Barcode]))
			none = ttl - length(unique(tmp1, tmp2))
			dd = rbind(dd, list(genes[i], genes[j], length(tmp1), length(tmp2), length(co), length(setdiff(tmp1, tmp2)), length(setdiff(tmp2, tmp1)), none, ttl))
		}
	}
	dd

	## add BRAF
	#genes = c('BRAF', 'NF1', 'TP53', 'NRAS', 'PTEN', 'CDKN2A', 'MAP2K1', 'KRAS', 'HRAS')
	i = 1
	for(j in (i+1):length(genes)){
		tmp1 = unique(tcga.maf.skcm[Hugo_Symbol == genes[i] & oncogenic %in% onco.sel, Tumor_Sample_Barcode])
		tmp2 = unique(tcga.maf.skcm[Hugo_Symbol == genes[j] & oncogenic %in% onco.sel, Tumor_Sample_Barcode])
		tmp1 = substr(tmp1, 1, 12)
		tmp2 = substr(tmp2, 1, 12)
		co = intersect(tmp1, tmp2)
		ttl = length(unique(tcga.maf.skcm$Tumor_Sample_Barcode))
		none = length(unique(tcga.maf.skcm[Hugo_Symbol %in% c(genes[i], genes[j]) & oncogenic %in% onco.sel, Tumor_Sample_Barcode]))
		dd = rbind(dd, list(genes[i], genes[j], length(tmp1), length(tmp2), length(co), length(setdiff(tmp1, tmp2)), length(setdiff(tmp2, tmp1)), none, ttl))
	}
	dd

	# (tab <- matrix(c(38, 25, 162, 75), nrow=2))
	#      [,1] [,2]
	# [1,]   38  162
	# [2,]   25   75
	#fisher.test(tab)

	dd[, or := 0]
	dd[, pvalue := 0]
	for(i in 1:nrow(dd)){
		(tab <- matrix(c(dd[i, both.cnt], dd[i, gene1.only], dd[i, gene2.only], dd[i, none.cnt]), nrow=2))
		dd[i, pvalue := fisher.test(tab)[[1]]]
		dd[i, or := fisher.test(tab)[[3]]]
	}
	dd[, padj := p.adjust(pvalue, 'fdr')]
	dd[, lbl := paste0(gene1, ' | ', gene2)]
	dd
	#dd[, lbl2 := ifelse(gene1 == 'NF1', gene2, gene1)]
	dd[, lbl2 := ifelse(gene1 == 'TP53', gene2, gene1)]
	dd[, clr := 'black']
	dd[padj > 0.05, clr := 'grey']

	ddd = dd[gene1 == 'TP53' | gene2 == 'TP53',]
	gg = ggplot(ddd, aes(log(or), -log(padj), label=lbl, size=both.cnt, color=clr)) + geom_point()
	gg = gg + geom_hline(yintercept = -log(0.05), color = 'grey')
	gg = gg + geom_text_repel(data=ddd, aes(label=lbl2), size=3)
	gg = gg + xlab('log Odd Ratio')
	gg = gg + ylab('-log(p value)')
	#ggsave(gg, file='nf1/impact_skcm_odd_ratio_oncogenic_mutations.pdf', width=4, height=4)
	ggsave(gg, file='nf1/tcga_impact_skcm_P53_odd_ratio.pdf', width=4.5, height=4)
	sync('nf1')

}

##  NF1 cases
## oncoplot
{
	tcga.maf.skcm = tcga.skcm
	tcga.maf.skcm = tcga.maf.skcm[, !duplicated(colnames(tcga.maf.skcm)), with=F]
	source('~/program/fun/cbioMaf.R')
	source('~/program/fun/param.r')
	i468 <- scan("/ifs/depot/resources/dmp/data/mskdata/interval-lists/VERSIONS/cv6/genelist", character())
	i468
	genes.ex = c('TP53', 'BRAF', 'NRAS', 'NF1', 'HRAS', 'KRAS', 'PTEN', 'CDKN2A') 
	genes.ex = c('BRAF', 'NRAS')
	sample.ex = tcga.maf.skcm[Hugo_Symbol %in% genes.ex, Tumor_Sample_Barcode]
	maf.ex = tcga.maf.skcm[!(Tumor_Sample_Barcode %in% sample.ex) & oncogenic %in% onco.sel, ] 
	tmp = cbioMaf(maf.ex)
	tmp.o = read.maf(tmp$maf, vc_nonSyn = tmp$vc.key)
	pdf('nf1/tcga_wt_BRAF_NRAS.pdf', width=12, height=6)
	oncoplot(tmp.o, keepGeneOrder=T, fontSize=6, color = tmp$color.table, GeneOrderSort = F, drawRowBar=F, drawColBar=F, sortByMutation=T, writeMatrix=T, showTumorSampleBarcodes=F, SampleNamefontSize=6, top=30, title='BRAF and NRAS WT cases')
	dev.off()
	sync('nf1')

	g = 'PTPRC'
	tmp = maf.skcm[!(Hugo_Symbol %in% g), ]
	length(unique(tmp))

	t1 = as.data.table(as.data.frame(table(maf.ex$Tumor_Sample_Barcode)))
	t2 = as.data.table(as.data.frame(table(maf.skcm$Tumor_Sample_Barcode)))
	median(t1$Freq)
	median(t2$Freq)
	summary(t2$Freq)
	nf1.id = maf.skcm[Hugo_Symbol == 'NF1' & oncogenic != '', Tumor_Sample_Barcode]
	maf.skcm[ Tumor_Sample_Barcode %in% nf1.id & Hugo_Symbol %in% i468, ] -> tmp

}

## pancancer melanoma odd ratio of two gene mutations
## save(maf.skcm, file='impact_maf_skcm.RData')
## and survival curve
{
	maffile = 'msk-impact_20180625/data_mutations_merged_maf2maf_facets_clin_oncokb_prepped.txt'
	fread(maffile) -> maf
	maf[, vaf := t_alt_count / t_depth]
	maf = maf[vaf > 0.05 & t_alt_count > 3, ]
	maf.skcm = maf[CANCER_TYPE == 'Melanoma', ] 
	save(maf.skcm, file='impact_maf_skcm.RData')
	tmp.maf = copy(maf.skcm)
	tmp.maf = tmp.maf[Hugo_Symbol == 'BRAF' & HGVSp_Short == 'p.V600E', Hugo_Symbol := 'V600E']
	tmp.maf = tmp.maf[Hugo_Symbol == 'BRAF', Hugo_Symbol := 'BRAF:NOS']

	genes = c('TP53', 'NF1', 'BRAF:NOS', 'V600E', 'NRAS', 'PTEN', 'CDKN2A') 
	i = 'NF1'; j = 'BRAF'; x = 'TP53'
	i = 'NF1'; j = 'V600E'; v = 'BRAF:NOS'; x = 'TP53'
	i = 'NRAS'; j = 'HRAS';
	tmp1 = tmp.maf[Hugo_Symbol == i,    Tumor_Sample_Barcode]; tmp1
	tmp2 = tmp.maf[Hugo_Symbol == j,    Tumor_Sample_Barcode]; tmp2
	tmp3 = tmp.maf[Hugo_Symbol == v,    Tumor_Sample_Barcode]
	tmp4 = tmp.maf[Hugo_Symbol == x,    Tumor_Sample_Barcode]
	tmp1 = substr(tmp1, 1, 12); tmp1
	tmp2 = substr(tmp2, 1, 12); tmp2
	tmp3 = substr(tmp3, 1, 12); tmp3
	tmp4 = substr(tmp4, 1, 12); tmp4
	tcga.skcm.clin[, mut := NULL]
	tcga.skcm.clin[, mut := 'NOS']
	tcga.skcm.clin[bcr_patient_barcode %in% tmp1, mut := i]
	tcga.skcm.clin[bcr_patient_barcode %in% tmp2, mut := j]
	tcga.skcm.clin[bcr_patient_barcode %in% tmp3, mut := v]
	tcga.skcm.clin[bcr_patient_barcode %in% tmp4, mut := x]
	tmp = intersect(tmp1, tmp3)
	tcga.skcm.clin[bcr_patient_barcode %in% tmp, mut := paste0(i, '&', v)]
	tmp = intersect(tmp1, tmp2)
	tcga.skcm.clin[bcr_patient_barcode %in% tmp, mut := paste0(i, '&', j)]
	table(tcga.skcm.clin$mut)
	tmp = intersect(tmp1, tmp2)
	tcga.skcm.clin[bcr_patient_barcode %in% tmp, mut := paste0(i, '&', j)]
	#tcga.skcm.clin[bcr_patient_barcode %in% tmp, mut := paste0(genes[i], ':', genes[j])]
	fname = paste0('surv/tcga_surv_', i, '_', j, '.pdf') 
	tmp = tcga.skcm.clin[mut != 'NOS', ]
	TCGAanalyze_survival(as.data.frame(tmp), risk.table = F, "mut", main = paste0('Survival curve w/ ', genes[i], ' & ', genes[j], ' mutations'), height = 5, width=8, filename = fname)
	TCGAanalyze_survival(as.data.frame(tcga.skcm.clin), risk.table = F, "mut", main = paste0('Survival curve w/ ', genes[i], ' & ', genes[j], ' mutations'), height = 5, width=8, filename = fname)
	sync('surv')

	tcga.maf.skcm[HGVSp_Short == 'p.D594E', .(Hugo_Symbol, HGVSp_Short) ]
	tcga.maf.skcm[grep('D594', HGVSp_Short), .(Hugo_Symbol, HGVSp_Short) ]
	tcga.maf.skcm[Hugo_Symbol == 'BRAF' & HGVSp_Short != 'p.V600E', .(Hugo_Symbol, HGVSp_Short) ]
	table(tcga.maf.skcm[Hugo_Symbol == 'BRAF' & HGVSp_Short != 'p.V600E', HGVSp_Short])
	table(tcga.maf.skcm[Hugo_Symbol == 'BRAF', HGVSp_Short])
	tmp1
	tmp2
	tmp = intersect(tmp1, tmp2); tmp
	tmp
	intersect(tmp, tmp4)
	colnames(tcga.maf.skcm)
	tcga.skcm.clin
	length(unique(tmp))
	length(unique(substr(tmp, 1, 12)))
	(tab <- matrix(c(23, 162, 37, 110), nrow=2))
	#      [,1] [,2]
	# [1,]   38  162
	# [2,]   25   75
	fisher.test(tab)

	i = 'NF1'; j = 'BRAF'; x = 'TP53'
	tmp1 = unique(tcga.maf.skcm[Hugo_Symbol == i,    Tumor_Sample_Barcode]); tmp1
	tmp2 = unique(tcga.maf.skcm[Hugo_Symbol == j,    Tumor_Sample_Barcode]); tmp2
	tmp = intersect(tmp1, tmp2)
	tmp4
	intersect(tmp, tmp4)
	setdiff(tmp1, tmp2)
	setdiff(tmp2, tmp1)
	length(unique(tcga.maf.skcm$Tumor_Sample_Barcode))
	length(unique(tcga.maf.skcm[(Hugo_Symbol %in% c(i, j)), Tumor_Sample_Barcode]))
	length(unique(tcga.maf.skcm[, Tumor_Sample_Barcode]))
}

## tcga NF1 LOH
{
	tmp = dbGetQuery(con, 'select * from tcga_cnv where Hugo_Symbol ="NF1" and  [.id] = "SKCM"')
	tmp = as.data.table(tmp)
	tmp = tmp[FACETS_CALL != 'DIPLOID', ]
	tmp[, cnv:='GAIN']
	tmp[grep('LOSS', FACETS_CALL), cnv := 'LOSS']
	tmp
	tmp2 = dbGetQuery(con, 'select * from tcga where Hugo_Symbol ="NF1" and DISEASE = "SKCM"')
	tmp2 = as.data.table(tmp2)
	tmp2 = tmp2[oncogenic %in% onco.sel, ]
	tmp2

	tcga.skcm.clin <- as.data.table(GDCquery_clinic("TCGA-SKCM","clinical"))
	dbWriteTable(con, 'tcga_skcm_clin', tcga.skcm.clin)

	# add cnv
	tmp3 = merge(tcga.skcm.clin, tmp[, .(patient, cnv)], by.x = 'bcr_patient_barcode', by.y = 'patient', all.x = T)
	# add mutation
	tmp4 = merge(tmp3, tmp2[, .(PATIENT_ID, Hugo_Symbol)] , by.x = 'bcr_patient_barcode', by.y = 'PATIENT_ID', all.x = T)
	tmp4[is.na(cnv), cnv := 'DIPLOID']
	tmp4[cnv != 'DIPLOID', mut := cnv]
	tmp4[!is.na(Hugo_Symbol),  mut := paste0('NF1 Mut')]
	tmp4[is.na(mut), mut := 'WT']
	tmp4[mut == 'GAIN', c('days_to_last_known_disease_status', 'days_to_death', 'days_to_last_follow_up')]
	as.character(tmp4[mut == 'GAIN', 'bcr_patient_barcode'])

	ss = tmp4$bcr_patient_barcode[tmp4$mut=='GAIN']
	xx = dbGetQuery(con, 'select * from tcga where DISEASE = "SKCM"')
	xx = as.data.table(xx)
	xx[, nf1gain := F]
	xx[PATIENT_ID %in% ss, nf1gain := T]
	ss
	unique(xx$PATIENT_ID)
	intersect(xx$PATIENT_ID, ss)
	ss
	table(xx$nf1gain)

	TCGAanalyze_survival(as.data.frame(tmp4), clusterCol='mut', file='res/nf1_survial.pdf', conf.int = F)
	sync('res')

}

## impact skcm co-occurring and mutual exclusive analysis
{
	library(discover)
	# https://github.com/NKI-CCB/DISCOVER
	onco.sel = c('Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic')
	cc.sel = c('Hugo_Symbol', 'Tumor_Sample_Barcode', 'HGVSp_Short', 'oncogenic', 'ccf_Mcopies', 'vaf')

	tmp = maf.skcm[oncogenic %in% onco.sel,]
	mtx = as.data.frame.matrix(table(tmp[oncogenic %in% onco.sel, .(Tumor_Sample_Barcode, Hugo_Symbol)])) # oncogenic
	mtx[mtx >1] = 1
	mtx = t(mtx)
	events = discover.matrix(mtx)
	genes = c('TP53', 'NF1', 'BRAF', 'NRAS', 'PTEN', 'CDKN2A', 'KRAS', 'HRAS', 'MEF1')
	genes.a = c('TP53')
	genes.b = c('NF1', 'BRAF', 'NRAS', 'HRAS', 'KRAS', 'BRAF.V600E')
	pairwise.discover.test(events[genes,], alternative='greater') -> events.gene.co
	pairwise.discover.test(events[genes,], alternative='less') -> events.gene.ex
	sum.co = as.data.frame(events.gene.co, 1)
	sum.ex = as.data.frame(events.gene.ex, 1)
	sum.co$alt = 'co'
	sum.ex$alt = 'ex'
	rbind(sum.co, sum.ex) -> sum.alt
	sum.alt
	write.table(sum.alt, file='impact_discover_co_ex.xls', sep="\t", quote=F)
	sum.alt$db = 'impact'
	sum.alt.impat = sum.alt

	sum.alt$case.gene1 = 0
	sum.alt$case.gene2 = 0
	sum.alt$case.gene12 = 0
	sum.alt$case.gene1o = 0
	sum.alt$case.gene2o = 0
	for(i in 1:nrow(sum.alt)){
		gene1 = sum.alt$gene1[i]
		gene2 = sum.alt$gene2[i]
		#t1 = tcga.maf.skcm[Hugo_Symbol == gene1 & oncogenic %in% onco.sel, Tumor_Sample_Barcode]
		#t2 = tcga.maf.skcm[Hugo_Symbol == gene2 & oncogenic %in% onco.sel, Tumor_Sample_Barcode]
		t1 = tmp[Hugo_Symbol == gene1, Tumor_Sample_Barcode]
		t2 = tmp[Hugo_Symbol == gene2, Tumor_Sample_Barcode]
		sum.alt[i, 'case.gene1'] = length(t1)
		sum.alt[i, 'case.gene2'] = length(t2)
		sum.alt[i, 'case.gene12'] = length(intersect(t1, t2))
		sum.alt[i, 'case.gene1o'] = length(setdiff(t1, t2))
		sum.alt[i, 'case.gene2o'] = length(setdiff(t2, t1))
	}
	sum.alt$case.total = length(unique(maf$Tumor_Sample_Barcode))
	sum.alt

	sum.alt.onco.all = sum.alt
	sum.alt.onco.sel = sum.alt
	sum.alt.onco.all[1:3,]
	sum.alt.onco.sel[1:3,]

	## discover tcga.skcm, see below
	tmp = tcga.skcm[oncogenic %in% onco.sel,]
	mtx = as.data.frame.matrix(table(tmp[oncogenic %in% onco.sel, .(Tumor_Sample_Barcode, Hugo_Symbol)])) # oncogenic
	mtx[mtx >1] = 1
	mtx = t(mtx)
	events = discover.matrix(mtx)
	genes = c('TP53', 'NF1', 'BRAF', 'NRAS', 'PTEN', 'CDKN2A', 'KRAS', 'HRAS')
	genes.a = c('TP53', 'NF1')
	genes.b = c('NF1', 'BRAF', 'NRAS', 'HRAS', 'KRAS', 'BRAF.V600E')
	pairwise.discover.test(events[genes,], alternative='greater') -> events.gene.co
	pairwise.discover.test(events[genes,], alternative='less') -> events.gene.ex
	sum.co = as.data.frame(events.gene.co, 1)
	sum.ex = as.data.frame(events.gene.ex, 1)
	sum.co$alt = 'co'
	sum.ex$alt = 'ex'
	rbind(sum.co, sum.ex) -> sum.alt.tcga
	sum.alt.tcga$db = 'tcga'
	rbind(sum.alt.impat, sum.alt.tcga) -> sum.alt
	sum.alt

	write.table(sum.alt, file='res/impact_tcga_discover_co_ex.xls', sep="\t", quote=F)
	sync('res')
}

## impact skcm co-occurring and mutual exclusive analysis
## Melanoma only
## NF1 mutated and their other accompannied mutations, 
## impact_skcm_nf1_percentage.pdf
{
	source('~/program/fun/param.r')
	maffile = 'msk_impact_2017_oncokb/data_mutations_merged_maf2maf_facets_clin_oncokb_cleaned.txt'
	fread(maffile) -> maf
	maf
	maf[, vaf := t_alt_count / t_depth]
	maf = maf[vaf > 0.05 & t_alt_count > 3, ]
	onco.sel = c('Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic')
	cc.sel = c('Hugo_Symbol', 'Tumor_Sample_Barcode', 'HGVSp_Short', 'oncogenic', 'ccf_Mcopies', 'vaf')
	maf.skcm = maf[CANCER_TYPE == 'Melanoma', ] 
	maf.skcm[grep('MAP2K1', Hugo_Symbol), cc.1, with=F]

	maf.skcm.nf1 = maf.skcm[Hugo_Symbol == 'NF1' & oncogenic %in% onco.sel, ]
	genes = c('NF1',  'TP53', 'BRAF','CDKN2A',  'NRAS', 'PTEN', 'HRAS', 'KRAS', 'MAP2K1')
	genes.2 = c('TP53', 'BRAF','CDKN2A',  'NRAS', 'PTEN', 'HRAS', 'KRAS', 'MAP2K1')
	genes.3 = c('NF1',  'TP53', 'BRAF', 'BRAF V600E', 'CDKN2A',  'NRAS', 'PTEN', 'HRAS', 'KRAS', 'MAP2K1')
	maf.skcm.nf1.genes = rbind(maf.skcm.nf1, maf.skcm[Tumor_Sample_Barcode %in% maf.skcm.nf1$Tumor_Sample_Barcode & Hugo_Symbol %in% genes.2 & oncogenic %in% onco.sel, ])
	maf.skcm.nf1.genes[Hugo_Symbol =='BRAF', cc.1, with=F]
	maf.skcm
	maf.skcm.nf1.genes[, Hugo_Symbol := factor(Hugo_Symbol, levels=genes)]
	## impact barplot
	maf.skcm[Hugo_Symbol == 'NF1' & oncogenic %in% onco.sel, .(Hugo_Symbol, Tumor_Sample_Barcode, HGVSp_Short, oncogenic)]
	maf.skcm.nf1.genes[, Hugo_Symbol := factor(Hugo_Symbol, levels=genes.3)]
	maf.skcm.nf1.genes
	gg = ggplot(maf.skcm.nf1.genes, aes(x = Hugo_Symbol)) + 
		geom_bar(aes(y = (..count..)/176), position = position_dodge(width = 0.8)) +
		scale_y_continuous(labels=percent) +  ylab('Accompaned mutation rate') +
		theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + xlab('Accompanied Mutations with NF1')
	ggsave(gg, file='res/impact_skcm_nf1_precentage.pdf', width=3.5, height=3)
	sync('res')
}

## impact skcm co-occurring and mutual exclusive analysis
## impact P53 mutation rate within each mutations 
## res/impact_skcm_tp53_precentage.pdf
{
	source('~/program/fun/param.r')
	maffile = 'msk_impact_2017_oncokb/data_mutations_merged_maf2maf_facets_clin_oncokb_cleaned.txt'
	fread(maffile) -> maf
	maf[, vaf := t_alt_count / t_depth]
	maf = maf[vaf > 0.05 & t_alt_count > 3, ]
	onco.sel = c('Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic')
	cc.sel = c('Hugo_Symbol', 'Tumor_Sample_Barcode', 'HGVSp_Short', 'oncogenic', 'ccf_Mcopies', 'vaf')

	maf.skcm = maf[CANCER_TYPE == 'Melanoma', ] 
	genes = c('NF1',  'BRAF', 'CDKN2A',  'NRAS', 'PTEN', 'HRAS', 'KRAS')
	genes = c('NF1',  'TP53', 'BRAF','CDKN2A',  'NRAS', 'PTEN', 'HRAS', 'KRAS', 'MAP2K1')
	tmp = data.table(genes = genes, genemut=rep(0, length(genes)), p53mut=rep(0, length(genes)))
	for(i in 1:nrow(tmp)){
		bcr1 = unique(maf.skcm[Hugo_Symbol == tmp[i, genes] & oncogenic %in% onco.sel, Tumor_Sample_Barcode])
		bcr2 = unique(maf.skcm[Tumor_Sample_Barcode %in% bcr1 & Hugo_Symbol == 'TP53' & oncogenic %in% onco.sel, Tumor_Sample_Barcode])
		tmp[i, genemut := length(bcr1)]
		tmp[i, p53mut := length(bcr2)]
	}
	tmp[, genes := paste0(genes, "(", genemut, ")")]
	tmp[, r := p53mut / genemut] 
	tmp[, genes := factor(genes, levels=genes)]
	gg = ggplot(tmp, aes(x = genes, y = r)) + geom_bar(stat='identity', position = position_dodge(width = 0.8)) + scale_y_continuous(labels=percent) +  ylab('Companied TP53 mutations') +
		theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + xlab('First mutations (total cases)')
	ggsave(gg, file='res/impact_skcm_tp53_precentage.pdf', width=3.5, height=3)
	sync('res')

	unique(maf[CANCER_TYPE == 'Melanoma', CANCER_TYPE_DETAILED] )
	unique(maf.skcm$Tumor_Sample_Barcode)

}

## impact skcm co-occurring and mutual exclusive analysis
## impact oncoplot, use all NF1 mutations with oncogenic unknown
{
	maf.skcm.nf1 = maf.skcm[Hugo_Symbol == 'NF1', ]
	maf.skcm.nf1.genes = rbind(maf.skcm.nf1, maf.skcm[Tumor_Sample_Barcode %in% maf.skcm.nf1$Tumor_Sample_Barcode & Hugo_Symbol %in% genes.2 & oncogenic %in% onco.sel,])
	maf.skcm.nf1.genes[, Hugo_Symbol := factor(Hugo_Symbol, levels=genes)]
	maf.skcm.nf1.genes[Hugo_Symbol == 'BRAF', cc.2, with=F]
	maf.skcm.nf1.genes[Hugo_Symbol == 'BRAF' & HGVSp_Short == 'p.V600E', Hugo_Symbol := 'BRAF V600E'] 
	maf.skcm.nf1.genes[Hugo_Symbol == 'BRAF V600E', cc.1, with=F]
	anno = data.table(Tumor_Sample_Barcode = unique(maf.skcm.nf1.genes$Tumor_Sample_Barcode))
	tmp = maf.skcm[Hugo_Symbol == 'NF1' & oncogenic %in% onco.sel, Tumor_Sample_Barcode] # onco
	anno[, NF1 := 'passage']
	anno[Tumor_Sample_Barcode %in% tmp, NF1 := 'oncogenic']
	tmp = maf.skcm[Hugo_Symbol == 'BRAF' & HGVSp_Short == 'p.V600E', Tumor_Sample_Barcode]
	anno[, BRAF := 'others']
	anno[Tumor_Sample_Barcode %in% tmp, BRAF := 'V600E']
	maf.skcm.nf1.genes[, Hugo_Symbol := factor(Hugo_Symbol, levels=genes.3)]
	source('~/program/fun/cbioMaf.R')
	tmp = cbioMaf(maf.skcm.nf1.genes)
	tmp$maf$Variant_Classification
	tmp.o = read.maf(tmp$maf, vc_nonSyn = tmp$vc.key, clinicalData=anno)
	pdf('res/impact_skcm_nf1_oncoplot.pdf', width=20, height=5)
	oncoplot(tmp.o, genes = genes.3, keepGeneOrder=T, color = tmp$color.table, GeneOrderSort = F, drawRowBar=F, drawColBar=F, sortByMutation=T, writeMatrix=T, sortByAnnotation=T, annotationDat = anno, clinicalFeatures = c('NF1', 'BRAF'), showTumorSampleBarcodes=T, SampleNamefontSize=6)
	dev.off()
	sync('res')

	maf.skcm.nf1.genes -> tmp
	xx = read.table('onco_matrix.txt', header=T, sep="\t")
	xx
	colnames(xx)
	ss = 'P-0019697-T01-IM6'
	maf.skcm[Tumor_Sample_Barcode == ss, cc.2, with=F]
	mean(table(maf.skcm$Tumor_Sample_Barcode))
	maf.skcm.nf1.genes[grep('196', Tumor_Sample_Barcode), c(cc.2, 'vaf'), with=F]
	maf.skcm.nf1.genes[Hugo_Symbol == 'NF1' & HGVSp_Short == 'p.R440*', cc.1, with=F]
}

## impact skcm NF1 related co-occurring and mutual exclusive analysis
## impact lollipop plot of NF1
{
	unique(maf.skcm$oncogenic)
	cc.sel = unique(names(maf.skcm))
	tmp = maf.skcm[Hugo_Symbol %in% c('NF1') & oncogenic == '', cc.sel, with=F]
	tmp[, cc.2, with=F]
	tmp.o = read.maf(tmp)
	pdf('res/impact_skcm_nf1_Nonco_lollipop.pdf', width=6, height=3)
	lollipopPlot(tmp.o, gene='NF1', cBioPortal=T)
	dev.off()
	tmp = maf.skcm[Hugo_Symbol %in% c('NF1') & oncogenic %in% onco.sel, ]
	tmp[, cc.2, with=F]
	tmp.o = read.maf(tmp)
	pdf('res/impact_skcm_nf1_onco_lollipop.pdf', width=6, height=3)
	lollipopPlot(tmp.o, gene='NF1', cBioPortal=T, repel=T)
	dev.off()
	sync('res')

	table(tmp$HGVSp_Short)

}

## still skcm, impact
## impact lollipop for BRAF
{
	aa.length = fread('~/program/cdd/Homo_sapiens.GRCh38.pep.all.aa.len', header=F)
	colnames(aa.length) = c('ENSP', 'aa.length')
	setkey(aa.length, 'ENSP')
	tmp1 = maf.skcm[Hugo_Symbol %in% c('NF1') & oncogenic == '', Tumor_Sample_Barcode]
	tmp2 = maf.skcm[Hugo_Symbol %in% c('NF1') & oncogenic %in% onco.sel, Tumor_Sample_Barcode]
	tmp2
	tmp1 = setdiff(tmp1, tmp2)
	tmp1.maf = maf.skcm[Hugo_Symbol == 'BRAF' & Tumor_Sample_Barcode %in% tmp1 & oncogenic %in% onco.sel, cc.sel, with=F] ## N onco
	tmp2.maf = maf.skcm[Hugo_Symbol == 'BRAF' & Tumor_Sample_Barcode %in% tmp2 & oncogenic %in% onco.sel, cc.sel, with=F] ##  onco
	tmp1.maf[, loc := sub(".*?([0-9]+).*", "\\1", HGVSp_Short)]
	tmp2.maf[, loc := sub(".*?([0-9]+).*", "\\1", HGVSp_Short)]
	setkey(tmp1.maf, 'ENSP')
	setkey(tmp2.maf, 'ENSP')
	tmp1.maf = merge(tmp1.maf, aa.length, by='ENSP')
	tmp2.maf = merge(tmp2.maf, aa.length, by='ENSP')
	tmp1.o = read.maf(tmp1.maf)
	tmp2.o = read.maf(tmp2.maf)
	fn = paste0('res/impact_skcm_Nonco_nf1_braf_lollipop.pdf')
	pdf(fn, width=7, height=4)
	lollipopPlot(tmp1.o, gene='BRAF', cBioPortal=T, labelPos=T)
	dev.off()
	fn = paste0('res/impact_skcm_onco_nf1_braf_lollipop.pdf')
	pdf(fn, width=7, height=4)
	lollipopPlot(tmp2.o, gene='BRAF', cBioPortal=T)
	dev.off()
	sync('res')

	table(tmp1.maf$HGVSp_Short)
	table(tmp2.maf$HGVSp_Short)

	# D594N
	t1 = tmp2.maf[HGVSp_Short == 'p.D594N', Tumor_Sample_Barcode]
	maf.skcm[Tumor_Sample_Barcode %in% t1 & Hugo_Symbol == 'TP53' & oncogenic %in% onco.sel, cc.1, with=F]

}

## tcga barplot
## pancancer melanoma
## pan cancer import data
{
	source('~/program/fun/param.r')
	maffile = '/home/ang46/luna/ext/mafs/mc3/old/mc3/tcga_mc3_20171217_oncokb_facets.tsv'
	tcga.maf = fread(maffile)
	cmd = paste0('grep SKCM ', maffile, ' > tcga_skcm.maf ')
	system(cmd)
	fread('tcga_skcm.maf') -> tcga.maf
	tmp = fread('hd.txt')
	colnames(tcga.maf) = colnames(tmp)
	head(tcga.maf)
	tcga.maf[, vaf := t_alt_count / t_depth]
	tcga.skcm = tcga.maf[DISEASE == 'SKCM' & vaf > 0.05 & t_alt_count > 3, ]
	save(tcga.skcm, file="tcga_skcm.RData")
	load('tcga_skcm.RData')
	rm(tcga.maf)
	rm(tmp)

	tcga.skcm = tcga.skcm[, !duplicated(colnames(tcga.skcm)), with=F]
}

## skcm tcga 
## choose oncogenic for barplot
## res/tcga_skcm_nf1_precentage_v2.pdf
{
	tcga.skcm[Hugo_Symbol == 'MAP2K1', cc.sel, with=F]
	genes = c('NF1',  'TP53', 'BRAF','CDKN2A',  'NRAS', 'PTEN', 'HRAS', 'KRAS', 'MAP2K1')
	genes.2 = c('TP53', 'BRAF','CDKN2A',  'NRAS', 'PTEN', 'HRAS', 'KRAS', 'NF1', 'MAP2K1')
	genes.3 = c('NF1',  'TP53', 'BRAF', 'BRAF V600E', 'CDKN2A',  'NRAS', 'PTEN', 'HRAS', 'KRAS', 'MAP2K1')
	tcga.skcm.nf1 = tcga.skcm[Hugo_Symbol == 'NF1' & oncogenic %in% onco.sel, ]
	tcga.skcm.nf1.genes = tcga.skcm[Tumor_Sample_Barcode %in% tcga.skcm.nf1$Tumor_Sample_Barcode & Hugo_Symbol %in% genes.2 & oncogenic %in% onco.sel, ]
	tcga.skcm.nf1.genes[Hugo_Symbol =='BRAF', cc.1, with=F]
	tcga.skcm.nf1.genes[, Hugo_Symbol := factor(Hugo_Symbol, levels=genes)]
	## tcga barplot
	length(unique(tcga.skcm.nf1.genes$Tumor_Sample_Barcode))
	length(unique(tcga.skcm.nf1$Tumor_Sample_Barcode))
	length(tcga.skcm.nf1$Hugo_Symbol)
	tcga.skcm[Hugo_Symbol == 'NF1' & oncogenic %in% onco.sel, .(Hugo_Symbol, Tumor_Sample_Barcode, HGVSp_Short, oncogenic)]
	tcga.skcm.nf1.genes[, Hugo_Symbol := factor(Hugo_Symbol, levels=genes.3)]
	tcga.skcm.nf1.genes
	gg = ggplot(tcga.skcm.nf1.genes, aes(x = Hugo_Symbol)) + 
		geom_bar(aes(y = (..count..)/47), position = position_dodge(width = 0.8)) +
		scale_y_continuous(labels=percent) +  ylab('Accompanied mutation rate') +
		theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + xlab('Accompanied Mutations with NF1')
	ggsave(gg, file='res/tcga_skcm_nf1_precentage_v2.pdf', width=3.5, height=3)
	sync('res')

	tcga.skcm[Hugo_Symbol == 'NF1' & oncogenic %in% onco.sel, .(Hugo_Symbol, Tumor_Sample_Barcode, HGVSp_Short, oncogenic)]
	tcga.skcm.nf1.genes[, Hugo_Symbol := factor(Hugo_Symbol, levels=genes.3)]
	gg = ggplot(tcga.skcm.nf1.genes, aes(x = Hugo_Symbol)) + 
		geom_bar(aes(y = (..count..)/176), position = position_dodge(width = 0.8)) +
		scale_y_continuous(labels=percent) +  ylab('Accompaned mutation rate') +
		theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + xlab('Accompanied Mutations with NF1')
	ggsave(gg, file='res/tcga_skcm_nf1_precentage.pdf', width=3.5, height=3)
	sync('res')
}

## tcga for p53 mutant cases, look them for other gene mutations
## res/tcga_skcm_tp53_precentage.pdf
{
	source('~/program/fun/param.r')
	genes = c('NF1',  'BRAF', 'CDKN2A',  'NRAS', 'PTEN', 'HRAS', 'KRAS', 'MAP2K1')
	tmp = data.table(genes = genes, genemut=rep(0, length(genes)), p53mut=rep(0, length(genes)))
	for(i in 1:nrow(tmp)){
		bcr1 = unique(tcga.skcm[Hugo_Symbol == tmp[i, genes] & oncogenic %in% onco.sel, Tumor_Sample_Barcode])
		bcr2 = unique(tcga.skcm[Tumor_Sample_Barcode %in% bcr1 & Hugo_Symbol == 'TP53' & oncogenic %in% onco.sel, Tumor_Sample_Barcode])
		tmp[i, genemut := length(bcr1)]
		tmp[i, p53mut := length(bcr2)]
	}
	tmp
	tmp[, genes := paste0(genes, "(", genemut, ")")]
	tmp[, r := p53mut / genemut] 
	tmp[, genes := factor(genes, levels=genes)]
	gg = ggplot(tmp, aes(x = genes, y = r)) + geom_bar(stat='identity', position = position_dodge(width = 0.8)) + scale_y_continuous(labels=percent) +  ylab('Companied TP53 mutations') +
		theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + xlab('First mutations (total cases)')
	ggsave(gg, file='res/tcga_skcm_tp53_precentage.pdf', width=3.5, height=3)
	sync('res')

	genes = c('NF1',  'TP53', 'BRAF','CDKN2A',  'NRAS', 'PTEN', 'HRAS', 'KRAS')
	genes.2 = c('NF1', 'TP53', 'BRAF','CDKN2A',  'NRAS', 'PTEN', 'HRAS', 'KRAS')
	genes.3 = c('TP53', 'NF1',  'BRAF', 'BRAF V600E', 'CDKN2A',  'NRAS', 'PTEN', 'HRAS', 'KRAS')
	tcga.skcm.nf1 = tcga.skcm[Hugo_Symbol == 'TP53' & oncogenic %in% onco.sel, ]
	tcga.skcm.nf1$Tumor_Sample_Barcode
	tcga.skcm.nf1.genes = tcga.skcm[Tumor_Sample_Barcode %in% tcga.skcm.nf1$Tumor_Sample_Barcode & Hugo_Symbol %in% genes.2 & oncogenic %in% onco.sel, ]
	tcga.skcm.nf1.genes
	tcga.skcm.nf1.genes[Hugo_Symbol =='BRAF', cc.1, with=F]
	tcga.skcm.nf1.genes[, Hugo_Symbol := factor(Hugo_Symbol, levels=genes)]
	## tcga barplot
	length(unique(tcga.skcm.nf1.genes$Tumor_Sample_Barcode))
	length(unique(tcga.skcm.nf1$Tumor_Sample_Barcode))
	length(tcga.skcm.nf1$Hugo_Symbol)
	tcga.skcm.nf1.genes[, Hugo_Symbol := factor(Hugo_Symbol, levels=genes.3)]
	tcga.skcm.nf1.genes$Hugo_Symbol
	colnames(tcga.skcm.nf1.genes)
	gg = ggplot(tcga.skcm.nf1.genes[, cc.1, with=F], aes(x = Hugo_Symbol)) + 
		geom_bar(aes(y = (..count..)/65), position = position_dodge(width = 0.8)) +
		scale_y_continuous(labels=percent) +  ylab('Case percentage') +
		theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) + xlab('')
	ggsave(gg, file='res/tcga_skcm_tp53_precentage.pdf', width=3.5, height=3)
	tcga.skcm.nf1.genes
	sync('res')
}

## tcga oncoplot
{
	tcga.skcm.nf1 = tcga.skcm[Hugo_Symbol == 'NF1', ]
	tcga.skcm.nf1.genes = rbind(tcga.skcm.nf1, tcga.skcm[Tumor_Sample_Barcode %in% tcga.skcm.nf1$Tumor_Sample_Barcode & Hugo_Symbol %in% genes.2 & oncogenic %in% onco.sel,])
	tcga.skcm.nf1.genes[, Hugo_Symbol := factor(Hugo_Symbol, levels=genes)]
	tcga.skcm.nf1.genes[Hugo_Symbol == 'BRAF', cc.2, with=F]
	tcga.skcm.nf1.genes[Hugo_Symbol == 'BRAF' & HGVSp_Short == 'p.V600E', Hugo_Symbol := 'BRAF V600E'] 
	tcga.skcm.nf1.genes[Hugo_Symbol == 'BRAF V600E', cc.1, with=F]
	anno = data.table(Tumor_Sample_Barcode = unique(tcga.skcm.nf1.genes$Tumor_Sample_Barcode))
	tmp = tcga.skcm[Hugo_Symbol == 'NF1' & oncogenic %in% onco.sel, Tumor_Sample_Barcode]
	anno[, NF1 := 'passage']
	anno[Tumor_Sample_Barcode %in% tmp, NF1 := 'oncogenic']
	tmp = tcga.skcm[Hugo_Symbol == 'BRAF' & HGVSp_Short == 'p.V600E', Tumor_Sample_Barcode]
	anno[, BRAF := 'others']
	anno[Tumor_Sample_Barcode %in% tmp, BRAF := 'V600E']
	tcga.skcm.nf1.genes[, Hugo_Symbol := factor(Hugo_Symbol, levels=genes.3)]
	source('~/program/fun/cbioMaf.R')
	tcga.skcm.nf1.genes = tcga.skcm.nf1.genes[, unique(names(tcga.skcm.nf1.genes)), with=F] 
	tmp = cbioMaf(tcga.skcm.nf1.genes)
	tmp$maf$Variant_Classification
	tmp.o = read.maf(tmp$maf, vc_nonSyn = tmp$vc.key, clinicalData=anno)
	pdf('res/tcga_skcm_nf1_oncoplot.pdf', width=10, height=6)
	oncoplot(tmp.o, genes = genes.3, keepGeneOrder=T, color = tmp$color.table, GeneOrderSort = F, drawRowBar=F, drawColBar=F, sortByMutation=T, writeMatrix=T, sortByAnnotation=T, annotationDat = anno, clinicalFeatures = c('NF1', 'BRAF'), showTumorSampleBarcodes=T, SampleNamefontSize=6)
	dev.off()
	sync('res')
}

## tcga lollipop plot
{
	unique(tcga.skcm$oncogenic)
	cc.sel = unique(names(tcga.skcm))
	tmp = tcga.skcm[Hugo_Symbol %in% c('NF1') & oncogenic == '', cc.sel, with=F]
	tmp[, cc.2, with=F]
	tmp.o = read.maf(tmp)
	pdf('res/tcga_skcm_nf1_Nonco_lollipop.pdf', width=6, height=3)
	lollipopPlot(tmp.o, gene='NF1', cBioPortal=T)
	dev.off()
	tmp = tcga.skcm[Hugo_Symbol %in% c('NF1') & oncogenic %in% onco.sel, cc.sel, with=F]
	tmp[, cc.2, with=F]
	tmp.o = read.maf(tmp)
	pdf('res/tcga_skcm_nf1_onco_lollipop.pdf', width=6, height=3)
	lollipopPlot(tmp.o, gene='NF1', cBioPortal=T)
	dev.off()
	sync('res')

	source('~/program/maftools/R/lollipopPlot_ngenes_vaf.R')
	aa.length = fread('~/program/cdd/Homo_sapiens.GRCh38.pep.all.aa.len', header=F)
	colnames(aa.length) = c('ENSP', 'aa.length')
	setkey(aa.length, 'ENSP')
	tmp1 = tcga.skcm[Hugo_Symbol %in% c('NF1') & oncogenic == '', Tumor_Sample_Barcode]
	tmp2 = tcga.skcm[Hugo_Symbol %in% c('NF1') & oncogenic %in% onco.sel, Tumor_Sample_Barcode]
	tmp1 = setdiff(tmp1, tmp2)
	tmp1.maf = tcga.skcm[Hugo_Symbol == 'BRAF' & Tumor_Sample_Barcode %in% tmp1, cc.sel, with=F] ## N onco
	tmp2.maf = tcga.skcm[Hugo_Symbol == 'BRAF' & Tumor_Sample_Barcode %in% tmp2, cc.sel, with=F] ##  onco
	tmp1.maf[, loc := sub(".*?([0-9]+).*", "\\1", HGVSp_Short)]
	tmp2.maf[, loc := sub(".*?([0-9]+).*", "\\1", HGVSp_Short)]
	setkey(tmp1.maf, 'ENSP')
	setkey(tmp2.maf, 'ENSP')
	tmp1.maf = merge(tmp1.maf, aa.length, by='ENSP')
	tmp2.maf = merge(tmp2.maf, aa.length, by='ENSP')
	lollipopPlot_ngenes_vaf(tmp1.maf, fn = paste0('res/tcga_skcm_Nonco_nf1_braf_lollipop.pdf'), ww=7, hh=4, header='Non oncogenic NF1 mutated cases')
	lollipopPlot_ngenes_vaf(tmp2.maf, fn = paste0('res/tcga_skcm_onco_nf1_braf_lollipop.pdf'), ww=7, hh=4, header='Oncogenic NF1 mutated cases')
	sync('res')
}

## impact co-occurring and mutual exclusive analysis
{
	maf.skcm.nf1.genes -> tmp
	xx = read.table('onco_matrix.txt', header=T, sep="\t")
	xx
	colnames(xx)
	ss = 'P-0019697-T01-IM6'
	maf.skcm[Tumor_Sample_Barcode == ss, cc.2, with=F]
	mean(table(maf.skcm$Tumor_Sample_Barcode))
	maf.skcm.nf1.genes[grep('196', Tumor_Sample_Barcode), c(cc.2, 'vaf'), with=F]
	maf.skcm.nf1.genes[Hugo_Symbol == 'NF1' & HGVSp_Short == 'p.R440*', cc.1, with=F]


	pairwise.discover.test(events[genes,], alternative='greater') -> events.gene.co
	pairwise.discover.test(events[genes,], alternative='less') -> events.gene.ex
	sum.co = as.data.frame(events.gene.co, 1)
	sum.ex = as.data.frame(events.gene.ex, 1)
	sum.co$alt = 'co'
	sum.ex$alt = 'ex'
	rbind(sum.co, sum.ex) -> sum.alt
	as.data.table(sum.alt) -> sum.alt

	sum.alt$case.gene1 = 0
	sum.alt$case.gene2 = 0
	sum.alt$case.gene12 = 0
	sum.alt$case.gene1o = 0
	sum.alt$case.gene2o = 0
	for(i in 1:nrow(sum.alt)){
		gene1 = sum.alt$gene1[i]
		gene2 = sum.alt$gene2[i]
		#t1 = tcga.maf.skcm[Hugo_Symbol == gene1 & oncogenic %in% onco.sel, Tumor_Sample_Barcode]
		#t2 = tcga.maf.skcm[Hugo_Symbol == gene2 & oncogenic %in% onco.sel, Tumor_Sample_Barcode]
		t1 = tmp[Hugo_Symbol == gene1, Tumor_Sample_Barcode]
		t2 = tmp[Hugo_Symbol == gene2, Tumor_Sample_Barcode]
		sum.alt[i, 'case.gene1'] = length(t1)
		sum.alt[i, 'case.gene2'] = length(t2)
		sum.alt[i, 'case.gene12'] = length(intersect(t1, t2))
		sum.alt[i, 'case.gene1o'] = length(setdiff(t1, t2))
		sum.alt[i, 'case.gene2o'] = length(setdiff(t2, t1))
	}
	sum.alt$case.total = length(unique(maf.skcm$Tumor_Sample_Barcode))
	sum.alt = as.data.table(sum.alt)
	sum.alt
	sum.alt[grep("53", gene2),]

	write.table(sum.alt[grep("53", gene2),] , file='P53_vs_others_sum_impact.xls', sep="\t", quote=F, append=T, row.names=F, col.names=T)
	sync()

	tmp = sum.alt[grep("53", gene2),]; tmp
	tmp = tmp[, .(gene1, gene2)]
	ggplot(tmp, aes(y=gene1

			sum.alt[c(grep('NF1', gene1), grep('NF1', gene2)),]

			tmp = sum.alt[, .(gene1, gene2, case.gene12)]
			tmp = rbind(tmp, tmp[, c(2,1,3)])
			tmp
			sum.alt2 = dcast.data.table(sum.alt, gene1 ~ gene2, value.var = 'case.gene12', fun.aggregate=function(x){mean(x)}, fill=0)
			sum.alt2 = sum.alt2[gene1 %in% genes.a, c('gene1', genes.b), with=F]
			sum.alt2
			write.table(sum.alt2, file='impact_nf1.xls')
			sync()

			maf.skcm[Hugo_Symbol == 'NF1', .(Hugo_Symbol, HGVSp_Short, Tumor_Sample_Barcode)]
			maf.skcm[Hugo_Symbol == 'BRAF' & oncogenic %in% onco.sel, .(Hugo_Symbol, HGVSp_Short, Tumor_Sample_Barcode)]
			xx1 = maf.skcm[Hugo_Symbol == 'BRAF' & oncogenic %in% onco.sel, .(HGVSp_Short, Tumor_Sample_Barcode)]
			table(maf.skcm[Hugo_Symbol == 'BRAF' & oncogenic %in% onco.sel, .(HGVSp_Short)])

			sum.alt.onco.all = sum.alt
			sum.alt.onco.sel = sum.alt
			sum.alt.onco.all[1:3,]
			sum.alt.onco.sel[1:3,]

			## apply dentrix
			## mutation file
			ids = unique(maf.skcm$Tumor_Sample_Barcode)
			for(i in ids){
				tmp = maf.skcm[Tumor_Sample_Barcode == i & oncogenic %in% onco.sel, Hugo_Symbol]
				if(length(tmp) > 0){
					cat(i, unique(tmp), "\n", sep="\t", file='impact_skcm_detrix_input.tsv', append=T)
				}
			}

			# python Dendrix.py mutations <- file K minFreqGene number_interations analyzed_genes_file num_exper step_length
			cmd = "~/local/bin/python ~/program/Dendrix/Dendrix.py impact_skcm_detrix_input.tsv 3 3 1000000 analyzed_genes_dentrix.txt 1 10"
			system(cmd)

			#row.sel = rowSums(mtx) > 30
			events.sel = events[row.sel,]
			events.sel[1:10,1:3]
			events[genes,]
			events.sel[genes,]
			dim(events.sel)
			pairwise.discover.test(events.sel, alternative='less') -> events.pair.ex
			pairwise.discover.test(events.sel, alternative='greater') -> events.pair.co
			events.pair = events.pair.ex
			print(events.pair, fdr.threshold=0.05)
			as.data.frame(events.pair, 0.2)
			events.pair = events.pair.co
			print(events.pair, fdr.threshold=0.2)
			as.data.frame(events.pair, 0.8)
			events.pair

			x1 = maf[Hugo_Symbol == 'NF1' & oncogenic %in% onco.sel, Tumor_Sample_Barcode]
			x2 = maf[Hugo_Symbol == 'TP53' & oncogenic %in% onco.sel, Tumor_Sample_Barcode]
			x3 = maf[Hugo_Symbol == 'BRAF'  & oncogenic %in% onco.sel, Tumor_Sample_Barcode]
			intersect(x1, intersect(x2, x3))

			xx1 = maf[Hugo_Symbol == 'TP53' & oncogenic %in% onco.sel, Tumor_Sample_Barcode]
			xx2 = maf[Hugo_Symbol == 'BRAF' & HGVSp_Short == 'p.V600E' & oncogenic %in% onco.sel, Tumor_Sample_Barcode]
			xx3 = maf[Hugo_Symbol == 'BRAF' & oncogenic %in% onco.sel, Tumor_Sample_Barcode]
			intersect(xx1, xx2)
			xx1
			xx1

			maf.skcm = maf[CANCER_TYPE == 'Melanoma', ] 
			x1 = maf.skcm[Hugo_Symbol == 'NF1' , Tumor_Sample_Barcode]
			x2 = maf.skcm[Hugo_Symbol == 'TP53', Tumor_Sample_Barcode]
			x3 = maf.skcm[Hugo_Symbol == 'BRAF' , Tumor_Sample_Barcode]
			triple.id = intersect(x1, intersect(x2, x3))
			P53wNf1 = intersect(x1, x2)
			P53wBraf = intersect(x2, x3)
			id = data.table(sampleid = unique(c(triple.id, P53wNf1, P53wBraf)))
			id[, NF1 := 'No mutations']
			id[, P53 := 'No mutations']
			id[, BRAF := 'No mutations']
			id[ sampleid %in% maf.skcm[Hugo_Symbol == 'NF1' & oncogenic %in% onco.sel, Tumor_Sample_Barcode], NF1 := 'oncogenic']
			id[ sampleid %in% maf.skcm[Hugo_Symbol == 'NF1' & !(oncogenic %in% onco.sel), Tumor_Sample_Barcode], NF1 := 'Not oncogenic']
			id[ sampleid %in% maf.skcm[Hugo_Symbol == 'TP53' & oncogenic %in% onco.sel, Tumor_Sample_Barcode], P53 := 'oncogenic']
			id[ sampleid %in% maf.skcm[Hugo_Symbol == 'TP53' & !(oncogenic %in% onco.sel), Tumor_Sample_Barcode], P53 := 'Not oncogenic']
			id[ sampleid %in% maf.skcm[Hugo_Symbol == 'BRAF' & oncogenic %in% onco.sel, Tumor_Sample_Barcode], BRAF := 'oncogenic']
			id[ sampleid %in% maf.skcm[Hugo_Symbol == 'BRAF' & !(oncogenic %in% onco.sel), Tumor_Sample_Barcode], BRAF := 'Not oncogenic']
			tickle = maf.skcm[ Hugo_Symbol == 'NF1', .(Tumor_Sample_Barcode, HGVSp_Short)]
			id = merge(id, tickle, by.x = 'sampleid', by.y = 'Tumor_Sample_Barcode', all.x=T)
			setnames(id, 'HGVSp_Short', 'NF1.HGVSp_Short')
			tickle = maf.skcm[ Hugo_Symbol == 'TP53', .(Tumor_Sample_Barcode, HGVSp_Short)]
			id = merge(id, tickle, by.x = 'sampleid', by.y = 'Tumor_Sample_Barcode', all.x=T)
			setnames(id, 'HGVSp_Short', 'TP53.HGVSp_Short')
			tickle = maf.skcm[ Hugo_Symbol == 'BRAF', .(Tumor_Sample_Barcode, HGVSp_Short)]
			id = merge(id, tickle, by.x = 'sampleid', by.y = 'Tumor_Sample_Barcode', all.x=T)
			setnames(id, 'HGVSp_Short', 'BRAF.HGVSp_Short')
			id

			write.table(id, file="sample_id.xls", sep="\t", row.names=F)
			sync()

			cat("triple mutation", "\n" file="ids.xls", sep="\t")

			maf.skcm = maf[CANCER_TYPE == 'Melanoma' & oncogenic %in% onco.sel, ] 
			x1 = maf.skcm[Hugo_Symbol == 'NF1', Tumor_Sample_Barcode]
			x2 = maf.skcm[Hugo_Symbol == 'TP53', Tumor_Sample_Barcode]
			x3 = maf.skcm[Hugo_Symbol == 'BRAF' , Tumor_Sample_Barcode]
			intersect(x1, intersect(x2, x3))
			ids = intersect(x1, intersect(x2, x3))
			ids
			cat("triple oncogenic mutations", ids, file="ids.xls", sep="\t", append=T)

			sync()

			write.table(tcga.maf.skcm, file='tcga_skcm.maf', sep="\t", quote=F)
			getwd()



			## cometExactTest to examine the mutual exclusive of NF1, TP53, BRAF
			#	      	           _____m3_______
			#			no		yes
			#      		 _____m1____		_____m1____
			#		no	yes		no	yes
			# m2 no		x000	x001		x100	x101
			# m2 yes	x010	x011		x110	x111
			#
			# c(x000, x001, x010, x011, x100, x101, x110, x111)

			x000 = 363 - length(unique(c(x1, x2, x3))) # 292
			x001 = length(setdiff(setdiff(x1, x2), x3)) # NF1 specific 26
			x010 = length(setdiff(setdiff(x2, x1), x3)) # TP53 specific 20
			x011 = length(setdiff(c(x1, x2), x3)) # not BRAF 140
			x100 = length(setdiff(unique(c(x3, x2)), x1)) # BRAF specific 
			x101 = length(setdiff(unique(c(x3, x1)), x2)) # not TP53 
			x110 = length(setdiff(unique(c(x3, x2)), x1)) # not NF1 
			x111 = length(unique(c(x1, x2, x3))) # either of three 243
			xx = c(x000, x001, x010, x011, x100, x101, x110, x111)
			xx

			library(cometExactTest)
			# pvalthresh=1.1, mutmatplot=T
			comet_exact_test(xx)

			ids = unique(maf.skcm[Hugo_Symbol == 'NF1' & oncogenic %in% onco.sel, Tumor_Sample_Barcode])
			ids
			maf.skcm.nf1 = as.data.table(table(maf.skcm[Tumor_Sample_Barcode %in% ids & oncogenic %in% onco.sel, Hugo_Symbol]))
			maf.skcm.nf1 = maf.skcm.nf1[order(N, decreasing=T),]
			maf.skcm.nf1.top = maf.skcm.nf1[N > 5,]
			maf.skcm.nf1.top

			write.table(maf.skcm.nf1.top, file='res/maf.skcm.nf1.top.txt')
			sync('res')

			sel1 = maf.skcm[Tumor_Sample_Barcode %in% ids & oncogenic %in% onco.sel, ] 
			sel2 = sel1[Hugo_Symbol %in% maf.skcm.nf1.top$V1,]
			sel2$Tumor_Sample_Barcode
			sel2$Hugo_Symbol

			source('~/program/maftools/R/lollipopPlot_ngenes_vaf.R')
			source('~/program/fun/maf_add_loc_aalen.r')
			source('~/program/fun/sync.r')

			aa.length = fread('~/program/cdd/Homo_sapiens.GRCh38.pep.all.aa.len', header=F)
			colnames(aa.length) = c('ENSP', 'aa.length')
			setkey(aa.length, 'ENSP')

			ids2 = unique(sel2$Tumor_Sample_Barcode)
			for(i in 1:length(ids2)){
				sub.maf = sel2[Tumor_Sample_Barcode == ids2[i],]
				sub.maf[, loc := sub(".*?([0-9]+).*", "\\1", HGVSp_Short)]
				setkey(sub.maf, 'ENSP')
				sub.maf = merge(sub.maf, aa.length, by='ENSP')
				lollipopPlot_ngenes_vaf(sub.maf[oncogenic %in% onco.sel & aa.length > 0 & loc > 0, ], fn = paste0('res/impact_skcm_', ids2[i], '_lollipop_bycase.pdf'))
			}
			sync('res')

			ids3 = intersect(intersect(x1, x2), x3)
			for(i in 1:length(ids3)){
				sub.maf = maf.skcm[Tumor_Sample_Barcode == ids3[i] & oncogenic %in% onco.sel,]
				sub.maf[, loc := sub(".*?([0-9]+).*", "\\1", HGVSp_Short)]
				setkey(sub.maf, 'ENSP')
				sub.maf = merge(sub.maf, aa.length, by='ENSP')
				lollipopPlot_ngenes_vaf(sub.maf[oncogenic %in% onco.sel & aa.length > 0 & loc > 0, ], fn = paste0('res/tri/impact_skcm_', ids3[i], '_lollipop_bycase.pdf'))
			}
			sync('res/tri/')


			## for using Dentrix
			## output format:
			## ID	genes 
			gg = c('NF1', 'BRAF', 'TP53')
			ids = unique(tcga.maf.skcm$Tumor_Sample_Barcode)
			for(i in ids){
				tt = tcga.maf.skcm[Tumor_Sample_Barcode == i & Hugo_Symbol %in%  gg, Hugo_Symbol]
				if(length(tt) > 0){
					cat(i, unique(tt), "\n", sep="\t", file='dentrix.txt', append=T)
				}
			}


}

## analysis akt mutations in impact dataset
{
	fread("msk_impact_2017_clinical_data.tsv", header = T) -> clin
	fread("data_mutations_annotated.txt", header = T) -> maf
	fread("data_CNA_width.txt") -> cnv
	colnames(cnv) = c("Hugo_Symbol", "Tumor_Sample_Barcode", "ploidy")
	cnv$Variant_Type = "CNV"
	cnv$Variant_Classification = "Del"
	cnv$Variant_Classification[cnv$ploidy > 0] = "Dup"

	maf_akt = maf[Hugo_Symbol == 'AKT1',]
	cnv_akt = cnv[ Hugo_Symbol == 'AKT1',]

	tmp = read.maf(maf_akt)
	pdf("lollipop_akt1_impact.pdf", width=6, height=3)
	lollipopPlot(tmp, 'AKT1', cBioPortal=T)
	dev.off()
	system("open lollipop_akt1_impact.pdf")

	## ATK1 mut, tumor type (oncotree code)
	merge(maf_akt, clin, by.x = 'Tumor_Sample_Barcode', by.y = 'Sample ID', all.x = T, all.y = F) -> maf_akt2
	tmp2 =table(maf_akt2$"Oncotree Code")
	as.data.frame(tmp2, stringsAsFactors = F) -> tmp2
	tmp2 = tmp2[order(tmp2$Freq, decreasing = T),]
	colnames(tmp2) = c('oncotreeCode', 'akt1_mut')
	head(tmp2)

	tmp =table(clin$"Oncotree Code")
	as.data.frame(tmp) -> tmp
	tmp = tmp[order(tmp$Freq, decreasing = T),]
	row.names(tmp) = tmp$Var1
	tmp2$totalcase = tmp[tmp2$oncotreeCode, 'Freq']
	head(tmp2)
	tmp2

	## ATK1 mut, cancer type
	tmp = unique(clin[,.(`Oncotree Code`, `Cancer Type`)])
	as.data.frame(tmp) -> tmp
	tmp = tmp[!is.na(tmp$"Oncotree Code"),]
	row.names(tmp) = tmp$"Oncotree Code"
	tmp2$CancerType = tmp[tmp2$oncotreeCode, "Cancer Type"]
	head(tmp2)

	## AKT1 E17K
	maf_akt_e17k = maf2[Hugo_Symbol == 'AKT1' & HGVSp_Short == 'p.E17K',]
	tmp =table(maf_akt_e17k$"Oncotree Code")
	as.data.frame(tmp, stringsAsFactors = F) -> tmp
	tmp = tmp[order(tmp$Freq, decreasing = T),]
	row.names(tmp) = tmp$Var1
	colnames(tmp) = c('Var1', 'akt1_e17k')
	tmp2$akt1_e17k = tmp[tmp2$oncotreeCode, 'akt1_e17k']
	head(tmp2)

	tmp2$r = tmp2$akt1_e17k / tmp2$totalcase
	tmp2$r2 = tmp2$akt1_mut / tmp2$totalcase
	tmp2[is.na(tmp2$r), 'r'] = 0
	tmp = tmp2[1:20,]

	akt_sum = tmp2
	row.names(akt_sum) = akt_sum$oncotreeCode
	akt_sum


	## one way 
	pdf("impact_e17k_sum.pdf", width=8, height=7);
	xx = barplot(tmp$akt1_mut, ylim = c(0, 80), width=0.5, col=c('green'), 
		     ylab = 'Case number with AKT1 mutations', main="AKT1 mutations (total case on top)")
	rect(xx - 0.25, 0, xx + 0.25, tmp$akt1_e17k, col='red')
	legend("topright", fill = c('green', 'red'), legend=c('total AKT1 mutations', 'AKT1 E17K'))
	tt = sub("TCGA-", "", tmp$oncotreeCode)
	text(xx, 0, tt, xpd = T, adj=c(1.2, 0.5), srt = 90, font = 7)
	text(xx, tmp$akt1_mut, tmp$totalcase, xpd = T, adj=c(0.5, 0), srt = 0, font = 7)
	text(xx, tmp$akt1_mut + 3, tmp$CancerType, xpd=T, adj=c(0, 0.5), srt = 90, font = 6)
	dev.off()
	system("open impact_e17k_sum.pdf")

	pdf("impact_e17k_sum_percentage.pdf", width=8, height=7)
	xx = barplot(tmp$r2 * 100, ylim = c(0, 40), width=0.5, col='green', ylab = '% AKT1 mutations', main="AKT1 mutations (total case on top)")
	rect(xx - 0.25, 0, xx + 0.25, tmp$r * 100, col='red')
	legend("topright", fill = c('green', 'red'), legend=c('total AKT1 mutations', 'AKT1 E17K'))
	tt = sub("TCGA-", "", tmp$oncotreeCode)
	text(xx, 0, tt, xpd = T, adj=c(1.2, 0.5), srt = 90, font = 7)
	text(xx, tmp$r2 * 100, tmp$totalcase, xpd = T, adj=c(0.5, 0), srt = 0, font = 7)
	aa = as.numeric(apply(tmp[, c('r', 'r2')], 1, max))
	text(xx,  aa*100 + 2, tmp$CancerType, xpd=T, adj=c(0, 0.5), srt = 90, font = 6)
	dev.off()
	system("open impact_e17k_sum_percentage.pdf")


	## another way
	pdf("impact_e17k_sum_f2.pdf", width=14, height=7);
	d = t(as.matrix(tmp[, c('akt1_mut', 'akt1_e17k')]))
	colnames(d) = tmp$oncotreeCode
	d[is.na(d)] = 0
	xx = barplot(d, ylim = c(0, 80), width=0.5, col=c('green', 'red'), las = 3,
		     ylab = 'Case number with AKT1 mutations', main="AKT1 mutations (total case on top)", beside = T)
	legend("topright", fill = c('green', 'red'), legend=c('total AKT1 mutations', 'AKT1 E17K'))
	apply(xx, 2, mean) -> xx
	text(xx, tmp$akt1_mut, tmp$totalcase, xpd = T, adj=c(0.5, -1), srt = 0, font = 7)
	text(xx, tmp$akt1_mut + 9, tmp$CancerType, xpd=T, adj=c(0, 0.5), srt = 90, font = 6)
	dev.off()
	system("open impact_e17k_sum_f2.pdf")

	pdf("impact_e17k_sum_percentage_f2.pdf", width=14, height=7)
	d = t(as.matrix(tmp[, c('r2', 'r')]))
	d = 100 * d
	colnames(d) = tmp$oncotreeCode
	d[is.na(d)] = 0
	xx = barplot(d, ylim = c(0, 25), width=0.5, col=c('green', 'red'), las = 3,
		     ylab = '% AKT1 mutations', main="AKT1 mutations (total case on top)", beside = T)
	legend("topright", fill = c('green', 'red'), legend=c('total AKT1 mutations', 'AKT1 E17K'))
	apply(xx, 2, mean) -> xx
	text(xx, tmp$r2 * 100, tmp$totalcase, xpd = T, adj=c(0.5, -1), srt = 0, font = 7)
	aa = as.numeric(apply(tmp[, c('r', 'r2')], 1, max))
	text(xx, aa * 100 + 2, tmp$CancerType, xpd=T, adj=c(0, 0.5), srt = 90, font = 6)
	dev.off()
	system("open impact_e17k_sum_percentage_f2.pdf")


	tmp = maf[Hugo_Symbol == 'AKT1',]
	read.maf(tmp) -> tmp
	pdf("lollipop_akt1_impact.pdf", width=6, height=3)
	lollipopPlot(tmp, 'AKT1', cBioPortal=T)
	dev.off()
	system("open lollipop_akt1_impact.pdf")

	for( i in 1:20 ){
		id = tmp2$oncotreeCode[i]
		tmp = maf2[Hugo_Symbol == 'AKT1' & `Oncotree Code` == id,]
		read.maf(tmp) -> tmp
		pdffile = paste0("lollipop_akt1_impact_", id, ".pdf")
		pdf(pdffile, width=6, height=3)
		lollipopPlot(tmp, 'AKT1', cBioPortal=T)
		dev.off()
		system(paste0("open ", pdffile))

	}



	cnv_akt_clin = merge(cnv_akt, clin, by.x = "Tumor_Sample_Barcode", by.y = "Sample ID", all.x = T, all.y = F)
	tmp = table(cnv_akt_clin$"Oncotree Code", cnv_akt_clin$Variant_Classification)
	tmp = as.data.frame.matrix(tmp)
	tmp = tmp[order(tmp$Dup, decreasing = T),]
	tmp = tmp[,2:1]
	pdf("impact_akt1_cnv.pdf", width=8, height=5)
	barplot(t(tmp), beside=T, las = 3, col = 2:3, main='AKT1 CNV', ylab = "# of case")
	legend("topright", fill = c('red', 'green'), legend=c('Duplication', 'Deletion'))
	dev.off()
	system("open impact_akt1_cnv.pdf")

	## there are total 34 case of AKT CNV, 4 case deletion, 30 case of duplication
	## there are 6 aptients have E17K and also CNV duplication in UEC (3 case), OPHSC, ILC, LUSC
	intersect(maf_akt$Tumor_Sample_Barcode, cnv_akt$Tumor_Sample_Barcode)
	tmp = maf_akt[maf_akt$Tumor_Sample_Barcode %in% cnv_akt$Tumor_Sample_Barcode[cnv_akt$ploidy > 0],] 
	tmp = maf_akt[maf_akt$Tumor_Sample_Barcode %in% cnv_akt$Tumor_Sample_Barcode,] 
	dim(tmp)
	dim(cnv_akt)


	table(maf_akt2$HGVSp_Short, maf_akt2$"Sample Type")
	#		       Metastasis	Primary
	# p.E17K               69      		42
	tmp = maf_akt2[HGVSp_Short == 'p.E17K',]
	table(tmp$"Sample Type", tmp$"Oncotree Code")

	dim(clin[clin$"Primary Tumor Site" == "Breast", ])
	# 1351 22
	dim(maf2[`Primary Tumor Site` == "Breast", ]) # 6547
	dim(maf2[`Primary Tumor Site` == "Breast" & `HGVSp_Short` == 'p.E17K', ]) #  59
	dim(maf2[`HGVSp_Short` == 'p.E17K', ]) #  119

	clin$"Metastatic Site"[ clin$"Metastatic Site" == "" ] = "NONE"
	table(clin$"Primary Tumor Site", clin$"Metastatic Site") -> t
	as.data.frame.matrix(t) -> t 
	t = t[, order(t['Breast',], decreasing = T)]
	formattable(t)
	t = t['Breast',]

	maf_akt2$"Metastatic Site"[maf_akt2$"Metastatic Site" == "" ] = "NONE"
	table(maf_akt2$"Primary Tumor Site", maf_akt2$"Metastatic Site") -> tt
	as.data.frame.matrix(tt) -> tt
	tt = tt[, order(tt['Breast',], decreasing = T)]
	formattable(tt)
	tt = tt['Breast',]
	wilcox.test(as.numeric(t), as.numeric(tt), paired=F)
	tt

	tn = names(tt)
	df = data.frame(row.names=tn, all_case = as.numeric(t[tn]), akt_case = as.numeric(tt[tn]))
	df = df[-1,]
	df = df[order(df$all_case, decreasing = T),]
	df$site = row.names(df)
	tmp = apply(df[,1:2], 1, sum)
	df = df[tmp > 0,]
	df

	df.glm = glm(cbind(all_case, akt_case) ~ site, data =df, family = binomial())
	anova(df.glm)
	lsmeans(df.glm, pairwise ~ site)
}

## UCEC
{
	unique(maf$CANCER_TYPE)
	tmp = maf[CANCER_TYPE == 'Cervical Cancer',]
	dim(tmp)
	onco.sel = c('Likely Oncogenic', 'Oncogenic', 'Predicted Oncogenic')
	cc.sel = c('Hugo_Symbol', 'Tumor_Sample_Barcode', 'HGVSp_Short', 'oncogenic', 'ccf_Mcopies', 'vaf')
	maf.skcm = maf[CANCER_TYPE == 'Melanoma'] 
	maf.skcm[, onco := 'non']
	maf.skcm[oncogenic %in% onco.sel, onco := 'onco'] 
	maf.skcm[Hugo_Symbol == 'BRAF' & HGVSp_Short == 'p.V600E', Hugo_Symbol := paste0(Hugo_Symbol, ':', HGVSp_Short)]
	maf.skcm[, Hugo_Symbol := paste0(Hugo_Symbol, onco)] 

	maf.skcm[Hugo_Symbol == 'NF1' & oncogenic %in% onco.sel, cc.sel, with=F]
	maf.skcm[Hugo_Symbol == 'TP53' & oncogenic %in% onco.sel, cc.sel, with=F]
	maf.skcm[Hugo_Symbol == 'BRAV' & oncogenic %in% onco.sel & HGVSp_Short == 'p.V600E', cc.sel, with=F]

	mtx = as.data.frame.matrix(table(maf.skcm[oncogenic %in% onco.sel, .(Tumor_Sample_Barcode, Hugo_Symbol)])) # oncogenic
	mtx = as.data.frame.matrix(table(maf.skcm[, .(Tumor_Sample_Barcode, Hugo_Symbol)])) # oncogenic
	#mtx = as.data.frame.matrix(table(maf[, .(Tumor_Sample_Barcode, Hugo_Symbol)])) # all mutations
	mtx[mtx >1] = 1
	mtx = t(mtx)
	events = discover.matrix(mtx)

	genes = c('TP53', 'NF1', 'BRAF', 'NRAS', 'PTEN', 'CDKN2A', 'BRAF:p.V600E')
	on = c('onco', 'non')
	genes = paste0(c(genes, genes), rep(on, each=length(genes)))
	genes = intersect(genes, maf.skcm$Hugo_Symbol)
	genes = sort(genes)
}


## prostate cancer
{
	save.image()
	unique(maf$CANCER_TYPE)
	maf.prad = maf[CANCER_TYPE == 'Prostate Cancer', ] 
	dim(maf.prad)
	source('~/program/fun/param.r')
	maf.prad[Hugo_Symbol == 'AR', cc.sel, with=F]
	unique(maf.prad$Tumor_Sample_Barcode)
	83/1215
}


## impact dscrt
{
	source('~/program/fun/param.r')
	maffile = 'msk_impact_2017_oncokb/data_mutations_merged_maf2maf_facets_clin_oncokb_cleaned.txt'
	fread(maffile) -> maf
	maf[, vaf := t_alt_count / t_depth]
	maf = maf[vaf > 0.05 & t_alt_count > 3, ]
	maf.dscrt = maf[CANCER_TYPE == 'Soft Tissue Sarcoma' & CANCER_TYPE_DETAILED =='Desmoplastic Small-Round-Cell Tumor' , ] 
	maf.dscrt = maf.dscrt[, cn, with=F]
	maf.dscrt
	unique(maf.dscrt$Tumor_Sample_Barcode)
	write.table(unique(maf.dscrt$Tumor_Sample_Barcode), file='aa', sep="\t", quote=F)
	maf.dscrt[Hugo_Symbol == 'PGBD5', ]

	source('~/program/fun/cbioMaf.R')
	tmp = cbioMaf(maf.dscrt)
	tmp.o = read.maf(tmp$maf, vc_nonSyn = tmp$vc.key, clinicalData=anno)
	pdf('res/impact_dscrt.pdf', width=6, height=7)
	oncoplot(tmp.o, keepGeneOrder=T, color = tmp$color.table, GeneOrderSort = F, drawRowBar=F, drawColBar=F, sortByMutation=T, writeMatrix=T, showTumorSampleBarcodes=T, SampleNamefontSize=6)
	dev.off()
	sync('res')

	tmp.o = read.maf(maf.dscrt)
	pdf('res/impact_dscrt_all.pdf', width=6, height=7)
	oncoplot(tmp.o, keepGeneOrder=T, GeneOrderSort = F, drawRowBar=F, drawColBar=F, sortByMutation=T, writeMatrix=T, showTumorSampleBarcodes=T, SampleNamefontSize=6)
	dev.off()
	sync('res')


}

## msisensor
{
	system('msisensor scan -d /ifs/depot/assemblies/H.sapiens/b37/b37.fasta -o ./msisensor/b37')
	system('msisensor scan -d /ifs/depot/assemblies/H.sapiens/b37_dmp/b37.fasta -o ./msisensor/dmp_b37', wait=F)

}
