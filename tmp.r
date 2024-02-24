	clin
	clin[, status := vital_status]
	clin[, surv.day := days_to_death]
	clin[, surv.mon := surv.day/30]
	clin[, surv.year := surv.mon/12]
	clin[, pfs.day := days_to_last_follow_up]
	clin[, pfs.mon := pfs/30]
	clin[, pfs.year := pfs.mon/12]

	load('../utuc_rnaseq/oct3/blca__cca.RData') # blca.cca
	

	tmp = clin
	tmp[, grp := get('kinase')]
	fit.os  <- survfit(Surv(surv.day,  status == 'dead') ~ grp, data = tmp)
	fit.pfs <- survfit(Surv(pfs.day,   status == 'dead') ~ grp, data = tmp)
	ov.plot = ggsurvplot(fit.os, pval=T, surv.median.line='v', risk.table = T, ylab='Overall Survival', xlab='Time(years)', 
			     xscale='d_y', break.x.by=365.25, xlim=c(0, 365.25*10))
	pfs.plot = ggsurvplot(fit.pfs,pval=T, risk.table = T,  surv.median.line='v', ylab = 'Progression Free Survival', xlab='Time(years)',  
			      xscale='d_y', break.x.by=365.25, xlim=c(0, 365.25*10))
	splots = list(ov = ov.plot, pfs = pfs.plot)
	fname = paste0('../utuc_rnaseq/oct3/survival_blca_v1.png')
	png(fname, width=3200, height=1800, res=300)
	arrange_ggsurvplots(splots, print = T, ncol = 2, nrow = 1, risk.table.height = 0.3, font=6)
	dev.off()
	scp(fname)



	#legend.title = legend, legend.labs = labels, palette = color, ...)

	bcr = unique(tcga.sel[oncogenic %in% onco.sel, ][Hugo_Symbol %in% kinase.gs, Tumor_Sample_Barcode])
	kinase.gs = c('FGFR3', 'PIK3CA', 'HRAS', 'KRAS', 'NRAS', 'BRAF', 'ERBB2')
	clin[bcr_patient_barcode  %in% bcr, kinase := 'Mut']
	clin[is.na(kinase),  kinase := 'Wt']
	clin[is.na(kinase),  clusterCol := kinase]
	fname = paste0('../utuc_rnaseq/oct3/tcga_blca_kinase_survival.png')
	source('TCGAanalyze_survival.R')
	TCGAanalyze_survival(clin, clusterCol="clusterCol", height = 6, width=6, filename = fname)
	TCGAanalyze_survival(clin, risk.table = F, clusterCol="gender", height = 6, width=6, filename = fname)
	scp(fname)

}
	clin.oo[, status := `Vital Status at Last Follow-up`]
	clin.oo[, surv.day := `date.of.last.followup` - date.of.surgery]
	clin.oo[, surv.mon := surv.day/30]
	clin.oo[, surv.year := surv.mon/12]
	clin.oo[, pfs.day := `date.of.recurrence.or.last.followup` - date.of.surgery]
	clin.oo[, pfs.mon := pfs/30]
	clin.oo[, pfs.year := pfs.mon/12]
	clin.oo[surv.day < 0, .(User.ID2, `Date of Surgery`,  date.of.surgery, date.of.last.followup)]
	clin.oo$surv.day
	clin.oo[, .(User.ID2, surv.mon, pfs.mon)]

	clin.oo[, local.recurrence := `Local Recurrence`]
'kinase')
	for(ii in genes){
		grp = ii
		tmp = clin[sample.type == 'UTUC', ];
		tmp[, grp := get(ii)]
		tmp[is.na(grp), grp := 'WT']
		tmp[grp != 'WT', grp := 'Mt']
		tmp[, grp := paste0(cluster, ':', grp)]
		fit.os  <- survfit(Surv(surv.day,  status == 'Dead') ~ grp, data = tmp)
		fit.pfs <- survfit(Surv(pfs.day,   status == 'Dead') ~ grp, data = tmp)
		ov.plot = ggsurvplot(fit.os, pval=T, surv.median.line='v', risk.table = T, ylab='Overall Survival', xlab='Time(years)', 
				      xscale='d_y', break.x.by=365.25, xlim=c(0, 365.25*10))
		pfs.plot = ggsurvplot(fit.pfs,pval=T, risk.table = T,  surv.median.line='v', ylab = 'Progression Free Survival', xlab='Time(years)',  
				      xscale='d_y', break.x.by=365.25, xlim=c(0, 365.25*10))
		splots = list(ov = ov.plot, pfs = pfs.plot)
		fname = paste0('survival_3cluster_cluster_', ii, '_v2.png')
		png(fname, width=3200, height=1800, res=300)
		arrange_ggsurvplots(splots, print = T, ncol = 2, nrow = 1, risk.table.height = 0.3, font=6)
		dev.off()
		scp(fname)
	}
