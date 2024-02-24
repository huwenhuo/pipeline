## load fun/simple_vc.tsv in param.r
## add simple.vc col here
## then add simple.vc.onco col here
## Notice: maftools are installed from program/maftools with R/oncomatrix.R modified by removing Multi_Hit
## Notice: a file lst.rdata will saved under the same folder
## https://docs.gdc.cancer.gov/Encyclopedia/pages/Mutation_Annotation_Format_TCGAv2/


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

#source('~/program/fun/sync.r')
#hmaf1 = fread("hmaf_fil_oncokb.hmaf")
#hmaf1
#hmaf2 = cbioMaf(hmaf1)
#hmaf2
#pdf('test.pdf')
#oncoplot(hmaf2)
#dev.off()
#sync()
#
#read.hmaf(hmaf2, vc_nonSyn = vc.key) -> hmaf.o
