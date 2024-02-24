require(data.table)
require(maftools)
library(ggplot2)

## a custom table of Variant Classification
vc.key = c("Missense_Mutation Oncogenic", "Inframe Oncogenic", 
	   "Inframe Unknown", "Missense_Mutation Unknown", 
	   "Truncating Unknown", "Truncating Oncogenic",  
	   "AML", "DEL", 'Multi_Hit')
vc.key

vc.non.syn = c('Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', 
	       'Missense_Mutation','Nonsense_Mutation', 'Nonstop_Mutation', 'Translation_Start_Site')

vc.table = data.table(
			  vc.value = factor(c('Truncating', 'Truncating', 'Inframe', 'Inframe', 
						  'Missense_Mutation', 'Truncating', 'Inframe', "Inframe",
						  'Inframe', 'Inframe', 'Inframe', 'Inframe'), 
						levels=c('Truncating', 'Missense_Mutation', 'Inframe')),
			  Variant.Classification = c('Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', 
						     'Missense_Mutation','Nonsense_Mutation', 'Nonstop_Mutation', "5'Flank", 
						     'Silent', 'Splice_Region', 'Splice_Site', 'Translation_Start_Site')
)
setkey(vc.table, 'Variant.Classification')

## a custom table of oncogenic 
oncogenic.table = data.table(oncogenic = c('Unknown', 'Likely Neutral', 'Inconclusive', 'Predicted Oncogenic', 'Likely Oncogenic', 'Oncogenic'),
			     oncogenic.value = factor(c('Unknown', 'Unknown', 'Unknown', 'Oncogenic', 'Oncogenic', 'Oncogenic'), levels=c('Oncogenic', 'Unknown')))
setkey(oncogenic.table, 'oncogenic')
colnames(oncogenic.table) = gsub('_', '.',  colnames(oncogenic.table))
oncogenic.table

onco.sel = c('Predicted Oncogenic', 'Likely Oncogenic', 'Oncogenic')

## a custom table of colors based on the VC and the oncogenic or not
color.table  = data.table( color.index = c("Missense_Mutation Oncogenic",
					   "Missense_Mutation Unknown",
					   "Truncating Oncogenic",
					   "Truncating Unknown",
					   "Inframe Oncogenic",
					   "Inframe Unknown",
					   "AML",
					   "DEL",
					   "Multi_Hit"),
			  color.value = c("forestgreen", 
					  adjustcolor("forestgreen", alpha.f=0.6), 
					  "black", 
					  adjustcolor("black", alpha.f=0.6),
					  "brown4", 
					  adjustcolor("brown4", alpha.f=0.6),
					  "red",
					  "blue",
					  "magenta3") )
color.table[, color.index := factor(color.index, levels=c("Missense_Mutation Oncogenic", "Truncating Oncogenic", 
							  "Inframe Oncogenic", "Missense_Mutation Unknown", 
							  "Truncating Unknown", "Inframe Unknown", 
							  "AML", "DEL", 'Multi_Hit') )]
setkey(color.table, 'color.index')
colnames(color.table) = gsub('_', '.', colnames(color.table))
color.table

color.vector = color.table$color.value
names(color.vector) = color.table$color.index
color.vector

## I have included an example of maf var
## this one need oncogenic column and others

## make a copy of original maf 
hmaf = copy(maf)

## convert Variant Classification to a custom form
setkey(hmaf, 'Variant_Classification')
hmaf = merge(hmaf, vc.table, all.x=T, by.x = 'Variant_Classification', by.y = 'Variant.Classification')

## convert oncogenic with a custom form
setkey(hmaf, 'oncogenic')
hmaf[oncogenic == '', oncogenic := 'Unknown']
hmaf = merge(hmaf, oncogenic.table, all.x=T, by = 'oncogenic')

## attach a color based on the custom forms of Variant Classification and oncogenic values
## and the above color.table
hmaf[, id := paste(vc.value, oncogenic.value, sep=' ')]
setkey(hmaf, 'id')
hmaf = merge(hmaf, color.table, by.x = 'id', by.y = 'color.index', all.x=T)
hmaf

hmaf[, Variant_Classification_old := Variant_Classification]
hmaf[, Variant_Classification := id]


genes = c('KMT2D', 'TP53', 'PIK3CA', 'ATM', 'ARID1A', 'EP300', 'FGFR3', 'STAG2', 'KDM6A', 'KMT2C', 'RB1', 'CUL1', 'ELF3', 'ERBB3', 'KMT2A', 'NFE2L2', 'RHOA', 'ACTB', 'CREBBP', 'SPTAN1')

hmaf.o = read.maf(hmaf, vc_nonSyn = vc.key)
fname = 'test_oncogenic.pdf'
pdf(fname, width=6, height=6)
oncoplot(hmaf.o, drawRowBar=F, genes = genes, drawColBar=F, showTumorSampleBarcodes=T, colors = color.vector)
dev.off()

save.image(file='oncoplot_test.RData')
