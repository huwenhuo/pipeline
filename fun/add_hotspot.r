add.hotspot.mmu = function(input.maf){
	require(data.table)
	#input.maf = maf
	input.maf[, loc := sub(".*?([0-9]+).*", "\\1", HGVSp_Short)]
	input.maf[, tag := paste0(Hugo_Symbol, ':', loc)]
	load('~/program/cdd/hotspots_v2_mmu.Rda')
	source('~/program/fun/param.r')
	input.maf[, hotspot := F]
	input.maf[tag %in% hotspots$mmu.symbol.loc, hotspot := T]
	input.maf[toupper(Hugo_Symbol) %in% c(mut.2017, 'TRP53'), hotspot := T]
	input.maf[toupper(Hugo_Symbol) %in% c(mut.2017, 'TRP53') & Variant_Classification %in% trunc.mut, mut2017 := T]
	input.maf[toupper(Hugo_Symbol) %in% hotspots[tsg == T, toupper(mmu.symbol)], tsg := T]
	input.maf[hotspot==T & FILTER == 'PASS' & Variant_Classification != 'Intron', c(cc.1, 'vaf', 'FILTER'), with=F]
	input.maf
}
