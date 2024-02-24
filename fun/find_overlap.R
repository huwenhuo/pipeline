find_overlape = function(gr1, gr2, gap=0){
	ov = findOverlaps(gr1, gr2, ignore.strand = T)
	tmp = as.data.frame(gr_peak[ov@from,])
	tmp2 = as.data.frame(gr_hg19_cpg[ov@to,])
	colnames(tmp2) = paste0("cpg_", colnames(tmp2))
	cbind(tmp, tmp2) -> tmp
	tmp = tmp[!duplicated(tmp$id),]
	row.names(tmp) = tmp$id
	anno3 = cbind(anno2, tmp[row.names(anno2),8:13])
	anno3$iscpg = !is.na(anno3$cpg_id)
	anno3
}


}
