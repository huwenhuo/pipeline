granges.2.bed = function(gr, bed.file){
	df <- data.frame(chr=seqnames(gr),
			 start=start(gr)-1,
			 end=end(gr),
			 name=c(rep(".", length(gr))),
			 score=c(rep(".", length(gr))),
			 strand=strand(gr))
	write.table(df, file=bed.file, quote=F, sep="\t", row.names=F, col.names=T)
}
