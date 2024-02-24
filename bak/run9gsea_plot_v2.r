xls_files = list.files("gsea_report.*neg.*.xls", path="gsea", recursive=T, full.names=T)
for(i in 1:length(xls_files)){
	type_file = xls_files[i];
	cat(type_file, " now\n")
	gsea_report = read.table(type_file, sep="\t", header=T, stringsAsFactors = F);
	if(nrow(gsea_report) < 1 ){ next;}
	ll = nrow(gsea_report);
	if(ll > 20) { ll =20 }
	names = gsea_report$NAME[1:ll];
	names = gsub("_", " ", names, perl=T);
	pvalues = gsea_report$NOM.p.val[1:ll];
	pvalues[pvalues == 0] = 0.0001;
	pvalues = -log10(pvalues);

	names = names[order(pvalues)];
	pvalues = pvalues[order(pvalues)];

	title = sub("gsea/gsea.", "", type_file)
	title = sub(".gmt.GseaPreranked.*", "", title)
	file_pre = title;
	title = sub("\\.rnk\\.", " ", title)

	pdf(file=paste0("gsea/", file_pre, "_neg.pdf"), width=8, height=5);
	par(mai=c(1, 3.5, .2, .1));
	xx = barplot(pvalues, horiz = T, main=paste0(toupper(title), " Down-regulatd"), col=3, xlab="-Log10(pvalues)");
	text(0, xx, names, adj=c(1,0.5), cex=.7, xpd=T);
	dev.off();
}

xls_files = list.files("gsea_report.*pos.*.xls", path="gsea", recursive=T, full.names=T)
for(i in 1:length(xls_files)){
	type_file = xls_files[i];
	cat(type_file, " now\n")
	gsea_report = read.table(type_file, sep="\t", header=T, stringsAsFactors = F);
	if(nrow(gsea_report) < 1 ){ next;}
	ll = nrow(gsea_report);
	if(ll > 20) { ll = 20 }
	names = gsea_report$NAME[1:ll];
	names = gsub("_", " ", names, perl=T);
	pvalues = gsea_report$NOM.p.val[1:ll];
	pvalues[pvalues == 0] = 0.0001;
	pvalues = -log10(pvalues);

	names = names[order(pvalues)];
	pvalues = pvalues[order(pvalues)];

	title = sub("gsea/gsea.", "", type_file)
	title = sub(".gmt.GseaPreranked.*", "", title)
	file_pre = title;
	title = sub("\\.rnk\\.", " ", title)

	pdf(file=paste0("gsea/", file_pre, "_pos.pdf"), width=8, height=5);
	par(mai=c(1, 3.5, .2, .1));
	xx = barplot(pvalues, horiz = T, main=paste0(toupper(title), "Up-regulated"), col=2, xlab="-Log10(pvalues)");
	text(0, xx, names, adj=c(1,0.5), cex=.7, xpd=T);
	dev.off();
}
