igraph.plot = function(maf, fn, ww = 9, hh=9){
	gs.oncogenic = unique(maf[oncogenic != '', Hugo_Symbol])
	cat(NULL, file='two_genes_nTp53.tsv')
	for(i in 1:length(gs.oncogenic)){
		for(j in (i+1):length(gs.oncogenic)){
			#i=8;j=52
			bcr.i = unique(maf2[Hugo_Symbol == gs.oncogenic[i] & oncogenic != '', Tumor_Sample_Barcode])
			bcr.j = unique(maf2[Hugo_Symbol == gs.oncogenic[j] & oncogenic != '', Tumor_Sample_Barcode])
			bcr.ij = intersect(bcr.i, bcr.j)
			if(length(bcr.ij) > 0){
				cat(gs.oncogenic[i], "\t", gs.oncogenic[j], "\t", length(bcr.ij),  "\t", paste(bcr.ij, collapse=";"), "\n", file='two_genes_nTp53.tsv', append=T)
			}else{
				cat(gs.oncogenic[i], "\t", gs.oncogenic[j], "\tno patients have both mutations\n")
			}
		}
	}

	two.genes = fread('two_genes_nTp53.tsv', header=F)
	two.genes = setNames(two.genes, c('gene1', 'gene2', 'nn',  'bcrs'))
	head(two.genes)
	two.genes = two.genes[order(nn, decreasing=T),]
	two.genes
	two.genes$nn
	two.genes[nn == 5, .(gene1, gene2)]
	two.genes.2 = two.genes[nn > 4,]
	two.genes.2.nodes = data.frame(id=unique(c(two.genes.2$gene1, two.genes.2$gene2)))
	two.genes.2.nodes
	two.genes.2.links = two.genes.2[, .(gene1, gene2, nn)]
	two.genes.2.links = setNames(two.genes.2.links, c('from', 'to', 'weight'))
	two.genes.2.links = two.genes.2.links[sample(1:nrow(two.genes.2.links)),]
	two.genes.2.links

	library(igraph)
	net = graph_from_data_frame(d=two.genes.2.links, vertices=two.genes.2.nodes, directed=F) 
	l <- layout_in_circle(net)
	E(net)$width <- E(net)$weight
	pdf('res/net.sample_nTp53.pdf', width=10, height=10)
	plot(net, layout=l, edge.label=two.genes.2.links$weight)
	dev.off()
}
