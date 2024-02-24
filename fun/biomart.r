getid = function(my.values = '', my.key = 'entrezgene', mart = 'human', my.attributes=c('entrezgene', 'hgnc_symbol')){
	require(biomaRt)
	ensembl=useMart("ensembl")  # using ensembl database data
	if(mart == 'human'){
		ensembl=useDataset("hsapiens_gene_ensembl",mart=ensembl)   # from ensembl using homosapien gene data
	}
	#my.values = row.sel; my.key = 'ensembl_gene_id'; mart = 'human'; my.attributes=c('entrezgene', 'hgnc_symbol','ensembl_gene_id')
	# refseq_mrna
	tf = getBM(attributes = my.attributes, filters = my.key, values = my.values, mart= ensembl) 
	tf = getBM(attributes = c('entrezgene', 'hgnc_symbol', 'refseq_mrna'), filters='refseq_mrna', values = my.key, mart=ensembl)
	aa = getBM(attributes = c('hgnc_symbol', 'refseq_mrna'), filters='refseq_mrna', values =g$V4, mart=ensembl)
	row.names(tf) = tf$entrezgene
	tf
}
