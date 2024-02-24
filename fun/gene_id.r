library(data.table)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79)

genes = c('ARID1A', 'ARID1B', 'ARID2', 'CREBBP', 'EP300', 'KDM6A', 'KMT2A', 'KMT2B', 'KMT2C', 'KMT2D', 'SETD2', 'SETDB1', 'SMARCA4', 'STAG2', 'TET1', 'TET2', 'TP53', 'TSC1') 
genes
rs.hsa = genes(EnsDb.Hsapiens.v86, filter=list(GenenameFilter(genes),GeneIdFilter("ENSG", "startsWith")), return.type="data.frame", columns=c("gene_id"))
rs.hsa

fwrite(rs.hsa, file = '~/program/epiCrop/rs.hsa')
write(rs.hsa$gene_id, file = '~/program/epiCrop/rs.hsa.engs')
getwd()

genes = c('Arid1a', 'Arid1b', 'Arid2', 'Crebbp', 'Ep300', 'Kdm6a', 'Kmt2a', 'Kmt2b', 'Kmt2c', 'Kmt2d', 'Setd2', 'Setdb1', 'Smarca4', 'Stag2', 'Tet1', 'Tet2', 'Tp53', 'Tsc1') 
rs.mmu = genes(EnsDb.Mmusculus.v79, filter=list(GeneNameFilter(genes),GeneIdFilter("ENSG", "startsWith")), return.type="data.frame", columns=c("gene_id"))
rs.mmu

rs.hsa
rs.mmu



##
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Homo.sapiens)
library(universalmotif)

## get promoter region
gene = c('FOXA1')
gs = select(Homo.sapiens, columns=c('SYMBOL', 'ENTREZID'), key=gene, keytype='SYMBOL'); gs
gs = gs[!is.na(gs$ENTREZID),]
gr = transcriptsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by = "gene")[gs$ENTREZ]
gr.seq = getPromoterSeq (gr, Hsapiens, upstream=2000, downstream=1000)
names(gr.seq) = gene
gr.seq

