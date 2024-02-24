
save.image()

options(bitmapType="cairo") 

source("https://bioconductor.org/biocLite.R")
biocLite("TCGAbiolinks")

library(DESeq2)
library(maftools)
library(SummarizedExperiment)
#library(TCGAbiolinks)
library(BiocParallel)
register(MulticoreParam(10))
library("AnnotationDbi");
library(foreach);
library("genefilter");
library("ggplot2");
library(gplots);
library("RColorBrewer");
library("org.Hs.eg.db");
library(ChIPQC)
library(tximport)
library("pheatmap");
library(ConsensusClusterPlus)

#RNA-Seq, old one
filename = c( "rt4_kdm6a_2_1", "rt4_kdm6a_2_2", "rt4_kdm6a_3_1", "rt4_kdm6a_3_2", "rt4_par_1", "rt4_par_2");
target = data.frame(filename=filename)
target$countfile = paste0("rnaseq_kdm6a/", filename, "/", filename, ".htseq.count")
target$sampleid = filename;
target$group = c( "rt4_kdm6a_2", "rt4_kdm6a_2", "rt4_kdm6a_3", "rt4_kdm6a_3", "rt4_par", "rt4_par");
target$group2 = c( "rt4_kdm6a", "rt4_kdm6a", "rt4_kdm6a", "rt4_kdm6a", "rt4_par", "rt4_par");
target$path = ".";
target$libtype = "PE";
if (exists("counttable")) {rm(counttable)}
for (i in 1:nrow(target)) {
	    sample <- read.table(as.character(target$countfile[i]));
    colnames(sample) = c("gene", as.character(target$sampleid[i]));
        if (exists("counttable")) {counttable <- merge(counttable, sample, by= "gene")} else {counttable <- sample}
        rm(sample)
}
row.names(counttable)= counttable$gene;
counttable$gene <- NULL;
head(counttable);
condition = target$group;

design = data.frame(
		    row.names       =colnames(counttable),
		    condition       =condition,
		    libType =target$libtype);

ddsmat = DESeqDataSetFromMatrix(countData = counttable,
				colData = design,
				design = ~ condition);
dds.ds <- estimateSizeFactors(ddsmat);

dds.ds$condition <- relevel(dds.ds$condition, "rt4_par");
dds <- DESeq(dds.ds, parallel=T);

kdm6a2<- results(dds, contrast=c("condition","rt4_kdm6a_2", "rt4_par"));
kdm6a3<- results(dds, contrast=c("condition","rt4_kdm6a_3", "rt4_par"));
head(kdm6a3);

kdm6a2 = as.data.frame(kdm6a2);
kdm6a3 = as.data.frame(kdm6a3);

vv = c("kdm6a2", 'kdm6a3');
vv.sel = paste0(vv, ".sel");

for(i in 1:length(vv)){
	    tmp = get(vv[i]);
    tmp = tmp[order(tmp$pvalue),];
        tmp$symbol = convertIDs(row.names(tmp), from = "ENSEMBL", to = "SYMBOL", db = org.Hs.eg.db);
        assign(vv[i], tmp);
	    assign(vv.sel[i], tmp[!is.na(tmp$padj) & tmp$padj < 0.05, ]);
}
kdm6a2 = kdm6a2[!is.na(kdm6a2$log2FoldChange), ]
kdm6a3 = kdm6a3[!is.na(kdm6a3$log2FoldChange), ]
write.table(kdm6a2, file='kdm6a2_vs_wt.txt', quote=F, sep="\t");
write.table(kdm6a3, file='kdm6a3_vs_wt.txt', quote=F, sep="\t");
kdm6a2.sel = kdm6a2[!is.na(kdm6a2$padj) & kdm6a2$padj < 0.05 & !is.na(kdm6a2$log2FoldChange) & abs(kdm6a2$log2FoldChange) > 1, ]
kdm6a3.sel = kdm6a3[!is.na(kdm6a3$padj) & kdm6a3$padj < 0.05 & !is.na(kdm6a3$log2FoldChange) & abs(kdm6a3$log2FoldChange) > 1, ]

source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/vennDia.R")
qq = venndiagram(x = kdm6a2.sel$symbol, y = kdm6a3.sel$symbol, type='2', labels=c("KDM6A 2 vs WT", "KDM6A 3 vs WT"), title = "Gene symbol");
qq = venndiagram(x = kdm6a2.sel$symbol[kdm6a2$log2FoldChange > 0], y = kdm6a3.sel$symbol[kdm6a3$log2FoldChange > 0], type='2', title="Gene symbol",
		                  labels=c("KDM6A 2 vs WT: UP", "KDM6A 3 vs WT: UP"));
qq = venndiagram(x = kdm6a2.sel$symbol[kdm6a2$log2FoldChange < 0], y = kdm6a3.sel$symbol[kdm6a3$log2FoldChange < 0], type='2', title="Gene symbol",
		                  labels=c("KDM6A 2 vs WT: DN", "KDM6A 3 vs WT: DN"));

cc = counts(dds, normalized=T);
sel.id = unique(c(row.names(kdm6a2.sel), row.names(kdm6a3.sel)));
cc[sel.id,] -> cc.sel;
dim(cc.sel)

col <- anno = data.frame(row.names = colnames(cc.sel), Cells = c("KDM6A 2", "KDM6A 2", "KDM6A 3", "KDM6A 3", "Parental", "Parental"), stringsAsFactors = F)
col <- anno
tmp = convertIDs(row.names(cc.sel), from = "ENSEMBL", to = "SYMBOL", db = org.Hs.eg.db);
cc <- sel <- gs = cc.sel
row.names(cc <- sel <- gs) = tmp

genelist <- vv = load("/Volumes/LaCie/huw/solit/study/memo/genelist.RData")
genelist <- vv
genelist <- uc = load("/Volumes/LaCie/huw/solit/study/hiseq/bcg/uc_basal_luminal_marks.RData")
genelist <- uc
pheatmap(cc.sel, show <- rownames = F, scale="row", clustering <- method = 'complete', cluster <- rows = TRUE, cluster <- cols = TRUE, color=greenred(75),
	         annotation <- col = col <- anno)
write.table(cc, file="rnaseq_kdm6a_normalized_reads_deseq2.txt")

read.table("rnaseq_kdm6a_normalized_reads_deseq2.txt")

read.table("rnaseq_kdm6a/kdm6a2_vs_wt.txt", header=T, stringsAsFactors = F) -> kdm6a2
read.table("rnaseq_kdm6a/kdm6a3_vs_wt.txt", header=T, stringsAsFactors = F) -> kdm6a3
read.table("rnaseq_kdm6a/rnaseq_kdm6a_normalized_reads_deseq2.txt") -> cc

## rna-seq of recent batch
## import RSEM results and combine it with tcga
target=data.frame(bases=c('K2R1', 'K2R2', 'K2R3', 'PR1', 'PR2', 'PR3'))
filenames = paste0("Project_7373_B/rsem/", target$bases, ".genes.results")
filenames

getwd()
cmd = paste0("/opt/common/CentOS_6/rsem/RSEM-1.2.25/rsem-generate-data-matrix ", paste(filenames, collapse = " "), " > rnaseq_new_rsem_matrix.txt")
cmd
system(cmd)

cmd = paste0("/opt/common/CentOS_6/rsem/RSEM-1.2.25/rsem-run-ebseq rnaseq_new_rsem_matrix.txt 3,3 rnaseq_new_DE.results")
system(cmd)

cmd = paste0("/opt/common/CentOS_6/rsem/RSEM-1.2.25/rsem-control-fdr rnaseq_new_DE.results 0.05 rnaseq_new_DE.res")
system(cmd)
rnaseq = fread("rnaseq_new_DE.res")
##setnames(rnaseq, c("symbol", "PPEE", "PPDE", "PostFC", "RealFC","C1Mean","C2Mean"))
head(rnaseq)

names(filenames) = target$bases
txi = tximport(filenames, type = "rsem")
colnames(txi$counts)
str(txi)

sel = apply(txi$length, 1, function(x){all(x>0)})
abundance = txi[["abundance"]][sel,]
counts = txi[["counts"]][sel,]
length= txi[["length"]][sel,]
countsFromAbundance = 'no'
txi = list(abundance = abundance, counts = counts, length = length, countsFromAbundance = 'no')
names(txi)
mode(txi[['counts']]) = "integer"
head(txi[['counts']])

condition = c('K2', 'K2', 'K2', 'PR', 'PR', 'PR')
condition
design = data.frame(
		    row.names       =c('K2R1', 'K2R2', 'K2R3', 'PR1', 'PR2', 'PR3'),
		    condition       =condition,
		    libType =rep("PE", 6));

ddsmat = DESeqDataSetFromTximport(txi, design, ~ condition)
ddsmat
dds.ds <- estimateSizeFactors(ddsmat);
dds <- DESeq(dds.ds, parallel=T);
res = results(dds, contrast=c("condition", "K2", "PR"))
res = res[order(res$pvalue),]
head(res)
save.image()

source("~/program/fun/write_rnk.r")
symbol = strsplit(row.names(res), "_")
symbol = unlist(lapply(symbol, "[[", 2))
res$symbol = symbol
res = as.data.frame(res)
res$symbol = sub(".*_", "", row.names(res))
head(res)

write_rnk(res, file="Project_7373_B/res/res.rnk")
write_tsv(res, path="Project_7373_B/res/res.tsv")

rld = rlog(dds)

pca = plotPCA(rld)
pca = pca + geom_text(aes(label=colnames(rld), vjust=0, hjust=0) )
pdf(file="res/rnaseq_new_pca.pdf")
print(pca)
dev.off()

rld = rlog(dds)

pca = plotPCA(rld)
pca = pca + geom_text(aes(label=colnames(rld), vjust=0, hjust=0) )
pdf(file="res/rnaseq_new_pca2.pdf")
print(pca)
dev.off()

pdf(file="pca2.pdf")
se <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1), colData=colData(dds))
print(se)
plotPCA( DESeqTransform( se ) )
dev.off()

rnaseq = read.table("Project_7373_B/res/res.tsv", header = T)
head(rnaseq)

##############################
### gsea

source("~/program/fun/run_gsea.R")
targets = c( "c2.cgp.v5.0.symbols.gmt", "c2.cp.biocarta.v5.0.symbols.gmt", "c2.cp.kegg.v5.0.symbols.gmt", "c2.cp.reactome.v5.0.symbols.gmt", "c2.cp.v5.0.symbols.gmt", "c5.bp.v5.0.symbols.gmt", "c5.cc.v5.0.symbols.gmt", "c5.mf.v5.0.symbols.gmt")
rnks = c("res.rnk")
run_gsea(rnkfile = rnks, gmtfile = targets)

############################################################
# ATAC-Seq
# dir ATAC-seq
# ATAC-Seq/run_loca.sh
getwd()
setwd("ATAC_seq")

atac_ko1_bam="KO1/Chr2306-KO1-AT-A11-AGGTTGGG_AGGTTGGG_BHG255BCXY_L002_001.R1_val.hg19.sorted.RmDup.bam"
atac_ko2_bam="KO2/Chr2307-KO2-AT-A11-CCGTTTGT_CCGTTTGT_BHG255BCXY_L002_001.R1_val.hg19.sorted.RmDup.bam"
atac_wt1_bam="WT1/Chr2304-WT1-AT-A11-TGCTGGGT_TGCTGGGT_BHG255BCXY_L002_001.R1_val.hg19.sorted.RmDup.bam"
atac_wt2_bam="WT2/Chr2305-WT2-AT-A11-GAGGGGTT_GAGGGGTT_BHG255BCXY_L002_001.R1_val.hg19.sorted.RmDup.bam"

cmd = past0("macs2 callpeak -t ", atac_ko1_bam, " ", atac_ko2_bam, " -c ", atac_wt1_bam, " ", atac_wt2_bam, " -n KOvsWT --outdir KO_WT -f BAM -g hs -B --nomodel --extsize 200  --SPMR ")
system(cmd)
cmd = past0("macs2 callpeak -c ", atac_ko1_bam, " ", atac_ko2_bam, " -t ", atac_wt1_bam, " ", atac_wt2_bam, " -n WTvsKO --outdir WT_KO -f BAM -g hs -B --nomodel --extsize 200  --SPMR ")
system(cmd)

peakfile = "KO_WT/KOvsWT_peaks.narrowPeak"
outfile = paste0(peakfile, ".annotate")
cmd = paste0("annotatePeaks.pl ", peakfile, " hg19 ", " > ", outfile); system(cmd)
source("process_homer_annotation.R")
ko.wt = process_homer_annotation(outfile, peakfile, "../hg19_CpG.bed")
head(ko.wt)
dim(ko.wt)
write.table(ko.wt, file=paste0(outfile, ".processed"))

## export genes that has cpg and peak in KO of KDM6A
table(anno_KOvsWT$anno)
dim(anno_KOvsWT[anno_KOvsWT$anno == 'promoter-TSS',])
gn = anno_KOvsWT[anno_KOvsWT$anno == 'promoter-TSS' & anno_KOvsWT$iscpg == T, 'Gene.Name']
gn = anno_KOvsWT[anno_KOvsWT$anno == 'promoter-TSS', 'Gene.Name']
cat("ATAC_K2_100", "", gn[1:100], file="~/program/atac_seq_top100_enriched_k2.gmt", sep="\t")

source("~/program/fun/run_gsea.R")
run_gsea(rnkfile = "../rnaseq_kdm6a/rnk_pvalue_kdm6a2.rnk", gmtfile = "atac_seq_top100_enriched_k2.gmt")

## if overlap with impact genes
impact = impact$V1
head(impact)
intersect(impact, gn[1:1000])

peakfile = "WT_KO/WTvsKO_peaks.narrowPeak"
outfile = paste0(peakfile, ".annotate")
cmd = paste0("annotatePeaks.pl ", peakfile, " hg19 ", " > ", outfile); print(cmd); system(cmd)
wt.ko = process_homer_annotation(outfile, peakfile, "../hg19_CpG.bed")
write.table(wt.ko, file=paste0(outfile, ".processed"))
nrow(wt.ko)
## there are no peaks enriched in WT vs KO. 

library(data.table)
wt.ko = as.data.table(wt.ko)
ko.wt = as.data.table(ko.wt)
sel.col = c('Gene.Name', 'FC', 'logpvalue', 'logqvalue', 'peakline')
ko.wt[anno == 'promoter-TSS', sel.col, with = F]
wt.ko[anno == 'promoter-TSS', sel.col, with = F]
ko.wt$vs = "koVSwt"
wt.ko$vs = "wtVSko"
rbind(wt.ko, ko.wt) -> atac
head(atac)
nrow(atac)

## add the RNA-Seq information of the older batch
head(kdm6a2)
tmp = kdm6a2
colnames(tmp) = paste0("kdm6a2.", colnames(kdm6a2))
as.data.table(tmp) -> tmp
setkey(tmp, 'kdm6a2.symbol')
setkey(atac, 'Gene.Name')
atac2 = atac[tmp, nomatch=0]
head(atac2)
atac2[,.(Gene.Name, kdm6a2.log2FoldChange)]
atac2[,.(Gene.Name, kdm6a2.log2FoldChange)]
atac2

head(rnaseq)
tmp = as.data.table(rnaseq)
head(tmp)
colnames(tmp) = paste0("new.", colnames(tmp))
setkey(tmp, 'new.symbol')
atac3 = atac2[tmp, nomatch=0]
atac3[,.(Gene.Name, FC, logpvalue, kdm6a2.log2FoldChange, new.log2FoldChange)]
atac3 = atac3[order(FC, decreasing = T)]
atac3[anno == 'promoter-TSS',.(Gene.Name, pos,Distance.to.TSS, FC, logpvalue, kdm6a2.log2FoldChange, new.log2FoldChange, anno, vs)]
nrow(ko.wt)
nrow(wt.ko)
head(atac3[anno == 'promoter-TSS' & kdm6a2.padj < 0.01 & new.padj < 0.01 & abs(kdm6a2.log2FoldChange) > 1, .(Gene.Name, pos,Distance.to.TSS, FC, logpvalue, kdm6a2.log2FoldChange, new.log2FoldChange, anno, vs)], 20)

## add 
colnames(atac3)
atac4 = atac3[anno == 'promoter-TSS', ]
atac4[,.(Chr, Start, End, cpg_start, cpg_end, iscpg)]

## redefine the cpg island
source("~/program/fun/bed_to_granges.r")
cpgfile = "../hg19_CpG.bed"
bed_to_granges(cpgfile) -> gr_hg19_cpg  # 28691
gr_hg19_cpg
atac5 = makeGRangesFromDataFrame(atac4, keep.extra.columns = T)
ov = findOverlaps(atac5, gr_hg19_cpg, ignore.strand=T, maxgap = 500)
atac4[ov@from, "iscpg"] = T

atac6 = as.data.table(atac4)

gg = ggplot(data = atac6, aes(x = kdm6a2.log2FoldChange, y = -log(kdm6a2.pvalue), color = iscpg)) + geom_point() 
ggsave(gg, file="../res/test2.pdf", width=15)

gg = ggplot(data = atac6, aes(x = kdm6a2.log2FoldChange, y = -log(kdm6a2.pvalue), color = iscpg)) + 
	geom_point()  + facet_wrap(~iscpg) + 
	geom_text(data = atac6[-log(kdm6a2.pvalue) > 100, ], aes(x = kdm6a2.log2FoldChange, y = -log(kdm6a2.pvalue), label=Gene.Name), nudge_x=0.18, nudge_y=5, size = 2)
ggsave(gg, file="../res/test3.pdf", width=15)

gg = ggplot(data = atac6[vs == "wtVSko",], aes(x = kdm6a2.log2FoldChange, y = -log(kdm6a2.pvalue), color = iscpg)) + 
	geom_point()  + facet_wrap(~iscpg) + 
	geom_text(data = atac6[-log(kdm6a2.pvalue) > 100 & vs == 'wtVSko', ], aes(x = kdm6a2.log2FoldChange, y = -log(kdm6a2.pvalue), label=Gene.Name), nudge_x=0.18, nudge_y=5, size = 2)
ggsave(gg, file="../res/test3.pdf", width=15)

colnames(atac6)
atac6$vs
table(atac6$iscpg)
nrow(atac3) ## total ATAC peaks
nrow(atac6) ## promoter TSS region peaks
nrow(atac7) ## promoter with CpG peaks

atac3 = atac3[order(logpvalue, decreasing = T),]
head(atac3[anno == 'promoter-TSS' & kdm6a2.padj < 0.01 & new.padj < 0.01 & abs(kdm6a2.log2FoldChange) > 1, .(Gene.Name, pos,Distance.to.TSS, FC, logpvalue, kdm6a2.log2FoldChange, new.log2FoldChange, anno, vs, iscpg)], 20)

## homer find motif for ATAC in the koVSwt peaks
"findMotifsGenome.pl <peak/BED file> <genome> <output directory> -size # [options]"
peafile = "KO_WT/KOvsWT_peaks.narrowPeak"

cmd = paste0("~/program/homer/bin/findMotifsGenome.pl ", peakfile, " hg19 KO_WT/homer -size 200")
cmd

Sys.setenv("PATH"=paste0("/home/huw/loca/weblog:", Sys.getenv("PATH")))
system(cmd)

motif = c('PAX6', 'ZNF415', 'GATA3', 'ZNF189', 'EKLF', 'KLF10')
motif = c('NFYA', 'NFYB', 'NFYC', 'ELF5', 'SP1', 'GATA4', "ETS1", "PITX1", "TFE3", "SOX4")
kdm6a2[kdm6a2$symbol %in% motif,] 
head(kdm6a2)


###### end here
## check the mRNA that has CpG island
source("~/program/fun/bed_to_granges.r")
bed_to_granges("../refGene_hg19_TSS.bed") -> gr_hg19_tss # 41088, 23635 unique genes
source("~/program/fun/extend_granges.R")
gr_hg19_promoter = extend_granges(gr_hg19_tss, 3000, 3000)
gr_hg19_promoter
bed_to_granges("../hg19_CpG.bed") -> gr_hg19_cpg  # 28691
ov = findOverlaps(gr_hg19_cpg, gr_hg19_promoter, ignore.strand=T)
subsetByOverlaps
ov
# all CpG islands are promoter TSS region 

## KO_WT/KOvsWT_peaks.narrowPeak peaks identified
bed_to_granges("KO_WT/KOvsWT_peaks.narrowPeak") -> gr_atac_peak # 15038
ov = findOverlaps(gr_hg19_promoter, gr_atac_peak, ignore.strand = T) #
length(unique(ov@to)) # 3878 of atac peak and promoter region
length(unique(ov@from)) # 13034 of atac peak and promoter region
gr_atac_peak_promoter = gr_hg19_promoter[unique(ov@to),] 
gr_atac_peak_promoter # 6810
ov = findOverlaps(gr_hg19_cpg, gr_atac_peak_promoter, ignore.strand = T) 
gr_atac_peak_promoter_cpg = gr_atac_peak_promoter[unique(ov@to),]
gr_atac_peak_promoter_cpg # 4707

## 
tt = table(anno_KOvsWT$anno)
tt = as.data.frame(tt)
colnames(tt) = c('Type', 'Peaks')
tt

pdf(file="../res/atac_anno_pie.pdf")
sp <- ggplot(tt, aes(x="", y=Peaks, fill=Type)) +
	geom_bar(width=1, stat="identity", color="black",) +
	scale_fill_brewer(palette="Pastel2") +
	coord_polar(theta="y")  + 
	geom_text(aes(x= 1.2, label = percent(percentage)), position = position_stack(vjust=0.5), size = 6)
print(sp)
dev.off()

tmp = table(anno_KOvsWT$anno, anno_KOvsWT$iscpg)
tmp = tmp[,2:1]
tmp = tmp[order(tmp[,1], decreasing = T),]
tmp

pdf(file="../res/atac_class_barplot.pdf", width=12, height=5)
barplot(t(tmp), beside=T, main='ATAC-Seq peaks in KDM6a KO RT4 cells', col=2:3)
legend("topright", c('Overlap w/ CpG island', 'No overlap'), fill = 2:3)
dev.off()

############################################################
#ChIP-Seq
# dir ATAC-seq
setwd("../ChIP-seq")
# run_chipseq.r
getwd()
target_chipseq = read.table("targets", sep="\t", header = T)
chipqc_chipseq = ChIPQC(target_chipseq, annotaiton="hg19")

## 

## merged
## new
## 10/2/2017
library(yaml)
library(data.table)

source('~/program/configure.R')
samtools
targetfile = 'targets.yaml'
source("~/program/fun/yaml2df.R")
targets2 = yaml2df(targetfile)
targets2 = as.data.table(targets2)
targets2[, path := dirname(oldbam)]
targets2[, batch := 'old']
colnames(targets2)
targets2$fastq1

pp = '/ifs/work/solitlab/huw/solit/study/hiseq/kdm6a_Mar2017/ATAC_seq'
atac.target = data.table(path2 = rep(pp,4))
atac.target[, batch := 'atac']
atac.target[, fastq1 := targets2$fastq1[9:12]]
atac.target[, fastq2 := targets2$fastq2[9:12]]
atac.target

#load the project Project_6852_E
tmp = fread("Project_6852_E/filelist")
dcast(tmp, samplename ~ R12, value.var='fastq', function(x){length(x)})
tmp = dcast(tmp, samplename ~ R12, value.var='fastq', function(x){paste(x, collapse=' ')})
tmp$batch = '6582'
as.data.table(tmp) -> tmp
tmp[, ID := samplename]
tmp[, cell := 'RT4']
tmp[, rep := 1]
tmp[, chip := samplename]
setnames(tmp, c('R1', 'R2'), c('fastq1', 'fastq2'))
tmp[, path := dirname(fastq1)]
tmp
colnames(tmp)

tt = targets[batch == '6582',]
tt

colnames(targets)

ov = intersect(colnames(tmp), colnames(targets2))
targets = rbind(targets2[, ov, with=F], tmp[, ov, with=F])
colnames(targets)
targets$path2[9:12]

# 
targets$fastq1
hg38bt2db = "/ifs/work/solitlab/huw/study/db/hsa/hg38_bowtie2/hg38"
hg38fai   = "/ifs/work/solitlab/huw/study/db/hsa/hg38_bowtie2/hg38.fai"
targets[, genome := hg38bt2db ]
targets$bamfile
file.exists(targets$bamfile)
colnames(targets)
targets$fastq2[1]
targets$fastq2n[9:12] = targets$fastq2[9:12]
targets[, bamfile := paste0(path2, "/", ID, '_hg38', ".bam") ]
targets[, samfile := paste0(path2, "/", ID, '_hg38.sam')]
targets[9:12, bamfile]
targets[, bam.posi.sorted.file := paste0(path2, "/", ID, '_hg38', "_posi_sorted.bam") ]
targets[, bam.name.sorted.file := paste0(path2, "/", ID, '_hg38', "_name_sorted.bam") ]
targets[, bam.posi.sorted      := paste0(path2, "/", ID, '_hg38', "_posi_sorted") ]
targets[, bam.posi.sorted.rmdup.file := paste0(path2, "/", ID, '_hg38', "_posi_sorted_rmdup.bam") ]
targets[, bam.posi.sorted.rmdup.filtered.file := paste0(path2, "/", ID, '_hg38', "_posi_sorted_rmdup_filtered.bam") ]
targets[, bam.sorted.rmdup.file := paste0(path2, "/", ID, '_hg38', "_sorted_rmdup.bam") ]
targets[, bases := ID]
targets$ID

targets[, path2 := paste0("/ifs/work/solitlab/huw/solit/study/hiseq/kdm6a_Mar2017/Project_6852_E/", ID),]
targets[, samfile := paste0(path2, "/", ID, '_hg38.sam')]
file.exists(targets$samfile)

targets[, fastq1n := paste0(path2, '/', basename(fastq1))]
targets[, fastq2n := paste0(path2, '/', basename(fastq2))]

# bowtie2
targets$ID[9:12] = c('RT41ATAC', 'RT42ATAC', 'K21ATAC', 'K22ATAC')
targets[9:12, fastq1 := paste0(path2, '/', basename(atac.target$fastq1))]
targets[9:12, fastq2 := paste0(path2, '/', basename(atac.target$fastq2))]
file.exists(targets[9:12, fastq1])
file.exists(targets[9:12, fastq2])
targets[9:12, fastq2]
targets$ID[9:12]
targets$path2[11:12] = sub('ATACK2', 'ATAC', targets$path2[11:12])
targets$path2[11:12]

targets[, bt.jobname := paste0("bt2.", ID)]
targets[, bt.cmd := paste0(BSUB, " -J ", bt.jobname, " -e ", bt.jobname, ".err -o ", bt.jobname, ".std ")]
targets[, bt.cmd := paste0(bt.cmd, "-cwd ", cwd, ' -We 2:59 -R "rusage[mem=30]" -R "rusage[iounits=0]" -n 10 ')]
#targets[, bt.cmd := paste(bt.cmd, '-w "post_done(', trim.jobname, ')"')]
targets[, bt.cmd := paste(bt.cmd, '"', bowtie2, '--local -p 9 -x', genome, '-1', fastq1, '-2', fastq2, '-k 1 -S', samfile, '"')]

sel$ID
sel = targets[ batch == 'atac', ]
sel$bt.cmd[3]
sel$bt.cmd
for(i in 1:nrow(sel)){system(sel$bt.cmd[i])}

## convert to bam
targets[, sam2bam.jobname := paste0("sam2bam.", ID)]
targets[, sam2bam.cmd := paste0(BSUB, " -J ", sam2bam.jobname, " -e ", sam2bam.jobname, ".err -o ", sam2bam.jobname, ".std ")]
targets[, sam2bam.cmd := paste0(sam2bam.cmd, "-cwd ", cwd, ' -We 2:29 -R "rusage[mem=8]" -R "rusage[iounits=0]" -n 2')]
targets[, sam2bam.cmd := paste(sam2bam.cmd, '-w "post_done(', bt.jobname, ')"')]
targets[, sam2bam.cmd := paste(sam2bam.cmd, '"', samtools, 'view -bt', hg38fai, samfile, ' >', bamfile, '"')]
sel = targets[9:12, ]
sel$sam2bam.cmd[1]

for(i in 1:nrow(sel)){system(sel$sam2bam.cmd[i])}

# coordinate sort bam for samtools view -f 2
targets[, posisortbam.jobname := paste0("posisortbam.", ID)]
targets[, posisortbam.cmd := paste0(BSUB, " -J ", posisortbam.jobname, " -e ", posisortbam.jobname, ".err -o ", posisortbam.jobname, ".std ")]
targets[, posisortbam.cmd := paste0(posisortbam.cmd, "-cwd ", cwd, ' -We 2:29 -R "rusage[mem=28]" -R "rusage[iounits=0]" -n 11 ')]
targets[, posisortbam.cmd := paste(posisortbam.cmd, '-w "post_done(', sam2bam.jobname, ')"')]
targets[, posisortbam.cmd := paste0(posisortbam.cmd, ' "', samtools, ' sort -@ 10 -m 2G ', bamfile, ' -o ', bam.posi.sorted.file, '"')]

sel = targets[ batch == 'atac', ]
sel$posisortbam.cmd[3]

for(i in 1:nrow(sel)){system(sel$posisortbam.cmd[i])}
system(sel$posisortbam.cmd[3])
system(sel$posisortbam.cmd[4])
file.exists(sel$bam.posi.sorted.file)

# index
targets[, index.jobname := paste0("index.", ID)]
targets[, index.cmd := paste0(BSUB, " -J ", index.jobname, " -e ", index.jobname, ".err -o ", index.jobname, ".std ")]
targets[, index.cmd := paste0(index.cmd, "-cwd ", cwd, ' -We 2:29 -R "rusage[mem=4]" -R "rusage[iounits=0]" -n 2 ')]
targets[, index.cmd := paste(index.cmd, '-w "post_done(', posisortbam.jobname, ')"')]
targets[, index.cmd := paste(index.cmd, ' "', samtools, 'index', bam.posi.sorted.file, '"')]
sel = targets[ batch == 'atac', ]
sel$index.cmd[4]

for(i in 1:nrow(sel)){system(sel$index.cmd[i])}
system(sel$index.cmd[4])

## picard mark dup
targets[, dupfile := paste0(path2, "/", ID, '_hg38.matrix') ]
targets[, bam.posi.sorted.rmdup.file := paste0(path2, "/", ID, '_hg38', "_posi_sorted_rmdup.bam") ]
targets[, rmdup.jobname := paste0("rmdup.", ID)]
targets[, rmdup.cmd := paste0(BSUB, " -J ", rmdup.jobname, " -e ", rmdup.jobname, ".err -o ", rmdup.jobname, ".std -cwd ", cwd, ' -We 2:59 -R "rusage[mem=50]" -R "rusage[iounits=0]" -n 12 ')]
targets[, rmdup.cmd := paste(rmdup.cmd, '-w "post_done(', index.jobname, ')"')]
targets[, rmdup.cmd := paste0(rmdup.cmd, ' "', java, ' -Xmx45g -Djava.io.tmpdir=', path2, ' -jar ', picard, ' MarkDuplicates REMOVE_DUPLICATES=TRUE METRICS_FILE=', dupfile, ' AS=TRUE INPUT=', bam.posi.sorted.file, ' OUTPUT=', bam.posi.sorted.rmdup.file, ' TMP_DIR=', path2, '"')]

sel = targets[batch == 'atac',] 
sel$rmdup.cmd[3]
system(sel$rmdup.cmd[3])
system(sel$rmdup.cmd[4])

write.table(sel$rmdup.cmd[3], file='a', quote=F, col.names=F, row.names=F)
write.table(sel$rmdup.cmd[4], file='a', append=T, quote=F, col.names=F, row.names=F)

for(i in 1:nrow(sel)){system(sel$rmdup.cmd[i])}

## filter bam file, based on fragment length
targets[, bam.posi.sorted.rmdup.fraglen.file := paste0(path2, "/", ID, '_hg38', "_posi_sorted_rmdup_fraglen.bam") ]
targets[, sambamba.jobname := paste0("sambamba.", ID)]
targets[, sambamba.cmd := paste0(BSUB, " -J ", sambamba.jobname, " -e ", sambamba.jobname, ".err -o ", sambamba.jobname, ".std ")]
targets[, sambamba.cmd := paste0(sambamba.cmd, "-cwd ", cwd, ' -We 2:59 -R "rusage[mem=15]" -R "rusage[iounits=0]" -n 12 ')]
targets[, sambamba.cmd := paste0(sambamba.cmd, ' -w "post_done(', rmdup.jobname, ')" ')]
targets[, sambamba.cmd := paste0(sambamba.cmd, '"', sambamba, ' view -t 9 -F \'template_length > -300 and template_length < 300 and proper_pair\' -f bam ')]
targets[, sambamba.cmd := paste0(sambamba.cmd, bam.posi.sorted.rmdup.file, ' -o ', bam.posi.sorted.rmdup.fraglen.file, '"')]

sel = targets[ batch == 'atac', ]
system(sel$sambamba.cmd[3])
system(sel$sambamba.cmd[4])
sel$sambamba.cmd
for(i in 1:nrow(sel)){system(sel$sambamba.cmd[i])}

# coordinate sort bam after sambamba
targets[, bam.posi.sorted.rmdup.sorted.file := paste0(path2, '/', ID, '_hg38_posi_sorted_rmdup_fraglen_sorted.bam')]
targets[, fraglensort.jobname := paste0("fraglensort.", ID)]
targets[, fraglensort.cmd := paste0(BSUB, " -J ", fraglensort.jobname, " -e ", fraglensort.jobname, ".err -o ", fraglensort.jobname, ".std ")]
targets[, fraglensort.cmd := paste0(fraglensort.cmd, "-cwd ", cwd, ' -We 2:29 -R "rusage[mem=28]" -R "rusage[iounits=0]" -n 11 ')]
targets[, fraglensort.cmd := paste(fraglensort.cmd, '-w "post_done(', sambamba.jobname, ')"')]
targets[, fraglensort.cmd := paste0(fraglensort.cmd, ' "', samtools, ' sort -@ 10 -m 2G ', bam.posi.sorted.rmdup.fraglen.file, ' -o ', bam.posi.sorted.rmdup.sorted.file, '"')]

sel = targets[ batch == 'atac', ]
sel$fraglensort.cmd

for(i in 1:nrow(sel)){system(sel$fraglensort.cmd[i])}
system(sel$fraglensort.cmd[3])
system(sel$fraglensort.cmd[4])
file.exists(sel$bam.posi.sorted.file)

# index2
targets[, index2.jobname := paste0("index2.", ID)]
targets[, index2.cmd := paste0(BSUB, " -J ", index2.jobname, " -e ", index2.jobname, ".err -o ", index2.jobname, ".std ")]
targets[, index2.cmd := paste0(index2.cmd, "-cwd ", cwd, ' -We 2:29 -R "rusage[mem=4]" -R "rusage[iounits=0]" -n 2 ')]
targets[, index2.cmd := paste(index2.cmd, '-w "post_done(', fraglensort.jobname, ')"')]
targets[, index2.cmd := paste(index2.cmd, ' "', samtools, 'index', bam.posi.sorted.rmdup.fraglen.sorted.file, '"')]
sel = targets[ batch == 'atac', ]
sel$index2.cmd[4]

for(i in 1:nrow(sel)){system(sel$index2.cmd[i])}
system(sel$index2.cmd[4])

## delete file after accident
sel = targets[ batch == '6582']
sel
tmp = paste0('rm -f ', sel$samfile)
tmp
for(i in 1:length(tmp)){system(tmp[i])}

fn = grep('bam', colnames(targets), value=T)
fn = grep('file', fn, value=T)
fn
for(i in fn){
	tmp = paste0('rm -f ', sel[, fn, with=F) 
	for(i in 1:length(tmp)){system(tmp[i])}
}

# call for peak
g = 'hs'
targets[, macs2.outname := paste0(path2, "/macs2_", ID, '_hg38') ]
targets[, macs2.jobname := paste0("macs2.", ID)]
targets[, macs2.cmd := paste0(BSUB, " -J ", macs2.jobname, " -e ", macs2.jobname, ".err -o ", macs2.jobname, ".std ")]
targets[, macs2.cmd := paste0(macs2.cmd, "-cwd ", cwd, ' -We 2:59 -R "rusage[mem=10]" -R "rusage[iounits=0]" -n 5 ')]
targets[, macs2.cmd := paste(macs2.cmd, '-w "post_done(', index2.jobname, ')"')]
targets[, macs2.cmd := paste(macs2.cmd, '"', macs2, 'callpeak -t', bam.posi.sorted.rmdup.fraglen.sorted.file, '-f BAMPE -g', g, '-n', macs2.outname, '-B -q 0.001"')]

sel = targets[ batch == 'atac']
for(i in 1:nrow(sel)){system(sel$macs2.cmd[i])}

# change the peak names in each file from the long file name to the ID
## here
targets[, npeakfile := paste0(path2, '/', 'macs2_', ID, '_hg38_summits_new.bed')]
sel = targets[batch == 'atac', ]
for(i in 1:nrow(sel)){
	pf = sel$peakfile[i]
	tmp = fread(pf)
	tmp[, V4 := paste0(sel$ID[i], '_', 1:nrow(tmp))]
	fwrite(tmp, file=sel$npeakfile[i], quote=F, col.names=F, row.names = F, sep="\t")
}

## merge peaks
## from homer
## lasttime
targets[, peakfile := paste0(path2, '/macs2_', ID, '_hg38_summits.bed')]
sel = targets[ batch == '6582' | batch == 'atac',]
nrow(sel)
peakfiles = paste(sel$npeakfile, collapse = " ")
file.exists(sel$npeakfile)
IDs = paste(sel$ID, collapse = " ")
peakfiles
IDs

if(all(file.exists(sel$npeakfile)){
	   system('mkdir -p merged')
	   cmd = paste0("mergePeaks -d 300 ", peakfiles, " -venn merged/merged.venn -matrix merged/merged.matrix -code > merged/merged.txt 2> mergepeaks.err")
	   system(cmd)
		    }

library(data.table)
fread("merged/merged.txt", colClasses = c("Parent files" = 'character')) -> merged.peaks
merged.peaks
setnames(merged.peaks, 1, 'PeakID')
setnames(merged.peaks, 7, 'Code')
setnames(merged.peaks, 9:ncol(merged.peaks), sel$ID)
#merged.peaks = merged.peaks[start > 0, ]
fwrite(merged.peaks, file='merged/merged.peaks.v1.txt', quote=F, col.names=F, row.names = F, sep="\t")
tmp = merged.peaks[,c(2,3,4,1,6,5,7,8)]
fwrite(tmp, file='merged/merged.peaks.v3.txt', quote=F, col.names=F, row.names = F, sep="\t") #chr start end name stat strand code subpeaks
merged.peaks[, chr := paste0('chr', chr)]
colnames(merged.peaks)
unique(merged.peaks$chr)
fwrite(merged.peaks, file='merged/merged.peaks.v2.txt', quote=F, col.names=F, row.names = F, sep="\t")
# v1 for peak score by bedgraph, since the chromosome name are the same  as '1'
# v2 for annotatePeaks.pl by which the hg38 version annotation use 'chr1'

## below for deeptools
## generate bigwig file
## bigwig -> wig -> feed to deeptools
## use bigwig file normalized data and use the generated wig file for annotatePeak.pl program
## to call the peak score
gsize = 2451960000
targets[, bigwigfile := paste0(path2, "/", ID, '_hg38.bigwig') ]
targets[, bw.jobname := paste0("bw.", ID)]
targets[, bw.cmd := paste0(BSUB, " -J ", bw.jobname, " -e ", bw.jobname, ".err -o ", bw.jobname, ".std ")]
targets[, bw.cmd := paste0(bw.cmd, "-cwd ", cwd, ' -We 2:59 -R "rusage[mem=10]" -R "rusage[iounits=0]" -n 5 ')]
targets[, bw.cmd := paste(bw.cmd, '-w "post_done(', index2.jobname, ')"')] 
targets[, bw.cmd := paste(bw.cmd, '"/home/huw/anaconda3/bin/bamCoverage -b', bam.posi.sorted.rmdup.fraglen.sorted.file)]
targets[, bw.cmd := paste(bw.cmd, ' --ignoreForNormalization chrX --normalizeTo1x', gsize, '--binSize 10 -o', bigwigfile, '"')]
sel = targets[batch == '6582' | batch == 'atac',]
sel$bw.cmd[1]
for( i in 1:nrow(sel)){system(sel$bw.cmd[i])}
nrow(sel)
targets$batch

## bigwig to bedgraph
## for annotatePeak.pl
targets[, bedgraphfile := paste0(path2, "/", ID, '_hg38.bdg') ]
targets[, bdg.jobname := paste0("bdg.", ID)]
targets[, bdg.cmd := paste0(BSUB, " -J ", bdg.jobname, " -e ", bdg.jobname, ".err -o ", bdg.jobname, ".std ")]
targets[, bdg.cmd := paste0(bdg.cmd, "-cwd ", cwd, ' -We 2:59 -R "rusage[mem=10]" -R "rusage[iounits=0]" -n 5 ')]
#targets[, bdg.cmd := paste(bdg.cmd, '-w "post_done(', bw.jobname, ')"')] 
targets[, bdg.cmd := paste0(bdg.cmd, ' "/home/huw/local/bin/bigWigToBedGraph ', bigwigfile, ' ', bedgraphfile, ' "')]
sel = targets[batch == '6582' | batch == 'atac',]
sel$bdg.cmd[1:4]
for( i in 1:nrow(sel)){system(sel$bdg.cmd[i])}

## flatstat
targets[, flagStatOutfile := paste0(path2, "/", ID, '_', "flagstat.txt")]
targets[, flagStat.jobname := paste0("flagStat.", ID)]
targets[, flagStat.cmd := paste0(BSUB, " -J ", flagStat.jobname, " -e ", flagStat.jobname, ".err -o ", flagStat.jobname, ".std ")]
targets[, flagStat.cmd := paste0(flagStat.cmd, "-cwd ", cwd, ' -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -n 1 ')]
targets[, flagStat.cmd := paste0(flagStat.cmd, ' "samtools flagstat ', bam.posi.sorted.rmdup.fraglen.sorted.file, ' > ', flagStatOutfile, ' "')]
sel = targets[batch == '6582' | batch == 'atac',]
dim(sel)

for( i in 1:nrow(sel)){system(sel$flagStat.cmd[i])}

## annotatePeak for each indiviual samples
## add chr to the chromosome
targets[, peakfilechr := paste0(path2, "/macs_", ID, '_hg38_summits_chr.bed') ]
sel = targets[batch == '6582' | batch == 'atac',]
for(i in 1:nrow(sel)){
	tmp = fread(sel$peakfile[i], header=F)
	tmp[, V1 := paste0('chr', V1)]
	fwrite(tmp, file=sel$peakfilechr[i], quote=F, col.names=F, row.names = F, sep="\t")
}
file.exists(sel$peakfilechr)

colnames(targets)
targets[, annotatedfile := paste0(path2, "/macs_",  ID, '_hg38_peak_annotated.txt')]
targets[, annstatsfile := paste0(path2, "/macs_",  ID, '_hg38_peak_annstats.txt')]
targets[, annotate.jobname := paste0("annotate.", ID)]
targets[, annotate.cmd := paste0(BSUB, " -J ", annotate.jobname, " -e ", annotate.jobname, ".err -o ", annotate.jobname, ".std ")]
targets[, annotate.cmd := paste0(annotate.cmd, "-cwd ", cwd, ' -We 2:59 -R "rusage[mem=10]" -R "rusage[iounits=0]" -n 5 ')]
targets[, annotate.cmd := paste0(annotate.cmd, ' "annotatePeaks.pl ', peakfilechr, ' -cpu 4 -log -mask -annStats ', annstatsfile)]
targets[, annotate.cmd := paste0(annotate.cmd, ' -bedGraph ', bedgraphfile, ' > ', annotatedfile, ' "')]
sel = targets[batch == '6582' | batch == 'atac',]
for( i in 1:nrow(sel)){system(sel$annotate.cmd[i])}

# homer findmotifgenoms.pl
colnames(targets)
targets[, findmotifdir := paste0(path2, '/motif')]
targets[, tmp := paste0('mkdir -p ', findmotifdir)]
sel = targets[batch == '6582' | batch == 'atac',]
sel$tmp[1]
for( i in 1:nrow(sel)){system(sel$tmp[i])}

targets[, findmotif.jobname := paste0("findmotif.", ID)]
targets[, findmotif.cmd := paste0(BSUB, " -J ", findmotif.jobname, " -e ", findmotif.jobname, ".err -o ", findmotif.jobname, ".std ")]
targets[, findmotif.cmd := paste0(findmotif.cmd, "-cwd ", cwd, ' -We 2:59 -R "rusage[mem=10]" -R "rusage[iounits=0]" -n 5 ')]
targets[, findmotif.cmd := paste0(findmotif.cmd, ' "findMotisGenome.pl ', peakfilechr, ' hg38 ', findmotifdir, ' -cpu 4 -mask "')]
sel = targets[batch == '6582' | batch == 'atac',]
sel$findmotif.cmd[3]

for( i in 1:nrow(sel)){system(sel$findmotif.cmd[i])}

sel = targets[batch == '6582' | batch == 'atac',]
sel$annotate.cmd[1:4]
for( i in 1:nrow(sel)){system(sel$annotate.cmd[i])}
## annotatePeak
##targets$bedGraph = paste0(targets$path, "/macs2_", targets$ID, "_treat_pileup.bdg")
##targets$bedGraph
##file.exists(targets$bedGraph)
##bdgfiles = paste(targets$bedGraph[sel], collapse = " ")

targets$wigfile
file.exists(sel$wigfile)
file.exists(targets$wigfile)
file.exists(targets$bedgraphfile)
dim(sel)
sel
bdgfiles = paste(sel$bedgraphfile, collapse = " ")
bdgfiles
cat("annotatePeaks.pl merged.peaks.v2.txt hg38 -cpu 10 ", file='merged/log_peakanno')
cat(" -log -mask ", file='merged/log_peakanno', append = T)
cat(" -annStats merged_annotate_Stats.txt > merged.annotation.v2.txt\n", file='merged/log_peakanno', append = T)

cat("annotatePeaks.pl merged.peaks.v1.txt hg38 -cpu 10 ", file='merged/log_peakanno', append=T)
cat(" -log -mask ", file='merged/log_peakanno', append = T)
cat(" -annStats merged_annotate_Stats.txt ", file='merged/log_peakanno', append = T)
cat(" -bedGraph ", bdgfiles, " > merged.annotation.v1.txt ", file='merged/log_peakanno', append = T)

system("cd merged & bash log_peakanno & cd ..")

## evluate the peaks by the bam file
targets[, bamCoverageFile := paste0(path2, '/merged.peaks.bam.count.txt') ]
targets[, coverage.jobname := paste0("coverage.", ID)]
targets[, coverage.cmd := paste0(BSUB, " -J ", coverage.jobname, " -e ", coverage.jobname, ".err -o ", coverage.jobname, ".std ")]
targets[, coverage.cmd := paste0(coverage.cmd, "-cwd ", cwd, ' -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -n 1 ')]
targets[, coverage.cmd := paste0(coverage.cmd, ' "/home/huw/program/bedtools-2.17.0/bin/multiBamCov -bams ', bam.posi.sorted.rmdup.fraglen.sorted.file, ' -bed merged/merged.peaks.v3.txt > ', bamCoverageFile, ' "')]
sel = targets[batch == '6582' | batch == 'atac',]
sel$coverage.cmd[3]

for( i in 1:nrow(sel)){system(sel$coverage.cmd[i])}

##
fread("merged/merged.annotation.v2.txt", colClasses = c("Focus Ratio/Region Size" = 'character')) ->  vv
setnames(vv, colnames(vv)[1], 'PeakID')
colnames(vv)
setkey(vv, PeakID)
for(i in 1:nrow(sel)){
	tt = fread(sel$bamCoverageFile[i], header=F)
	setnames(tt, colnames(tt)[4], 'peakID')
	setkey(tt, peakID)
	vv[tt[,.(peakID, V9)]] -> vv
}

sel$ID
setnames(vv, 20:41, sel$ID)

# normalize by the total reads in each sample
sel$flagStatOutfile
read.size = sapply(sel$flagStatOutfile, function(x){
       x = readLines(x)
       x = x[1]
       x = unlist(strsplit(x, " "))
       x = x[1]
       x
})
read.size = as.numeric(read.size)
names(read.size) = sel$ID
read.size

merged.peaks.anno = vv
for(i in 20:41){
	merged.peaks.anno[, i] = merged.peaks.anno[, i, with=F] * 1000000 / read.size[i-19]
}
merged.peaks.anno

merged.peaks.anno[, anno := 'Enh']
merged.peaks.anno[abs(`Distance to TSS`) < 3000, anno := 'Promoter']
colnames(merged.peaks.anno)

for(j in 20:41){
	set(merged.peaks.anno, which(merged.peaks.anno[[j]] < 1), j, 1)
}

merged.peaks.anno[, atac :=log2((K21ATAC + K22ATAC) / (RT41ATAC + RT42ATAC))]
merged.peaks.anno[, h3k27ac :=log2((RT4_K2_H3K27ac + RT4_K316_H3K27ac) / (2 * RT4_P_H3K27ac) )]
merged.peaks.anno[, h3k27me3 :=log2((RT4_K2_H3K27me3 + RT4_K316_H3K27me3) /( 2 * RT4_P_H3K27me3) )]
merged.peaks.anno[, kdm6a :=log2((RT4_K2_KDM6A + RT4_K316_KDM6A) / (2 * RT4_P_KDM6A)) ]
merged.peaks.anno[, ezh2 :=log2((RT4_K2_EZH2 + RT4_K316_EZH2) / (2 * RT4_P_EZH2) )]
merged.peaks.anno[, h3k4me3 :=log2((RT4_K2_H3K4me3 + RT4_K316_H3K4me3 ) / (2 * RT4_P_H3K4me3) )]

merged.peaks.anno[, isatac 	:= abs(atac) > 1]
merged.peaks.anno[, ish3k27ac 	:= abs(h3k27ac) > 1]
merged.peaks.anno[, ish3k27me3 	:= abs(h3k27me3) > 1]
merged.peaks.anno[, iskdm6a 	:= abs(kdm6a) > 1]
merged.peaks.anno[, isezh2 	:= abs(ezh2) > 1]
merged.peaks.anno[, ish3k4me3 	:= abs(h3k4me3) > 1]

merged.peaks.anno[, maxATAC := apply(merged.peaks.anno[,grep('\\dATAC', colnames(merged.peaks.anno)), with=F], 1, max)]
merged.peaks.anno[, maxKDM6A := apply(merged.peaks.anno[,grep('_KDM6A', colnames(merged.peaks.anno)), with=F], 1, max)]
merged.peaks.anno[, maxEZH2 := apply(merged.peaks.anno[,grep('_EZH2', colnames(merged.peaks.anno)), with=F], 1, max)]
merged.peaks.anno[, maxH3K4me3 := apply(merged.peaks.anno[,grep('_H3K4me3', colnames(merged.peaks.anno)), with=F], 1, max)]
merged.peaks.anno[, maxH3K27me3 := apply(merged.peaks.anno[,grep('_H3K27me3', colnames(merged.peaks.anno)), with=F], 1, max)]
merged.peaks.anno[, maxH3K27ac := apply(merged.peaks.anno[,grep('_H3K27ac', colnames(merged.peaks.anno)), with=F], 1, max)]

##is.tests = paste0('is.', tests)
##merged.peaks.anno[, isATAC := get.ss('ATAC')]
##merged.peaks.anno[, isKDM6A := get.ss('KDM6A')]
##merged.peaks.anno[, isEZH2 := get.ss('EZH2')]
##merged.peaks.anno[, isH3K4me3 := get.ss('H3K4me3')]
##merged.peaks.anno[, isH3K27me3 := get.ss('H3K27me3')]
##merged.peaks.anno[, isH3K27ac := get.ss('H3K27ac')]

table(merged.peaks.anno[, .(anno, isatac)])
table(merged.peaks.anno[, .(anno, ish3k27ac)])
table(merged.peaks.anno[, .(anno, ish3k27me3)])
table(merged.peaks.anno[, .(anno, ish3k4me3)])
table(merged.peaks.anno[, .(anno, iskdm6a)])
table(merged.peaks.anno[, .(anno, isezh2)])

merged.peaks.anno[, tag := paste0(Chr, ':', Start, '-', End)]
colnames(merged.peaks.anno)[7] = 'binCode'
merged.peaks.anno[, peakWidth := End - Start]
fwrite(merged.peaks.anno, file='merged/merged.peaks.anno.txt', sep="\t")

getwd()
library(ggplot2)
library(GGally)
g = ggpairs(merged.peaks.anno[, .(atac, h3k27ac, h3k27me3, h3k4me3, ezh2, kdm6a)])
ggsave(g, file="res/ggpair.pdf")

tmp = merged.peaks.anno[, .(h3k27ac, h3k27me3, h3k4me3, ezh2, kdm6a)]

cn = c('PeakID', 'Gene Name', 'Distance to TSS', 'anno', 'atac')
merged.peaks.anno[, cn, with=F]
mean(merged.peaks.anno$atac)

## ngsplot
Sys.setenv(NGSPLOT = "/home/huw/program/ngsplot/")

targets[, ngsplot.jobname := paste0("ngsplot.", ID)]
targets[, ngsplot.cmd := paste0(BSUB, " -J ", ngsplot.jobname, " -e ", ngsplot.jobname, ".err -o ", ngsplot.jobname, ".std ")]
targets[, ngsplot.cmd := paste0(ngsplot.cmd, "-cwd ", cwd, ' -We 2:59 -R "rusage[mem=5]" -R "rusage[iounits=0]" -n 2 ')]
#targets[, ngsplot.cmd := paste(ngsplot.cmd, '-w "post_done(', trim.jobname, ')"')]
targets[, ngsplot.cmd := paste0(ngsplot.cmd, ' "/home/huw/program/ngsplot/bin/ngs.plot.r -G hg38 -R tss -C ', bam.posi.sorted.rmdup.fraglen.sorted.file, ' -O ', path2, '/ngsplot_', ID, '_tss_ "')]
sel = targets[ batch == '6582' | batch == 'atac',]
sel$ngsplot.cmd[1]
for(i in 1:nrow(sel)){system(sel$ngsplot.cmd[i])}

targets[, ngsplot.jobname := paste0("ngsplot.", ID)]
targets[, ngsplot.cmd := paste0(BSUB, " -J ", ngsplot.jobname, " -e ", ngsplot.jobname, ".err -o ", ngsplot.jobname, ".std ")]
targets[, ngsplot.cmd := paste0(ngsplot.cmd, "-cwd ", cwd, ' -We 2:59 -R "rusage[mem=5]" -R "rusage[iounits=0]" -n 2 ')]
#targets[, ngsplot.cmd := paste(ngsplot.cmd, '-w "post_done(', trim.jobname, ')"')]
targets[, ngsplot.cmd := paste0(ngsplot.cmd, ' "/home/huw/program/ngsplot/bin/ngs.plot.r -G hg38 -R genebogy -C ', bam.posi.sorted.rmdup.fraglen.sorted.file, ' -O ', path2, '/ngsplot_', ID, '_genebody_ "')]
sel = targets[ batch == '6582' | batch == 'atac',]
for(i in 1:nrow(sel)){system(sel$ngsplot.cmd[i])}

targets[, ngsplot.jobname := paste0("ngsplot.", ID)]
targets[, ngsplot.cmd := paste0(BSUB, " -J ", ngsplot.jobname, " -e ", ngsplot.jobname, ".err -o ", ngsplot.jobname, ".std ")]
targets[, ngsplot.cmd := paste0(ngsplot.cmd, "-cwd ", cwd, ' -We 2:59 -R "rusage[mem=5]" -R "rusage[iounits=0]" -n 2 ')]
#targets[, ngsplot.cmd := paste(ngsplot.cmd, '-w "post_done(', trim.jobname, ')"')]
targets[, ngsplot.cmd := paste0(ngsplot.cmd, ' "/home/huw/program/ngsplot/bin/ngs.plot.r -G hg38 -R cgi -C ', bam.posi.sorted.rmdup.fraglen.sorted.file, ' -O ', path2, '/ngsplot_', ID, '_cgi_ "')]
sel = targets[ batch == '6582' | batch == 'atac',]
for(i in 1:nrow(sel)){system(sel$ngsplot.cmd[i])}


#chipqc
#qc = copy(sel)
#setNames(qc, 'ID', 'SampleID')
#sel$Tissue = 'cells'
#colnames(sel)
#qc
#target_chipseq = read.table("targets", sep="\t", header = T)
#chipqc_chipseq = ChIPQC(target_chipseq, annotaiton="hg19")


## deeptools
bin.2.dec = function(x) sum(2^(which(rev(unlist(strsplit(as.character(x), '')) == 1)) -1 ))

merged.peaks.anno$bincode = sapply(merged.peaks.anno$Code, bin.2.dec)
merged.peaks.anno[, Chr := sub('chr', '', Chr)]
colnames(merged.peaks.anno)[7]
setNames(merged.peaks.anno, colnames(merged.peaks.anno)[7], 'binCode')
get.ss = function(x){
	cc = rep('0', 22)
	cc[grep(paste0('\\d', x), colnames(merged.peaks.anno))-19] = '1'
	cc = paste0(cc, collapse='')
	cc = bin.2.dec(cc)
	ss = sapply(merged.peaks.anno$bincode, function(x){bitwAnd(cc,x) > 1}) # 
	ss
}
get.ss('ATAC')
colnames(merged.peaks.anno)
merged.peaks.anno[,7]

colnames(merged.peaks.anno)

colorlist = '\"blue,white,red\"'
tests = c('ATAC', 'KDM6A', 'H3K4me3', 'H3K27me3', 'H3K27ac', 'EZH2')
for(i in tests){
	tmp = merged.peaks.anno[get.ss(i),  ]
	tmp
	bedfn = paste0(i, '.peaks.bed')
	fwrite(tmp[,c(2,3,4,1,5)], file=paste0('merged/', bedfn), quote=F, col.names=F, row.names = F, sep="\t")
	tt = sel[grep(i, sel$ID),]; bwf = tt$bigwigfile
	cal.matrix(bedfn, bwf, 'merged/cmd.sh', app=T)
	plotheatmap.wrap(paste0(i, '.peaks.gz'), 3, colorlist, tt$ID, 'merged/cmd.sh', T)
}

cal.matrix = function(bedfile, bwfiles, cmdfile, app=F){
	outfile = ''
	bwfiles = paste(bwfiles, sep=" ")
	if(length(bedfile) > 1){
		outfile=sub(".bed", paste0("_", length(bedfile)), bedfile[1])
	}else{
		outfile=sub(".bed", "", bedfile[1])
	}
	outfile.gz = paste0(outfile, ".gz")
	outfile.bed = paste0(outfile, "_matrix_norm.bed")
	cat('computeMatrix reference-point -R ', bedfile, file = cmdfile, append = app)
	cat(' -S ', bwfiles, file = cmdfile, append = T)
	cat(' -b 500 -a 500 -out ', outfile.gz, ' --referencePoint=center --missingDataAsZero --skipZeros', file = cmdfile, append = T)
	cat(' --outFileSortedRegions ', outfile.bed, file=cmdfile, append = T)
	#cat(' --blackListFileName ', blacklistfile, file=cmdfile, append = T)
	cat(" \n\n", file=cmdfile, append = T)
}
plotheatmap.wrap = function(matrixfile, kmean, colorlist, samplelabel, cmdfile, app){
	samplelabel = paste(samplelabel, collapse=" ")
	basename = sub(".gz", "", matrixfile)
	outfile = paste0(basename, "_kmeans", kmean, ".pdf")
	outmatrixfile = paste0(basename, "_kmeans", kmean, "_out.txt")
	outbedfile = paste0(basename, "_kmeans", kmean, "_out.bed")
	cat('plotHeatmap --matrixFile ', matrixfile, ' ', file = cmdfile, append = app, sep='')
	cat('--outFileSortedRegions ', outbedfile, ' ', file=cmdfile, append = T, sep='')
	cat('--outFileNameMatrix ', outmatrixfile,' ', file=cmdfile, append = T, sep='')
	cat('--refPointLabel "Center" ', file = cmdfile, append = T, sep='')
	cat('-out ', outfile, ' ', file = cmdfile, append = T, sep='')
	#cat('--zMin 0 --zMax 10 10 30 30 ', file=cmdfile, append = T, sep='')
	cat('--zMin 0 ', file=cmdfile, append = T, sep='')
	cat('--yMin 0 ', file=cmdfile, append = T, sep='')
	cat('--colorList ', colorlist, ' ', file=cmdfile, append = T, sep='')
	cat('--samplesLabel ', samplelabel, ' ', file=cmdfile, append = T, sep='')
	if(kmean > 0){
		cat(' --kmeans ', kmean, ' ', file = cmdfile, append = T, sep = '')
	}
	cat(" \n\n", file=cmdfile, append = T)
}


