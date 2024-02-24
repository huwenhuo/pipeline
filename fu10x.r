library(data.table)
library(openxlsx)
library(Seurat)
library(googlesheets4)
library(ggplot2)
library(dplyr)

options(width=222)

gs4_deauth()
dd = read_sheet('https://docs.google.com/spreadsheets/d/1lnZiGCWl2KfLpzVnNkJox0I9FKEEKwVJJ5ln7BHvRe0/edit?usp=sharing')
dd = as.data.table(dd)
dd = dd[4:20, .(patient, `Tube Code`, celltype, sample.dir, Gender, age, Dx, basedir)]
dd[, var.name := paste0(patient, '.', celltype)]
dd

dir.exists(dd$sample.dir)

# merge two pt7 chond samples
mtx.1 = Read10X(dd[13, sample.dir])
mtx.2 = Read10X(dd[14, sample.dir])
mtx.1 = CreateSeuratObject(counts = mtx.1, min.cells = 20, min.features=1000, project='pt7.1')
mtx.2 = CreateSeuratObject(counts = mtx.2, min.cells = 20, min.features=1000, project='pt7.2')
mtx.o <- merge(mtx.1, y = c(mtx.2), add.cell.ids = c("pt7.1", 'pt7.2'), project = "pt7")
fname = paste0('pt7.chond.RData')
save(mtx.o, file=fname)

# merge two pt11 chond samples
mtx.1 = Read10X(dd[16, sample.dir])
mtx.2 = Read10X(dd[17, sample.dir])
mtx.1 = CreateSeuratObject(counts = mtx.1, min.cells = 20, min.features=1000, project='mtx.1')
mtx.2 = CreateSeuratObject(counts = mtx.2, min.cells = 20, min.features=1000, project='mtx.2')
mtx.o <- merge(mtx.1, y = c(mtx.2), add.cell.ids = c("pt11.1", 'pt11..2'), project = "pt11")
fname = paste0('pt11.chond.RData')
save(mtx.o, file=fname)

d2 = dd[!duplicated(var.name), ]
d2

for(i in 1:nrow(d2)){
	var.name = d2[i, var.name]
	sample.dir = d2[i, sample.dir]
	if(!exists(var.name)){
		cat('Dealing with', var.name, '\n')
		mtx = Read10X(sample.dir)
		mtx.o = CreateSeuratObject(counts = mtx, min.cells = 20, min.features=1000, project=var.name)
		assign(var.name, mtx.o)
		fname = paste0(var.name, '.RData')
		save(mtx.o, file=fname)
	}
}

d2[, rdata := paste0(var.name, '.RData')]
file.exists(d2$rdata)

dir.pre = 'res'

d3 = d2[celltype == 'chond', ]
d3

d3[, {
	#var.name = 'pt11.chond'
	st = get(var.name)
	st[["percent.mt"]] <- PercentageFeatureSet(st, pattern = "^MT-")
	st.fil = subset(x= st, subset = percent.mt < 10)

	fname = paste0(dir.pre, '/', var.name, '_live_dead.pdf')
	pdf(fname)
	VlnPlot(st, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
	VlnPlot(st.fil, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
	dev.off()

	st.fil = SCTransform(st.fil, verbose=T)
	st.fil <- RunPCA(object = st.fil, verbose = FALSE)
	st.fil <- FindNeighbors(object = st.fil, dims = 1:20, verbose = FALSE)
	st.fil <- RunTSNE(object=st.fil, dims=1:5,  tsne.method='FIt-SNE', nthreads=50, max_iter=1000, fast_tsne_path='~/program/FIt-SNE/bin/fast_tsne')
	st.fil <- FindClusters(object = st.fil, verbose = FALSE)
	st.fil

	fname = paste0('rdata/', var.name, '.Rdata')
	save(st.fil, file=fname)
}, by=1:nrow(d3) ]

d3[, {
	fname = paste0('rdata/', var.name, '.Rdata')
	load(fname)

	gg <- FeatureScatter(st.fil, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
	fname = paste0(dir.pre, '/', var.name, '_nn.pdf'); fname
	ggsave(gg, file=fname, width=5, height=4)

	mm <- FindAllMarkers(object = st.fil, only.pos = T, min.pct = 0, logfc.threshold = log(2))
	fname = paste0('rdata/', var.name, '_findAllMarkers.Rdata'); fname
	save(mm, file=fname)

	topn = mm %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
	fname = paste0(dir.pre, '/', var.name, '_markers2.png'); fname
	png(fname, width=10*300, height=3*300)
	DoHeatmap(object = st.fil, features = topn$gene, raster=F)# + NoLegend()
	dev.off()

	fname = paste0(dir.pre, '/', var.name, 'st.fil_dimplot_orig.pdf'); fname
	gg = DimPlot(st.fil, group.by='orig.ident', reduction = 'tsne', label=T)
	ggsave(gg, file=fname, width=5, height=4)

	fname = paste0(dir.pre, '/', var.name, '_st.fil_dimplot.pdf'); fname
	gg = DimPlot(st.fil, reduction = 'tsne', label=T)
	ggsave(gg, file=fname, width=5, height=4)

	#ff = c('PTPRC', 'CD3G', 'GNLY', 'KRT19', 'EPCAM', 'CD14', 'MS4A2', 'COL1A1')
	#fname = paste0(dir.pre, '/', var.name, '_featureplot_', paste0(ff, collapse='_'), '_v2.pdf'); fname
	#gg = FeaturePlot(st.fil, features = ff, ncol = 3, label=T, reduction='tsne') + NoLegend()
	#ggsave(gg, file=fname, width=11, height=3)

	ff = c('COL2A1', 'SOX9', 'ACAN', 'DCN', 'VCAN')
	fname = paste0(dir.pre, '/', var.name, '_featureplot_', paste0(ff, collapse='_'), '_v2.pdf'); fname
	gg = FeaturePlot(st.fil, features = ff, ncol = 3, label=T, reduction='tsne') + NoLegend()
	ggsave(gg, file=fname, width=11, height=6)

	ff = c('GRN', 'YWHAH', 'TNFRSF1B', 'SCN9A')
	fname = paste0(dir.pre, '/', var.name, '_featureplot_', paste0(ff, collapse='_'), '_v2.pdf'); fname
	gg = FeaturePlot(st.fil, features = ff, ncol = 2, label=T, reduction='tsne') + NoLegend()
	ggsave(gg, file=fname, width=8, height=6)

	ff = c('CCL20', 'SPP1', 'CXCL8', 'CRYAB', 'CHI3L2', 'S100A2', 'HBB', 'COL2A1', 'IGFBP3')
	fname = paste0(dir.pre, '/', var.name, '_featureplot_', paste0(ff, collapse='_'), '_v2.png'); fname
	gg = FeaturePlot(st.fil, features = ff, ncol = 3, label=T, reduction='tsne') + NoLegend()
	ggsave(gg, file=fname, width=11, height=6) }, by=1:nrow(d3) ]
