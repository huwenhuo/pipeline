
## chip-seq from penn state
{


	chip.dir = 'penn'
	chip.dsn = data.table(fastq = system(paste0('find ', chip.dir, ' -name *.fastq.gz'), intern=T))
	chip.dsn[, sample.id := str_replace(fastq, ".*Penn_(.*)_hs.*", "\\1")]
	chip.dsn = chip.dsn[grep('G432', sample.id, invert=T), ]
	chip.dsn[, sample.id := sub('-', '_', sample.id)]
	chip.dsn = chip.dsn[!duplicated(basename(fastq)), ][order(sample.id), ]
	chip.dsn[, .(sample.id)]
	chip.dsn[, sample.name := c('Ctrl1', 'Ctrl2', 'FOXA1KO1', 'FOXA1KO2', 'Pooled2', 'Pooled1', 'UMUC1FOXA11', 'UMUC1FOXA12')]
	chip.dsn[, ControlID := c(rep('Pooled1', 6), 'Pooled2', 'Pooled2']
	chip.dsn[, .(sample.id, sample.name, ControlID)]
	chip.dsn


	chip.dsn[, bam.file := paste0('penn/data/bamfiles/', sample.name, '.bam')]
	chip.dsn[, sam.file := paste0('penn/data/bamfiles/', sample.name, '.sam')]
	chip.dsn[, sorted.bam.file := paste0('penn/data/bamfiles/', sample.name, '_sorted.bam')]
	file.exists(chip.dsn$sorted.bam.file)

	chip.dsn[, btw.jobname := paste0('bw.', sample.id)]
	chip.dsn[, bam.file := paste0('penn/', sample.id, '.bam')]
	chip.dsn[, btw.cmd := bsub.head(cpu=10, jobname=btw.jobname, mem=2, hmem=30), by=1:nrow(chip.dsn)]
	chip.dsn[, btw.cmd := paste0(btw.cmd,  ' " bowtie2 --local -U ', fastq, ' -x ', cfg$btw2db)]
	chip.dsn[, btw.cmd := paste0(btw.cmd,  ' -N 1 --norc -p 5 > ', sam.file, '"')]
	chip.dsn[, btw.cmd]

	write(chip.dsn$btw.cmd, file='btw.sh')
	exe.jobs(chip.dsn[1:6, btw.cmd], logger, if.wait=F)

	## sort sam
	chip.dsn[, sorted.bam.file := paste0('penn/data/bamfiles/', sample.name, '_sorted.bam')]
	chip.dsn[, bam.sort.jobname := paste0('bam.sort.', sample.name)]
	chip.dsn[, bam.sort.cmd := bsub.head(cpu=6, mem=5, hmem =32, jobname=bam.sort.jobname, We='9999', postdone=btw.jobname), by=1:nrow(chip.dsn)]
	chip.dsn[, bam.sort.cmd := paste0(bam.sort.cmd, " \"", cfg$JAVA, '/java -Xms256M -Xmx32g -Djava.io.tmpdir=', cfg$tmpdir, ' -jar ')]
	chip.dsn[, bam.sort.cmd := paste0(bam.sort.cmd, cfg$PICARD, ' SortSam I=', sam.file, " O=", sorted.bam.file)]
	chip.dsn[, bam.sort.cmd := paste0(bam.sort.cmd, " CREATE_INDEX=true SORT_ORDER=coordinate\"")]
	chip.dsn[, bam.sort.cmd]
	write(chip.dsn[, bam.sort.cmd], file='sortsam.sh')

	exe.jobs(chip.dsn[1:6, bam.sort.cmd], logger, if.wait=F)

	chip.dsn[, sample.dir := 'penn/data/']
	chip.dsn[, controlID := c('Pooled1', 'Pooled1', 'Pooled1', 'Pooled1', '', '', 'Pooled2', 'Pooled2')]
	chip.dsn[, bamControl := paste0(dirname(sorted.bam.file), '/', controlID, '_sorted.bam')]
	chip.dsn[, .(sample.name, sorted.bam.file, controlID, bamControl)]

	# call for peak
	# bedtools intersect -v -abam FILE.BAM -b BLACKLIST.BED > FILTERED.BAM
	chip.dsn[, peak.file := paste0(sample.dir, '/macs2_', sample.name, '_peaks.narrowPeak')]
	chip.dsn[, macs2.jobname := paste0("macs2.", sample.id)]
	#chip.dsn[, macs2.cmd := bsub.head(jobname = macs2.jobname, mem=20, hmem=20, cpu=5, We='30:30', postdone=bam.sort.jobname), by=1:nrow(chip.dsn)]
	chip.dsn[, macs2.cmd := paste0(' macs2 callpeak -f BAM --nomodel -g hs -t ', sorted.bam.file, ' -c ', bamControl, ' --outdir ', sample.dir, ' -n macs2_', sample.name, ' -B -q 0.001')]
	chip.dsn[, macs2.cmd]

	exe.jobs(chip.dsn[Tissue != '', macs2.cmd], logger, if.wait = F)

	chip.dsn[, macs2.cmd := paste0( ' macs2 callpeak -f BAM --nomodel -g hs -t ', sorted.bam.file, ' -c ', bamControl, ' --outdir ', sample.dir, ' -n macs2_', sample.name, ' -B -q 0.001')]
	write(chip.dsn[Tissue != '', macs2.cmd], file='peak.sh')

	# call for broad peak
	chip.dsn[, peak.brd.file := paste0(sample.dir, '/macs_broad_', sample.name, '_peaks.broadPeak')]
	chip.dsn[, macs2.broad.outname := paste0(sample.dir, "/", sample.name, "_macs2_hg19_broad") ]
	chip.dsn[, macs2.broad.jobname := paste0("macs2.brd.", sample.name)]
	chip.dsn[, macs2.broad.cmd := bsub.head(jobname = macs2.broad.jobname, mem=20, hmem=20, cpu=5, We='30:30', postdone = bam.sort.jobname), by=1:nrow(chip.dsn)]
	chip.dsn[, macs2.broad.cmd := paste0(' macs2 callpeak -f BAM --broad --nomodel -g hs -t ', sorted.bam.file, ' -c ', bamControl, ' --outdir ', sample.dir, ' -n macs2_broad_', sample.name, ' -B -q 0.001')]
	chip.dsn[, macs2.broad.cmd]

	exe.jobs(chip.dsn[Tissue != '', macs2.broad.cmd], logger)

	write(chip.dsn[Tissue != '', macs2.broad.cmd], file='peak_broad.sh')

	## ChIPQC
	chip.dsn[, Replicate := 1]
	chip.dsn[, peakcaller := 'narrow']
	chip.dsn[, Tissue := c('WT', 'WT', 'KO', 'KO', '', '', 'WT', 'WT')]
	chip.dsn[, Condition := NA]
	chip.dsn[, Factor := c('H3K27ac', 'H3K27ac', 'H3K27ac', 'H3K27ac', '', '', 'FOXA1', 'FOXA1')]

	save(chip.dsn, file='chip_dsn.rdata')
	load('chip_dsn.rdata')

	tmp = chip.dsn[Tissue != '', .(sample.name, Tissue, Replicate, sorted.bam.file, ControlID, bamControl, peak.file, peakcaller, Condition, Factor)]
	setnames(tmp, c('SampleID', 'Tissue', 'Replicate', 'bamReads', 'ControlID', 'bamControl', 'Peaks', 'PeakCaller', 'Condition', 'Factor'))
	tmp = setDF(tmp)
	tmp$SampleID
	tmp$ControlID
	dim(tmp)
	tmp 

	chipObj <- ChIPQC::ChIPQC(tmp, annotation="hg19") 
	ChIPQC::ChIPQCreport(chipObj, reportName="SCC", reportFolder="r_fang/ChIPQCreport")

	exe.jobs(chip.dsn[, bam.sort.cmd], logger)

	# call for peak
	chip.dsn[, macs2.outname := paste0("penn/macs2_", sample.id, '_hg19') ]
	chip.dsn[, macs2.jobname := paste0("macs2.", sample.id)]
	chip.dsn[, macs2.cmd := bsub.head(jobname = macs2.jobname, mem=20, hmem=20, cpu=5, We='30:30', postdone=index.jobname), by=1:nrow(chip.dsn)]
	chip.dsn[, macs2.cmd := paste0(macs2.cmd, '" macs2 callpeak -f BAM --nomodel -g hsa -t ', sorted.bam.file, ' --outdir penn/data -n ', sample.id, ' -B -q 0.001"')]
	chip.dsn[, macs2.cmd]

	exe.jobs(chip.dsn[, macs2.cmd], logger)

	# call for broad peak
	chip.dsn[, macs2.broad.outname := paste0("penn/data/macs2_", sample.id, '_hg19_broad') ]
	chip.dsn[, macs2.broad.jobname := paste0("macs2.", sample.id)]
	chip.dsn[, macs2.broad.cmd := bsub.head(jobname = macs2.broad.jobname, mem=20, hmem=20, cpu=5, We='30:30'), by=1:nrow(chip.dsn)]
	chip.dsn[, macs2.broad.cmd := paste0(macs2.broad.cmd, '" macs2 callpeak -f BAM --nomodel --broad -g hs -t ', sorted.bam.file, ' --outdir penn/data -n ', sample.id, '_broad -B -q 0.001"')]
	chip.dsn[, macs2.broad.cmd]

	exe.jobs(chip.dsn[, macs2.broad.cmd], logger)

	f1 = fread('/juno/work/solitlab/huw/solit/study/hiseq/tf3/r_fang/initFiles/FOXA1ChIP/macs2_FOXA1ChIP_peaks.narrowPeak')
	f2 = fread('/juno/work/solitlab/huw/solit/study/hiseq/blca_cmo_06155_2016/penn/data/macs2_UMUC1FOXA11_peaks.narrowPeak')
	setnames(f1, c('chromosome', 'start', 'end', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak'))
	setnames(f2, c('chromosome', 'start', 'end', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak'))
	f1.gr = GenomicRanges::makeGRangesFromDataFrame(f1)
	f2.gr = GenomicRanges::makeGRangesFromDataFrame(f2)
	f12 = ChIPpeakAnno::findOverlapsOfPeaks(f1.gr, f2.gr)
	f12
	



	## bigwig to bedgraph
	## for annotatePeak.pl

	# generate bigwig by deeptools
	use_condaenv('dt')

	gsize = 2451960000
	chip.dsn[, bigwig.file := paste0("penn/", sample.id, '_hg19_normCPM.bigwig') ]
	chip.dsn[, bw.jobname := paste0("bw.", sample.id)]
	chip.dsn[, bw.cmd := bsub.head(jobname=bw.jobname, cpu=5, mem=30, hmem=30)]
	chip.dsn[, bw.cmd := paste0(bw.cmd, '" bamCoverage -b ', sorted.bam.file)]
	chip.dsn[, bw.cmd := paste0(bw.cmd, ' --ignoreDuplicates --minMappingQuality 10 --numberOfProcessors 4 --binSize 10 -o ', bigwig.file, ' --extendReads 200 --normalizeUsing CPM "')]
	exe.jobs(chip.dsn[, bw.cmd], logger)
	write(chip.dsn[, bw.cmd], 'a.sh')

	chip.dsn[, bedgraph.file := paste0("penn/", sample.id, '_hg19.bdg') ]
	chip.dsn[, bdg.jobname := paste0("bdg.", sample.id)]
	chip.dsn[, bdg.cmd := bsub.head(jobname = bdg.jobname, mem=20, hmem=20, cpu=5)]
	chip.dsn[, bdg.cmd := paste0(bdg.cmd, ' "bigWigToBedGraph ', bigwig.file, ' ', bedgraph.file, ' "')]
	write(chip.dsn[, bdg.cmd], 'b.sh')

	exe.jobs(chip.dsn[, bdg.cmd], logger)

	## flatstat
	chip.dsn[, flagStatOutfile := paste0(path2, "/", ID, '_', "flagstat.txt")]
	chip.dsn[, flagStat.jobname := paste0("flagStat.", ID)]
	chip.dsn[, flagStat.cmd := paste0(BSUB, " -J ", flagStat.jobname, " -e ", flagStat.jobname, ".err -o ", flagStat.jobname, ".std ")]
	chip.dsn[, flagStat.cmd := paste0(flagStat.cmd, "-cwd ", cwd, ' -We 0:59 -R "rusage[mem=2]" -R "rusage[iounits=0]" -n 1 ')]
	chip.dsn[, flagStat.cmd := paste0(flagStat.cmd, ' "samtools flagstat ', bam.posi.sorted.rmdup.fraglen.sorted.file, ' > ', flagStatOutfile, ' "')]
	sel = chip.dsn[batch == '6582' | batch == 'atac',]
	dim(sel)

	for( i in 1:nrow(sel)){system(sel$flagStat.cmd[i])}

	## annotatePeak for each indiviual samples
	## add chr to the chromosome
	chip.dsn[, peakfilechr := paste0(path2, "/macs_", ID, '_hg38_summits_chr.bed') ]
	sel = chip.dsn[batch == '6582' | batch == 'atac',]
	for(i in 1:nrow(sel)){
		tmp = fread(sel$peakfile[i], header=F)
		tmp[, V1 := paste0('chr', V1)]
		fwrite(tmp, file=sel$peakfilechr[i], quote=F, col.names=F, row.names = F, sep="\t")
	}
	file.exists(sel$peakfilechr)


	chip.dsn[, peak.file := paste0("penn/data/", sample.id, '_peaks.narrowPeak')]
	chip.dsn[, annstats.file := paste0("penn/", sample.id, '_hg19_peak_annstats.txt')]
	chip.dsn[, annotated.file := paste0("penn/", sample.id, '_hg19_peak_annotated.txt')]
	chip.dsn[, annotate.jobname := paste0("annotate.", sample.id)]
	chip.dsn[, annotate.cmd := bsub.head(jobname = annotate.jobname,  mem=20, cpu=5, hmem=20, postdone = macs2.jobname)]
	chip.dsn[, annotate.cmd := paste0(annotate.cmd, ' "annotatePeaks.pl ', peak.file, ' hg19 -log -mask -annStats ', annstats.file)]
	chip.dsn[, annotate.cmd := paste0(annotate.cmd, ' -bedGraph ', bedgraph.file, ' -genomeOntology penn/genomeOntology > ', annotated.file, ' "')]
	exe.jobs(chip.dsn[, annotate.cmd], logger)

	## transform for DiffBind
	chip.dsn[, peak.file.2 := sub('$', '.2', peak.file)]
	chip.dsn[, {tmp = fread(peak.file)
		 tmp[, V1 := paste0('chr', V1)]
		 fwrite(tmp[, c(1,2,3,5)], file=peak.file.2, sep="\t", col.names=F)}, by=1:nrow(chip.dsn)]
	file.exists(chip.dsn[, peak.file.2])

	peak.all = lapply(chip.dsn$peak.file, fread)
	names(peak.all) = chip.dsn$sample.id
	peak.all = rbindlist(peak.all, idcol=T)
	peak.all
	
	chip.dsn[, {
		findmotif.dir = paste0('penn/data/motif_', sample.id)
		sync(findmotif.dir)
		 }, by=1:nrow(chip.dsn)]

	chip.dsn[, {
		findmotif.dir = paste0('penn/data/motif_', sample.id)
		dir.create(findmotif.dir)
		findmotif.jobname = paste0("findmotif.", sample.id)
		findmotif.cmd = bsub.head(jobname = findmotif.jobname, mem=20, hmem=20, cpu=6)
		findmotif.cmd = paste0(findmotif.cmd, ' "findMotifsGenome.pl ', peak.file.2, ' hg19 ', findmotif.dir, ' -cpu 4 -mask "')
		system(findmotif.cmd) }, by=1:nrow(chip.dsn)]

	library(DiffBind)

	db.dsn = chip.dsn[c(1,3,4,5), ]
	db.dsn
	db.dsn[, SampleID := sample.id]
	#db.dsn[, Tissue := 'umuc1]'
	#db.dsn[, Factor := '']
	db.dsn[, Condition := c('wt', 'wt', 'ko', 'ko')]
	db.dsn[, Treatment := c('wt', 'wt', 'ko', 'ko')]
	db.dsn[, Replicate := c(1, 2, 1, 2)]
	db.dsn[, bamReads := sorted.bam.file]
	db.dsn[, ControlID := 'wt']
	db.dsn[, bamControl := db.dsn[1, sorted.bam.file]]
	db.dsn[, Peaks := peak.file.2]
	db.dsn[, PeakCaller := 'MACS2']
	db.dsn

	dbi = dba(sampleSheet = db.dsn)

	dbi.2 = dba(sampleSheet = db.dsn)
		dba.blacklist(blacklist = DBA_BLACKLIST_HG19) %>%
		dba.count()     %>%
		dba.normalize() %>%
		dba.contrast()  %>%
		dba.analyze()


	fname = 'a.pdf'
	pdf(file=fname)
	plot(dbi)
	dev.off()
	scp(fname)

	dbi = dba.count(dbi)#  summits=250)

	fname = 'b.pdf'
	pdf(file=fname)
	plot(dbi)
	dev.off()
	scp(fname)

	dbi = dba.contrast(dbi, minMembers = 2, categories = DBA_CONDITION)

	dbi = dba.analyze(dbi)

	fname = 'c.pdf'
	pdf(file=fname)
	plot(dbi)
	dev.off()
	scp('c.pdf')

	fname = 'dbi.RData'
	load(fname)
	save(dbi, file=fname)
	dbi

	dbi.report = dba.report(dbi)
	dbi.report

	dbi.report.dt = as.data.table(dbi.report)
	dbi.report.dt[, logp := -log10(FDR)]
	dbi.report.dt = dbi.report.dt[logp > 4, ]
	dbi.report.dt[, chr := paste0('chr', seqnames)]
	dbi.report.dt[, tag := paste0(chr, '_', start, '_', end)]
	dbi.report.dt

	fname = 'dbi_report.tsv'
	fwrite(dbi.report.dt, file=fname, sep="\t")
	fread(fname)
	fname = 'dbi_report.bed'
	fwrite(dbi.report.dt[, .(chr, start, end, tag, logp, strand)], file=fname, sep="\t", col.names=F)

	## ** annotate peaks
	## deeptools clustering peaks
	## RNA-Seq to the top up/dn clusters
	## find motifs in the top up/dn clusters
	fname = 'dbi_report.bed'
	annstats.file = paste0("penn/summary_hg19_peak_annstats.txt")
	annotated.file = paste0("penn/dbi_report_annotated.tsv")
	annotate.jobname = paste0("annotate.sum")
	#annotate.cmd = bsub.head(jobname = annotate.jobname,  mem=20, cpu=5, hmem=20)
	annotate.cmd = paste0(' annotatePeaks.pl ', fname, ' hg19 -log -mask -annStats ', annstats.file)
	annotate.cmd = paste0(annotate.cmd, ' -genomeOntology penn/genomeOntology > ', annotated.file)
	annotated.cmd = bsub.head(jobname = annotate.jobname, mem=20, hmem=20, cpu=2)
	annotate.cmd = paste0(annotated.cmd, '"', annotate.cmd, '"')
	annotate.cmd
	exe.jobs(annotate.cmd, logger)

	dbi.anno = fread(annotated.file)
	dbi.anno
	setnames(dbi.anno, 1, 'id')
	dbi.anno[, anno := 'genebody']
	dbi.anno[grep('promoter', `Annotation`), anno := 'promoter']
	dbi.anno[grep('intron', `Annotation`), anno := 'intron']
	dbi.anno[grep('Intergenic', `Annotation`), anno := 'intergenic']
	dbi.anno[, tag := paste0(Chr, '_', Start-1, '_', End)]
	dbi.anno = merge(dbi.anno, dbi.report.dt, by='tag', all.x=T)
	table(dbi.anno$anno)
	nrow(dbi.anno)


	fname = 'penn/data/dbi_anno.RData'
	save(dbi.anno, file=fname)

	updns = c('up', 'dn', 'updn')
	for(updn in updns){
		findmotif.dir = paste0('penn/motif_h3k27ac_', updn, '_logp4')
		sync(findmotif.dir)
	}

	fname = 'penn/pie.pdf'
	pdf(fname)
	pie( table(dbi.anno$anno) )
	dev.off()
	scp(fname)

	## annotate peaks
	## ** deeptools clustering peaks
	## RNA-Seq to the top up/dn clusters
	## find motifs in the top up/dn clusters

	## deeptools
	## wait for the annotation
	annos = c('intron', 'intergenic', 'promoter', 'genebody', 'all')
	dd.dsn = data.table(annos = rep(annos, each=3), kk = rep(3:5, times=5))
	dd.dsn[, mtx.in.fname := paste0('penn/data/dbi_mtx_', annos, '.tsv')]
	dd.dsn[, mtx.file := paste0('penn/data/dbi_mtx_', annos, '.gz')]
	dd.dsn[, dp.heatmap.file := paste0('penn/', annos, '_deep_tools_heatmap_kmean', kk, '_v2.pdf')]
	dd.dsn[, dp.out.cluster.file := paste0('penn/data/', annos, '_out_cluster_mtx_kmean', kk, '.tsv')]
	dd.dsn[, bigwig.file := paste0("penn/", sample.id, '_hg19_normCPM.bigwig') ]

	dd.dsn[, {
		if(annos == 'all'){
			tmp = dbi.anno
		}else{
			tmp = dbi.anno[anno == annos, ]
		}
		fwrite(tmp[, .(chr, Start, End, Strand, id)], file=mtx.in.fname, sep="\t", col.names=F)
		dp.mtx.jobname = paste0(ann, '.dp.mtx')
		bigwig.file = paste0(chip.dsn$bigwig.file[c(1,3,4,5)], collapse=' ')
		dp.mtx.cmd = bsub.head(jobname = dp.mtx.jobname, cpu=5, mem=20, hmem=20)
		dp.mtx.cmd = paste0(dp.mtx.cmd, '"', ' computeMatrix scale-regions -S ', bigwig.file)
		dp.mtx.cmd = paste0(dp.mtx.cmd, ' -R ', mtx.in.fname, ' -b 500 -out ', mtx.file, '"')
		dp.mtx.cmd
		exe.jobs(dp.mtx.cmd, logger)
	}, by=1:nrow(dd.dsn)]

	dd.dsn[, dp.plot.h.jobname := paste0(annos, '.dp.heatmap.', kk)]
	dd.dsn[, dp.plot.cmd = bsub.head(jobname = dp.plot.h.jobname, cpu=5, mem=20, hmem=20, postdone=dp.mtx.jobname), by=1:nrow(dd.dsn)]
	dd.dsn[, dp.plot.cmd := paste0(dp.plot.cmd, '"', ' plotHeatmap -m ', mtx.file, ' --samplesLabel Ctrl1 Ctrl2 KO1 KO2')]
	dd.dsn[, dp.plot.cmd := paste0(dp.plot.cmd, ' --startLabel S --endLabel E -out ', dp.heatmap.file, ' --kmeans ', k)]
	dd.dsn[, dp.plot.cmd := paste0(dp.plot.cmd, ' --zMin 0 0 0 0 --zMax 6 6 6 6 --whatToShow \\"heatmap and colorbar\\"')]
	dd.dsn[, dp.plot.cmd := paste0(dp.plot.cmd, ' --heatmapWidth 3 --heatmapHeight 4')]
	dd.dsn[, dp.plot.cmd := paste0(dp.plot.cmd, ' --outFileSortedRegion ', dp.out.cluster.file, '"')]
	dd.dsn[, dp.plot.cmd]
	exe.jobs(dd.dsn[, dp.plot.cmd], logger)

	## promoter
	{

		annos = 'promoter'
		tmp = dbi.anno[anno %in% c(annos), ][symbol %in% row.names(res.cc.log2), ]
		table(sign(tmp$Fold))
		mtx.in.fname = paste0('penn/data/dbi_mtx_', annos, '.tsv')
		mtx.file = paste0('penn/data/dbi_mtx_', annos, '.gz')
		fwrite(tmp[, .(chr, Start, End, Strand, id)], file=mtx.in.fname, sep="\t", col.names=F)
		dp.mtx.jobname = paste0(annos, '.dp.mtx')
		dp.mtx.jobname
		bigwig.file = paste0(chip.dsn$bigwig.file[c(1,3,4,5)], collapse=' ')
		#dp.mtx.cmd = bsub.head(jobname = dp.mtx.jobname, cpu=5, mem=20, hmem=20)
		dp.mtx.cmd = paste0(' computeMatrix scale-regions -S ', bigwig.file)
		dp.mtx.cmd = paste0(dp.mtx.cmd, ' -R ', mtx.in.fname, ' -b 500 -a 500 -p 10 -out ', mtx.file)
		dp.mtx.cmd
		exe.jobs(dp.mtx.cmd, logger)

		kk = 2
		dp.heatmap.file = paste0('penn/cor/', annos, '_deep_tools_heatmap_kmean', kk, '2.pdf')
		dp.out.cluster.file = paste0('penn/data/', annos, '_out_cluster_mtx_kmean', kk, '.tsv')
		dp.plot.jobname = paste0(annos, '.dp.heatmap.', kk)
		#dp.plot.cmd = bsub.head(jobname = dp.plot.h.jobname, cpu=5, mem=20, hmem=20, postdone=dp.mtx.jobname)
		dp.plot.cmd = paste0(' plotHeatmap -m ', mtx.file, ' --samplesLabel Ctrl1 Ctrl2 KO1 KO2')
		dp.plot.cmd = paste0(dp.plot.cmd, ' --startLabel S --endLabel E -out ', dp.heatmap.file, ' --kmeans ', kk)
		dp.plot.cmd = paste0(dp.plot.cmd, ' --zMin 0 0 0 0 --zMax 6 6 6 6 --whatToShow "heatmap and colorbar"')
		dp.plot.cmd = paste0(dp.plot.cmd, ' --heatmapWidth 3 --heatmapHeight 5')
		dp.plot.cmd = paste0(dp.plot.cmd, ' --outFileSortedRegion ', dp.out.cluster.file)
		dp.plot.cmd
		exe.jobs(dp.plot.cmd, ifwait=T, logger)

		# gene expression of the genes next to these peaks for promoter

		kk = 2
		annos = 'promoter'
		dp.out.cluster.file = paste0('penn/data/', annos, '_out_cluster_mtx_kmean', kk, '.tsv')
		dp.out.cluster.file
		tmp = fread(dp.out.cluster.file)
		setnames(tmp, 1, 'chr')
		tmp[, id := paste0('chr', chr, '_', start-1, '_', end)]
		setkey(dbi.anno, 'id')
		tmp2 = dbi.anno[tmp$id, ]
		tmp2[, Fold := -Fold]
		hp.title = paste0(paste0(as.vector(table(sign(tmp2$Fold))), c(' up', ' dn'), sep=''), collapse=' and '); hp.title


		mtx = res.cc.log2[tmp2$symbol, c(1,2,5,6)]
		cn = colnames(mtx)
		mtx = t(apply(mtx, 1, scale))
		colnames(mtx) = cn
		mtx[mtx > ll] = ll
		mtx[mtx < -ll] = ll 

		ht.1 = Heatmap(mtx, col =col.fun(seq(-2, 2, length.out=20)),  show_column_dend=F, show_row_dend=F, show_column_names=T, show_heatmap_legend = T,  cluster_columns = T, cluster_rows=F, show_row_names=F)
		fname = paste0("penn/cor/res_H3K27_peak_gene_", annos, "_heatmap.pdf"); fname
		cat(fname, "\n")
		pdf(file=fname, width=3, height=3)
		draw(ht.1)
		dev.off()

	}

	## combine intronic and intergenic
	{

		tmp = dbi.anno[anno %in% c('intron', 'intergenic'), ][`Distance to TSS` < 10000, ][symbol %in% row.names(res.cc.log2), ]
		tmp
		table(sign(tmp$Fold))
		annos = 'introgenic'
		mtx.in.fname = paste0('penn/data/dbi_mtx_', annos, '.tsv')
		mtx.file = paste0('penn/data/dbi_mtx_', annos, '.gz')
		fwrite(tmp[, .(chr, Start, End, Strand, id)], file=mtx.in.fname, sep="\t", col.names=F)
		dp.mtx.jobname = paste0(annos, '.dp.mtx')
		bigwig.file = paste0(chip.dsn$bigwig.file[c(1,3,4,5)], collapse=' ')
		#dp.mtx.cmd = bsub.head(jobname = dp.mtx.jobname, cpu=5, mem=20, hmem=20)
		dp.mtx.cmd = paste0(' computeMatrix scale-regions -S ', bigwig.file)
		dp.mtx.cmd = paste0(dp.mtx.cmd, ' -R ', mtx.in.fname, ' -b 500 -a 500 -p 10 -out ', mtx.file)
		dp.mtx.cmd

		exe.jobs(dp.mtx.cmd, logger)

		for(kk in 2:6){
			dp.heatmap.file = paste0('penn/cor/', annos, '_deep_tools_heatmap_kmean', kk, '.pdf')
			dp.out.cluster.file = paste0('penn/data/', annos, '_out_cluster_mtx_kmean', kk, '.tsv')
			dp.plot.jobname = paste0(annos, '.dp.heatmap.', kk)
			#dp.plot.cmd = bsub.head(jobname = dp.plot.h.jobname, cpu=5, mem=20, hmem=20, postdone=dp.mtx.jobname)
			dp.plot.cmd = paste0(' plotHeatmap -m ', mtx.file, ' --samplesLabel Ctrl1 Ctrl2 KO1 KO2')
			dp.plot.cmd = paste0(dp.plot.cmd, ' --startLabel S --endLabel E -out ', dp.heatmap.file, ' --kmeans ', kk)
			dp.plot.cmd = paste0(dp.plot.cmd, ' --zMin 0 0 0 0 --zMax 6 6 6 6 --whatToShow "heatmap and colorbar"')
			dp.plot.cmd = paste0(dp.plot.cmd, ' --heatmapWidth 3 --heatmapHeight 5')
			dp.plot.cmd = paste0(dp.plot.cmd, ' --outFileSortedRegion ', dp.out.cluster.file)
			dp.plot.cmd
			exe.jobs(dp.plot.cmd, ifwait=F, logger)
		}

		# gene expression of the genes next to these peaks
		{

			kk = 2
			annos = 'introgenic'
			dp.out.cluster.file = paste0('penn/data/', annos, '_out_cluster_mtx_kmean', kk, '.tsv')
			dp.out.cluster.file
			tmp = fread(dp.out.cluster.file)
			setnames(tmp, 1, 'chr')
			dbi.anno[Start == 45987503, ]
			tmp[, id := paste0('chr', chr, '_', start-1, '_', end)]
			setkey(dbi.anno, 'id')
			tmp2 = dbi.anno[tmp$id, ]
			tmp2[, Fold := -Fold]
			hp.title = paste0(paste0(as.vector(table(sign(tmp2$Fold))), c(' up', ' dn'), sep=''), collapse=' and '); hp.title

			head(tmp2)
			mtx = res.cc.log2[tmp2$symbol, ]
			#mtx = mtx[!duplicated(row.names(mtx)), ]
			cn = colnames(mtx)
			mtx = t(apply(mtx, 1, scale))
			colnames(mtx) = cn
			mtx[mtx > ll] = ll
			mtx[mtx < -ll] = ll 

			mtx.rr = apply(mtx, 1, function(x){any(is.na(x))})
			mtx = mtx[!mtx.rr, ]
			rn = unique(row.names(mtx)); 
			an.row = tmp[symbol %in% rn, .(symbol, Fold, deepTools_group)]
			an.row = an.row[!is.na(Fold), ]
			an.row[, Fold.m := median(Fold), by='symbol']
			an.row = an.row[!duplicated(symbol), ]
			an.row = data.frame(row.names = an.row$symbol, Fold = an.row$Fold.m, deepcluster = an.row$deepTools_group) 
			an.row$Fold[an.row$Fold > 3] = 3
			an.row$Fold[an.row$Fold < -3] = -3
			an.row$updn[an.row$Fold < 0] = 'Par' 
			an.row$updn[an.row$Fold > 0] = 'KO' 
			head(an.row)

			col.fun.dp = colorRamp2(c(-3,0,3), colors=c('firebrick', 'white', 'midnightblue'))
			col.clu = c(cluster_1 = adjustcolor('midnightblue', 1), cluster_3 = adjustcolor('dodgerblue', 1), cluster_2 = adjustcolor('firebrick', 1), cluster_4 = adjustcolor('maroon1', 1))
			mtx = mtx[row.names(an.row), ]
			ha.r = rowAnnotation(df = an.row, col=list(updn=c('KO'='midnightblue', 'Par'='firebrick'), Fold=col.fun.dp, deepcluster=col.clu))

			ht.1 = Heatmap(mtx, col =col.fun(seq(-2, 2, length.out=20)),  show_column_dend=F, show_row_dend=F, show_column_names=T, show_heatmap_legend = T,  cluster_columns = T, cluster_rows=F, show_row_names=F)
			fname = paste0("penn/cor/res_H3K27_peak_gene_", annos, "_heatmap.pdf"); fname
			cat(fname, "\n")
			pdf(file=fname, width=3, height=3)
			draw(ht.1)
			dev.off()

			dp.out.cluster.file = paste0('penn/data/', annos, '_out_cluster_mtx_kmean', kk, '.tsv')
			dp.out.cluster.file
			tmp = fread(dp.out.cluster.file)
			setnames(tmp, 1, 'chr')
			tmp[, id := paste0('chr', chr, '_', start-1, '_', end)]
			tmp = merge(tmp, dbi.anno, all.x=T, by.x='id', by.y='id')
			tmp[, Fold := -Fold]
			tmp
			tmp2 = merge(tmp, res, by = 'symbol', all.x = T)
			head(tmp2)
			gg = ggplot(tmp2, aes(Fold, log2FoldChange, color=deepTools_group)) + geom_point(aes(shape=deepTools_group)) + xlab('Fold Change of H3K27ac modification (KO vs Parental)') +
				ylab('log2FC RNA-Seq (KO vs Parental)') + ggtitle('FOXA1 deletion w/ intergenic/intron H3K27ac modifcation and gene expression')
			fname = paste0("penn/cor/res_H3K27_peak_gene_introgenic_heatmap_cor.pdf"); fname
			ggsave(gg, file=fname, width=5, height=4)

		}

	}

	dd.dsn[, { scp(dp.heatmap.file) }, by=1:nrow(dd.dsn)]
	dd.dsn
	dir.create('penn/heatmap')
	dd.dsn[, {system(paste0('cp ', dp.heatmap.file, ' penn/heatmap'))}, by=1:nrow(dd.dsn)]

	dd.dsn
	aa = fread(dd.dsn$dp.out.cluster.file[7])
	table(aa$deepTools_group)

	## annotate peaks
	## deeptools clustering peaks
	## ** find motifs in the top up/dn clusters
	## RNA-Seq to the top up/dn clusters

	## motif for deeptool clusters
	dd.dsn[, tag := paste0(annos, '_', kk)]
	annos = c('intron', 'intergenic', 'promoter', 'genebody', 'all')
	annos.k = c(4, 4, 3, 5, 4)
	tag.sel = paste0(annos, '_', annos.k)
	tag.sel

	xx = dd.dsn[tag %in% tag.sel, ]
	xx[, up := c(1,1,1,1,1)]
	xx[, dn := c(2,2,3,2,2)]
	xx = xx[rep(1:nrow(xx), each=4), ]
	xx[, updn := rep(c('up', 'dn', 'updn', 'all'), times=5)]
	xx

	pp.dsn = xx
	pp.dsn[, .(annos, kk, up, dn)]

	pp.dsn[, findmotif.dir := paste0('penn/motif_', annos, '_kmean_', updn)]
	pp.dsn[, findmotif.in.fname := paste0(findmotif.dir, '/findmotif_input_mtx.tsv')]
	pp.dsn[, dp.out.cluster.full.fname := paste0(findmotif.dir, '/dp_out_cluster_full.tsv')]
	pp.dsn[, findmotif.jobname := paste0("findmotif.", annos, '.', updn)]

	pp.dsn[, {
		tmp.ann = fread(dp.out.cluster.file)
		setnames(tmp.ann, 1, 'chr')
		tmp.ann[, chr := paste0('chr', chr)]
		tmp.ann[, tag := paste0(chr, '_', start - 1, '_', end)] 
		tmp.ann = merge(tmp.ann[, .(tag, deepTools_group)], dbi.anno, by='tag', all.x=T)
		tmp.ann
		if(updn == 'up'){
			tmp = tmp.ann[deepTools_group %in% paste0('cluster_', up), ]
		}else if(updn == 'dn'){
			tmp = tmp.ann[deepTools_group %in% paste0('cluster_', dn), ]
		}else if(updn == 'updn'){
			tmp = tmp.ann[deepTools_group %in% paste0('cluster_', c(up, dn)), ]
		}else if(updn == 'all'){
			tmp = tmp.ann
		}
		fwrite(tmp[, .(chr, start, end)], file=findmotif.in.fname, sep="\t", col.names=F)
		fwrite(tmp, file=dp.out.cluster.full.fname, sep="\t", col.names=T)
	}, by=1:nrow(pp.dsn)]

	file.exists(pp.dsn[, findmotif.in.fname])

	pp.dsn[, findmotif.cmd.2 := bsub.head(jobname = findmotif.jobname, mem=20, hmem=20, cpu=6), by=1:nrow(pp.dsn)]
	pp.dsn[, findmotif.cmd.2 := paste0(findmotif.cmd.2, ' "findMotifsGenome.pl ', findmotif.in.fname, ' hg19 ', findmotif.dir, ' -cpu 4 -mask "')]
	pp.dsn[, findmotif.cmd.2]
	exe.jobs(pp.dsn[, findmotif.cmd.2], logger)

	pp.dsn[, {tmp = fread(findmotif.in.fname);
	       tmp = nrow(tmp)
	       cat(basename(findmotif.dir), '\t', tmp, '\n')
	}, by=1:nrow(pp.dsn)]

	pp.dsn[, motif.fname := paste0(findmotif.dir, '/knownResults.txt')]
	file.exists(pp.dsn$motif.fname)
	file.size(pp.dsn$motif.fname)
	pp.dsn$sample.name

	mm = lapply(pp.dsn$motif.fname, fread)
	names(mm) = pp.dsn$findmotif.dir
	mm = rbindlist(mm, idcol=T, use.names=F)
	setnames(mm, c('.id', 'motif.name', 'sequence', 'pvalue', 'logp', 'qvalue', 'target.hits', 'percentage.hits', 'background.hits', 'percentage.background.hits'))
	mm[, logp2 := -log(pvalue)]
	mm[, motif.name2 := sub('\\(.*', '', motif.name)]
	mm[, pp := as.numeric(sub('%', '', percentage.hits))]
	mm[, sig := 'no']
	mm[ qvalue < 0.05, sig := 'yes']
	mm = mm[order(qvalue), ]
	mm[, .id := basename(.id)]

	mm.l = mm %>% group_by(.id) %>% top_n(n=10, wt=qvalue)
	dim(mm.l)
	dim(mm)
	mm[grep('TEAD', motif.name), ]
	mm[grep('STAT', motif.name), ][grep('promoter', .id), ]
	mm[grep('STAT', motif.name), ][sig == 'yes', ]

	mm.1 = mm[grep('all_kmean_all', .id), ]
	mm.1[, yanse := 'others']
	mm.1[!is.na(ll), yanse := 'Sig']
	mm.1 = mm.1[order(qvalue), ]
	mm.1[1:10, ll := motif.name2]
	mm.1[grepl('IRF|ISRE', motif.name) & qvalue < 0.05, ll := motif.name2]
	mm.1[grep('TEAD', motif.name), ]

	gg = ggplot(mm.1[qvalue < 0.05, ], aes(x=pp, y=logp2, label=ll, color=yanse)) + geom_point(size=2)  + 
		geom_text_repel(size=2)  + xlab('H3K27ac peaks (%) with motif binding sites') + ylab('-log(p value)')  +
		theme(legend.position = 'none')
	fname = paste0('penn/dotplot_', basename(findmotif.dir), '_dotplot_v3.pdf'); fname
	ggsave(gg, file=fname, width=4, height=3)
	scp(fname)

	pp.dsn[, motif.sum.fname := paste0(findmotif.dir, 'mtx_sum.tsv')]
	file.exists(pp.dsn$motif.fname)
	file.size(pp.dsn$motif.fname)
	pp.dsn[ , {
		mtx = fread(motif.fname)
		setnames(mtx, c('motif.name', 'sequence', 'pvalue', 'logp', 'qvalue', 'target.hits', 'percentage.hits', 'background.hits', 'percentage.background.hits'))
		mtx[, logp2 := -log(pvalue)]
		mtx[, motif.name2 := sub('\\(.*', '', motif.name)]
		mtx[, pp := as.numeric(sub('%', '', percentage.hits))]
		mtx[, sig := 'no']
		mtx[ qvalue < 0.05, sig := 'yes']
		mtx = mtx[order(qvalue), ]
		mtx[1:10, ll := motif.name2]
		mtx[grepl('IRF|ISRE', motif.name) & qvalue < 0.05, ll := motif.name2]
		mtx[, yanse := 'others']
		mtx[!is.na(ll), yanse := 'Sig']
		mtx
		fwrite(mtx, file=motif.sum.fname)
		gg = ggplot(mtx[qvalue < 0.05, ], aes(x=pp, y=logp2, label=ll, color=yanse)) + geom_point(size=2)  + 
			geom_text_repel(size=2)  + xlab('H3K27ac peaks (%) with motif binding sites') + ylab('-log(p value)')  +
			theme(legend.position = 'none')
		fname = paste0('penn/dotplot_', basename(findmotif.dir), '_dotplot_v2.pdf'); fname
		ggsave(gg, file=fname, width=4, height=3)
		scp(fname)
	}, by=1:nrow(pp.dsn)]

	dbi.anno[, symbol := `Gene Name`]
	setkey(dbi.anno, symbol)
	dbi.anno

	pp.dsn[, .(annos, kk, updn)]
	dp.out.cluster.full.fname = pp.dsn$dp.out.cluster.full.fname[8]
	dp.out.cluster.full.fname = NULL

	pp.sel = pp.dsn[updn %in% c('all'), ]
	pp.sel[, .(annos, kk, updn)]
	for(ii in 1:nrow(pp.sel)){
		ii.dp.fname = pp.sel[ii, dp.out.cluster.full.fname]
		ii.annos = pp.sel[ii, annos]
		ii.kk = pp.sel[ii, kk]
		ii.updn = pp.sel[ii, updn]
		tmp = fread(ii.dp.fname)
		hp.title = paste0(paste0(as.vector(table(sign(tmp$Fold))), c(' up', ' dn'), sep=''), collapse=' and '); hp.title
		mtx = res.cc.log2[row.names(res.cc.log2) %in% tmp$`Gene Name`, c(1,2,5,6)]
		mtx = mtx[!duplicated(row.names(mtx)), ]
		cn = colnames(mtx)
		mtx = t(apply(mtx, 1, scale))
		colnames(mtx) = cn
		mtx[mtx > ll] = ll
		mtx[mtx < -ll] = ll 
		mtx.rr = apply(mtx, 1, function(x){any(is.na(x))})
		mtx = mtx[!mtx.rr, ]
		rn = unique(row.names(mtx)); 
		an.row = tmp[symbol %in% rn, .(symbol, Fold, deepTools_group)]
		an.row = an.row[!is.na(Fold), ]
		an.row[, Fold.m := median(Fold), by='symbol']
		an.row = an.row[!duplicated(symbol), ]
		an.row = data.frame(row.names = an.row$symbol, Fold = an.row$Fold.m, deepcluster = an.row$deepTools_group) 
		an.row$Fold[an.row$Fold > 3] = 3
		an.row$Fold[an.row$Fold < -3] = -3
		an.row$updn[an.row$Fold < 0] = 'Par' 
		an.row$updn[an.row$Fold > 0] = 'KO' 
		col.fun.dp = colorRamp2(c(-3,0,3), colors=c('red', 'white', 'blue'))
		col.fun.dp = colorRamp2(c(-3,0,3), colors=c('firebrick', 'white', 'midnightblue'))
		if(ii.kk == 3){ col.clu = c(cluster_1 = adjustcolor('midnightblue', .7), 
					 cluster_2 = adjustcolor('midnightblue', 0.5), 
					 cluster_3 = adjustcolor('firebrick', .4))}
		if(ii.kk == 4){ col.clu = c(cluster_1 = adjustcolor('midnightblue', .7), 
					 cluster_3 = adjustcolor('midnightblue', 0.5), 
					 cluster_2 = adjustcolor('firebrick', .7),
					 cluster_4 = adjustcolor('firebrick', .4))}
		if(ii.kk == 5){ col.clu = c(cluster_1 = adjustcolor('midnightblue', .7), 
					 cluster_3 = adjustcolor('midnightblue', 0.5), 
					 cluster_5 = adjustcolor('midnightblue', 0.2), 
					 cluster_2 = adjustcolor('firebrick', .7),
					 cluster_4 = adjustcolor('firebrick', .4))}
		mtx = mtx[row.names(an.row), ]
		ha.r = rowAnnotation(df = an.row, width=1, col=list(updn=c('KO'='midnightblue', 'Par'='firebrick'), Fold=col.fun.dp, deepcluster=col.clu))
		ht.1 = Heatmap(mtx, col =col.fun(seq(-2, 2, length.out=20)),  show_column_dend=F, show_row_dend=F, 
			       show_column_names=T, show_heatmap_legend = T,  cluster_columns = T, cluster_rows=T,
			       show_row_names=F, left_annotation=ha.r)
		fname = paste0("penn/cor/res_H3K27_peak_gene_", ii.annos, '_', ii.updn, "_heatmap_v8.pdf"); fname
		cat(fname, "\n")
		pdf(file=fname, width=3, height=3)
		draw(ht.1)
		dev.off()
	}

	pp.dsn[, .(annos, kk, updn)]
	pp.sel = pp.dsn[updn %in% c('all'), ]
	pp.sel[, {
		tmp = fread(dp.out.cluster.full.fname)
		tmp
		tmp2 = merge(tmp, res, by = 'symbol')
		tmp2
		gg = ggplot(tmp2, aes(Fold, log2FoldChange, color=deepTools_group)) + geom_point(aes(shape=deepTools_group))
		fname = paste0("penn/cor/res_H3K27_peak_gene_", annos, '_', updn, "_heatmap_cor.pdf"); fname
		ggsave(gg, file=fname, width=5.8, height=4)
	}, by=1:nrow(pp.sel)]
	
	tmp = peak.all[V1 == '9' &  V2 > 5440000 &  V3< 5460000, ]
	colnames(tmp) = c('sample', 'chr', 'start', 'end', 'id', 'score', 'strand', 'signalvalue', 'pvalue', 'qvalue', 'peakCenter')
	tmp

	pp.dsn[ , {
		tmp.cor = merge(tmp, res, by.x = 'Nearest Ensembl', by.y='rn', all=T)
		tmp.cor = tmp.cor[!is.na(Chr), ]
		tmp.cor = tmp.cor[!is.na(symbol), ]
		#tmp.cor = tmp.cor[abs(log2FoldChange) > 1, ]
		tmp.cor[, exp.log2 := -log(pvalue)]
		tmp.cor
		gg = ggplot(tmp.cor, aes(x=log2FoldChange, y = Fold)) + geom_point()
		#fname = paste0("penn/res_H3K27_peak_gene_", k.ann, "_heatmap_v6.pdf");
		fname = sub('.pdf',  '_cor.pdf', fname)
		ggsave(gg, file=fname)
		scp(fname)
	}, by=1:nrow(pp.dsn)]

	pp.lst = lapply(pp.dsn$motif.sum.fname, fread)
	names(pp.lst) = pp.dsn$findmotif.dir
	pp.dt = rbindlist(pp.lst, idcol=T)
	pp.dt

	gg = ggplot(pp.dt, aes(x=pp, y=logp2)) + geom_point(size=2)  + 
		geom_text_repel(size=2, aes(label=ll)) +
		facet_wrap(~ pp.dt$.id, ncol=4, scale='free')
	fname = paste0('penn/motisum_motifsum.pdf')
	ggsave(gg, file=fname, width=24, height=23)
	scp(fname)

	col.fun = colorRamp2(c(-2, 0, 2), colors=c('red', 'white', 'blue'))
	cols = col.fun(seq(-2, 2, length.out=31))
	cols

	ht.1 = Heatmap(mtx, col =col.fun(seq(-2, 2, length.out=20)),  show_column_dend=F, show_row_dend=T, 
		       show_column_names=F, show_heatmap_legend = F,  cluster_columns = T, cluster_rows=T,
		       show_row_names=F)
	fname = 'penn/res_heatmap.pdf'
	pdf(file=fname, width=7, height=8)
	draw(ht.1)
	dev.off()
	scp(fname)

	fname = 'penn/rnaseq_coverage_cd274.pdf'
	cmd = paste0('plotCoverage -b ', paste(sort(ko.dsn$bamfile), collapse=' '), ' -r chr9:5439265:5457265i -p 10 -o ', fname)
	system(cmd)
	scp(fname)

}


	### RNA Seq of FOXA1 knockout
{

	ko.dsn = data.table(r1 = system('find penn/pennfastq -name "*.R1.fastq.gz"', intern=T))
	ko.dsn[, r2 := sub('R1', 'R2', r1)]
	ko.dsn[, sample.name := sub('.*G432-UMUC1-(.*)_merged.*', '\\1', r1)]
	ko.dsn[, sample.dir := paste0('penn/', cfg$initDir, '/', sample.name)]
	ko.dsn[, {dir.create(sample.dir)}, by=1:nrow(ko.dsn)]
	ko.dsn

	ko.dsn[, rsem.jobname := paste0("rsem.", sample.name)]
	ko.dsn[, rsem.cmd := bsub.head(jobname = rsem.jobname, cpu=8, mem=70, hmem=7), by=1:nrow(ko.dsn)]
	ko.dsn[, rsem.cmd := paste0(rsem.cmd, ' "', cfg$RSEM, '/rsem-calculate-expression -p 8 --paired-end --star --star-path ', cfg$STAR_DIR)]
	ko.dsn[, rsem.cmd := paste0(rsem.cmd, ' --gzipped-read-file --append-names --estimate-rspd ', r1, ' ', r2, ' ', cfg$RSEM_REF, ' ')]
	ko.dsn[, rsem.cmd := paste0(rsem.cmd, sample.dir, '/', sample.name, '"')]
	ko.dsn[, rsem.cmd]

	exe.jobs(ko.dsn[, rsem.cmd], logger)


	ko.dsn[, bamfile := paste0(sample.dir, '/', sample.name, '.transcript.sorted.bam')]

	rsem = fread('penn/G432.txt')
	rsem = setDF(rsem[, 2:7], rownames=rsem$name)
	head(rsem)

	condition = c('wt', 'wt', 'ko1', 'ko1', 'ko2', 'ko2')
	condition = c('wt', 'wt', 'ko', 'ko', 'ko', 'ko')
	dsn = data.frame(
			 row.names       = colnames(rsem),
			 condition       = condition,
			 patientID 	 = colnames(rsem),
			 libType         = rep("PE", ncol(rsem)));
	dsn$condition = factor(dsn$condition, levels=c('wt', 'ko'))
	dsn

	ddsmat = DESeqDataSetFromMatrix(countData = rsem,
					colData =dsn,
					design = ~ condition);
	dds.ds <- estimateSizeFactors(ddsmat);
	dds <- DESeq(dds.ds, parallel=F);
	fname = 'dds_FOXA1KO.RData'
	save(dds, file=fname)

	resultsNames(dds)
	res.ko1 = results(dds, name = 'condition_ko1_vs_wt')
	res.ko2 = results(dds, name = 'condition_ko2_vs_wt')
	res.ko1 = as.data.table(res.ko1, keep.rownames=T)
	res.ko2 = as.data.table(res.ko2, keep.rownames=T)
	res.12 = rbindlist(list(res.ko1 = res.ko1, res.ko2 = res.ko2), idcol=T)
	res.12 = res.12[order(padj), ]
	res.12 = merge(res.12, id.dt, by.x = 'rn', by.y='ensg', all.x=T)
	table(res.12[padj < 0.05, paste0(.id, sign(log2FoldChange))])
	intersect(res.12[padj < 0.05 & .id == 'res.ko1', rn], res.12[padj < 0.05 & .id == 'res.ko1', rn])
	setdiff(res.12[padj < 0.05 & .id == 'res.ko1', rn], res.12[padj < 0.05 & .id == 'res.ko2', rn])
	setdiff(res.12[padj < 0.05 & .id == 'res.ko2', rn], res.12[padj < 0.05 & .id == 'res.ko1', rn])

	res = results(dds)
	res = as.data.table(as.data.frame(res), keep.rownames=T)
	res = res[order(padj),][!is.na(padj), ]
	res = merge(res, id.dt, by.x = 'rn', by.y='ensg', all.x=T)
	res = res[order(padj), ]
	head(res)
	fname = 'res_FOXA1KO.RData'
	load(fname)
	res
	save(res, file=fname)

	res[grep('IRF', symbol), ]

	table(res[padj < 0.05 & abs(log2FoldChange) > 1, sign(log2FoldChange)])

	setkey(id.dt, 'ensg')
	id.dt[row.names(res.cc), symbol]
	res.cc = counts(dds, normalized = T)
	res.cc = res.cc[row.names(res.cc) %in% res$rn, ]
	row.names(res.cc) = id.dt[row.names(res.cc), symbol]
	head(res.cc)
	fname = 'res_FOXA1KO_cc.RData'
	save(res.cc, file=fname)
	head(res.cc)
	load(fname)

	id.dt = data.table(id = row.names(res.tp.cc))
	id.dt[, symbol := sub('.*_', '', id)]
	id.dt[, ensg := sub('_.*', '', id)]
	id.dt
	fname = '../id_dt.RData'
	save(id.dt, file=fname)

	source("~/program/fun/write_rnk.r")
	write_rnk(res, file='penn/res.rnk', neg=T)
	write.csv(res, file='penn/res.csv')

	source("~/program/fun/run_gsea.R")
	run_gsea('penn/res.rnk', ifsvg=T)
	sync('penn/gsea')

	## sum_gsea.r

	res.cc.log2 = log2(res.cc + 1)
	head(res.cc.log2)

	## 
	rn.sig = res[padj < 0.01 & !is.na(symbol),]
	hp.title = paste0(paste0(as.vector(table(sign(rn.sig$log2FoldChange))), c(' up', ' dn'), sep=''), collapse=' and ')
	hp.title
	mtx = res.cc.log2[rn.sig$symbol, ]
	cn = colnames(mtx)
	mtx = t(apply(mtx, 1, scale))
	colnames(mtx) = cn
	#mtx[mtx > ll] = ll
	#mtx[mtx < -ll] = -ll 
	col.fun = colorRamp2(c(-2, 0, 2), colors=c('red', 'white', 'blue'))
	ht.1 = Heatmap(mtx, col =col.fun(seq(-2, 2, length.out=20)),  show_column_dend=F, show_row_dend=T, 
		       show_column_names=F, show_heatmap_legend = F,  cluster_columns = T, cluster_rows=T,
		       show_row_names=F)
	fname = 'penn/res_heatmap.pdf'
	pdf(file=fname, width=7, height=8)
	draw(ht.1)
	dev.off()
	scp(fname)

	ss.sel = c('ISG15', 'IFIT2', 'IFIT3', 'CD274', 'IFI44L', 'STAT2', 'IFI35', 'IRF9', 'PDCD1LG2')
	mtx = res.cc.log2[row.names(res.cc.log2) %in% ss.sel, ]
	cn = colnames(mtx)
	mtx = t(apply(mtx, 1, scale))
	colnames(mtx) = cn
	mtx 

	col.fun = colorRamp2(c(-2, 0, 2), colors=c('red', 'white', 'blue'))
	ht.1 = Heatmap(mtx, col =col.fun(seq(-2, 2, length.out=20)),  show_column_dend=F, show_row_dend=F, 
		       show_column_names=T, show_heatmap_legend = T,  cluster_columns = T, cluster_rows=T,
		       show_row_names=T)
	fname = paste0("penn/res_heatmap_ss_sel.pdf")
	pdf(fname, width=3, height=3)
	draw(ht.1)
	dev.off()
	scp(fname)

	## res.david maplot
	source('~/pipeline/fun/gg_maplot.r')

	table(res.david[abs(log2FoldChange) > 1 & padj < 0.05, sign(log2FoldChange)])
	res.david[abs(log2FoldChange) > 1 & padj < 0.05, ][grep('DSC', symbol), ]

	fname = 'penn/res_maplot.png'
	gg = gg_maplot(res.david[order(padj),], labx = 'log2 FC (KO vs WT)', filename='',  topn=0, wi=3, hi=2)
	gg = gg + xlim(-10, 15) + ylim(0, 35)
	ggsave(gg, file=fname, width=2.5, height=2, unit='in', dpi=300)
	scp(fname)

	## pca analysis
	rld = rlog(dds, blind = F)

	pcafile = paste0('penn/pca.pdf')
	pdf(pcafile)
	rd = plotPCA(rld, intgroup = 'condition', returnData = T)
	dev.off()
	scp(pcafile)

	rd
	pcafile = paste0('penn/res_ggplot.pdf')
	pdf(pcafile, width=8, height=6)
	ggplot(rd, aes(PC1, PC2, label = group,  color=group)) + geom_point() + theme(legend.position="right") 
	dev.off()

	pcafile = paste0('penn/pca_ggplot2.pdf')
	percentvar = round(100 * attr(rd, 'percentVar'))
	pdf(pcafile)
	ggplot(rd, aes(x = PC1, y = PC2, color = condition, shape = condition)) +
		geom_point(size=3) +
		xlab(paste0("PC1: ", percentvar[1], "% variance")) +
		ylab(paste0("PC2: ", percentvar[2], "% variance")) +
		coord_fixed()
	dev.off()
	scp(pcafile)

	tmp = melt(as.data.table(res.cc.log2[row.names(res.cc.log2) %in% c('CD274', 'FOXA1', 'IRF9', 'IRF3', 'IRF5','GATA3', 'PPARG', 'INFGR1', 'INFGR2', 'INFG'), ], keep.rownames=T))
	tmp[, gt := 'wt']
	tmp[grep('KO', variable), gt := 'ko']
	tmp[, gt := factor(gt, levels=c('wt', 'ko'))]
	tmp

	gg = ggplot(tmp, aes(x = gt, y = value)) + geom_jitter(width=.2) + facet_wrap(~rn, nrow=2, scale='free_y') +
		stat_summary(fun.data = "mean_cl_boot", geom='crossbar', width=.5, color='red')  +
		theme(axis.text.x=element_text(angle=90, vjust=.5, hjust=1)) + xlab('') + ylab('log2 # reads') 
	fname = 'penn/dotplot_genes_v2.pdf'
	ggsave(gg, file=fname, width=6, height=4)
	scp(fname)

	cmd = paste('bedtools coverage -a cd274.bed -d -o a -b', paste(ko.dsn$bamfile, collapse=','), sep=' ')
	cmd
	write(cmd, file='run.sh')
	system(cmd)

	res[symbol == 'ARID1A', ]
	table(scc.maf$bcr)

}
