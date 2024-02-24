bismark = function(tgt, cfg=cfg){

	if('fastq' %in%  colnames(tgt)){
		tgt[grep('R2', fastq), R12 := 'R2']
		tgt[grep('R1', fastq), R12 := 'R1']
		tgt[, samplename := sub(".*\\/(.*)_IGO.*", "\\1", fastq)]
		tgt[grep("^\\d+", samplename), samplename := paste0("s_", samplename)]
		tgt[, samplename := gsub('-', '_', samplename)]
		tgt = dcast(tgt, samplename ~ R12, value.var = 'fastq')
	}

	if(!('sampledir' %in% colnames(tgt))){
		tgt[, sampledir := paste0(cfg$initDir, '/', samplename)]
		if(any(dir.exists(tgt$sampledir))){stop(paste0(tgt$sampledir, 'sample dir already there, stopped!'))}
		tgt[, {system(paste('mkdir -p ', sampledir))}, by=1:nrow(tgt)]
	}

	ifelse('bamfile' %in% colnames(tgt), 1, 		{ tgt[, bamfile := paste0(sampledir, '/', samplename, '_pe.bam')] })
	ifelse('name.sorted.bam.file' %in% colnames(tgt), 1, 	{ tgt[, name.sorted.bam.file := paste0(sampledir, '/', samplename, '_pe_name_sorted.bam')] })
	ifelse('sorted.bam.file' %in% colnames(tgt), 1, 	{ tgt[, sorted.bam.file := paste0(sampledir, '/', samplename, '_pe_sorted.bam')] })
	ifelse('dedup.bam.file' %in% colnames(tgt), 1, 		{ tgt[, dedup.bam.file := paste0(sampledir, '/', samplename, '_pe_sorted.deduplicated.bam')] })
	ifelse('sorted.dedup.bam.file' %in% colnames(tgt), 1, 	{ tgt[, sorted.dedup.bam.file := paste0(sampledir, '/', samplename, '_pe_sorted_dedup_sorted.bam')] })
	ifelse('align.sum.file' %in% colnames(tgt), 1, 		{ tgt[, align.sum.file := paste0(sampledir, '/', samplename, '_align_sum.txt')] })
	ifelse('gc.bias.metrics.file' %in% colnames(tgt), 1, 	{ tgt[, gc.bias.metrics.file:= paste0(sampledir, '/', samplename, '_gc_bias_metrics.txt')] })
	ifelse('gc.bias.metrics.pdf.file' %in% colnames(tgt), 1,{ tgt[, gc.bias.metrics.pdf.file := paste0(sampledir, '/', samplename, '_gc_bias_metrics.pdf')] })
	ifelse('bismark.cov.file' %in% colnames(tgt), 1,{ tgt[, bismark.cov.file := paste0(sampledir, '/', samplename, '_pe.bismark.cov.gz')] })
	ifelse('gc.bias.metrics.sum.file' %in% colnames(tgt), 1,{ tgt[, gc.bias.metrics.sum.file := paste0(sampledir, '/', samplename, '_gc_bias_metrics_sum.txt')] })
	ifelse('hs.metrics.file' %in% colnames(tgt), 1,		{ tgt[, hs.metrics.file := paste0(sampledir, '/', samplename, '_hs_metrics.txt')] })
	ifelse('insert.size.metrics.file' %in% colnames(tgt), 1,		{ tgt[, insert.size.metrics.file := paste0(sampledir, '/', samplename, '_insert_size_metrics.txt')] })
	ifelse('insert.size.metrics.histgram.file' %in% colnames(tgt), 1,	{ tgt[, insert.size.metrics.histgram.file := paste0(sampledir, '/', samplename, '_insert_size_metrics_histgram.pdf')] })
	ifelse('targeted.pcr.metrics.file' %in% colnames(tgt), 1,		{ tgt[, targeted.pcr.metrics.file := paste0(sampledir, '/', samplename, '_targeted_pcr_metrics.txt')] })
	ifelse('targeted.pcr.metrics.histgram.file' %in% colnames(tgt), 1,	{ tgt[, targeted.pcr.metrics.histgram.file := paste0(sampledir, '/', samplename, '_targeted_pcr_metrics_histgram.txt')] })
	ifelse('depth.file' %in% colnames(tgt), 1,				{ tgt[, depth.file := paste0(sampledir, '/', samplename, '_bam.depth')] })
	ifelse('dedup.report.file' %in% colnames(tgt), 1,			{ tgt[, dedup.report.file := paste0(sampledir, '/', samplename, '_pe_dedup_report.txt')] })
	ifelse('coord.sorted.bam.file' %in% colnames(tgt), 1,			{ tgt[, coord.sorted.bam.file := paste0(sampledir, '/', samplename, '_pe_coord_sorted.bam')] })
	ifelse('bigwig.file' %in% colnames(tgt), 1,				{ tgt[, bigwig.file := paste0(sampledir, '/', samplename, '_pe_coord_sorted.bigwig')] })

	if(!(all(dir.exists(tgt$sampledir)))) {stop('not all sample dir exists!')}
	tgt$bismark.cov.file
	sum(file.size(tgt$methylkit.file))
	file.exists(tgt$sorted.dedup.bam.file)

	## make genome for bismark
#	bismark.db.jobname = 'bismark.db'
#	bismark.db.cmd = bsub.head(bismark.db.jobname, mem=20, cpu=10)
#	bismark.db.cmd = paste0(bismark.db.cmd, " \"", cfg$bismark, '/bismark_genome_preparation --path_to_bowtie ', cfg$bowtie2dir, ' --verbose /ifs/work/solitlab/huw/study/db/hsa/grch37_bismark \"')
#	bismark.db.cmd
#	exe.jobs(bismark.db.cmd, cfg$logger)

	## align
	tgt[, bismark.jobname := paste0(samplename, '.bismark')]
	tgt[, bismark.cmd := bsub.head(bismark.jobname, mem=30, hmem = 30, cpu=10, postdone = 'bismark.db', W = '290:00')]  
	tgt[, bismark.cmd := paste0(bismark.cmd, " \"", cfg$bismark, "/bismark --quiet  -B ", samplename, " -o ", sampledir)]
	#tgt[, bismark.cmd := paste0(bismark.cmd, " --samtools_path ", cfg$samtoolsdir, " -B ", samplename, " --nucleotide_coverage ")]
	tgt[, bismark.cmd := paste0(bismark.cmd, " --samtools_path ", cfg$conda3bin, " --path_to_bowtie2 ", cfg$conda3bin, " --nucleotide_coverage ")]
	tgt[, bismark.cmd := paste0(bismark.cmd, " --temp_dir . ")]
	tgt[, bismark.cmd := paste0(bismark.cmd, " ", cfg$bismarkGRCh37, " -1 ", R1, " -2 ", R2, " \"")]
	tgt[1, bismark.cmd]

	exe.jobs(tgt[!file.exists(tgt$bamfile), bismark.cmd], cfg$logger)
	file.exists(tgt$bamfile)

	## sort bam file by coord for IGV
	file.exists(tgt$bamfile)
	tgt[, bam.coord.sort.jobname := paste0(samplename, '.coord.sort.bam')]
	tgt[, bam.coord.sort.cmd := bsub.head(bam.coord.sort.jobname, mem=40, hmem = 40, cpu=10, postdone = bismark.jobname, We = '130:00'), by=1:nrow(tgt)]  
	tgt[, bam.coord.sort.cmd := paste0(bam.coord.sort.cmd, " \"", cfg$conda3bin, '/samtools sort -@ 10 -m 4G -o ', coord.sorted.bam.file, ' ', bamfile, "\"")]
	tgt[53, bam.coord.sort.cmd]

	exe.jobs(tgt$bam.coord.sort.cmd, cfg$logger)
	write(tgt$bam.coord.sort.cmd, file='run.sh')
	tgt[, fname := paste0(sampledir, '/', samplename, '_pe.deduplicated_sorted.bam')]

	sname.sel = c('s_5632', 's_639v', 'umuc3', 'ku1919', 'rt112', 'RT112_p3', 'RT4_p4', 'SCaBER_p2', 'SW1710_p3','SW780_p11')
	fnames = tgt[samplename %in% sname.sel, ][file.exists(coord.sorted.bam.file), coord.sorted.bam.file]
	fcmd = 'r_fang/initFiles/bamfile/cmd.txt'
	cat(NULL, file=fcmd, append=F)
	cat(paste0('scp selene:', fnames, ' /Volumes/SolitLab/huw/.solit/study/hiseq/Project_06000_DD/r_fang/initFiles/bamfile/'), file=fcmd, append=T, sep="\n")
	scp(fcmd)

	fnames = tgt[samplename %in% sname.sel, ][file.exists(fname),fname]
	fnames
	tgt[, idx := 1:nrow(tgt)]
	tgt[grep('RT4', samplename), ]
	scp('a.txt')

	tgt[, bam.index.jobname := paste0(samplename, '.index.bam')]
	tgt[, bam.index.cmd := bsub.head(bam.index.jobname, mem=10, hmem = 4, cpu=2, We = '130:00'), by=1:nrow(tgt)]  
	tgt[, bam.index.cmd := paste0(bam.index.cmd, " \"", cfg$JAVA, '/java -Xms256M -Xmx10g -Djava.io.tmpdir=', cfg$tmpdir, ' -jar ', cfg$PICARD, ' BuildBamIndex I=', coord.sorted.bam.file, '"')]
	tgt[63, bam.index.cmd]

	exe.jobs(tgt[, bam.index.cmd], cfg$logger)
	write(tgt$bam.index.cmd, file='index.sh')

	## bigwig file
	tgt[, bigwig.jobname := paste0(samplename, ".bigwig")]
	tgt[, bigwig.cmd := bsub.head(bigwig.jobname, mem=30, We ='30:00', postdone=bam.index.jobname, cpu=10), by=1:nrow(tgt)]  
	tgt[, bigwig.cmd := paste0(bigwig.cmd, " \"", cfg$conda3bin, "/bamCoverage -b ", coord.sorted.bam.file, " -o ", bigwig.file, ' --normalizeUsing CPM -p 10 "')]
	tgt[, bigwig.cmd]

	exe.jobs(tgt[, bigwig.cmd], cfg$logger)

	## bismark_methylation_extractor
	tgt[, bismark.ex.jobname := paste0(samplename, ".bismark.ex")]
	tgt[, bismark.ex.cmd := bsub.head(bismark.ex.jobname, mem=30, hmem = 30, We ='30:00', cpu=10), by=1:nrow(tgt)]  
	tgt[, bismark.ex.cmd := paste0(bismark.ex.cmd, " \"", cfg$bismark, "/bismark_methylation_extractor --samtools_path ", cfg$samtoolsdir, " -o ", sampledir, " --parallel 8 --cytosine_report ")]
	tgt[, bismark.ex.cmd := paste0(bismark.ex.cmd, "--genome_folder ", cfg$bismarkGRCh37, ' ')] 
	tgt[, bismark.ex.cmd := paste0(bismark.ex.cmd, name.sorted.bam.file, " \"")]
	tgt[, bismark.ex.cmd]

	write(tgt[1:62, bismark.ex.cmd], file='bismark.ex.cmd.sh')
	exe.jobs(tgt[1:62, bismark.ex.cmd], cfg$logger)

	## bismark2report
	tgt[, bismark.report.jobname := paste0(samplename, ".bismark.report")]
	tgt[, bismark.report.cmd := bsub.head(bismark.report.jobname, mem=2, We ='10:00', cpu=1, postdone = bismark.jobname), by=1:nrow(tgt)]  
	tgt[, bismark.report.cmd := paste0(bismark.report.cmd, " \"", cfg$bismark, "/bismark2report --dir ", sampledir, " --alignment_report ", dedup.report.file, "\"")]
	exe.jobs(tgt$bismark.report.cmd, cfg$logger)

	## skip
	{
		## bismark_deduplication is more recommended for whole genome seq
		## sort bam file by seq name for dedup
		## bam index files only works with coord sorted file!
		tgt[, bam.name.sort.jobname := paste0(samplename, '.name.sort.bam')]
		tgt[, bam.name.sort.cmd := bsub.head(bam.name.sort.jobname, mem=40, cpu=10, We = '130:00'), by=1:nrow(tgt)]  
		tgt[, bam.name.sort.cmd := paste0(bam.name.sort.cmd, " \"", cfg$JAVA, '/java -Xms256M -Xmx40g -Djava.io.tmpdir=', cfg$tmpdir, ' -jar ', cfg$PICARD, ' SortSam I=', coord.sorted.bam.file, " O=", name.sorted.bam.file, " CREATE_INDEX=true SORT_ORDER=queryname\"")]
		tgt[1, bam.name.sort.cmd]

		exe.jobs(tgt[63, bam.name.sort.cmd], cfg$logger)

		## sort bam file by seq name for dedup
		file.exists(tgt$bamfile)
		tgt[, bam.name.sort.jobname := paste0(samplename, '.name.sort.bam')]
		tgt[, bam.name.sort.cmd := bsub.head(bam.name.sort.jobname, mem=40, hmem = 4, cpu=10, postdone = bismark.jobname, We = '130:00'), by=1:nrow(tgt)]  
		tgt[, bam.name.sort.cmd := paste0(bam.name.sort.cmd, " \"", cfg$JAVA, '/java -Xms256M -Xmx40g -Djava.io.tmpdir=', cfg$tmpdir, ' -jar ', cfg$PICARD, ' SortSam I=', bamfile, " O=", name.sorted.bam.file, " CREATE_INDEX=true SORT_ORDER=queryname\"")]

		exe.jobs(tgt[, bam.name.sort.cmd], cfg$logger)
		file.exists(tgt$sorted.bam.file)

		## deduplicate need to use name sorted bam file
		tgt[, dedup.jobname := paste0(samplename, '.dedup')]
		tgt[, dedup.cmd := bsub.head(dedup.jobname, mem=10, hmem=5, cpu=2, postdone = bismark.jobname, We = '91:00'), by=1:nrow(tgt)]  
		tgt[, dedup.cmd := paste0(dedup.cmd, " \"", cfg$bismark, '/deduplicate_bismark -p --output_dir ', sampledir, ' --bam --outfile ', dedup.bam.file, ' --samtools_path ', cfg$conda3bin, ' ', sorted.bam.file, "\"")]
		tgt[44, dedup.cmd]
		file.exists(tgt$dedup.bam.file)

		exe.jobs(tgt[, dedup.cmd], cfg$logger)

		## sort coord sorted bam file
		tgt[, bam.sort2.jobname := paste0(samplename, '.sort2.bam')]
		tgt[, bam.sort2.cmd := bsub.head(bam.sort2.jobname, mem=20, cpu=10, postdone = dedup.jobname, We = '30:00'), by=1:nrow(tgt)]  
		tgt[, bam.sort2.cmd := paste0(bam.sort2.cmd, " \"", cfg$JAVA, '/java -Xms256M -Xmx10g -Djava.io.tmpdir=', cfg$tmpdir, ' -jar ', cfg$PICARD, ' SortSam I=', dedup.bam.file, " O=", sorted.dedup.bam.file, " CREATE_INDEX=true SORT_ORDER=queryname\"")]
		tgt[6, bam.sort2.cmd]

		write(tgt$bam.sort3.cmd, file='bam.sort2.cmd.sh')
		exe.jobs(tgt[, bam.sort2.cmd], cfg$logger)
		file.size(tgt$dedup.bam.file)
	}

	## use methylKit processBismarkAln to call the methylation that can import into methylKit
	## CpG methylation levels
	tgt[, methcall.jobname := paste0(samplename, '.methcall')]
	tgt[, methcall.cmd := bsub.head(methcall.jobname, mem=20, cpu=2, We = '30:00'), by=1:nrow(tgt)]  
	tgt[, methcall.cmd := paste0(methcall.cmd, " \" ", cfg$conda3bin, "/Rscript --vanilla ~/pipeline/fun/methcall.r --loc ", coord.sorted.bam.file)]
	tgt[, methcall.cmd := paste0(methcall.cmd, " --outdir ", sampledir, " --sample.name ", samplename, " --treatment 1")]
	tgt[, methcall.cmd := paste0(methcall.cmd, " \" ")]
	tgt[, methcall.cmd]

	exe.jobs(tgt[, methcall.cmd], cfg$logger)
	tgt[, cpg.file := paste0(tgt$sampledir, '/', samplename, '_CpG.txt')]

	tgt[grep('187', samplename), ][, {scp(coord.sorted.bam.file)}]
	tgt[grep('187', samplename), coord.sorted.bam.file]
	file.exists(tgt$coord.sorted.bam.file)
	file.size(tgt$coord.sorted.bam.file)
	## SKIP deduplication BEGIN

	## collect alignment summary metrics
	tgt[, align.sum.jobname := paste0(samplename, '.align.sum')]
	tgt[, align.sum.cmd := bsub.head(align.sum.jobname, mem=10, hmem = 20, cpu=3, We = '10:00'), by = 1:nrow(tgt)]  
	tgt[, align.sum.cmd := paste0(align.sum.cmd, " \"", cfg$JAVA, '/java -Xms256M -Xmx10g -Djava.io.tmpdir=', cfg$tmpdir, ' -jar ', cfg$PICARD, ' CollectAlignmentSummaryMetrics R=', cfg$genomeGRCh37Fasta, " I=", coord.sorted.bam.file, " O=", align.sum.file, " \"")]
	tgt[, align.sum.cmd]

	exe.jobs(tgt[, align.sum.cmd], cfg$logger)

	## collect gc bias metrics
	tgt[, gc.bias.jobname := paste0(samplename, '.gc.bias')]
	tgt[, gc.bias.cmd := bsub.head(gc.bias.jobname, mem=10, hmem=10, cpu=3, We = '10:00'), by = 1:nrow(tgt)]  
	tgt[, gc.bias.cmd := paste0(gc.bias.cmd, " \"", cfg$JAVA, '/java -Xms256M -Xmx10g -Djava.io.tmpdir=', cfg$tmpdir, ' -jar ', cfg$PICARD, ' CollectGcBiasMetrics R=', cfg$genomeGRCh37Fasta, " I=", sorted.dedup.bam.file, " O=", gc.bias.metrics.file, " CHART=", gc.bias.metrics.pdf.file, " S=", gc.bias.metrics.sum.file, "\"")]
	tgt[6, gc.bias.cmd]

	exe.jobs(tgt[, gc.bias.cmd], cfg$logger)

	## collect hs metrics
	# mmuMethylInterval	= '/home/huw/program/truseq-methyl-capture-epic-manifest-file.bed'
	tgt[, hs.jobname := paste0(samplename, '.hs')]
	tgt[, hs.cmd := bsub.head(hs.jobname, mem=10, cpu=3, We = '10:00'), by = 1:nrow(tgt)]  
	tgt[, hs.cmd := paste0(hs.cmd, " \"", cfg$JAVA, '/java -Xms256M -Xmx10g -Djava.io.tmpdir=', cfg$tmpdir, ' -jar ', cfg$PICARD, ' CollectHsMetrics R=', cfg$genomeGRCh37Fasta, " I=", sorted.dedup.bam.file, " O=", hs.metrics.file, " TARGET_INTERVALS=", cfg$hsaMethylTargetInterval, " BAIT_INTERVALS=", cfg$hsaMethylBaitInterval, "\"")]
	exe.jobs(tgt[, hs.cmd], cfg$logger)

	## collect insert size
	tgt[, insert.size.jobname := paste0(samplename, '.insert.size')]
	tgt[, insert.size.cmd := bsub.head(insert.size.jobname, mem=10, cpu=3,  We = '10:00'), by = 1:nrow(tgt)]  
	tgt[, insert.size.cmd := paste0(insert.size.cmd, " \"", cfg$JAVA, '/java -Xms256M -Xmx10g -Djava.io.tmpdir=', cfg$tmpdir, ' -jar ', cfg$PICARD, ' CollectInsertSizeMetrics R=', cfg$genomeGRCh37Fasta, " I=", sorted.dedup.bam.file, " O=", insert.size.metrics.file, " H=", insert.size.metrics.histgram.file, " M=0.5\"")]
	exe.jobs(tgt[, insert.size.cmd], cfg$logger)

	## collect targeted pcr metrics
	tgt[, targeted.pcr.jobname := paste0(samplename, '.targeted.pcr')]
	tgt[, targeted.pcr.cmd := bsub.head(targeted.pcr.jobname, mem=10, cpu=3, We = '10:00'), by = 1:nrow(tgt)]  
	tgt[, targeted.pcr.cmd := paste0(targeted.pcr.cmd, " \"", cfg$JAVA, '/java -Xms256M -Xmx10g -Djava.io.tmpdir=', cfg$tmpdir, ' -jar ', cfg$PICARD, ' CollectInsertSizeMetrics R=', cfg$genomeGRCh37Fasta, " I=", sorted.dedup.bam.file, " O=", targeted.pcr.metrics.file, " H=", targeted.pcr.metrics.histgram.file, " M=0.5\"")]
	exe.jobs(tgt[, targeted.pcr.cmd], cfg$logger)

	## depth per base
	tgt[, depth.jobname := paste0(samplename, '.depth')]
	tgt[, depth.cmd := bsub.head(depth.jobname, mem=10, We ='30:00', cpu=1), by=1:nrow(tgt)]  
	tgt[, depth.cmd := paste0(depth.cmd, " \"", cfg$samtools, " depth -a ", sorted.dedup.bam.file, '>', depth.file, '"')]
	exe.jobs(tgt[, depth.cmd], cfg$logger)


	#tgt[, {cmd = paste0('rsync -avur ', sorted.bam.file, ' mski1925:/Volumes/LaCie/huw/solit/study/hiseq/Project_06000_DD/r_fang/initFiles/bamfile/'); system(cmd)}, by=1:nrow(tgt)]
	#tgt[, sorted.bai.file := sub('bam', 'bai', sorted.bam.file)]
	#tgt[, {cmd = paste0('rsync -avur ', sorted.bai.file, ' mski1925:/Volumes/LaCie/huw/solit/study/hiseq/Project_06000_DD/r_fang/initFiles/bamfile/'); system(cmd)}, by=1:nrow(tgt)]

}
