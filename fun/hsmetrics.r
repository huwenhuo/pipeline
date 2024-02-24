hsmetrics = function(tg){
	# tg has col of bamfile
	# bases
	require(data.table)
	require(log4r)
	cfg=new.env()
	cfg$assay = ''
	cfg$species='none'
	cwd = getwd()
	source('~/pipeline/configure_env.R')
	tg[, hmsFile := paste0(cfg$tmpdir, '/', 'hms_', bases, '_HsMatrix.txt')]
	tg[, hms.jobname := paste0('hms.', bases)]
	tg[, hms.cmd := bsub.head(hms.jobname, mem=10, cpu=1, We='8:26', cwd=cwd, statsDir = '.')]
	tg[, hms.cmd := paste0(hms.cmd, ' "', cfg$JAVA, '/java -Xms256m -Xmx30g -Djava.io.tmpdir=', cfg$tmpdir, ' -jar ', cfg$PICARD, '/picard.jar CalculateHsMetrics I=', V1, ' TMP_DIR=', cfg$tmpdir, ' O=', hmsFile, ' REFERENCE_SEQUENCE=', cfg$genomeGRCh37Fasta, ' METRIC_ACCUMULATION_LEVEL=SAMPLE BAIT_INTERVALS=', cfg$AgilentExon51Mbv3BaitsIlist, ' BAIT_SET_NAME=', cfg$AgilentExon51Mbv3tgt, ' TARGET_INTERVALS=', cfg$agilentExon51v3tgtIlist, ' VALIDATION_STRINGENCY=LENIENT "')]
	tg$hms.cmd

	exe.jobs(targets.smpBam$hms.cmd, logger)
	tg
}
