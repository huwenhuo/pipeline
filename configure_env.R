#library(log4r)
#library(data.table)

#cfg$cwd = getwd()

source('~/pipeline//fun/param.r')

if(cfg$assay == 'g10x'){
	cfg$statsDir = paste0(cfg$resDir, '/stats')
	system(paste0('mkdir -p ', cfg$statsDir)) 
	cfg$initDir = paste0(cfg$resDir, '/initFiles')
	system(paste0('mkdir -p ', cfg$initDir)) 
	cfg$g10x= paste0(cfg$resDir, '/g10x')
	dir.create(cfg$g10x)
}

if(cfg$assay == '' | is.null(cfg$assay)){
	source('~/pipeline//fun/param.r')
	cfg$statsDir = paste0(cfg$resDir, '/stats')
	system(paste0('mkdir -p ', cfg$statsDir)) 
}

source(paste0('~/pipeline/configure_env_dir.R'))

if(cfg$assay == 'exonseq'){
	source('~/pipeline//fun/param.r')
	cfg$initDir = paste0(cfg$resDir, '/initFiles')
	cfg$matricsDir = paste0(cfg$resDir, '/matrics')
	#cfg$fastqcDir = paste0(cfg$resDir, '/fastqc')
	cfg$statsDir = paste0(cfg$resDir, '/stats')
	#cfg$rseqcDir = paste0(cfg$resDir, '/RSeqC')
	#cfg$geoDir= paste0(cfg$resDir, '/geo')
	#system(paste0('mkdir -p ', cfg$resDir)) ## r_fang
	#system(paste0('mkdir -p ', cfg$matricsDir))  ## collected metrics
	#system(paste0('mkdir -p ', cfg$initDir)) 
	#system(paste0('mkdir -p ', cfg$fastqcDir)) 
	system(paste0('mkdir -p ', cfg$statsDir)) 
	#system(paste0('mkdir -p ', cfg$rseqcDir)) 
	#system(paste0('mkdir -p ', cfg$geoDir)) 
}


if(cfg$assay == 'rnaseq'){
	source('~/pipeline/fun/param.r')
	cfg$initDir = paste0(cfg$resDir, '/initFiles')
	cfg$matricsDir = paste0(cfg$resDir, '/matrics')
	cfg$fastqcDir = paste0(cfg$resDir, '/fastqc')
	cfg$statsDir = paste0(cfg$resDir, '/stats')
	cfg$rseqcDir = paste0(cfg$resDir, '/RSeqC')
	cfg$geoDir= paste0(cfg$resDir, '/geo')
	#system(paste0('mkdir -p ', cfg$resDir)) ## r_fang
	system(paste0('mkdir -p ', cfg$matricsDir))  ## collected metrics
	system(paste0('mkdir -p ', cfg$initDir)) 
	system(paste0('mkdir -p ', cfg$fastqcDir)) 
	system(paste0('mkdir -p ', cfg$statsDir)) 
	system(paste0('mkdir -p ', cfg$rseqcDir)) 
	system(paste0('mkdir -p ', cfg$geoDir)) 
}

if(cfg$assay == 'methylation'){
	source('~/pipeline/fun/param.r')
	cfg$initDir = paste0(cfg$resDir, '/initFiles')
	cfg$progressDir = paste0(cfg$resDir, '/progress')
	cfg$matricsDir = paste0(cfg$resDir, '/matrics')
	cfg$fastqcDir = paste0(cfg$resDir, '/fastqc')
	cfg$statsDir = paste0(cfg$resDir, '/stats')
	cfg$methylDir = paste0(cfg$resDir, '/methyl')
	cfg$geoDir= paste0(cfg$resDir, '/geo')
	system(paste0('mkdir -p ', cfg$resDir)) ## r_fang
	system(paste0('mkdir -p ', cfg$progressDir)) ## jobname.done
	system(paste0('mkdir -p ', cfg$matricsDir))  ## collected metrics
	system(paste0('mkdir -p ', cfg$initDir)) 
	system(paste0('mkdir -p ', cfg$fastqcDir)) 
	system(paste0('mkdir -p ', cfg$statsDir)) 
	system(paste0('mkdir -p ', cfg$methylDir)) 
	system(paste0('mkdir -p ', cfg$geoDir)) 
}

cfg$logfile = paste0(cwd, '/', '_', gsub("/", "_", cfg$resDir), '_log')
create.logger() -> cfg$logger
logfile(cfg$logger) = cfg$logfile
level(cfg$logger) = 'INFO' 

## provide for everywhere

logger = cfg$logger
## this for set up @ENV
path			= '/home/huw/pipeline/fun//homer/bin:/home/huw/local/bin:/home/huw/perl5/CPAN/bin:/home/huw/program/cufflink221/:/home/huw/program/tophat213/:/home/huw/program/bowtie2/:/home/huw/program/bin:/home/huw/program/blat:/home/huw/program/weblogo:/home/huw/program/gs/bin:/home/huw/program/ngsplot:$PATH'
JAVA_HOME 		= '/home/huw/pipeline/fun//jdk7/bin/java'
JAVA_HOME 		= '/opt/common/CentOS_6-dev/java/jdk1.8.0_31/'
NGSPLOT 		= '/home/huw/pipeline/fun//ngsplot'
TMPDIR			= '/scratch/huw'

if(is.null(cfg$species)){ stop("species is not assigned")}

if(cfg$species == 'mm9') {
	cfg$genomeFasta = cfg$MM9_FASTA
	cfg$genomeBWA =  cfg$MM9_BWA_INDEX
	cfg$genomeFAI = cfg$MM9_FAI
	cfg$DB_SNP = "" 
	cfg$FP_INT = "" 
	cfg$FP_TG  = ""
}

if(cfg$species == 'mm10_custom') {
	cfg$genomeFasta = cfg$MM10_CUSTOM_FASTA
	cfg$genomeBWA =  cfg$MM10_BWA_INDEX
	cfg$DB_SNP = paste0(cfg$dataDir, "/mm10/mm10_snp142.vcf") 
	cfg$FP_INT = "" 
	cfg$FP_TG  = ""
}

if(cfg$species == 'mm10') {
	cfg$genomeFasta = cfg$MM10_FASTA
	cfg$genomeBWA =  cfg$MM10_BWA_INDEX
	cfg$DB_SNP = paste0(cfg$dataDir, "/mm10/mm10_snp142.vcf") 
	cfg$FP_INT = "" 
	cfg$FP_TG  = ""
}

if(cfg$species == 'hybrid') {
	cfg$genomeFasta = cfg$B37_MM10_HYBRID_FASTA
	#genomeBWA =  B37_BWA_INDEX
	cfg$DB_SNP = paste0(cfg$dataDir, "/b37/dbsnp_138.b37.vcf")
	cfg$FP_INT = paste0(cfg$dataDir, "/b37/Agilent51MBExome__b37__FP_intervals.list")
	cfg$FP_TG  = paste0(cfg$dataDir, "/b37/Agilent51MBExome__b37__FP_tiling_genotypes.txt")
}

if(grepl('hg19', cfg$species, ignore.case=T)){
	   cfg$genomeFasta = cfg$HG19_FASTA
	   #bismarkHg19
	   cfg$genomeFAI = cfg$HG19_FAI
}

## configure files
if(cfg$assay == 'dnaseq'){
	source('~/pipeline/fun//fun/param.r')
	cfg$pairfile = paste0(cfg$resDir, '/', cfg$pre, '_sample_pairing.txt')
	cfg$mapfile = paste0(cfg$resDir, '/', cfg$pre, '_sample_mapping.txt')
	cfg$groupfile = paste0(cfg$resDir, '/', cfg$pre, '_sample_grouping.txt')

	cfg$statsDir = paste0(cfg$resDir, '/stats')
	system(paste0('mkdir -p ', cfg$statsDir)) ## err std files

	if(!file.exists(cfg$groupfile)){
		ifgroupfile = F
		#stop('groupfile not exist.')
	}

	if(!file.exists(cfg$mapfile)){
		ifmapfile = F
		#stop('mapfile not exist.')
	}

	if(is.null(cfg$species)){
		cat('Species must provided\n\n')
		cat('Environment parameters are not imported!!\n\n')
		return
	}

	if(is.null(cfg$wesImpactTarget)){
		cat('wesImpactTarget must provided\n\n')
		cat('Environment parameters are not imported!!\n\n')
		return
	}

	cfg$initDir = paste0(cfg$resDir, '/initFiles')
	cfg$progressDir = paste0(cfg$resDir, '/progress')
	cfg$matricsDir = paste0(cfg$resDir, '/matrics')
	cfg$varDir = paste0(cfg$resDir, '/variation')
	cfg$haploDir = paste0(cfg$varDir, '/haplotypecaller')
	cfg$mutectDir = paste0(cfg$varDir, '/mutect')
	cfg$sniperDir = paste0(cfg$varDir, '/somaticsniper')
	cfg$facetsDir = paste0(cfg$varDir, '/facets')
	cfg$strvarDir = paste0(cfg$varDir, '/strvar')
	system(paste0('mkdir -p ', cfg$resDir)) ## r_fang
	system(paste0('mkdir -p ', cfg$progressDir)) ## jobname.done
	system(paste0('mkdir -p ', cfg$matricsDir))  ## collected metrics
	system(paste0('mkdir -p ', cfg$varDir))  
	system(paste0('mkdir -p ', cfg$haploDir)) 
	system(paste0('mkdir -p ', cfg$mutectDir)) 
	system(paste0('mkdir -p ', cfg$sniperDir)) 
	system(paste0('mkdir -p ', cfg$facetsDir)) 
	system(paste0('mkdir -p ', cfg$facetsDir, '/sample')) 
	system(paste0('mkdir -p ', cfg$strvarDir)) 
	system(paste0('mkdir -p ', cfg$initDir)) 

	cfg$wesImpactTarget = 'AgilentExon_51MB_b37_v3' # baits
	cfg$baits_ilist = paste0(cfg$targetsDir, '/', cfg$wesImpactTarget, '/', cfg$wesImpactTarget, '_baits.ilist')
	cfg$targets_ilist = paste0(cfg$targetsDir, '/', cfg$wesImpactTarget, '/', cfg$wesImpactTarget, '_targets.ilist')
	cfg$targets_bed = paste0(cfg$targetsDir, '/', cfg$wesImpactTarget, '/', cfg$wesImpactTarget, '_targets.bed')
	cfg$targets5bp_ilist = paste0(cfg$targetsDir, '/', cfg$wesImpactTarget, '/', cfg$wesImpactTarget, '_targets_plus5bp.ilist')
	cfg$targets5bp_bed = paste0(cfg$targetsDir, '/', cfg$wesImpactTarget, '/', cfg$wesImpactTarget, '_targets_plus5bp.bed')

}

## this for set up @ENV
path			= '/home/huw/pipeline/fun//homer/bin:/home/huw/local/bin:/home/huw/perl5/CPAN/bin:/home/huw/program/cufflink221/:/home/huw/program/tophat213/:/home/huw/program/bowtie2/:/home/huw/program/bin:/home/huw/program/blat:/home/huw/program/weblogo:/home/huw/program/gs/bin:/home/huw/program/ngsplot:$PATH'
JAVA_HOME 		= '/home/huw/pipeline/fun//jdk7/bin/java'
JAVA_HOME 		= '/opt/common/CentOS_6-dev/java/jdk1.8.0_31/'
NGSPLOT 		= '/home/huw/pipeline/fun//ngsplot'
TMPDIR			= '/scratch/huw'

if(cfg$species == 'mm9') {
	cfg$genomeFasta = cfg$MM9_FASTA
	cfg$genomeBWA =  cfg$MM9_BWA_INDEX
	cfg$genomeFAI = cfg$MM9_FAI
	cfg$DB_SNP = "" 
	cfg$FP_INT = "" 
	cfg$FP_TG  = ""
}

if(cfg$species == 'mm10_custom') {
	cfg$genomeFasta = cfg$MM10_CUSTOM_FASTA
	cfg$genomeBWA =  cfg$MM10_BWA_INDEX
	cfg$DB_SNP = paste0(cfg$dataDir, "/mm10/mm10_snp142.vcf") 
	cfg$FP_INT = "" 
	cfg$FP_TG  = ""
}

if(cfg$species == 'mm10') {
	cfg$genomeFasta = cfg$MM10_FASTA
	cfg$genomeBWA =  cfg$MM10_BWA_INDEX
	cfg$DB_SNP = paste0(cfg$dataDir, "/mm10/mm10_snp142.vcf") 
	cfg$FP_INT = "" 
	cfg$FP_TG  = ""
}

if(cfg$species == 'hybrid') {
	cfg$genomeFasta = cfg$B37_MM10_HYBRID_FASTA
	#genomeBWA =  B37_BWA_INDEX
	cfg$DB_SNP = paste0(cfg$dataDir, "/b37/dbsnp_138.b37.vcf")
	cfg$FP_INT = paste0(cfg$dataDir, "/b37/Agilent51MBExome__b37__FP_intervals.list")
	cfg$FP_TG  = paste0(cfg$dataDir, "/b37/Agilent51MBExome__b37__FP_tiling_genotypes.txt")
}

if(grepl('hg19', cfg$species, ignore.case=T)){
	   cfg$genomeFasta = cfg$HG19_FASTA
	   cfg$genomeBWA = cfg$HG19_BWA_INDEX
	   cfg$genomeFAI = cfg$HG19_FAI
	   cfg$DB_SNP = paste0(cfg$dataDir, "/hg19/dbsnp_138.hg19.vcf")
	   cfg$FP_INT = paste0(cfg$dataDir, "/hg19/Agilent51MBExome__hg19__FP_intervals.list")
	   cfg$FP_TG  = paste0(cfg$dataDir, "/hg19/Agilent51MBExome__hg19__FP_tiling_genotypes.txt")
	   cfg$ExAC_VCF = paste0(cfg$VEP, "/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz");
	   cfg$FACETS_DB_SNP = paste0(cfg$dataDir, "/hg19/dbsnp_137.hg19__RmDupsClean__plusPseudo50__DROP_SORT.vcf");
	   cfg$MILLS_1000G = paste0(cfg$dataDir, "/hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf");
	   cfg$HAPMAP = paste0(cfg$dataDir, "/hg19/hapmap_3.3.hg19.vcf");
	   cfg$OMNI_1000G = paste0(cfg$dataDir, "/hg19/1000G_omni2.5.hg19.vcf");
	   cfg$PHASE1_SNPS_1000G = paste0(cfg$dataDir, "/hg19/1000G_phase1.snps.high_confidence.hg19.vcf");
	   cfg$COSMIC = paste0(cfg$dataDir, "/hg19/CosmicCodingMuts_v67_20131024.vcf");
}

if(grepl('b37', cfg$species, ignore.case=T)){
	cfg$genomeFasta = cfg$B37_FASTA
	cfg$genomeBWA =  cfg$B37_BWA_INDEX
	cfg$genomeFAI = cfg$B37_FAI
	cfg$DB_SNP = paste0(cfg$dataDir, "/b37/dbsnp_138.b37.vcf")
	cfg$FP_INT = paste0(cfg$dataDir, "/b37/Agilent51MBExome__b37__FP_intervals.list")
	cfg$FP_TG  = paste0(cfg$dataDir, "/b37/Agilent51MBExome__b37__FP_tiling_genotypes.txt")
	cfg$ExAC_VCF = paste0(cfg$VEP, "/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz");
	cfg$FACETS_DB_SNP = paste0(cfg$dataDir, "/b37/dbsnp_137.b37__RmDupsClean__plusPseudo50__DROP_SORT.vcf");
	cfg$MILLS_1000G = paste0(cfg$dataDir, "/b37/Mills_and_1000G_gold_standard.indels.b37.vcf");
	cfg$HAPMAP = paste0(cfg$dataDir, "/b37/hapmap_3.3.b37.vcf");
	cfg$btw2db = paste0('/juno/work/solitlab/huw/study/db/hsa/hg19-bowtie2/hg19')
	cfg$OMNI_1000G = paste0(cfg$dataDir, "/b37/1000G_omni2.5.b37.vcf");
	cfg$PHASE1_SNPS_1000G = paste0(cfg$dataDir, "/b37/1000G_phase1.snps.high_confidence.b37.vcf");
	cfg$COSMIC = paste0(cfg$dataDir, "/b37/CosmicCodingMuts_v67_b37_20131024__NDS.vcf");
	cfg$COSMIC_HOTSPOTS = paste0(cfg$dataDir, "/b37/dmp_cosmic_for_hotspots.vcf");
	cfg$vep.species = 'homo_sapiens'
}


#>Parameters
isPE			= '1'
clusters		= 'LSF'
#>Versions

#>Convention
##name1fastq		= 'base_xx_R1.fastq.gz'
##name2fastq		= 'base_xx_R2.fastq.gz'
##name1mrgedfq		= 'base_R1.fq.gz'
##name2mrgedfq		= 'base_R2.fq.gz'
##name1trimmedfq		= 'base_ptrimmed_R1.fq.gz'
##name2trimmedfq		= 'base_ptrimmed_R2.fq.gz'
##sai1			= 'base_R1.sai'
##sai2			= 'base_R2.sai'
##sam			= 'base_genome.sam'
##bam			= 'base_genome'
##bamext			= 'base_genome.bam'
##bamsorted		= 'base_genome_sorted'
##bamsortedext		= 'base_genome_sorted.bam'
##bamsortedrmdup		= 'base_genome_sorted_rmdup.bam'
##htseqcount		= 'base.htseq.count'
##bigwig			= 'base_genome_sorted_rmdup_10m.bw'
##MACS14_peaks_file 	= 'base_genome_macs14_peaks.bed'
##MACS2_peaks_file	= 'base_genome_macs2_peaks.xls'
##MACS2_broad_peaks_file  = 'base_genome_macs2broad_peaks.broadPeak'
##peakannoout		= 'base_genome_macs2_peaks_anno.txt'
##mcspeaks		= 'base_genome_macs2_peaks.narrowPeak'
##mcsbroadpeaks		= 'base_genome_macs2broad_peaks.broadPeak'
##peakannout		= 'base_genome_macs2_peaks_anno.txt'
##peakannout2		= 'base_genome_macs2broad_peaks_anno.txt'
##peakannoFiltered	= 'base_genome_macs2_peaks_anno_filtered.txt'
##peakannoFiltered2	= 'base_genome_macs2broad_peaks_anno_filtered.txt'
##ngsplot_file		= 'ngsplot_base.tss.avgprof.pdf'
##ngsplot_file_heatmap	= 'ngsplot_base.tss.heatmap.pdf'
##TSSCounts		= 'base_genome_tss2kb.bed'
##tssout			= 'ngsplot.base.tss'
##genebodyout		= 'ngsplot.base.genebody'
##enhancerout		= 'ngsplot.base.enhancer'
##motifout		= 'motif'
##starbam        		= 'base.Aligned.sortedByCoord.out.bam'
##
##>Checkpoints
##mergefastq		= 'name1mrgedfq'
##trimfq			= 'name1trimmedfq'
##alignStarAlign		= 'bamsortedext'
##alignBowtie2		= 'bamsrotedext'
##alignBwaAln		= 'sai1'
##alignBwaSampe		= 'sam'
##sortSam			= 'bamsortedext'
##sam2bam			= 'bamsortedext'
##removeDup		= 'bamsortedrmdup'
##bigwig			= 'bigwig'
##indexBam		= 'bamInput_bai'
##MACS2			= 'mcspeaks'
##MACS2BroadPeak		= 'mcsbroadpeaks'
##annotatePeaksHomer	= 'peakannoout_peakannoout2'
##annotateFilter		= 'peakannoFiltered_peakannoFiltered2'
##macs14			= 'mcspeaks'
##TSSCounts		= 'TSSCounts'
##Ngsplot			= 'tssout.avgprof.pdf'
##motifHomerFinder	= 'motifout'
##htseqcount		= 'htseqcount'

## below could import from ~/pipeline//fun/param.r
simple_table = data.table(
			  simple_val = c('Truncating', 'Truncating', 'Inframe', 'Inframe', 
					 'Missense_Mutation', 'Truncating', 'Inframe', 
					 'Inframe', 'Inframe', 'Inframe', 'Inframe'),
			  simple_name = c('Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', 
					  'Missense_Mutation','Nonsense_Mutation', 'Nonstop_Mutation', 
					  'Silent', 'Splice_Region', 'Splice_Site', 'Translation_Start_Site')
			  )
simple_table
setkey(simple_table, 'simple_name')
color_table  = data.table(
			 color_index = c("Missense_Mutation Oncogenic",
					 "Truncating Oncogenic",
					 "Inframe Oncogenic",
					 "Missense_Mutation Passage",
					 "Truncating Passage",
					 "Inframe Passage") ,
			 color_value = c("#bebebe", 
					 "#000000", 
					 "#993404", 
					 adjustcolor("#bebebe", alpha.f=0.6), 
					 adjustcolor("#000000", alpha.f=0.6),
					 adjustcolor("#993404", alpha.f=0.6)) 
			 )
color_table
setkey(color_table, 'color_index')

rsem.fpkm = function(rsem){
	rsem$fpkm = rsem$counts * 10^9 / rsem$length
	rsem$fpkm = sweep(rsem$fpkm, 2, colSums(rsem$counts), '/')
	rsem
}

del.files = function(filesL){
	if(all(file.exists(files))){
		for(i in 1:nrow(files)){
			info(logger, paste('deleting file:', files[i]))
			system(paste0('rm ', files[i]))
		}
	}else{
		info(logger, paste0('Not all files can be found. Nothing are deleted!'))
	}
}

source('~/pipeline/fun/jobs.r')
