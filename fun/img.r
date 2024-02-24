

img = system(paste0('find /juno/work/taylorlab/cmopipeline/singularity_images/ -name *.img'), intern=T)
img.name = sub('.img', '', basename(img))
img.name = gsub('cmopipeline-', '', img.name)
img.name = gsub('-', '.', img.name)
names(img) = img.name
img = as.list(img)
img

gnome = list(
	     acLoci           = "/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/smallGRCh37/1000G_phase3_20130502_SNP_maf0.3.small.loci",
	     acLociGC         = "/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/smallGRCh37/1000G_phase3_20130502_SNP_maf0.3.small.loci.gc",
	     dbsnp            = "/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/smallGRCh37/dbsnp_138.b37.small.vcf",
	     dbsnpIndex       = "/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/smallGRCh37/dbsnp_138.b37.small.vcf.idx",
	     genomeDict       = "/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/smallGRCh37/human_g1k_v37_decoy.small.dict",
	     genomeFile       = "/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/smallGRCh37/human_g1k_v37_decoy.small.fasta",
	     genomeIndex      = "/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/smallGRCh37/human_g1k_v37_decoy.small.fasta.fai",
	     intervals        = "/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/smallGRCh37/small.intervals",
	     knownIndels      = "/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/smallGRCh37/1000G_phase1.indels.b37.small.vcf",
	     knownIndels      = "/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/smallGRCh37/Mills_and_1000G_gold_standard.indels.b37.small.vcf",
	     msiSensorList    = "/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/smallGRCh37/small.msi.list",
	     snpeffDb         = "GRCh37.75",
	     vepCacheVersion  = "95",
	     svCallingExcludeRegions 	= "/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/smallGRCh37/human.hg19.excl.tsv",
	     svCallingIncludeRegions 	= "/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/smallGRCh37/b37.test.bed.gz",
	     idtTargets 		= "/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/smallGRCh37/bedfile.idt.bed.gz",
	     agilentTargets 		= "/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/smallGRCh37/bedfile.agilent.bed.gz",
	     wgsTargets 		= "/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/smallGRCh37/bedfile.wgs.bed.gz",
	     exomePoN 			= "/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/smallGRCh37/pon_test.vcf.gz",
	     wgsPoN 			= "/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/smallGRCh37/pon_test.1.vcf.gz",
	     idtTargetsList 		= "/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/smallGRCh37/IDT_Exome_v1_FP_b37_targets.small.interval_list",
	     idtBaitsList 		= "/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/smallGRCh37/IDT_Exome_v1_FP_b37_baits.small.interval_list",
	     agilentTargetsList 	= "/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/smallGRCh37/AgilentExon_51MB_b37_v3_targets.small.interval_list",
	     agilentBaitsList 		= "/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/smallGRCh37/AgilentExon_51MB_b37_v3_baits.small.interval_list" )

