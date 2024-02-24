cfg$blacklist.hg19		= '/home/huw/program/Blacklist/lists/Blacklist_v1/hg19-blacklist.bed.gz'
cfg$rsem.sq			= '/juno/work/solitlab/huw/solit/study/hiseq/blca_cmo_06155_2016/rnaseq_rsem_raw_counts.RData'
cfg$rsem.mcp			= '/juno/work/solitlab/huw/solit/study/hiseq/micropapillary/rsem.RData'
cfg$ssgsea 			= '/home/huw/program/ssGSEA2.0/src/ssGSEA2.0.R'
cfg$grch38.gtf.file 		= '/juno/work/solitlab/huw/igenome/Homo_sapiens.GRCh38.88.gtf'
cfg$p27 			= '/home/huw/anaconda3/envs/p27/bin/python2.7'
cfg$tmp.dir 			= "/scracth/huw"
cfg$genomeFile 			= "/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta"
cfg$genomeFileFai		= '/juno/work/taylorlab/cmopipeline/mskcc-igenomes/igenomes/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.fasta.fai' 
cfg$singularity 		= '/opt/local/singularity/3.3.0/bin/singularity'
cfg$vep.dir			= '/juno/work/taylorlab/cmopipeline/mskcc-igenomes/grch37/vep/'
cfg$vep.data			= '/juno/work/taylorlab/cmopipeline/mskcc-igenomes/grch37/vep/'
cfg$vep.cache			= '/juno/work/taylorlab/cmopipeline/mskcc-igenomes/grch37/vep/cache/'
cfg$simg 			= '/juno/work/taylorlab/cmopipeline/singularity_images/'
cfg$simg.vcf2maf		= paste0(cfg$simg, '/cmopipeline-vcf2maf-1.6.17.img')
cfg$polysolver_bin		= '/juno/res/taylorlab/bandlamc/neoantigens/pipelines/polysolver/shell_call_hla_type_singularity'
cfg$neoantigen 			= '/home/huw/program/neoantigen/neoantigen-dev/neoantigen.py'
cfg$neoantigen.config		= '/home/huw/program/neoantigen/neoantigen-dev/neoantigen-luna.config'
cfg$oncokb			= '~/program/oncokb-annotator/'
cfg$hg19.picard.list		= '~/program/refgene_hg19.txt'
cfg$cmo				='~/program/cmo/cmo'
cfg$ngsfilter			='~/program/ngs-filters'
cfg$vcf2maf			='~/program/vcf2maf/'
cfg$vcf2maf			='/opt/common/CentOS_6-dev/vcf2maf/v1.6.17/'
cfg$vep.dir			='/opt/common/CentOS_6-dev/vep'
cfg$vep.cache			='/opt/common/CentOS_6-dev/vep/cache'
cfg$vep.ver			='95'
cfg$omni			= '/juno/work/solitlab/huw/program/socciCode/Pipelines/CBE/Variant/variants_pipeline/data/b37/1000G_omni2.5.b37.vcf'
cfg$hapmap 			= '/juno/work/solitlab/huw/program/socciCode/Pipelines/CBE/Variant/variants_pipeline/data/b37/hapmap_3.3.b37.vcf'
cfg$dbsnp			= '/juno/work/solitlab/huw/program/socciCode/Pipelines/CBE/Variant/variants_pipeline/data/b37/dbsnp_138.b37.vcf'
cfg$G1000G			= '/juno/work/solitlab/huw/program/socciCode/Pipelines/CBE/Variant/variants_pipeline/data/b37/1000G_phase1.snps.high_confidence.b37.vcf'
cfg$cranger.ref    		= '/home/huw/program/cellranger/cellranger-3.0.2/refdata-cellranger-GRCh38-3.0.0/'
cfg$cranger			= '/home/huw/program/cellranger/cellranger-3.0.2/cellranger'
cfg$funcotatorData.somatic	= '/home/huw/program/funcotator/funcotator_dataSources.v1.6.20190124s'
cfg$bic				= '/home/huw/program/BIC-variants_pipeline'
cfg$program 			= '/home/huw/program'
cfg$pipeline			= '/home/huw/pipeline'
cfg$cibersort			= paste0(cfg$program, '/cibersort')
cfg$anaconda2 			= paste0(cfg$program, '/anaconda2')
cfg$anaconda3 			= paste0('/home/huw/anaconda3')
cfg$conda3bin 			= paste0('/home/huw/conda3/envs/r36/bin')
cfg$hsa.housekeeping.bed	= paste0(cfg$program, '/hg19_RefSeq.bed')
cfg$hsa.housekeeping.bed.ii	= paste0(cfg$program, '/hg19_RefSeq_ii.bed')
cfg$hsaMethylBaitIntervalOrig	= paste0(cfg$program, '/truseq-methyl-capture-epic-manifest-file.bed')
cfg$hsaMethylBaitInterval	= paste0(cfg$program, '/truseq-methyl-capture-epic-manifest-file-bait.interval')
cfg$hsaMethylTargetInterval	= paste0(cfg$program, '/truseq-methyl-capture-epic-manifest-file-plusminus5bp.interval')
cfg$hsaMethylAnno		= paste0(cfg$program, '/truseq-methyl-capture-epic-manifest-file.annotated.islands.txt')
cfg$hg19.cpgi			= paste0(cfg$program, '/hg19_cpg_island.bed.txt')
cfg$hg38.cpgi			= paste0(cfg$program, '/cpgIslandExt.hg38.bed')
cfg$lncGTFGRCh38		= paste0(cfg$program, '/lncipedia_5_2_hg38.gtf')
cfg$hg38.refseq			= paste0(cfg$program, '/hg38_RefSeq.bed')
cfg$hg38.refseq.n		= paste0(cfg$program, '/hg38_refseq_n.bed')
cfg$deeptools			= paste0(cfg$program, '/anaconda3/bin/')
cfg$getGenomeBuild		= '/home/socci/Code/FillOut/FillOut/GenomeData/getGenomeBuildBAM.sh' 
cfg$mm10.refseq.bed		= '~/program/mm10_RefSeq.bed'

cfg$fastqc 			= paste0(cfg$program, '/FastQC/fastqc')
#options(scipen=999)
#system(paste0('cp /home/huw/program/truseq-methyl-capture-epic-header.txt ', hsaMethylBaitInterval))
#tmp = fread(hsaMethylBaitIntervalOrig)
#tmp = tmp[order(V1),]
#library(tibble)
#tmp = add_column(tmp, string = "+", .after = 3)
#tmp[, V1:=sub("chr", "", V1)]
#fwrite(tmp, file=hsaMethylBaitInterval, quote=F, sep="\t", col.names=F, append=T)
#write.table(as.data.frame(tmp), file=hsaMethylBaitInterval, quote=F, sep="\t", col.names=F, append=T, row.names=F)
#tmp[, V2:=V2-5]
#tmp[, V3:=V3+5]
#head(tmp)
#tmp[, V1 := as.character(V1)]
#tmp[, V2 := as.character(V2)]
#tmp[, V3 := as.character(V3)]
#tmp[, V4:=paste0("chr", V1, "_", V2, "-", V3)]
#system(paste0('cp /home/huw/program/truseq-methyl-capture-epic-header.txt ', hsaMethylTargetInterval))
#write.table(as.data.frame(tmp), file=hsaMethylTargetInterval, quote=F, sep="\t", col.names=F, append=T, row.names=F)
#rm(tmp)

cfg$bismarkGRCh37		= '/juno/work/solitlab/huw/study/db/hsa/grch37_bismark'
cfg$bismark			= paste0(cfg$program, '/Bismark-0.23.0')
cfg$methylTarget		= paste0(cfg$program, '/truseq-methyl-capture-epic-manifest-file.bed')
cfg$SOAPHLA			= paste0(cfg$program, '/SOAP-HLA')
cfg$tmpdir			= '/scratch/huw'
cfg$mutsig			= '/opt/common/CentOS_6-dev/mutsig/cv-1.4/'

cfg$belly2			= paste0(cfg$program, '/delly_v0.7.7_linux_x86_64bit')
cfg$bamReadCount		= paste0(cfg$program, 'am-readcount')
cfg$IEDB			= paste0(cfg$program, '/iedb')
cfg$IEDBMHC1			= paste0(cfg$program, '/iedb/mhc_i')
cfg$IEDBMHC2			= paste0(cfg$program, '/iedb/mhc_ii')
cfg$pvacseq			= paste0(cfg$program, '/anaconda3/bin/pvacseq')
cfg$polysolver			= paste0(cfg$program, '/polysolver')
cfg$ponfile			= paste0(cfg$program, '/pon56Sample.vcf.gz')
cfg$make.trinuc			= paste0(cfg$program, '/make_trinuc_maf.py')
cfg$ConvertQualityScore    	= paste0(cfg$program, '/ConvertQualityScore')
cfg$STAR			= paste0(cfg$program, '/STAR-STAR_2.5.0a/bin/Linux_x86_64_static/STAR')
cfg$STAR_DIR			= paste0(cfg$program, '/STAR-STAR_2.5.0a/bin/Linux_x86_64_static')
cfg$sambamba			= paste0(cfg$program, '/sambamba' )
cfg$java			= paste0(cfg$program, '/jdk7/bin/java')
cfg$picard			= paste0(cfg$program, '/picard2.20.jar')
cfg$NGSPLOTexe			= paste0(cfg$program, '/ngsplot/bin/ngs.plot.r')
cfg$NGSFilter			= paste0(cfg$program, '/ngs-filters/')
cfg$trim_galore			= paste0(cfg$program, '/trim_galore')
cfg$bowtie2			= paste0(cfg$program, '/bowtie2/bowtie2')
cfg$bowtie2dir                  = paste0(cfg$program, '/bowtie2')
cfg$homerFindMotif              = paste0(cfg$program, '/homer/bin/findMotifsGenome.pl')
cfg$homerAnnotatePeaks		= paste0(cfg$program, '/homer/bin/annotatePeaks.pl')
cfg$tssCounts			= paste0(cfg$program, '/tsscounts.pl')
cfg$bwa				= paste0(cfg$program, '/bwa-0.7.12/bwa ')
cfg$bedtools			= paste0(cfg$program, '/bedtools/bin')
cfg$bedtools2			= paste0(cfg$program, '/bedtools2/bin')
cfg$bam2bigwig			= paste0(cfg$program, '/bam2bigwig.sh')
cfg$annotationFilter		= paste0(cfg$program, '/annotationFilter.pl')
cfg$FACETS_SUITE		= paste0(cfg$program, '/facets-suite/facets-suite-1.5.5_MW_MOD_huw/')
cfg$RLIB_PATH= paste0(cfg$program, '/facets_lib/facets-0.5.6_huw/')

cfg$enhancerMm9Bed          	= '/home/program/PLtest/Enhancers_mm9_CreyghtonPNAS/PutativeEnhancers_mm9_5types_SortMerge300_noPromOL_GT50bp.txt'
cfg$exeDir			= '/juno/work/solitlab/huw/program/socciCode/Pipelines/CBE/Variant/variants_pipeline/'
cfg$targetsDir			= '/juno/work/solitlab/huw/program/socciCode/Pipelines/CBE/Variant/variants_pipeline/targets/'
cfg$facets			= '/juno/work/solitlab/huw/program/socciCode/Pipelines/CBE/Variant/variants_pipeline/facets'
cfg$vepDir			= '/juno/work/solitlab/huw/program/socciCode/Pipelines/CBE/Variant/variants_pipeline/vep'
cfg$dataDir			=  '/juno/work/solitlab/huw/program/socciCode/Pipelines/CBE/Variant/variants_pipeline/data'
cfg$RSEM			= '/opt/common/CentOS_6/rsem/RSEM-1.2.25/'
cfg$KALLISTO			= '/home/huw/local/bin/kallisto '
##cfg$STAR			= '/opt/common/CentOS_6/star/STAR-STAR_2.5.0a/bin/Linux_x86_64/STAR'
cfg$fg$singularity 		= '/opt/local/singularity/3.3.0/bin/singularity'
cfg$simg 			= '/juno/work/taylorlab/cmopipeline/singularity_images/'
cfg$simg.vcf2maf		= paste0(cfg$simg, '/cmopipeline-vcf2maf-1.6.17.img')
cfg$polysolver_bin= '/juno/res/taylorlab/bandlamc/neoantigens/pipelines/polysolver/shell_call_hla_type_singularity'
cfg$neoantigen = '/home/huw/program/neoantigen/neoantigen-dev/neoantigen.py'
cfg$neoantigen.config = '/home/huw/program/neoantigen/neoantigen-dev/neoantigen-luna.config'
cfg$oncokb= '~/program/oncokb-annotator/'
cfg$hg19.picard.list = '~/program/refgene_hg19.txt'
cfg$cmo='~/program/cmo/cmo'
cfg$ngsfilter='~/program/ngs-filters'
cfg$vcf2maf='~/program/vcf2maf/'
cfg$vcf2maf='/opt/common/CentOS_6-dev/vcf2maf/v1.6.17/'
cfg$vep.dir='/opt/common/CentOS_6-dev/vep'
cfg$vep.cache='/opt/common/CentOS_6-dev/vep/cache'
cfg$vep.ver='95'
cfg$omni			= '/juno/work/solitlab/huw/program/socciCode/Pipelines/CBE/Variant/variants_pipeline/data/b37/1000G_omni2.5.b37.vcf'
cfg$hapmap 			= '/juno/work/solitlab/huw/program/socciCode/Pipelines/CBE/Variant/variants_pipeline/data/b37/hapmap_3.3.b37.vcf'
cfg$dbsnp			= '/juno/work/solitlab/huw/program/socciCode/Pipelines/CBE/Variant/variants_pipeline/data/b37/dbsnp_138.b37.vcf'
cfg$G1000G			= '/juno/work/solitlab/huw/program/socciCode/Pipelines/CBE/Variant/variants_pipeline/data/b37/1000G_phase1.snps.high_confidence.b37.vcf'
cfg$cranger.ref    		= '/home/huw/program/cellranger/cellranger-3.0.2/refdata-cellranger-GRCh38-3.0.0/'
cfg$cranger			= '/home/huw/program/cellranger/cellranger-3.0.2/cellranger'
cfg$funcotatorData.somatic	= '/home/huw/program/funcotator/funcotator_dataSources.v1.6.20190124s'
cfg$bic				= '/home/huw/program/BIC-variants_pipeline'
cfg$program 			= '/home/huw/program'
cfg$pipeline			= '/home/huw/pipeline'
cfg$cibersort			= paste0(cfg$program, '/cibersort')
cfg$anaconda2 			= paste0(cfg$program, '/anaconda2')
cfg$anaconda3 			= paste0(cfg$program, '/anaconda3')
cfg$conda3bin 			= paste0('/home/huw/conda3/envs/r36/bin')
cfg$hsa.housekeeping.bed	= paste0(cfg$program, '/hg19_RefSeq.bed')
cfg$hsaMethylBaitIntervalOrig	= paste0(cfg$program, '/truseq-methyl-capture-epic-manifest-file.bed')
cfg$hsaMethylBaitInterval	= paste0(cfg$program, '/truseq-methyl-capture-epic-manifest-file-bait.interval')
cfg$hsaMethylTargetInterval	= paste0(cfg$program, '/truseq-methyl-capture-epic-manifest-file-plusminus5bp.interval')
cfg$hsaMethylAnno		= paste0(cfg$program, '/truseq-methyl-capture-epic-manifest-file.annotated.islands.txt')
cfg$hg19.cpgi			= paste0(cfg$program, '/hg19_cpg_island.bed.txt')
cfg$hg38.cpgi			= paste0(cfg$program, '/cpgIslandExt.hg38.bed')
cfg$lncGTFGRCh38		= paste0(cfg$program, '/lncipedia_5_2_hg38.gtf')
cfg$hg38.refseq			= paste0(cfg$program, '/hg38_RefSeq.bed')
cfg$deeptools			= paste0(cfg$program, '/anaconda3/bin/')
cfg$getGenomeBuild		= '/home/socci/Code/FillOut/FillOut/GenomeData/getGenomeBuildBAM.sh' 

cfg$fastqc 			= paste0(cfg$program, '/FastQC/fastqc')
#options(scipen=999)
#system(paste0('cp /home/huw/program/truseq-methyl-capture-epic-header.txt ', hsaMethylBaitInterval))
#tmp = fread(hsaMethylBaitIntervalOrig)
#tmp = tmp[order(V1),]
#library(tibble)
#tmp = add_column(tmp, string = "+", .after = 3)
#tmp[, V1:=sub("chr", "", V1)]
#fwrite(tmp, file=hsaMethylBaitInterval, quote=F, sep="\t", col.names=F, append=T)
#write.table(as.data.frame(tmp), file=hsaMethylBaitInterval, quote=F, sep="\t", col.names=F, append=T, row.names=F)
#tmp[, V2:=V2-5]
#tmp[, V3:=V3+5]
#head(tmp)
#tmp[, V1 := as.character(V1)]
#tmp[, V2 := as.character(V2)]
#tmp[, V3 := as.character(V3)]
#tmp[, V4:=paste0("chr", V1, "_", V2, "-", V3)]
#system(paste0('cp /home/huw/program/truseq-methyl-capture-epic-header.txt ', hsaMethylTargetInterval))
#write.table(as.data.frame(tmp), file=hsaMethylTargetInterval, quote=F, sep="\t", col.names=F, append=T, row.names=F)
#rm(tmp)

#cfg$bismarkGRCh37		= '/juno/work/solitlab/huw/study/db/hsa/grch37_bismark'
cfg$btw2db			= '/juno/depot/assemblies/H.sapiens/hg19/index/bowtie/2.2.4/hg19_bowtie2'
cfg$bismark			= paste0(cfg$program, '/Bismark-0.22.3')
cfg$methylTarget		= paste0(cfg$program, '/truseq-methyl-capture-epic-manifest-file.bed')
cfg$SOAPHLA			= paste0(cfg$program, '/SOAP-HLA')
cfg$tmpdir			= '/scratch/huw'
cfg$mutsig			= '/opt/common/CentOS_6-dev/mutsig/cv-1.4/'

cfg$belly2			= paste0(cfg$program, '/delly_v0.7.7_linux_x86_64bit')
cfg$bamReadCount		= paste0(cfg$program, 'am-readcount')
cfg$IEDB			= paste0(cfg$program, '/iedb')
cfg$IEDBMHC1			= paste0(cfg$program, '/iedb/mhc_i')
cfg$IEDBMHC2			= paste0(cfg$program, '/iedb/mhc_ii')
cfg$pvacseq			= paste0(cfg$program, '/anaconda3/bin/pvacseq')
cfg$polysolver			= paste0(cfg$program, '/polysolver')
cfg$ponfile			= paste0(cfg$program, '/pon56Sample.vcf.gz')
cfg$make.trinuc			= paste0(cfg$program, '/make_trinuc_maf.py')
cfg$ConvertQualityScore    	= paste0(cfg$program, '/ConvertQualityScore')
cfg$STAR			= paste0(cfg$program, '/STAR-STAR_2.5.0a/bin/Linux_x86_64_static/STAR')
cfg$STAR_DIR			= paste0(cfg$program, '/STAR-STAR_2.5.0a/bin/Linux_x86_64_static')
cfg$sambamba			= paste0(cfg$program, '/sambamba' )
cfg$java			= paste0(cfg$program, '/jdk7/bin/java')
cfg$picard			= paste0(cfg$program, '/picard2.20.jar')
cfg$NGSPLOTexe			= paste0(cfg$program, '/ngsplot/bin/ngs.plot.r')
cfg$NGSFilter			= paste0(cfg$program, '/ngs-filters/')
cfg$trim_galore			= paste0(cfg$program, '/trim_galore')
cfg$bowtie2			= paste0(cfg$program, '/bowtie2/bowtie2')
cfg$bowtie2dir                  = paste0(cfg$program, '/bowtie2')
cfg$homerFindMotif              = paste0(cfg$program, '/homer/bin/findMotifsGenome.pl')
cfg$homerAnnotatePeaks		= paste0(cfg$program, '/homer/bin/annotatePeaks.pl')
cfg$tssCounts			= paste0(cfg$program, '/tsscounts.pl')
cfg$bwa				= paste0(cfg$program, '/bwa-0.7.12/bwa ')
cfg$bedtools			= paste0(cfg$program, '/bedtools/bin')
cfg$bedtools2			= paste0(cfg$program, '/bedtools2/bin')
cfg$bam2bigwig			= paste0(cfg$program, '/bam2bigwig.sh')
cfg$annotationFilter		= paste0(cfg$program, '/annotationFilter.pl')
cfg$FACETS_SUITE		= paste0(cfg$program, '/facets-suite/facets-suite-1.5.5_MW_MOD_huw/')
cfg$RLIB_PATH= paste0(cfg$program, '/facets_lib/facets-0.5.6_huw/')

cfg$enhancerMm9Bed          	= '/home/program/PLtest/Enhancers_mm9_CreyghtonPNAS/PutativeEnhancers_mm9_5types_SortMerge300_noPromOL_GT50bp.txt'
cfg$exeDir			= '/juno/work/solitlab/huw/program/socciCode/Pipelines/CBE/Variant/variants_pipeline/'
cfg$targetsDir			= '/juno/work/solitlab/huw/program/socciCode/Pipelines/CBE/Variant/variants_pipeline/targets/'
cfg$facets			= '/juno/work/solitlab/huw/program/socciCode/Pipelines/CBE/Variant/variants_pipeline/facets'
cfg$vepDir			= '/juno/work/solitlab/huw/program/socciCode/Pipelines/CBE/Variant/variants_pipeline/vep'
cfg$dataDir			=  '/juno/work/solitlab/huw/program/socciCode/Pipelines/CBE/Variant/variants_pipeline/data'
cfg$RSEM			= '/opt/common/CentOS_6/rsem/RSEM-1.2.25/'
cfg$KALLISTO			= '/home/huw/local/bin/kallisto '
##cfg$STAR			= '/opt/common/CentOS_6/star/STAR-STAR_2.5.0a/bin/Linux_x86_64/STAR'
cfg$samtools			= '/opt/common/CentOS_6-dev/bin/current/samtools'
cfg$samtoolsdir 		= '/opt/common/CentOS_6-dev/samtools/samtools-1.3.1'
cfg$java			= '/opt/common/CentOS_6-dev/java/jdk1.8.0_31/bin/java'
cfg$htseqcountexe		= '/home/huw/local/bin/htseq-count'
cfg$QSUB			= '/common/sge/bin/lx24-amd64/qsub'
cfg$BSUB			= '/admin/lsfjuno/lsf/10.1/linux3.10-glibc2.17-x86_64/bin/bsub'
cfg$macs14                      = '/home/huw/local/bin/macs14'
cfg$macs2                  	= '/home/huw/local/bin/macs2'

cfg$enhancerMm10Bed		= '/juno/work/solitlab//huw/study/db/mmu/dna/Mus_musculus.GRCm38.79.enhancer.bed'
cfg$tssMm9Bed              	= '/home/huw/program/PLtest/refFlat_mm9_TSSpm2kbSORTED.txt'
cfg$tssMm10Bed			= '/juno/work/solitlab//huw/study/db/mmu/dna/Mus_musculus.GRCm38.79.chr.symbol.tss2kb.bed'
cfg$genomeMm9Bowtie2		= '/juno/work/solitlab//huw/study/db/mmu/mm9-bowtie2/mm9'
cfg$genomeMm9StarDir		= '/juno/work/solitlab//huw/study/db/mmu/mm9_star'
cfg$genomeMm9GTF		= '/juno/work/solitlab//huw/study/db/mmu/dna/Mus_musculus.NCBIM37.67.chr.gtf'
cfg$genomeMm10Bowtie2		= '/juno/work/solitlab//huw/study/db/mmu/mm10-bowtie2/mm10'
cfg$genomeMm10StarDir		= '/juno/work/solitlab//huw/study/db/mmu/mm10star'
cfg$genomeMm10GTF		= '/juno/work/solitlab//huw/study/db/mmu/dna/Mus_musculus.GRCm38.79.chr.gtf'
cfg$mm10chromsize		= '/juno/work/solitlab//huw/study/db/mmu/dna/Mus_musculus.GRCm38.79.chrom.size'

cfg$ABRA			= '/opt/common/CentOS_6/abra2/abra2-2.07'
cfg$BCFTOOLS			= '/opt/common/CentOS_6/bcftools/bcftools-1.2/bin'
cfg$BEDTOOLS			= '/opt/common/CentOS_6/bedtools/bedtools-2.22.0/bin'
cfg$BWA				= '/opt/common/CentOS_6/bwa/bwa-0.7.12/bwa'
cfg$CUTADAPT			= '/opt/common/CentOS_6/cutadapt/cutadapt-1.9.1/bin/cutadapt'
cfg$DELLY			= '/opt/common/CentOS_6/delly/delly_v0.6.1'
cfg$dRANGER			= '/opt/common/CentOS_6/dranger/dRanger_annotate_v2.0'
cfg$FACETS_LIB			= '/opt/common/CentOS_6-dev/facets_lib/facets-0.5.6_huw'
cfg$FACETS_LIB			= '/home/huw/program/facets_lib/facets-0.5.6_huw'
cfg$FIXMULTIINDEL		= '/opt/common/CentOS_6/FixMultiInDel/FixMultiInDel-2.0.1'
cfg$gatk413			= '/home/huw/program/gatk-4.1.3.0/gatk-package-4.1.3.0-local.jar'

cfg$gatk4			= '/home/huw/program/gatk-4.1.2.0/gatk'
cfg$gnomadfile.gz		= '/home/huw/program/af-only-gnomad.raw.sites.b37.vcf.gz'
cfg$gnomadfile			= '/home/huw/program/af-only-gnomad.raw.sites.b37.vcf'
cfg$gatk36 			= '/home/huw/program/GenomeAnalysisTK36.jar'


cfg$JAVA			= '/opt/common/CentOS_6/java/jdk1.8.0_31/bin'
cfg$JAVA7_MUTECT            	= '/opt/common/CentOS_6/java/jdk1.7.0_75/bin'
cfg$MCR				= '/opt/common/CentOS_6/matlab/v81'
cfg$MUTECT			= '/opt/common/CentOS_6/muTect/muTect-1.1.7'
cfg$PERL			= '/opt/common/CentOS_6/perl/perl-5.22.0/bin'
cfg$PERL			= '/home/huw/local/perl/perl-5.26.0/bin'
cfg$GATK			= paste0(cfg$program, '/GenomeAnalysisTK-3.8-0-ge9d806836')
cfg$GATK4			= paste0(cfg$program, '/gatk-4.0.2.1/gatk')
cfg$PICARD			= paste0(cfg$program, '/picard2.20.jar')
cfg$PYTHON			= '/opt/common/CentOS_6/python/python-2.7.8/bin/python'
cfg$R				= '/opt/common/CentOS_6/R/R-3.1.2/bin'
cfg$SAMTOOLS                	= '/opt/common/CentOS_6/samtools/samtools-1.2'
cfg$SCALPEL			= '/opt/common/CentOS_6/scalpel/scalpel-0.2.2'
cfg$SOMATICSNIPER		= '/opt/common/CentOS_6/somaticsniper/somatic-sniper-1.0.4'
cfg$STRELKA			= '/opt/common/CentOS_6/strelka/strelka_1.0.11'
cfg$TABIX			= '/opt/common/CentOS_6/samtools/samtools-1.2/htslib-1.2.1'
cfg$VARSCAN			= '/opt/common/CentOS_6/varscan/v2.3.7'
cfg$VCF2MAF			= paste0(cfg$program, '/vcf2maf-1.6.14')
cfg$VEP_plugins             	= paste0(cfg$program, '/VEP_plugins')
cfg$VCFTOOLS			= '/opt/common/CentOS_6-dev/vcftools/v0.1.14/bin'
cfg$VEP                     	= '/opt/common/CentOS_6-dev/vep/v86'
cfg$VEPv88                  	= '/opt/common/CentOS_6-dev/vep/v88/variant_effect_predictor.pl'
cfg$VIRMID			= '/opt/common/CentOS_6/virmid/Virmid-1.1.1'
cfg$WES_FILTER			= '/opt/common/CentOS_6/wes-filter/wes-filters-1.1.3'
cfg$CMOBIN			= '/opt/common/CentOS_6-dev/python/python-2.7.10/bin/'


cfg$B37_BWA_INDEX		= '/juno/depot/assemblies/H.sapiens/b37/index/bwa/0.7.12/b37.fasta'
cfg$B37_FASTA			= '/juno/depot/assemblies/H.sapiens/b37/b37.fasta'
cfg$B37_FAI			= '/juno/depot/assemblies/H.sapiens/b37/b37.fasta.fai'
cfg$B37_MM10_HYBRID_BWA_INDEX	= '/juno/depot/assemblies/hybrids/H.sapiens_M.musculus/b37_mm10/index/bwa/0.7.12/b37_mm10.fasta'
cfg$B37_MM10_HYBRID_FASTA	= '/juno/depot/assemblies/hybrids/H.sapiens_M.musculus/b37_mm10/b37_mm10.fasta'
cfg$B37_MM10_HYBRID_FAI		= '/juno/depot/assemblies/hybrids/H.sapiens_M.musculus/b37_mm10/b37_mm10.fasta.fai'
cfg$HG19_BWA_INDEX		= '/juno/depot/assemblies/H.sapiens/hg19/index/bwa/0.7.4-r385/hg19.fasta'
cfg$HG19_FASTA			= '/juno/depot/assemblies/H.sapiens/hg19/hg19.fasta'
cfg$HG19_FAI			= '/juno/depot/assemblies/H.sapiens/hg19/hg19.fasta.fai'
cfg$HG19_MM10_HYBRID_BWA_INDEX	= '/juno/depot/assemblies/hybrids/H.sapiens_M.musculus/hg19_mm10/index/bwa/0.7.8/hg19_mm10.fasta'
cfg$HG19_MM10_HYBRID_FASTA	= '/juno/depot/assemblies/hybrids/H.sapiens_M.musculus/hg19_mm10/hg19_mm10.fasta'
cfg$HG19_MM10_HYBRID_FAI	= '/juno/depot/assemblies/hybrids/H.sapiens_M.musculus/hg19_mm10/hg19_mm10.fasta.fai'
cfg$MM9_BWA_INDEX		= '/juno/depot/assemblies/M.musculus/mm9/index/bwa/0.7.4-r385/mm9.fasta'
cfg$MM9_FASTA			= '/juno/depot/assemblies/M.musculus/mm9/mm9.fasta'
cfg$MM9_FAI			= '/juno/depot/assemblies/M.musculus/mm9/mm9.fasta.fai'
cfg$MM10_BWA_INDEX		= '/juno/depot/assemblies/M.musculus/mm10/index/bwa/0.7.8/mm10.fasta'
cfg$MM10_FASTA			= '/juno/depot/assemblies/M.musculus/mm10/mm10.fasta'
cfg$MM10_FAI			= '/juno/depot/assemblies/M.musculus/mm10/mm10.fasta.fai'
cfg$clipR1			= 'AGATCGGAAGAGCACACGTCT'
cfg$clipR2			= 'AGATCGGAAGAGCACACGTCT'
cfg$bqtrim			= '3'

## for homer annotatePeaks.pl mm10
cfg$homerMm10AnnotationFile	= '/juno/work/solitlab//huw/study/db/mmu/dna/Mus_musculus.GRCm38.79.chr.gtf.annotations.final.txt'

cfg$genomeHg19Bwa		= '/juno/work/solitlab//huw/study/db/hsa/hg37-bwa/hg37 '
cfg$genomeHg19StarDir		= '/juno/work/solitlab//huw/study/db/hsa/hg37star'
cfg$genomeHg19Bowtie2		= '/juno/work/solitlab//huw/study/db/hsa/hg37bowtie2/hg37'
cfg$genomeHg38Bowtie2		= '/juno/work/solitlab//huw/study/db/hsa/hg38_bowtie2/hg38'
cfg$genomeHg38Fai   		= "/juno/work/solitlab/huw/study/db/hsa/hg38_bowtie2/hg38.fai"
cfg$genomeHg19GTF		= '/juno/work/solitlab//huw/study/db/hsa/dna/Homo_sapiens.GRCh37.75.gtf'
cfg$genomeHg19GTFv2		= '/juno/work/solitlab//huw/study/db/hsa/dna/Homo_sapiens.GRCh37.75.noSelencysteine.gtf'
cfg$genomeHg19chromsize 	= '/juno/work/solitlab//huw/study/db/hsa/dna/hg19.chrom.sizes'
cfg$genomeHg19Fasta		= '/juno/work/solitlab/huw/study/db/hsa/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa'


cfg$genomeGRCh37Fasta		= '/juno/work/solitlab/huw/study/db/hsa/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa'
cfg$genomeGRCh37StarDir		= '/juno/depot/assemblies/H.sapiens/b37/index/star/2.4.1d/gencode/v18/overhang49 /ifs/work/solitlab/huw/park/study/hiseq/chipseq/template.conf'
cfg$genomeGRCh37chromsize 	= '/juno/work/solitlab//huw/study/db/hsa/dna/Homo_sapiens.GRCh37.75.chrom.size'

cfg$genomeGRCh38GTF		= '/juno/work/solitlab/huw/study/db/hsa/dna/Homo_sapiens.GRCh38.88.gtf'
cfg$genomeGRCh38Fasta		= '/juno/work/solitlab/huw/study/db/hsa/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
cfg$genomeGRCh38StarDir		= '/juno/work/solitlab/huw/study/db/hsa/hg38_star'
cfg$genomeGRCh38Rsem		= '/juno/work/solitlab/huw/program/RSEM_tutorial/ref/GRCh38/GRCh38'
cfg$genomeGRCh38chromsize	= '/juno/work/solitlab/huw/study/db/hsa/dna/Homo_sapiens.GRCh38.dna.primary_assembly_chrom_size.txt'
cfg$RSEM_REF			= '/juno/work/solitlab/huw/program/RSEM_tutorial/ref/GRCh38/GRCh38'
cfg$RSEM_REF_MMU		= '/juno/work/solitlab/huw/program/RSEM_tutorial/ref/GRCm38/GRCm38'
cfg$lnc_RSEM_REF		= '/juno/work/solitlab/huw/program/RSEM_tutorial/ref/LncGRCh38/lncGRCh38'

cfg$genomeGRCh38Kallisto	= '/juno/work/solitlab/huw/study/db/hsa/hg38_kallisto/kallisto_GRCh38'

cfg$genomeGRCm38Fasta		= '/juno/work/solitlab/huw/study/db/mmu/dna/Mus_musculus.GRCm38.79.dna.primary_assembly.chr.fa'
cfg$genomeGRCm38GTF		= '/juno/work/solitlab/huw/study/db/mmu/dna/Mus_musculus.GRCm38.79.chr.gtf'
cfg$genomeGRCm38Rsem		= '/juno/work/solitlab/huw/program/RSEM_tutorial/ref/GRCm38/GRCm38'

cfg$AgilentExon51Mbv3BaitsIlist	= '/juno/depot/resources/targets/AgilentExon_51MB_b37_v3/AgilentExon_51MB_b37_v3_baits.ilist'
cfg$AgilentExon51Mbv3tgtIlist	= '/juno/depot/resources/targets/AgilentExon_51MB_b37_v3/AgilentExon_51MB_b37_v3_targets.ilist'
cfg$AgilentExon51Mbv3tgt	= '/juno/depot/resources/targets/AgilentExon_51MB_b37_v3/AgilentExon_51MB_b37_v3_targets.bed'

