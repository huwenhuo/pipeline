index:
	#grep '>' ~/program/cellranger/cellranger-3.0.2/refdata-cellranger-GRCh38Crop-3.0.0/fasta/genome.fa |cut -d ' ' -f 1 > /juno/work/solitlab/huw/program/salmon/decoys.txt
	#sed -i.bak -e 's/>//g' /juno/work/solitlab/huw/program/salmon/decoys.txt
	cat /home/huw/program/cellranger/cellranger-3.0.2/cropseq /juno/work/solitlab/huw/program/salmon/gencode.v35.transcripts.fa  >> /juno/work/solitlab/huw/program/salmon/gencode.v35.crop.fa
	salmon index -t /juno/work/solitlab/huw/program/salmon/gencode.v35.crop.fa -p 12 -i salmon_index --gencode #-d /juno/work/solitlab/huw/program/salmon/decoys.txt 
txmap:
	grep '>' /juno/work/solitlab/huw/program/salmon/gencode.v35.crop.fa | sed "s/|/\t/g"  |  sed "s/>//" | awk 'BEGIN{OFS="\t"}{print $$1, $$2}'> /juno/work/solitlab/huw/program/salmon/txp2gene.tsv
alevin:
	salmon alevin -l ISR -1 /juno/work/solitlab/huw/solit/study/hiseq/10xgenomic/Project_11218/JAX_0480/Sample_639V_CROP_IGO_11218_1/639V_CROP_IGO_11218_1_S66_R1_001.fastq.gz /juno/work/solitlab/huw/solit/study/hiseq/10xgenomic/Project_11218/PITT_0509/Sample_639V_CROP_IGO_11218_1/639V_CROP_IGO_11218_1_S30_R1_001.fastq.gz -2 /juno/work/solitlab/huw/solit/study/hiseq/10xgenomic/Project_11218/JAX_0480/Sample_639V_CROP_IGO_11218_1/639V_CROP_IGO_11218_1_S66_R2_001.fastq.gz /juno/work/solitlab/huw/solit/study/hiseq/10xgenomic/Project_11218/PITT_0509/Sample_639V_CROP_IGO_11218_1/639V_CROP_IGO_11218_1_S30_R2_001.fastq.gz --chromiumV3 -i /juno/work/solitlab/huw/program/salmon/salmon_index -p 8 -o /juno/work/solitlab/huw/solit/study/hiseq/10xgenomic/cranger/Crop639V_Alevin --tgMap  /juno/work/solitlab/huw/program/salmon/txp2gene.tsv

	#bash /home/huw/program/SalmonTools/scripts/generateDecoyTranscriptome.sh -g ~/program/cellranger/cellranger-3.0.2/refdata-cellranger-GRCh38Crop-3.0.0/fasta/genome.fa -a ~/program/cellranger/cellranger-3.0.2/refdata-cellranger-GRCh38Crop-3.0.0/genes/genes.gtf -t /juno/work/solitlab/huw/program/gencode.v35.transcripts.fa -o /juno/work/solitlab/huw/program/salmon/
