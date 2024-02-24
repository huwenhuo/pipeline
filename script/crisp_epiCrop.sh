#/opt/local/singularity/3.3.0/bin/singularity pull docker//boutroslab/cld_docker
#folder for CLD CRISPR Library Designer
#cd /juno/work/solitlab/huw/program/cld/
#/admin/lsfjuno/lsf/10.1/linux3.10-glibc2.17-x86_64/bin/bsub -J cld.makedatabase -W 99:99 -R "mem>30" -M 15 -n 5 "/opt/local/singularity/3.3.0/bin/singularity exec ~/program/img/cld_docker_latest.sif cld --task=make_database --organism=homo_sapiens"
#/opt/local/singularity/3.3.0/bin/singularity exec ~/program/img/cld_docker_latest.sif cld --task=make_database --organism=homo_sapiens
#cd /juno/work/solitlab/huw/solit/study/hiseq/epiCrop
#/admin/lsfjuno/lsf/10.1/linux3.10-glibc2.17-x86_64/bin/bsub -J cld.makedatabase -W 99:99 -R  "mem>30" -M 30 -n 20 "/opt/local/singularity/3.3.0/bin/singularity exec ~/program/img/cld_docker_latest.sif cld --task=target_ident --output-dir=/juno/work/solitlab/huw/solit/study/hiseq/epiCrop --parameter-file=/home/huw/program/cld/params.txt --gene-list=/juno/work/solitlab/huw/solit/study/hiseq/epiCrop/rs.hsa.engs"
#/opt/local/singularity/3.3.0/bin/singularity exec -B /juno/work/solitlab/huw/program/cld/homo_sapiens:/data/homo_sapiens ~/program/img/cld_docker_latest.sif cld --task=target_ident --output-dir=/juno/work/solitlab/huw/solit/study/hiseq/epiCrop --parameter-file=/juno/work/solitlab/huw/solit/study/hiseq/epiCrop/params.txt --gene-list=/juno/work/solitlab/huw/solit/study/hiseq/epiCrop/rs.hsa.engs
#/opt/local/singularity/3.3.0/bin/singularity exec ~/program/img/cld_docker_latest.sif cld --help


