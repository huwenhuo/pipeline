#!/bin/bash
#$ -cwd
#$ -l h_vmem=16G

########################################################## ######################
### COPYRIGHT ########################################################## ########
# New York Genome Center
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright (2014) by the New York
# Genome Center. All rights are reserved. This software is supplied without
# any warranty or guaranteed support whatsoever. The New York Genome Center 
# cannot be responsible for its use, misuse, or functionality.
# Version: 0.1
# Author: Avinash Abhyankar
########################################################## ####### /COPYRIGHT ###
########################################################## ######################

if [ $# -ne 2 ]
then   
        echo "Usage: <path_to_bam><sampleID>"
        exit 65
fi

bamdir=$(dirname ${1})
fastqdir=$(dirname ${bamdir})/fastq
mkdir -p ${fastqdir}/tmpdir
tmpdir=${fastqdir}/tmpdir

#my_java="~/local/bin/java -Djava.io.tmpdir=${tmpdir} -Xmx12g"
#picard="${my_java} -jar ~/program/picard-tools-1.128/picard.jar"
my_java="/opt/common/CentOS_6/java/jdk1.8.0_25/bin/java -Djava.io.tmpdir=${tmpdir} -Xmx80g"
picard="${my_java} -jar /home/huw/program/picard-tools-1.128/picard.jar"

## Convert BAM to FASTQ
${picard} SamToFastq \
VALIDATION_STRINGENCY=SILENT \
INPUT=${1} \
OUTPUT_PER_RG=true \
RG_TAG=ID \
OUTPUT_DIR=${fastqdir}
