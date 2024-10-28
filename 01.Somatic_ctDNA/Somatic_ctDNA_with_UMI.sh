#!/bin/sh

# Copyright (c) 2016-2020 Sentieon Inc. All rights reserved

# *******************************************
# Script to perform TNscope variant calling on
# UMI-tagged samples
# *******************************************

set -eo pipefail

sample=$1
fq1=$2
fq2=$3
Out=$4/${sample}

check() { $1 ; if [ $? -ne 0 ]; then exit -1; fi }

path_out="$Out/00.QC"
mkdir -p $path_out && cd $path_out

fq1_clean="$path_out/${sample}.clean.1.fq.gz"
fq2_clean="$path_out/${sample}.clean.2.fq.gz"
json_p="$path_out/${sample}.fastp_p.json"
html_p="$path_out/${sample}.fastp_p.html"

run_QC()
{
fastp -i $fq1 -I $fq2 \
  --out1 $fq1_clean --out2 $fq2_clean \
  -w 10 --detect_adapter_for_pe -j $json_p -h $html_p
}

path_out="$Out/01.TNscope"
mkdir -p $path_out && cd $path_out

if [[ -f ../00.QC/${sample}.clean.1.fq.gz ]];then
        sleep 1
else
        check run_QC
fi

# Update with the fullpath location of your sample fastq
TUMOR_SM=${sample} #sample name
TUMOR_RGID="rg_$TUMOR_SM" #read group ID
PL="ILLUMINA" #or other sequencing platform
TUMOR_FASTQ_1=${fq1_clean}
TUMOR_FASTQ_2=${fq2_clean} #If using Illumina paired data

# Update with the location of the reference data files
FASTA_DIR="~/DataBase/hg19"
FASTA="$FASTA_DIR/hg19"
KNOWN_DBSNP="~/DataBase/gatk-bundles/dbsnp_138.hg19.vcf"
INTERVAL_FILE="~/DataBase/1345-fusion.sorted.bed"

# Update with the location of the Sentieon software package and license file
SENTIEON_INSTALL_DIR=~/Soft/Sentieon
TNSCOPE_FILTER=~/Soft/sentieon-scripts/tnscope_filter/tnscope_filter.py
export SENTIEON_LICENSE=XXX.XXX.XXX.XXX:8000 #or using licsrvr: c1n11.sentieon.com:5443

#UMI information
READ_STRUCTURE="1S3M3S+T,1S3M3S+T" #an example duplex: "3M2S+T,3M2S+T" where duplex UMI extraction requires an identical read structure for both strands
DUPLEX_UMI="false" #set to "true" if duplex

# Other settings
NT=4 #number of threads to use in computation, set to number of cores in the server
START_DIR="$PWD" #Determine where the output files will be stored

# ******************************************
# 0. Setup
# ******************************************
WORKDIR="$START_DIR"
mkdir -p $WORKDIR
LOGFILE=$WORKDIR/run.log
# exec >$LOGFILE 2>&1
cd $WORKDIR

# ******************************************
# 1. Pre-processing of FASTQ containing UMIs
# ******************************************
if [ "$DUPLEX_UMI" = "true" ] ; then
    READ_STRUCTURE="-d $READ_STRUCTURE"
fi

( $SENTIEON_INSTALL_DIR/bin/sentieon umi extract $READ_STRUCTURE $TUMOR_FASTQ_1 $TUMOR_FASTQ_2 || \
    { echo -n 'Extract error' >&2; exit -1; } ) | \
( $SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -p -C \
  -R "@RG\tID:$TUMOR_RGID\tSM:$TUMOR_SM\tPL:$PL" -t $NT \
  -K 10000000 $FASTA - || { echo -n 'BWA error'; exit 1; } ) | \
  $SENTIEON_INSTALL_DIR/bin/sentieon umi consensus -t $NT -o umi_consensus.fastq.gz || \
  { echo "Alignment/Consensus failed"; exit 1; }

( $SENTIEON_INSTALL_DIR/bin/sentieon bwa mem -p -C \
    -R "@RG\tID:$TUMOR_RGID\tSM:$TUMOR_SM\tPL:$PL" -t $NT -K 10000000 \
    $FASTA umi_consensus.fastq.gz || { echo -n 'BWA error'; exit 1; } ) | \
    $SENTIEON_INSTALL_DIR/bin/sentieon util sort -t $NT --umi_post_process --sam2bam -i - \
    -o umi_consensus.bam || { echo "Consensus alignment failed"; exit 1; }

# ******************************************
# 2. Somatic and Structural variant calling
# ******************************************
# Consider adding `--disable_detector sv --trim_soft_clip` if not interested in SV calling
FASTA=${FASTA}.fa
$SENTIEON_INSTALL_DIR/bin/sentieon driver -r $FASTA -t $NT -i umi_consensus.bam \
    ${INTERVAL_FILE:+--interval $INTERVAL_FILE} \
    --algo TNscope \
    --tumor_sample $TUMOR_SM \
    --dbsnp $KNOWN_DBSNP \
    --min_tumor_allele_frac 0.001 \
    --min_tumor_lod 3.0 \
    --min_init_tumor_lod 3.0 \
    --pcr_indel_model NONE \
    --min_base_qual 40 \
    --resample_depth 100000 \
    --assemble_mode 4 \
    ${sample}_tnscope.pre_filter.vcf.gz || \
    { echo "TNscope failed"; exit 1; }

#### The output of TNscope requires filtering to remove false positives.
#### Filter design depends on the specific sample and user needs to modify the following accordingly.
$SENTIEON_INSTALL_DIR/bin/sentieon pyexec $TNSCOPE_FILTER \
    -v ${sample}_tnscope.pre_filter.vcf.gz \
    --tumor_sample $TUMOR_SM \
    -x ctdna_umi --min_tumor_af 0.001 --min_depth 100 \
    ${sample}_tnscope.filter.vcf.gz || \
    { echo "TNscope filter failed"; exit 1; }