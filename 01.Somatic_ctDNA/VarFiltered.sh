#!/bin/sh
set -eo pipefail

sample=$1
Out=$2

path_out="$Out/01.TNscope"
cd $path_out

ref=~/DataBase/hg19/hg19.fa

JAVA_MEM=-Xmx1G
gatk --java-options ${JAVA_MEM} SelectVariants \
 -R ${ref} \
 -V ${sample}_tnscope.filter.vcf.gz \
 -select 'vc.isNotFiltered()' \
 -O ${sample}_tnscope.filtered.vcf.gz