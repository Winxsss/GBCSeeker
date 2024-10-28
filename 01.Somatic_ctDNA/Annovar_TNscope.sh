#!/bin/sh
set -eo pipefail

sample=$1
Out=$2

path_out="$Out/02.Anno_TNscope"
mkdir -p $path_out && cd $path_out

convert2annovar.pl -format vcf4 -allsample -withfreq ../01.TNscope/${sample}_tnscope.filtered.vcf.gz > ${sample}.avinput

table_annovar.pl ${sample}.avinput ~/DataBase/humandb/ -buildver hg19 -out ${sample}.anno \
    -remove -nastring NA \
    --thread 10 \
    -protocol "refGene,cytoBand,knownGene,ensGene,avsift,dbnsfp42a,ljb26_all,exac03,clinvar_20220320,1000g2015aug_EAS" \
    -operation "g,r,g,g,f,f,f,f,f,f"

table_annovar.pl ../01.TNscope/${sample}_tnscope.filtered.vcf.gz --vcfinput ~/DataBase/humandb/ -buildver hg19 -out ${sample}.anno \
    -remove -nastring . \
    --thread 10 \
    -protocol "refGene,cytoBand,knownGene,ensGene,avsift,dbnsfp42a,ljb26_all,exac03,clinvar_20220320,1000g2015aug_EAS" \
    -operation "g,r,g,g,f,f,f,f,f,f"
