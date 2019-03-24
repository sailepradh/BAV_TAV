#!/bin/bash -l


echo "infile: "$1
orig_name=$1
substring=${orig_name:0:7}
echo "oufile: "$substring"_snp_indel.recalibrated.vcf"
stringtie \
$i \
-p 20 \
-G gencode.v26lift37.annotation.gff3 \
-o ballgown/$sample_name/$sample_name.stringtie.gtf \
-B \
-e \
-A $sample_name.gene_abund.tab \
-x chrM,chrY
