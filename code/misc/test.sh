#!/bin/bash -l


echo "infile: "$1
orig_name=$1
substring=${orig_name:0:7}
echo "oufile: "$substring"_snp_indel.recalibrated.vcf"
