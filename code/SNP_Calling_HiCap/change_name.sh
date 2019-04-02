#!/bin/bash -l

#SBATCH -A snic2016-7-108
#SBATCH -p devel -n 8
#SBATCH -t 1:00:00
#SBATCH -J BAV_TAV.genotypeGVCFs
#SBATCH -e BAV_TAV.genotypeGVCFs.err
#SBATCH -o BAV_TAV.genotypeGVCFs.out
#SBATCH --mail-type=All
#SBATCH --mail-user=sailendra.pradhananga@scilifelab.se

echo "$(date) Running on: $(hostname)"
echo "$(date) Running on: $(hostname)" >&2



module load bioinfo-tools
module load samtools/1.3
module load picard
module load GATK/3.8-0
module load vcftools
module load vcflib
module load tabix

picard=/sw/apps/bioinfo/picard/2.10.3/rackham/picard.jar
GATK=/sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar
GATKbundle=/crex/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/hg19/

time awk -v OFS="\t" '{if ($1 == "#CHROMPOS") print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","TAV2431"; else print $0}' TAV2431.gvcf > tmp ; mv tmp TAV2431.gvcf
time awk -v OFS="\t" '{if ($1 == "#CHROMPOS") print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","TAV2709"; else print $0}' TAV2709.gvcf > tmp ; mv tmp TAV2709.gvcf
time awk -v OFS="\t" '{if ($1 == "#CHROMPOS") print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","TAV2515"; else print $0}' TAV2515.gvcf > tmp ; mv tmp TAV2515.gvcf
time awk -v OFS="\t" '{if ($1 == "#CHROMPOS") print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","BAV2424"; else print $0}' BAV2424.gvcf > tmp ; mv tmp BAV2424.gvcf
time awk -v OFS="\t" '{if ($1 == "#CHROMPOS") print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","BAV2375"; else print $0}' BAV2375.gvcf > tmp ; mv tmp BAV2375.gvcf
time awk -v OFS="\t" '{if ($1 == "#CHROMPOS") print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","TAV2515"; else print $0}' TAV2515.gvcf > tmp ; mv tmp TAV2515.gvcf



echo "$(date) Running on: $(hostname)"
echo "GenotypeGVCFs"


wait

echo "Done. $(date) Running on: $(hostname)"
wait
echo "$(date) AllDone"
echo "$(date) AllDone" >&2