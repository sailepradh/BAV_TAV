#!/bin/bash -l

#SBATCH -A snic2016-7-108
#SBATCH -n 16 -p node
#SBATCH -t 16:00:00
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
input=/crex/proj/sllstore2017025/nobackup/private/Comparison_BAV_TAV/SNP_calling/


echo "$(date) Running on: $(hostname)"
echo "GenotypeGVCFs"

java -Xmx100g -jar $GATK -nt 16 -T GenotypeGVCFs  \
-R $GATKbundle/ucsc.hg19.fasta \
--variant /crex/proj/snic2016-7-108/sail/tmp_file_vcf_BAV_TAV/TAV2431.gvcf \
--variant /crex/proj/snic2016-7-108/sail/tmp_file_vcf_BAV_TAV/TAV2515.gvcf \
--variant TAV2709.gvcf \
--variant BAV2375.gvcf \
--variant BAV2424.gvcf \
--variant BAV2714.gvcf \
--dbsnp $GATKbundle/dbsnp_138.hg19.vcf \
-L /domus/h1/sail/SNP_calling_BAV/sorted_PD_uniq.bed \
-o $input/raw_variants.vcf.gz

wait

echo "Done. $(date) Running on: $(hostname)"
wait
echo "$(date) AllDone"
echo "$(date) AllDone" >&2
