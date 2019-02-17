#!/bin/bash -l

#SBATCH -A snic2016-7-108
#SBATCH -n 8 -p node
#SBATCH -t 10:00:00
#SBATCH -J BAV_TAV.genotypeGVCFs
#SBATCH -e BAV_TAV.genotypeGVCFs.err
#SBATCH -o BAV_TAV.genotypeGVCFs.out
#SBATCH --mail-type=All
#SBATCH --mail-user=sailendra.pradhananga@scilifelab.se

echo "$(date) Running on: $(hostname)"
echo "$(date) Running on: $(hostname)" >&2

"ava -Xmx4g -jar  /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar -T CombineVariants -R /crex/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta --variant TAV2709.vcf  --variant TAV2431.vcf --variant TAV2515.vcf  --variant  BAV2714.vcf --variant BAV2424.vcf  --variant BAV2375.vcf -o combined.vcf -genotypeMergeOptions UNIQUIFY"

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

java -Xmx50g -jar $GATK -nt 8 -T GenotypeGVCFs \
-R $GATKbundle/ucsc.hg19.fasta \
--variant $input/TAV2431.vcf \
--variant $input/TAV2515.vcf \
--variant $input/TAV2709.vcf \
--variant $input/BAV2375.vcf \
--variant $input/BAV2424.vcf \
--variant $input/BAV2714.vcf \
--dbsnp $GATKbundle/dbsnp_138.hg19.vcf \
-A StrandBiasBySample \
-A FisherStrand \
-A MappingQualityZero \
-A MappingQualityZeroBySample \
-A QualByDepth \
-o $input/raw_variants.vcf.gz

wait

echo "Done. $(date) Running on: $(hostname)"
wait
echo "$(date) AllDone"
echo "$(date) AllDone" >&2
