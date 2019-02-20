#!/bin/bash -l
#SBATCH -A snic2016-7-108
#SBATCH -t 2:00:00
#SBATCH -p core -n 8
#SBATCH -o ReCalibrateSNP.out
#SBATCH -e ReCalibrateSNP.err
#SBATCH -J ReCalibrate_snps_all_TAVs
#SBATCH --mail-user sailendra.pradhananga@scilifelab.se
#SBATCH --mail-type=ALL

module load bioinfo-tools
module load samtools/1.3
module load picard
module load GATK/3.8-0
module load vcftools
module load vcflib
module load tabix

echo "infile: "$1
orig_name=$1
substring=${orig_name:0:7}
echo "oufile: "$substring"_snp_indel.recalibrated.vcf"


java -Xmx60g -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /crex/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta \
--filterExpression  "QUAL < 50.0" --filterName "LowQual" --filterExpression "QD<1.0" --filterName "QualByDepth" \
--variant $1 \
--out filtered.vcf

java -Xmx60g -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar  \
-T VariantRecalibrator \
-R  /crex/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta \
-input filtered.vcf \
-resource:hapmap,known=false,training=true,truth=true,prior=15.0 \
/crex/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/hg19/hapmap_3.3.hg19.vcf    -resource:omni,known=false,training=true,truth=false,prior=12.0  \
/crex/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/hg19/1000G_omni2.5.hg19.vcf -resource:1000G,known=false,training=true,truth=false,prior=10.0 \
/crex/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/hg19/1000G_phase1.snps.high_confidence.hg19.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \
/crex/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/hg19/dbsnp_138.hg19.vcf \
-minNumBad 1000 -recalFile recalibrate_SNP.recal -tranchesFile recalibrate_SNP.tranches \
-rscriptFile recalibrate_SNP_plots.R \
-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
-mode SNP

java -Xmx60g -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R  /crex/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta \
-input filtered.vcf -recalFile recalibrate_SNP.recal \
-tranchesFile recalibrate_SNP.tranches \
-o recalibrated_snps_raw_indels.vcf  \
--ts_filter_level 99.0 \
-mode SNP

java -Xmx24g -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar \
-T VariantRecalibrator \
-R /crex/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta \
-input recalibrated_snps_raw_indels.vcf -recalFile indels.raw.recal  \
-tranchesFile indels.raw.tranches \
-resource:mills,known=true,training=true,truth=true,prior=12.0 /crex/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf \
-an DP -an FS -an MQRankSum -an ReadPosRankSum  \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0  \
--maxGaussians 4 -rscriptFile recalibrate_INDEL_plots.R \
-mode INDEL

java -Xmx72g -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar \
-T  ApplyRecalibration \
-R  /crex/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta \
-input recalibrated_snps_raw_indels.vcf  \
-tranchesFile indels.raw.tranches -recalFile indels.raw.recal \
-o $substring"_snp_indel.recalibrated.vcf"  \
--ts_filter_level 99.0 \
-mode INDEL
