#!/bin/bash -l

#SBATCH -A snic2016-7-108
#SBATCH -n 16 -p node
#SBATCH -t 12:00:00
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

#picard=/sw/apps/bioinfo/picard/2.10.3/rackham/picard.jar
#GATK=/sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar
#GATKbundle=/crex/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/hg19/
#input=/crex/proj/sllstore2017025/nobackup/private/Comparison_BAV_TAV/SNP_calling/


echo "$(date) Running on: $(hostname)"
echo "HaplotypeCaller"

java -Xmx50g -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-nct 16  \
-R /crex/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta \
-I BAV2714EC.trunc.dedup.sortedbycoord.withreadgroups.recal.bam \
-ERC GVCF \
-G Standard \
-G AS_Standard \
--dbsnp /crex/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/hg19/dbsnp_138.hg19.vcf \
-variant_index_type LINEAR -variant_index_parameter 128000 \
-o BAV2714.gvcf

wait

echo "Done. $(date) Running on: $(hostname)"
wait
echo "$(date) AllDone"
echo "$(date) AllDone" >&2
