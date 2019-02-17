#! /bin/bash -l
#SBATCH -A b2012058
#SBATCH -n 8 -p node
#SBATCH -t 72:00:00
#SBATCH -J BWA.VQSR.INDELs
#SBATCH -e /proj/bils2015003/nobackup/BWA_Exome/scr/stderr.vqsr.INDELs.txt
#SBATCH -o /proj/bils2015003/nobackup/BWA_Exome/scr/stdout.vqsr.INDELs.txt
#SBATCH --mail-type=All
#SBATCH --mail-user=sailendra.pradhananga@scilifelab.se
echo "$(date) Running on: $(hostname)"
echo "$(date) Running on: $(hostname)" >&2
cd /proj/bils2015003/nobackup/BWA_Exome/

reference=/proj/b2012058/private/benjamin/eb/GATKbundle/human_g1k_v37.fasta
picard=/sw/apps/bioinfo/picard/1.127/milou/picard.jar
GATK=/sw/apps/bioinfo/GATK/3.3.0/GenomeAnalysisTK.jar
GATKbundle=/proj/b2012058/private/benjamin/eb/GATKbundle/
echo "$(date) Running on: $(hostname)"
echo "VQSR" 
java -Xmx24g -jar $GATK \
   -T VariantRecalibrator \
   -nt 8 \
   -R $reference \
   -input recalibrated_snps_raw_indels.vcf \
   -recalFile indels.raw.recal \
   -tranchesFile indels.raw.tranches \
   -resource:mills,known=true,training=true,truth=true,prior=12.0 $GATKbundle/Mills_and_1000G_gold_standard.indels.b37.vcf \
   -an DP -an FS -an MQRankSum -an ReadPosRankSum \
   -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
   --maxGaussians 4 \
   -rscriptFile recalibrate_INDEL_plots.R \
   -mode INDEL
echo "Done. $(date) Running on: $(hostname)"
wait
echo "$(date) Running on: $(hostname)"
echo "---> apply recal <---"
java -Xmx72g -jar $GATK \
   -nt 16 \
   -T ApplyRecalibration \
   -R $reference \
   -input recalibrated_snps_raw_indels.vcf \
   -tranchesFile indels.raw.tranches \
   -recalFile indels.raw.recal \
   -o indels.recalibrated.vcf \
   --ts_filter_level 99.0 \
   -mode INDEL

echo "Done. $(date) Running on: $(hostname)"
wait
echo $file;
echo "$(date) AllDone"
echo "$(date) AllDone" >&2
