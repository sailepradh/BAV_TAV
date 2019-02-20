#! /bin/bash -l
#SBATCH -A b2012058
#SBATCH -n 8 -p node
#SBATCH -t 20:00:00
#SBATCH -J BWA.VQSR.SNPs
#SBATCH -e /proj/bils2015003/nobackup/BWA_Exome/scr/stderr.vqsr.SNPs.txt
#SBATCH -o /proj/bils2015003/nobackup/BWA_Exome/scr/stdout.vqsr.SNPs.txt
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


echo "---> Variant filtration <---"
echo "---> Variant filtration <---" 1>&2

java -Xmx60g -jar $GATK -T \
VariantFiltration -R \
$reference --filterExpression  "QUAL < 50.0" --filterName "LowQual" \
--filterExpression "FS > 60" --filterName "FisherSB" --filterExpression "QD<1.0" --filterName "QualByDepth" \
--filterExpression "(MQ0 >= 4 && ((MQ0/(1.0 * DP)) > 0.1))" --filterName "FUBAR" \
--variant raw_variants.vcf \
--out filtered_variants.vcf

echo "Variant filtration done. $(date) Running on: $(hostname)"
echo "---> VQSR <---"
echo "---> VQSR <---" 1>&2

java -Xmx60g -jar $GATK -T VariantRecalibrator -R $reference -input filtered_variants.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $GATKbundle/hapmap_3.3.b37.vcf -resource:omni,known=false,training=true,truth=false,prior=12.0 $GATKbundle/1000G_omni2.5.b37.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=6.0 $GATKbundle/dbsnp_138.b37.vcf -an QD -an MQRankSum -an ReadPosRankSum -an FS -an DP -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -minNumBad 1000 -recalFile recalibrate_SNP.recal  -tranchesFile recalibrate_SNP.tranches  -rscriptFile recalibrate_SNP_plots.R

echo "Done. $(date) Running on: $(hostname)"
wait
echo "-----"

echo "$(date) Running on: $(hostname)"


echo "---> apply recal <---"
echo "---> apply recal <---" 1>&2

java -Xmx60g -jar $GATK -T ApplyRecalibration -R $reference  -input filtered_variants.vcf  -recalFile recalibrate_SNP.recal  -tranchesFile recalibrate_SNP.tranches  -o recalibrated_snps_raw_indels.vcf  --ts_filter_level 99.0 -mode SNP


echo "Done. $(date) Running on: $(hostname)"
wait
echo "$(date) AllDone"
echo "$(date) AllDone" >&2
