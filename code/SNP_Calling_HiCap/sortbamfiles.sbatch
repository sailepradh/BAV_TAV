#! /bin/bash -l

#SBATCH -A a2014205
#SBATCH -J mapBowtie2.job
#SBATCH -p node
#SBATCH -t 1-00:00:00
#SBATCH --mail-user pelin.akan@scilifelab.se
#SBATCH --mail-type=ALL

base_file_name=$1
sample_name=$2
echo sample_name $sample_name
echo base_file_name $base_file_name

#SBATCH -o $1.out
#SBATCH -e $1.err

module load bioinfo-tools
module load bowtie2/2.2.6
module load samtools

#Get current directory
curdir=`readlink -f .`
echo "Current directory is " $curdir
#mkdir -p $sample_name
#echo "Created a directory within current directory " $sample_name

#Copy the file needs to be sorted to scratch
cp $sample_name/$base_file_name\.trunc.dedup.sortedbyreadname.bam $SNIC_TMP/

#To load to IGV, sort and index BAM
samtools sort $SNIC_TMP/$base_file_name\.trunc.dedup.sortedbyreadname.bam $SNIC_TMP/$base_file_name\.trunc.dedup.sortedbycoord
samtools index $SNIC_TMP/$base_file_name\.trunc.dedup.sortedbycoord.bam $SNIC_TMP/$base_file_name\.trunc.dedup.sortedbycoord.bai

rm $SNIC_TMP/*.sortedbyreadname.bam

echo -ne '\n\nCopying back the resulting analysis files at '
date
rsync -rptoDLv $SNIC_TMP// $curdir/$sample_name//


#awk -v OFS="\t" '{if (NR==194) print $1,$2,$3,$4,$5,$6,$7,$8,$9,"TAV2515"; else print $0}' test.vcf > tmp ; mv tmp test.vcf
#java -Xmx50g -jar /sw/apps/bioinfo/GATK/3.8-0/GenomeAnalysisTK.jar -T GenotypeGVCFs  -R /crex/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/hg19/ucsc.hg19.fasta --variant test.vcf --variant test2.vcf --dbsnp  /crex/data/uppnex/reference/biodata/GATK/ftp.broadinstitute.org/bundle/2.8/hg19/dbsnp_138.hg19.vcf -G Standard -G AS_Standard -o raw_variants.vcf.gz
