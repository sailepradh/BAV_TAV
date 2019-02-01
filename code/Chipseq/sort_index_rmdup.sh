#! /bin/bash -l

#SBATCH -A snic2016-7-108
#SBATCH -t 9:00:00
#SBATCH -p core -n 20
#SBATCH -J rmdup
#SBATCH -o rmdup.out
#SBATCH -e rmdup.err
#SBATCH --mail-user sailendra.pradhananga@scilifelab.se
#SBATCH --mail-type=ALL
#SBATCH --tmp=20480

module load bioinfo-tools
module load picard
module load java/sun_jdk1.8.0_40
module load samtools


input=/crex/proj/sllstore2017025/nobackup/private/BAV/Chipseq/data/BAM/
cd ${input}

wait
# INPUT BAM: filename (no extension e.g sample NOT sample.bam)
samtools sort -@ 10 $1".bam" -o $1".sort.bam"

wait

samtools index $1".sort.bam"

wait

#input filename - ".bam"
java -Xmx124G -jar $PICARD_HOME/picard.jar MarkDuplicates \
       I=$1".sort.bam" \
       O=$1".sort.dup.bam" \
       M=$1".sort.markduplicates.txt" VALIDATION_STRINGENCY=LENIENT \
       REMOVE_DUPLICATES=false ASSUME_SORTED=true

wait

module load samtools/1.6

samtools view -@ 5 -h -b -F 1804 -o $1".sort.rmdup.bam" $1".sort.dup.bam"
