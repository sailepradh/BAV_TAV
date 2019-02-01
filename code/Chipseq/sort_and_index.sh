#! /bin/bash -l

#SBATCH -A snic2016-7-108
#SBATCH -t 3:00:00
#SBATCH -p core -n 10
#SBATCH -J sort_index
#SBATCH -o sort_index.out
#SBATCH -e sort_index.err
#SBATCH --mail-user pontus.hojer@scilifelab.se
#SBATCH --mail-type=ALL
#SBATCH --tmp=20480

module load bioinfo-tools
module load samtools

wait 
# INPUT BAM: filename (no extension e.g sample NOT sample.bam)
samtools sort -@ 10 $1".bam" -o $1".sort.bam"

wait

samtools index $1".sort.bam"

