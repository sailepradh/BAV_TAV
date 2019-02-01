#! /bin/bash -l

#SBATCH -A snic2016-7-108
#SBATCH -t 1:00:00
#SBATCH -p core -n 5
#SBATCH -J bamcoverage
#SBATCH -o bamcoverage.out
#SBATCH -e bamcoverage.err
#SBATCH --mail-user pontus.hojer@scilifelab.se
#SBATCH --mail-type=FAIL
#SBATCH --tmp=20480

module load bioinfo-tools
module load deepTools

wait 

echo "bamCoverage for file: " $1
# Effective genome size from http://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
# Command from: https://deeptools.readthedocs.io/en/develop/content/tools/bamCoverage.html 

bamCoverage -b $1".sort.rmdup.bam" -o $1".sort.rmdup.seqdepthnorm_mapq20.bw" \
	--normalizeTo1x 2736124973 \
	--minMappingQuality 20 \
	--binSize 10 \
	--extendReads 200 \
	-p 1
	

echo "Done!"

