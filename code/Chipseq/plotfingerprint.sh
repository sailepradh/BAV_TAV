#! /bin/bash -l

#SBATCH -A snic2016-7-108
#SBATCH -t 3:00:00
#SBATCH -p core -n 10
#SBATCH -J plotfingerprint
#SBATCH -o plotfingerprint.out
#SBATCH -e plotfingerprint.err
#SBATCH --mail-user pontus.hojer@scilifelab.se
#SBATCH --mail-type=FAIL
#SBATCH --tmp=20480

module load bioinfo-tools
module load deepTools
cd /proj/snic2016-7-108/private/pontus_temp_fastq

wait 

echo "Plot all files"
plotFingerprint -b *.rmdup.bam \
	--labels Ker0r2Ctrl Ker0r2K27Ac Ker0r2ZNF365 Ker3r2Ctrl Ker3r2K27Ac Ker3r2ZNF365 Ker7r2Ctrl Ker7r2K27Ac Ker7r2ZNF365 \
	--minMappingQuality 20 \
	--skipZeros \
	-p 10 \
	-T Ker_all \
	-plot plotFingerprint_all.png \
	--extendReads 200

wait

echo "Plot d0"
plotFingerprint -b Ker0r2Ctrl_S3_R1_001.fastq.gz.sort.rmdup.bam Ker0r2K27Ac_S7_R1_001.fastq.gz.sort.rmdup.bam Ker0r2ZNF365_S12_R1_001.fastq.gz.sort.rmdup.bam \
	--labels Ctrl K27Ac ZNF365 --minMappingQuality 20 --skipZeros -p 10 -T Ker_d0 -plot plotFingerprint_d0.png --extendReads 200

wait

echo "Plot d3"
plotFingerprint -b Ker3r2Ctrl_S5_R1_001.fastq.gz.sort.rmdup.bam Ker3r2K27Ac_S6_R1_001.fastq.gz.sort.rmdup.bam Ker3r2ZNF365_S11_R1_001.fastq.gz.sort.rmdup.bam \
	--labels Ctrl K27Ac ZNF365 --minMappingQuality 20 --skipZeros -p 10 -T Ker_d3 -plot plotFingerprint_d3.png --extendReads 200

wait

echo "Plot d7"
plotFingerprint -b Ker7r2Ctrl_S2_R1_001.fastq.gz.sort.rmdup.bam Ker7r2K27Ac_S4_R1_001.fastq.gz.sort.rmdup.bam Ker7r2ZNF365_S9_R1_001.fastq.gz.sort.rmdup.bam \
	--labels Ctrl K27Ac ZNF365 --minMappingQuality 20 --skipZeros -p 10 -T Ker_d7 -plot plotFingerprint_d7.png --extendReads 200

echo "Done!"
