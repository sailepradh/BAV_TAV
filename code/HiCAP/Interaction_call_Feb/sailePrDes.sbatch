#!/bin/bash -l
#SBATCH -A snic2016-7-108
#SBATCH -p core -n 4
#SBATCH -t 24:00:00
#SBATCH -J runBAVTAVFeb042019
#SBATCH --mail-user sailendra.pradhananga@scilifelab.se
#SBATCH --mail-type=ALL
#SBATCH --tmp=20480


#SBATCH -o $1.out
#SBATCH -e $1.err

export LD_LIBRARY_PATH=/crex2/proj/sllstore2017025/private/inhouse_software/internal/PrDe/bamtools/:$LD_LIBRARY_PATH

#gdb -batch -ex "run ProximityDetector -c chr22 -m PrintProximities -p Both" -ex "bt" HiCapTools 2>&1 | grep -v ^"No stack."$ > gdblogT1.out

#./HiCapTools ProximityDetector -m ComputeStatsOnly > ComputeStats.out

chr=$1

./bin/HiCapTools ProximityDetector -c $1 -m PrintProximities -p Both > PrDet.$1.out
