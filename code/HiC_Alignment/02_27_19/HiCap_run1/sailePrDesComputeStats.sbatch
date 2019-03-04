#!/bin/bash -l
#SBATCH -A snic2016-7-108
#SBATCH -p core -n 4
#SBATCH -t 10:00:00
#SBATCH -J compute_stats_02_27
#SBATCH --mail-user sailendra.pradhananga@scilifelab.se
#SBATCH --mail-type=ALL
#SBATCH --tmp=20480

export LD_LIBRARY_PATH=/crex2/proj/sllstore2017025/private/inhouse_software/internal/PrDe/bamtools/:$LD_LIBRARY_PATH

# bin/HiCapTools ProximityDetector -m ComputeStatsOnly > PrDet.out
# gdb -batch -ex "run ProximityDetector -m ComputeStatsOnly" -ex "bt" HiCapTools 2>&1 a| grep -v ^"No stack."$ > gdblogT1.out

chr=$1
./bin/HiCapTools ProximityDetector -c $1 -m ComputeStatsOnly > PrDet.$1.out
