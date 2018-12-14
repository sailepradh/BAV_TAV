#!/bin/bash -l
#SBATCH -A snic2017-7-367
#SBATCH -p core -n 4
#SBATCH -t 12:00:00
#SBATCH -J runBAVTAVreal
#SBATCH --mail-user sailendra.pradhananga@scilifelab.se
#SBATCH --mail-type=ALL
#SBATCH --tmp=20480

chr=$1

#SBATCH -o $1.out
#SBATCH -e $1.err

export LD_LIBRARY_PATH=/domus/h1/sail/HiCapTools/bamtools/:$LD_LIBRARY_PATH

#gdb -batch -ex "run ProximityDetector -c chr22 -m PrintProximities -p Both" -ex "bt" HiCapTools 2>&1 | grep -v ^"No stack."$ > gdblogT1.out

#./HiCapTools ProximityDetector -m ComputeStatsOnly > ComputeStats.out

./bin/HiCapTools ProximityDetector -c $1 -m PrintProximities -p NonNeg > PrDet.$1.out
