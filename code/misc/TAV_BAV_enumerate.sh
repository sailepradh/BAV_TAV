#!/bin/bash -l
#SBATCH -A snic2016-7-108
#SBATCH -p core -n 8
#SBATCH -t 2:00:00
#SBATCH -J Filtered_Subset_18_April
#SBATCH --mail-user sailendra.pradhananga@scilifelab.se
#SBATCH --mail-type=ALL
#SBATCH --tmp=20480

awk -F '\t' -v OFS="\t" '{if (NR ==1) print $0 ; if (($12 != -1 && $14 >= 0.01 && $17 >= 0.01 && $20 >= 0.01  && $22 >= 4 && $23 < 0.01 && $26 >= 0.01 && $29 >= 0.01)) print $0}' BAVTAV_mergedDigestPeak.hg19.Proximities.Probe_Distal.txt > tmp1.txt

echo "BAV1 significant interaction only"
wc -l tmp1.txt

awk -F '\t' -v OFS="\t" '{if (NR ==1) print $0 ; if (($12 != -1 && $14 >= 0.01 && $17 >= 0.01 && $20 >= 0.01  && $23 >= 0.01 && $25 >= 4 && $26 < 0.01 && $29 >= 0.01)) print $0}' BAVTAV_mergedDigestPeak.hg19.Proximities.Probe_Distal.txt > tmp1.txt

echo "BAV2 significant interaction only"
wc -l tmp1.txt

awk -F '\t' -v OFS="\t" '{if (NR ==1) print $0 ; if (($12 != -1 && $14 >= 0.01 &&  $17 >= 0.01 && $20 >= 0.01  && $23 >= 0.01 && $26 >= 0.01 && $28 >=4 && $29 < 0.01)) print $0}' BAVTAV_mergedDigestPeak.hg19.Proximities.Probe_Distal.txt > tmp1.txt

echo "BAV3 significant interaction only"
wc -l tmp1.txt

### two interaction samples

awk -F '\t' -v OFS="\t" '{if (NR ==1) print $0 ; if (($12 != -1 && $14 >= 0.01 && $17 >= 0.01 && $20 >= 0.01  && $22 >= 4 && $23 < 0.01 && $25 >=4 && $26 < 0.01 && $29 >= 0.01)) print $0}' BAVTAV_mergedDigestPeak.hg19.Proximities.Probe_Distal.txt > tmp1.txt

echo "BAV1_BAV2 significant interaction only"
wc -l tmp1.txt

awk -F '\t' -v OFS="\t" '{if (NR ==1) print $0 ; if (($12 != -1 && $14 >= 0.01 && $17 >= 0.01 && $20 >= 0.01  && $23 >= 0.01 && $25 >= 4 && $26 < 0.01 && $28 >=4 && $29 < 0.01)) print $0}' BAVTAV_mergedDigestPeak.hg19.Proximities.Probe_Distal.txt > tmp1.txt

echo "BAV2_BAV3 significant interaction only"
wc -l tmp1.txt

awk -F '\t' -v OFS="\t" '{if (NR ==1) print $0 ; if (($12 != -1 && $14 >= 0.01 &&  $17 >= 0.01 && $20 >= 0.01  && $22 >= 4 && $23 < 0.01 && $26 >= 0.01 && $28 >=4 && $29 < 0.01)) print $0}' BAVTAV_mergedDigestPeak.hg19.Proximities.Probe_Distal.txt > tmp1.txt

echo "BAV1_BAV3 significant interaction only"
wc -l tmp1.txt



# awk -F '\t' -v OFS="\t" '{if (NR ==1) print $0 ; if (($12 != -1 && $14 >= 0.01 && $17 >= 0.01 && $20 >= 0.01  && $22 >=4  && $23 < 0.01 && $25 >= 4 && $26 < 0.01 && $28 >= 4 && $29 < 0.01)) print $0}' BAVTAV_mergedDigestPeak.hg19.Proximities.Probe_Distal.txt  > TAV.txt
#
# echo "BAV1_BAV2_BAV3"
# wc -l TAV.txt
#
# awk -F '\t' -v OFS="\t" '{if (NR ==1) print $0 ; if (($12 != -1 && $13 >= 4 && $14 < 0.01 && $17 >= 0.01 && $20 >= 0.01  && $22 >=4  && $23 < 0.01 && $25 >= 4 && $26 < 0.01 && $29 >= 0.01)) print $0}' BAVTAV_mergedDigestPeak.hg19.Proximities.Probe_Distal.txt  > TAVRandom.txt
#
# echo "BAV1_BAV2_TAV1"
# wc -l TAVRandom.txt
#
# awk -F '\t' -v OFS="\t" '{if (NR ==1) print $0 ; if (($12 != -1  && $14 >= 0.01 && $16 >=4  && $17 < 0.01 && $20 >= 0.01  && $22 >=4  && $23 < 0.01 && $25 >= 4 && $26 < 0.01 && $29 >= 0.01)) print $0}' BAVTAV_mergedDigestPeak.hg19.Proximities.Probe_Distal.txt  > TAVRandom.txt
#
# echo "BAV1_BAV2_TAV2"
# wc -l TAVRandom.txt
#
# awk -F '\t' -v OFS="\t" '{if (NR ==1) print $0 ; if (($12 != -1  && $14 >= 0.01 && $17 >= 0.01 && $19 >= 4 && $20 < 0.01  && $22 >=4  && $23 < 0.01 && $25 >= 4 && $26 < 0.01 && $29 >= 0.01)) print $0}' BAVTAV_mergedDigestPeak.hg19.Proximities.Probe_Distal.txt  > TAVRandom.txt
#
# echo "BAV1_BAV2_TAV3"
# wc -l TAVRandom.txt

#
# awk -F '\t' -v OFS="\t" '{if (NR ==1) print $0 ; if (($12 != -1 && $13 >= 4 && $14 < 0.01 && $17 >= 0.01 && $20 >= 0.01  && $22 >=4  && $23 < 0.01 && $26 >= 0.01 && $28 >= 4 && $29 < 0.01)) print $0}' BAVTAV_mergedDigestPeak.hg19.Proximities.Probe_Distal.txt  > TAVRandom.txt
#
# echo "BAV1_TAV1_BAV3"
# wc -l TAVRandom.txt
#
# awk -F '\t' -v OFS="\t" '{if (NR ==1) print $0 ; if (($12 != -1  && $14 >= 0.01 && $16 >=4  && $17 < 0.01 && $20 >= 0.01  && $22 >=4  && $23 < 0.01 &&  $26 >= 0.01 && $28 >=4  && $29 < 0.01)) print $0}' BAVTAV_mergedDigestPeak.hg19.Proximities.Probe_Distal.txt  > TAVRandom.txt
#
# echo "BAV1_TAV2_BAV3"
# wc -l TAVRandom.txt
#
# awk -F '\t' -v OFS="\t" '{if (NR ==1) print $0 ; if (($12 != -1  && $14 >= 0.01 && $17 >= 0.01 && $19 >= 4 && $20 < 0.01  && $22 >=4  && $23 < 0.01 && $26 >= 0.01 && $28 >= 4 && $29 >= 0.01)) print $0}' BAVTAV_mergedDigestPeak.hg19.Proximities.Probe_Distal.txt  > TAVRandom.txt
#
# echo "BAV1_TAV3_BAV3"
# wc -l TAVRandom.txt
#
# awk -F '\t' -v OFS="\t" '{if (NR ==1) print $0 ; if (($12 != -1 && $13 >= 4 && $14 < 0.01 && $17 >= 0.01 && $20 >= 0.01  && $23 >= 0.01 && $25 >=4 && $26 >= 0.01 && $28 >= 4 && $29 < 0.01)) print $0}' BAVTAV_mergedDigestPeak.hg19.Proximities.Probe_Distal.txt  > TAVRandom.txt
#
# echo "TAV1_BAV2_BAV3"
# wc -l TAVRandom.txt
#
# awk -F '\t' -v OFS="\t" '{if (NR ==1) print $0 ; if (($12 != -1  && $14 >= 0.01 && $16 >=4  && $17 < 0.01 && $20 >= 0.01   && $23 >= 0.01 && $25 >=4 && $26 >= 0.01 && $28 >= 4 && $29 < 0.01)) print $0}' BAVTAV_mergedDigestPeak.hg19.Proximities.Probe_Distal.txt  > TAVRandom.txt
#
# echo "TAV2_BAV2_BAV3"
# wc -l TAVRandom.txt
#
# awk -F '\t' -v OFS="\t" '{if (NR ==1) print $0 ; if (($12 != -1  && $14 >= 0.01 && $17 >= 0.01 && $19 >= 4 && $20 < 0.01   && $23 >= 0.01 && $25 >=4 && $26 >= 0.01 && $28 >= 4 && $29 < 0.01)) print $0}' BAVTAV_mergedDigestPeak.hg19.Proximities.Probe_Distal.txt  > TAVRandom.txt
#
# echo "TAV3_BAV2_BAV3"
# wc -l TAVRandom.txt
