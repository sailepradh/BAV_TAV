#!/bin/bash -l
#SBATCH -A snic2016-7-108
#SBATCH -p devcore -n 8
#SBATCH -t 1:00:00
#SBATCH -J Filtered_Subset_18_April
#SBATCH --mail-user sailendra.pradhananga@scilifelab.se
#SBATCH --mail-type=ALL
#SBATCH --tmp=20480


#awk -F '\t' '{ if (($12 != -1 && $13 >= 4 && $14 < 0.01 && $16 >= 4 && $17 < 0.01) || ($12 != -1 && $19 >= 4 && $20 < 0.001 && $22 >= 4 && $23 < 0.001) ) print $0}' THP1run.hg19.Proximities.Probe_Distal.txt >  THP1run.hg19.Proximities.Probe_Distal_SP4_p0.001_filtered.txt


#awk -F '\t' '{ if (($17 != -1 && $18 >= 4 && $19 < 0.001 && $21 >= 4 && $22 < 0.001) || ($17 != -1 && $24 >= 4 && $25 < 0.001 && $27 >= 4 && $28 < 0.001) ) print $0}' THP1run.hg19.Proximities.Probe_Probe.txt > THP1run.hg19.Proximities.Probe_Probe_SP4_p0.001_filtered.txt

#awk -F '\t' '{ if (($12 != -1 && $13 >= 4 && $14 < 0.05 && $16 >= 4 && $17 < 0.05) || ($12 != -1 && $19 >= 4 && $20 < 0.05 && $22 >= 4 && $23 < 0.05) ) print $0}' THP1run.hg19.Proximities.Probe_Distal.txt > THP1run.hg19.Proximities.Probe_Distal_SP4_p0.05_filtered.txt


#awk -F '\t' '{ if (($17 != -1 && $18 >= 4 && $19 < 0.05 && $21 >= 4 && $22 < 0.05) || ($17 != -1 && $24 >= 4 && $25 < 0.05 && $27 >= 4 && $28 < 0.05) ) print $0}' THP1run.hg19.Proximities.Probe_Probe.txt > THP1run.hg19.Proximities.Probe_Probe_SP4_p0.05_filtered.txt

#awk -F '\t' '{ if ($12 != -1 ) print $0}' THP1run.hg19.Proximities.Probe_Distal.txt >  THP1run.hg19.Proximities.Probe_Distal_intrachrom_filtered.txt

#wc -l THP1run.hg19.Proximities.Probe_Distal_intrachrom_filtered.txt

#awk -F '\t' -v OFS="\t" '{ if (($12 != -1 && $13 >= 1  && sqrt(($10-$11)*($10-$11)) < 6500) || ($12 != -1 && $16 >= 1  && sqrt(($10-$11)*($10-$11)) < 6500) ||($12 != -1 && $19 >= 1 && sqrt(($10-$11)*($10-$11)) < 6500))  print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21}' BAVTAV_mergedDigestPeak.hg19.Proximities.Probe_Distal.txt >  TAV_PD.txt

#wc -l TAV_PD.txt

#awk -F '\t' -v OFS="\t" '{ if (($12 != -1 && $22 >= 1  && sqrt(($10-$11)*($10-$11)) < 6500 ) || ($12 != -1 && $25 >= 1  && sqrt(($10-$11)*($10-$11)) < 6500) ||  ($12 != -1 && $28 >= 1  && sqrt(($10-$11)*($10-$11)) < 6500))  print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$22,$23,$24,$25,$26,$27,$28,$29,$30}' BAVTAV_mergedDigestPeak.hg19.Proximities.Probe_Distal.txt  > BAV_PD.txt

#wc -l BAV_PD.txt

#awk -F '\t' -v OFS="\t" '{ if (($17 != -1 && $18 >= 1) || ($17 != -1 && $21 >= 1) ||  ($17 != -1 && $24 >= 1)) print $1,$2,$3,$4,$5,$6,$7,$8,$12,$13,$14,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26}' BAVTAV_mergedDigestPeak.hg19.Proximities.Probe_Probe.txt | sed 1d > tmp3

#wc -l tmp3

#awk -F '\t' -v OFS="\t" '{ if (($17 != -1 && $27 >= 1) || ($17 != -1 && $30 >= 1) || ($17 != -1 && $33 >= 1)) print $1,$2,$3,$4,$5,$6,$7,$8,$12,$13,$14,$17,$27,$28,$29,$30,$31,$32,$33,$34,$35}' BAVTAV_mergedDigestPeak.hg19.Proximities.Probe_Probe.txt | sed 1d  > tmp4

#wc -l tmp4


#awk -F '\t' -v OFS="\t"  '{if (($17 != -1 && $18 >= 1 && $16=="-") || ($17 != -1 && $21 >= 1 && $16=="-") || ($17 != -1 && $24 >= 1 && $16=="-"))  print $1,$2,$3,$4,$5,$6,$7,$8,$12,$14-2500,$14+2500,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26}' BAVTAV_mergedDigestPeak.hg19.Proximities.Probe_Probe.txt > tmp1.txt

#awk -F '\t' -v OFS="\t"  '{if (($17 != -1 && $18 >= 1 && $16=="+") || ($17 != -1 && $21 >= 1 && $16=="+") || ($17 != -1 && $24 >= 1 && $16=="+"))  print $1,$2,$3,$4,$5,$6,$17,$8,$12,$13-2500,$13+2500,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26}' BAVTAV_mergedDigestPeak.hg19.Proximities.Probe_Probe.txt > tmp2.txt

#wc -l tmp1.txt tmp2.txt
#cat tmp1.txt tmp2.txt > TAV_PP.txt

#awk -F '\t' -v OFS="\t" '{if (($17 != -1 && $27 >= 1 && $16 == "-") || ($17 != -1 && $30 >= 1 && $16 == "-") || ($17 != -1 && $33 >= 1 && $16 == "-")) print $1,$2,$3,$4,$5,$6,$7,$8,$12,$14-2500,$14+2500,$17,$27,$28,$29,$30,$31,$32,$33,$34,$35}' BAVTAV_mergedDigestPeak.hg19.Proximities.Probe_Probe.txt > tmp1.txt

#awk -F '\t' -v OFS="\t" '{if (($17 != -1 && $27 >= 1 && $16 == "+") || ($17 != -1 && $30 >= 1 && $16 == "+") || ($17 != -1 && $33 >= 1 && $16 == "+")) print $1,$2,$3,$4,$5,$6,$7,$8,$12,$13-2500,$13+2500,$17,$27,$28,$29,$30,$31,$32,$33,$34,$35}' BAVTAV_mergedDigestPeak.hg19.Proximities.Probe_Probe.txt > tmp2.txt

#wc -l tmp1.txt tmp2.txt
#cat tmp1.txt tmp2.txt > BAV_PP.txt

#rm -r tmp1.txt tmp2.txt

echo "temp files created"

cat TAV_PD.txt TAV_PP.txt > TAV_all_inter_19_apl.txt

wc -l TAV_all_inter_19_apl.txt

cat BAV_PD.txt BAV_PP.txt > BAV_all_inter_19_apl.txt
wc -l BAV_all_inter_19_apl.txt
echo "all interaction files made"
#cat tmp1 tmp3 > TAV_All_Filtered_three_replicates.txt
#rm -r tmp1 tmp3
#cat tmp2 tmp4 > BAV_All_Filtered_three_replicates.txt
#rm -r tmp2 tmp4

#awk '{print $1}' THP1run.hg19.Proximities.both_nLPS.txt | sort | uniq > features_nLPS.txt
#awk '{print $1}' THP1run.hg19.Proximities.both_wLPS.txt | sort | uniq > features_wLPS.txt

#echo "all done"

#time python3 hicap_expr_PP_raw_18_Sep_2018.py SMC_filtered_16_11_18.txt  Common_Expressed_gene_18_Nov.txt  --o SMC_CPM_Exprs_18_Nov.txt
#time python3 hicap_expr_PP_raw_18_Sep_2018.py SMC_filtered_16_11_18.txt  Common_UnExpressed_gene_18_Nov.txt --o SMC_CPM_Unexprs_18_Nov.txt

#hicap_expr_PP_raw_SMC_21_Nov.py
#time python3 hicap_expr_PP_raw_SMC_21_Nov.py SMC_filtered_16_11_18.txt  Genes_target_SMC  --o SMC_All_21_Nov.txt
