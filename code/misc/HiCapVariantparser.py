#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
An updated script on variant manipulator from the HiCap data.
Rare variant are 0.5% within a population at the current criteria.

current use:
python3 ./HiCapVariantparser.py /Volumes/Work_drive/prj/BAV_TAV/data/raw_internal/Interaction_march_corrected/BAVTAV.Proximities.Probe_Distal_SP4_p01_either_sample_dist_width.txt /Volumes/Work_drive/prj/BAV_TAV/data/raw_internal/SNP_calling/Individual_and_combined_VCFs/April_05/Annotated_snp_indel.vcf.gz
'''


import sys
import argparse
import re
from pysam import VariantFile

def Main ():
    parser = argparse.ArgumentParser(description="loading vcf and interaction files")
    parser.add_argument("interactionfile", help = "Interaction calls from HiCap method")
    parser.add_argument("vcfile", help = "Variant calls from either HiCap or sequencing samples")
    parser.add_argument("-o", "--output", help ="output of interaction files", action='store', default=None)
    args = parser.parse_args()
    Vcfin = VariantFile(args.vcfile)


    with open (args.interactionfile, 'r') as f:
        next(f)

        for line in f:
            line = line.strip().split("\t")
            gene = line [0]
            gene_tss = line [2]
            print (gene, gene_tss)
            chr = ((line[8])[3:],line[9], line[10])

            TAV2431 = [line [12],line[13]]
            TAV2515 = [line [15],line[16]]
            TAV2709 = [line [18],line[19]]
            BAV2375 = [line [21],line[22]]
            BAV2424 = [line [24],line[25]]
            BAV2714 = [line [27],line[28]]

            interaction_sample = [TAV2431, TAV2515, TAV2709, BAV2375, BAV2424, BAV2714]
            interaction_binary = []

            for int_sta in interaction_sample:
                if int(int_sta[0]) >= 4 and float(int_sta[1]) < 0.01:
                    interaction_binary.append("1")
                else:
                    interaction_binary.append("0")

            print (interaction_binary)

            sample_list = [3,4,5,0,1,2]
            for rec in Vcfin.fetch (chr[0], int(chr[1]), int(chr[2])):
                genotype_binary = []
                for test in rec.samples.values():
                    genotype = test["GT"]
                    if genotype == (0,1) or genotype == (1,1):
                        genotype_binary.append("1")
                    else:
                        genotype_binary.append("0")

                sorted_genotype = [x for _,x in sorted(zip(sample_list,genotype_binary))]
                print (sorted_genotype)

                zip_array = list(zip (interaction_binary, sorted_genotype))
                count = 0
                for a,b in zip_array:
                    if a == b:
                        count =count+1
                print (count)
                if count == 6:
                    print (zip_array)


            break

if __name__=="__main__":
    Main()
