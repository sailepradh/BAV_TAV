#!/usr/bin/python

import os
os.chdir ("/crex/proj/sllstore2017025/nobackup/private/Comparison_BAV_TAV/run2/data/script/")

file1 = "probe_transcript_in_both.txt"
file2 = "BAV_TAV.Proximities.Probe_Distal_SP1_filtered_in_one.txt"
file3 = "BAV_TAV.Proximities_common_probes_SP1_filtered.txt"
fh = open (file3, "w")

common_probe = []

with open (file1, "r") as probes:
    for lines in probes:
        line = lines.strip().split("\t")
        common = line[3]+line[6]+line[0]+line[1]+line[2]
        common_probe.append(common)

with open (file2, "r") as interaction_files:
    for lines in interaction_files:
        line = lines.strip().split("\t")
        gene_id = line[0]+line[1]+line[3]+line[4]+line[5]
        if gene_id in common_probe:
            fh.write ("\t".join(line))
            fh.write("\n")
fh.close()
