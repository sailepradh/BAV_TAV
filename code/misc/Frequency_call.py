#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import re
from pysam import VariantFile


os.chdir ("/Volumes/Work_drive/prj/Rare_variants_atherosclerosis/data/swegen_20161223")
VCFin = VariantFile("SNP_VCF.gz")
for rec in VCFin.fetch("1",8964100,8964500):
    print (rec.allele)
