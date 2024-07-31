import sys
import random
from os.path import join as pjoin
from os.path import isfile
import gzip


FA = config["fa"]
GTF = config["gtf"]
VCF = config["vcf"]
GENES = config["glist"]
ODIR = config["odir"]

software_folder = config["software"]

S1 = config["s1"]
S2 = config["s2"]

L = config["L"]

Ws = [1, 3, 5, 10]

FQs = {}
C1 = {}
for fq in S1:
    bn, fq1, fq2 = fq.strip("\n").split(",")
    C1[bn] = f"{fq1},{fq2}"
    FQs[bn] = (fq1, fq2)

C2 = {}
for fq in S2:
    bn, fq1, fq2 = fq.strip("\n").split(",")
    C2[bn] = f"{fq1},{fq2}"
    FQs[bn] = (fq1, fq2)
