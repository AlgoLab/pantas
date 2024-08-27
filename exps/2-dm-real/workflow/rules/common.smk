import sys
import random
from os.path import join as pjoin
from os.path import isfile
import gzip


software_folder = config["software"]
ENVS = pjoin(os.getcwd(), "..", "envs")
ODIR = config["odir"]

FA = config["fa"]
GTF = config["gtf"]
VCF = config["vcf"]

L = config["L"]
S1 = config["s1"]
S2 = config["s2"]

Ws = [5]
p_value = 0.05
min_dpsi = 0.05
min_prob = 0.9
relax = 0
min_coverage = 5

chroms = []
for line in open(FA + ".fai"):
    chroms.append(line.split("\t")[0])

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

rule index_fa:
    input:
        FA
    output:
        FA + ".fai"
    conda: pjoin(ENVS, "samtools.yaml")
    shell:
        "samtools faidx {input}"
