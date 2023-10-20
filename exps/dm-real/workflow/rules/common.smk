import sys
import random
from os.path import join as pjoin
from os.path import isfile
import gzip


software_folder = config["software"]

ODIR = config["odir"]

FA = config["fa"]
GTF = config["gtf"]
VCF = config["vcf"]

L = config["L"]
S1 = config["s1"]
S2 = config["s2"]

Ws = [-1, 5, 10]

p_value = 10
min_dpsi = -1
min_prob = -1
relax = 4
min_coverage = 10

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
    conda: "../envs/samtools.yaml"
    shell:
        "samtools faidx {input}"
