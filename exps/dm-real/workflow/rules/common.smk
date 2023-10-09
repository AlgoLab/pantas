import sys
import random
from os.path import join as pjoin
from os.path import isfile
import gzip

## params
seed = 23
iterations = 1

random.seed(seed)

software_folder = config["software"]
#input_folder = config["input"]
results_folder = config["res"]

ODIR = config["odir"]

FA = config["fa"]
GTF = config["gtf"]
VCF = config["vcf"]


S1 = config["s1"]
S2 = config["s2"]

L = config["L"]
N = config["N"]

Ks = [31]
#Ws = [-1, 1, 3, 5, 10]
Ws = [-1,2]
# if not isfile(FA + ".fai"):
#     print("\n\nInput reference not indexed, please index with samtools faidx\n\n")
#     sys.exit(1)

chroms = []
for line in open(FA + ".fai"):
    chroms.append(line.split("\t")[0])

samples = []
with gzip.open(VCF, "rt") as fin:
    for line in fin:
        if line.startswith("#C"):
            samples = random.sample(line.split("\t")[9:], iterations)
            break
assert samples != []

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
