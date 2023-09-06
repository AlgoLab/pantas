#!/usr/bin/env python3

import sys
import os
from Bio import SeqIO


def main():
    fq1_path = sys.argv[1]
    fq2_path = sys.argv[2]

    new_fq1_path = os.path.splitext(fq1_path)[0] + ".clean.fq"
    new_fq2_path = os.path.splitext(fq2_path)[0] + ".clean.fq"
    new_fq1 = open(new_fq1_path, "w")
    new_fq2 = open(new_fq2_path, "w")

    reads_to_remove = set()
    for record in SeqIO.parse(fq1_path, "fastq"):
        if "mate1Start:1;mate2Start:1" in record.id:
            reads_to_remove.add(record.id.split("/")[0])
    for record in SeqIO.parse(fq2_path, "fastq"):
        if "mate1Start:1;mate2Start:1" in record.id:
            reads_to_remove.add(record.id.split("/")[0])

    for record in SeqIO.parse(fq1_path, "fastq"):
        if record.id.split("/")[0] in reads_to_remove:
            continue
        SeqIO.write(record, new_fq1, "fastq")
    for record in SeqIO.parse(fq2_path, "fastq"):
        if record.id.split("/")[0] in reads_to_remove:
            continue
        SeqIO.write(record, new_fq2, "fastq")

    new_fq1.close()
    new_fq2.close()


if __name__ == "__main__":
    main()
