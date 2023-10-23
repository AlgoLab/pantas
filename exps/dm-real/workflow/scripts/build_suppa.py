#!/usr/bin/env python3

import sys
import numpy as np
import statistics

ETYPES = ["ES", "IR", "A3", "A5"]
EMAP = {"SE": "ES", "RI": "IR", "A3SS": "A3", "A5SS": "A5"}





def get_interval_s(c, s, e):
    return f"{c}:{s}-{e}"




def parse_suppa(suppa_dpsi, pvalue=0.05):
    EVENTS = {"ES": [], "A3": [], "A5": [], "IR": []}
    for i, line in enumerate(open(suppa_dpsi)):
        if i == 0:
            continue
        idx, dpsi, pv = line.strip("\n").split("\t")
        dpsi, pv = float(dpsi), float(pv)

        if pv > pvalue:
            continue
        gene, rest = idx.split(";")
        etype, chrom, *positions, strand = rest.split(":")
        #print(etype, chrom, positions, strand)

        if etype == "SE":
            ab, cd = positions
            intron1 = tuple(int(x) for x in ab.split("-"))
            intron1 = (intron1[0], intron1[1] - 1)
            intron2 = tuple(int(x) for x in cd.split("-"))
            intron2 = (intron2[0], intron2[1] - 1)
            # k = f"{chrom}:{intron1[0]}-{intron1[1]}-{intron2[0]}-{intron2[1]}"
            k = [
                chrom,
                gene,
                strand,
                get_interval_s(chrom, intron1[0] + 1, intron2[1]),
                get_interval_s(chrom, intron1[0] + 1, intron1[1]),
                get_interval_s(chrom, intron2[0] + 1, intron2[1]),
                "W1",
                "w2",
                "NAN",
                "NAN",
                dpsi,
            ]

            EVENTS["ES"].append(k)
        elif (etype == "A5" and strand == "+") or (etype == "A3" and strand == "-"):
            ab, cd = positions
            shorter_intron = tuple(int(x) for x in ab.split("-"))
            shorter_intron = (shorter_intron[0], shorter_intron[1])
            longer_intron = tuple(int(x) for x in cd.split("-"))
            longer_intron = (longer_intron[0], longer_intron[1])
            k = [
                chrom,
                gene,
                strand,
                get_interval_s(chrom, longer_intron[0] + 1, longer_intron[1] - 1),
                get_interval_s(chrom, shorter_intron[0]+ 1, shorter_intron[1] - 1),
                ".",
                "W1",
                "w2",
                "NAN",
                "NAN",
                dpsi,
            ]
            EVENTS[etype].append(k)
        elif (etype == "A3" and strand == "+") or (etype == "A5" and strand == "-"):
            ab, cd = positions
            shorter_intron = tuple(int(x) for x in ab.split("-"))
            shorter_intron = (shorter_intron[0], shorter_intron[1] - 1)
            longer_intron = tuple(int(x) for x in cd.split("-"))
            longer_intron = (longer_intron[0], longer_intron[1] - 1)
            k = [
                chrom,
                gene,
                strand,
                get_interval_s(chrom, longer_intron[0] + 1, longer_intron[1]),
                get_interval_s(chrom, shorter_intron[0] + 1, shorter_intron[1]),
                ".",
                "W1",
                "w2",
                "NAN",
                "NAN",
                dpsi,
            ]
            EVENTS[etype].append(k)
        elif etype == "RI":
            a, bc, d = positions
            a = int(a)
            d = int(d)
            intron = tuple(int(x) for x in bc.split("-"))
            intron = (intron[0] + 1, intron[1] - 1)
            k = [
                chrom,
                gene,
                strand,
                get_interval_s(chrom, intron[0], intron[1]),
                ".",
                ".",
                "W1",
                "w2",
                "NAN",
                "NAN",
                dpsi,
            ]
            EVENTS["IR"].append(k)
    return EVENTS

def main(suppa_file, p_value, output):

    suppa_events = parse_suppa(suppa_file, p_value)
    # rmats = {x: set() for x in ETYPES}


    with open(output, "w") as f:
        for etype in suppa_events:
            for e in suppa_events[etype]:
                print(etype, "annotated", *map(lambda x: str(x).strip('"'), e), sep=",", file=f)



if __name__ == "__main__":
    suppa_path = snakemake.input.suppa
    min_pvalue = snakemake.params.p_value
    output = snakemake.output.csv
    # suppa_path = sys.argv[1]
    # min_pvalue = float(sys.argv[2])
    # output = sys.argv[3]
    main(suppa_path, min_pvalue, output)
