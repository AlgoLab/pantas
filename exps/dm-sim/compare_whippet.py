import sys
import numpy as np

ETYPES = ["ES", "IR", "A3", "A5"]
EMAP = {"CE": "ES", "RI": "IR", "AD": "A5", "AA": "A3"}


def get_interval(region):
    if region == ".":
        return region
    s, e = [int(x) for x in region.split(":")[1].split("-")]
    return s, e


def main():
    truth_path = sys.argv[1]
    whippet_path = sys.argv[2]
    min_dpsi = float(sys.argv[3])
    min_prob = float(sys.argv[4])
    relax = int(sys.argv[5])

    truth = {x: set() for x in ETYPES}
    for line in open(truth_path):
        etype, chrom, gene, strand, i1, i2, i3, W1, W2, psi1, psi2 = line.strip(
            "\n"
        ).split(",")
        if psi1 == "NaN" or psi2 == "NaN":
            continue

        i1 = get_interval(i1)
        i2 = get_interval(i2)
        i3 = get_interval(i3) if i3 != "." else "."
        if etype == "ES":
            e = (i1[1], i2[0])
        elif etype == "IR":
            e = (i1[0] + 1, i1[1] - 1)
        else:
            if i1[0] == i2[0]:
                e = (min(i1[1], i2[1]), max(i1[1], i2[1]))
            else:  # i1[1] = i2[1]
                e = (min(i1[0], i2[0]) + 1, max(i1[0], i2[0]))
        truth[etype].add(f"{chrom}:{e[0]}-{e[1]}")

    whippet = {x: set() for x in ETYPES}
    for line in open(whippet_path):
        if line.startswith("Gene"):
            continue
        (
            gene,
            _,
            region,
            strand,
            etype,
            psi1,
            psi2,
            dpsi,
            prob,
            compl,
            entr,
        ) = line.strip("\t \n").split("\t")
        dpsi = float(dpsi)
        prob = float(prob)
        if dpsi <= min_dpsi or prob < min_prob:
            continue
        etype = EMAP[etype]
        whippet[etype].add(region)

    print("Event", "TP", "FP", "FN", "P", "R", "F1", sep=",")
    for e in ETYPES:
        TP = len(whippet[e] & truth[e])
        FP = len(whippet[e] - truth[e])
        FN = len(truth[e] - whippet[e])
        if relax > 0:
            TP = 0
            FP = 0
            FN = 0
            for c in whippet[e]:
                c1, c2 = get_interval(c)
                hit = False
                for t in truth[e]:
                    t1, t2 = get_interval(t)
                    if abs(t1 - c1) <= relax and abs(t2 - c2) <= relax:
                        hit = True
                        break
                if hit:
                    TP += 1
            FP = len(whippet[e]) - TP
            FN = len(truth[e]) - TP
        P = TP / (TP + FP) if TP + FP != 0 else 0
        R = TP / (TP + FN) if TP + FN != 0 else 0
        F1 = 2 * (P * R) / (P + R) if (P + R) != 0 else 0
        print(e, TP, FP, FN, P, R, F1, sep=",")


if __name__ == "__main__":
    main()
