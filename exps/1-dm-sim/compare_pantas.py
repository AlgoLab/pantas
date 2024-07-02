import sys
import numpy as np

ETYPES = ["ES", "IR", "A3", "A5"]


def parse_truth(truth_path, novel):
    truth = {x: set() for x in ETYPES}
    truth_w = {x: {} for x in ETYPES}

    for line in open(truth_path):
        etype, chrom, gene, strand, i1, i2, i3, W1, W2, psi1, psi2 = line.strip(
            "\n"
        ).split(",")
        if psi1 == "NaN" or psi2 == "NaN":
            continue
        k = None
        if novel:
            if W1[-2:] == "/0" or W2[-2:] == "/0":
                continue
            if i3 == ".":
                k = (chrom, i1, i2)
            else:
                k = (chrom, i1, i2, i3)
        else:
            i1 = get_interval(i1)
            i2 = get_interval(i2)
            if etype == "ES":
                k = f"{chrom}:{i1[0]}-{i1[1]}-{i2[0]}-{i2[1]}"
            elif etype[0] == "A":
                if i1[0] == i2[0]:
                    k = f"{chrom}:{i1[0]}-{min(i1[1], i2[1])}-{max(i1[1], i2[1])}"
                else:
                    k = f"{chrom}:{min(i1[0], i2[0])}-{max(i1[0], i2[0])}-{i1[1]}"
            elif etype == "IR":
                k = f"{chrom}:{i1[0]}-{i1[1]}"
        truth[etype].add(k)
        truth_w[etype][k] = (W1, W2)
    return truth, truth_w


def get_interval(region):
    if region == ".":
        return region
    s, e = [int(x) if x != "?" else -1 for x in region.split(":")[1].split("-")]
    return s, e


def main_anno(truth_path, pantas_path):
    truth, truth_w = parse_truth(truth_path, False)

    pantas = {x: set() for x in ETYPES + [_ + "-all" for _ in ETYPES]}
    for line in open(pantas_path):
        (
            etype,
            novel,
            chrom,
            gene,
            strand,
            i1,
            i2,
            i3,
            W1,
            W2,
            psi1,
            psi2,
            dpsi,
        ) = line.strip("\n").split(",")
        if etype not in ETYPES:
            continue
        assert novel == "annotated"
        if psi1 == "NaN" or psi2 == "NaN":
            etype += "-all"

        k = ""
        if etype[:2] == "ES":
            i2 = get_interval(i2)
            i3 = get_interval(i3)
            k = f"{chrom}:{i2[0]-1}-{i2[1]+1}-{i3[0]-1}-{i3[1]+1}"
        elif etype[0] == "A":
            i1 = get_interval(i1)
            i2 = get_interval(i2)
            if i1[0] == i2[0]:
                k = f"{chrom}:{i1[0]-1}-{min(i1[1], i2[1])+1}-{max(i1[1], i2[1])+1}"
            else:
                k = f"{chrom}:{min(i1[0], i2[0])-1}-{max(i1[0], i2[0])-1}-{i1[1]+1}"
        elif etype[:2] == "IR":
            i1 = get_interval(i1)  # intron
            i2 = get_interval(i2)  # nothing
            # print(line)
            k = f"{chrom}:{i1[0]-1}-{i1[1]+1}"

        pantas[etype].add(k)
        if not etype.endswith("-all"):
            pantas[etype + "-all"].add(k)

    print("Event", "C", "T", "TP", "FP", "FN", "P", "R", "F1", sep=",")
    for etype in ETYPES:
        TP = len(pantas[etype] & truth[etype])
        FP = len(pantas[etype] - truth[etype])
        FN = len(truth[etype] - pantas[etype])
        # TODO: relax this
        P = TP / (TP + FP) if TP + FP != 0 else 0
        R = TP / (TP + FN) if TP + FN != 0 else 0
        F1 = 2 * (P * R) / (P + R) if (P + R) != 0 else 0
        P = round(P * 100, 2)
        R = round(R * 100, 2)
        F1 = round(F1 * 100, 2)
        print(
            etype, len(pantas[etype]), len(truth[etype]), TP, FP, FN, P, R, F1, sep=","
        )

        for e in pantas[etype] - truth[etype]:
            print("FP", etype, e, file=sys.stderr)
        for e in truth[etype] - pantas[etype]:
            print("FN", etype, e, file=sys.stderr)

        TP = len(pantas[etype + "-all"] & truth[etype])
        FP = len(pantas[etype + "-all"] - truth[etype])
        FN = len(truth[etype] - pantas[etype + "-all"])
        # TODO: relax this
        P = TP / (TP + FP) if TP + FP != 0 else 0
        R = TP / (TP + FN) if TP + FN != 0 else 0
        F1 = 2 * (P * R) / (P + R) if (P + R) != 0 else 0
        P = round(P * 100, 2)
        R = round(R * 100, 2)
        F1 = round(F1 * 100, 2)
        print(
            etype + "-all",
            len(pantas[etype + "-all"]),
            len(truth[etype]),
            TP,
            FP,
            FN,
            P,
            R,
            F1,
            sep=",",
        )


def relaxed_intersect(t1, t2, relax=3):
    if t1[0] != t2[0]:
        # Different chromosome
        return 0
    t1 = [get_interval(x) for x in t1[1:]]
    t2 = [get_interval(x) for x in t2[1:]]
    i = 0
    for x1 in t1:
        hit = False
        for x2 in t2:
            if abs(x2[0] - x1[0]) <= relax and abs(x2[1] - x1[1]) <= relax:
                hit = True
                break
        i += hit
    # print(t1, t2, i)
    return i


def main_novel(truth_path, pantas_path, relax=3):
    truth, truth_w = parse_truth(truth_path, True)

    pantas = {x: set() for x in ETYPES + [_ + "-all" for _ in ETYPES]}
    for line in open(pantas_path):
        print(line, file=sys.stderr)
        (
            etype,
            novel,
            chrom,
            gene,
            strand,
            i1,
            i2,
            i3,
            W1,
            W2,
            psi1,
            psi2,
            dpsi,
        ) = line.strip("\n").split(",")
        if etype == "CE":
            etype == "SE"
        if etype not in ETYPES:
            continue
        assert novel == "novel"
        if psi1 == "NaN" or psi2 == "NaN":
            etype += "-all"

        k = [chrom]
        for i in [i1, i2, i3]:
            if i != "?" and i != ".":
                k.append(i)
        pantas[etype].add(tuple(k))
        if not etype.endswith("-all"):
            pantas[etype + "-all"].add(tuple(k))

    print("Event", "C", "T", "TP", "FP", "FN", "P", "R", "F1", sep=",")
    for etype in ETYPES:
        TP = 0
        for e1 in pantas[etype]:
            hit = False
            for e2 in truth[etype]:
                if e1[0] != e2[0]:
                    # different chrom
                    continue
                if relaxed_intersect(e1, e2, relax) == len(e1) - 1:
                    hit = True
                    break
            if hit:
                TP += 1
            else:
                print("FP", etype, e1, file=sys.stderr)
        FP = len(pantas[etype]) - TP
        FN = len(truth[etype]) - TP
        P = TP / (TP + FP) if TP + FP != 0 else 0
        R = TP / (TP + FN) if TP + FN != 0 else 0
        F1 = 2 * (P * R) / (P + R) if (P + R) != 0 else 0
        P = round(P * 100, 2)
        R = round(R * 100, 2)
        F1 = round(F1 * 100, 2)
        print(
            etype, len(pantas[etype]), len(truth[etype]), TP, FP, FN, P, R, F1, sep=","
        )

        for e1 in truth[etype]:
            hit = False
            for e2 in pantas[etype]:
                if e1[0] != e2[0]:
                    # different chrom
                    continue
                if relaxed_intersect(e1, e2, relax) == len(e2) - 1:
                    hit = True
                    break
            if not hit:
                print("FN", etype, e1, truth_w[etype][e1], file=sys.stderr)

        TP = 0
        for e1 in pantas[etype + "-all"]:
            hit = False
            for e2 in truth[etype]:
                if e1[0] != e2[0]:
                    # different chrom
                    continue
                if relaxed_intersect(e1, e2, relax) == len(e1) - 1:
                    hit = True
                    break
            if hit:
                TP += 1
            else:
                pass  # print("FP", etype, e1, file=sys.stderr)
        FP = len(pantas[etype + "-all"]) - TP
        FN = len(truth[etype]) - TP
        P = TP / (TP + FP) if TP + FP != 0 else 0
        R = TP / (TP + FN) if TP + FN != 0 else 0
        F1 = 2 * (P * R) / (P + R) if (P + R) != 0 else 0
        P = round(P * 100, 2)
        R = round(R * 100, 2)
        F1 = round(F1 * 100, 2)
        print(
            etype + "-all",
            len(pantas[etype + "-all"]),
            len(truth[etype]),
            TP,
            FP,
            FN,
            P,
            R,
            F1,
            sep=",",
        )


if __name__ == "__main__":
    mode = sys.argv[1]
    truth_path = sys.argv[2]
    pantas_path = sys.argv[3]
    relax = int(sys.argv[4])
    if mode == "anno":
        main_anno(truth_path, pantas_path)
    elif mode == "novel":
        main_novel(truth_path, pantas_path, relax)
