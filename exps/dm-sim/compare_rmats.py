import sys
import numpy as np

ETYPES = ["ES", "IR", "A3", "A5"]
EMAP = {"SE": "ES", "RI": "IR", "A3SS": "A3", "A5SS": "A5"}


def get_interval(region):
    if region == ".":
        return region
    s, e = [int(x) for x in region.split(":")[1].split("-")]
    return s, e


def get_interval_s(c, s, e):
    return f"{c}:{s}-{e}"


def parse_truth(truth_path, novel):
    truth = {x: set() for x in ETYPES}
    truth_w = {x: {} for x in ETYPES}

    for line in open(truth_path):
        etype, chrom, gene, strand, i1, i2, i3, W1, W2, psi1, psi2 = line.strip(
            "\n"
        ).split(",")
        if psi1 == "NaN" and psi2 == "NaN":
            continue
        k = None
        if novel:
            # if any([int(x) == 0 for x in W1.split("/")]) or any(
            #     [int(x) == 0 for x in W2.split("/")]
            # ):
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
                k = f"{chrom}:{i2[0]}-{i1[0]}-{i1[1]}-{i2[1]}"
        truth[etype].add(k)
        truth_w[etype][k] = (W1, W2)
    return truth, truth_w


def parse_rmats_se(fpath, novel, pvalue=0.05):
    EVENTS = set()
    for line in open(fpath):
        if line.startswith("ID"):
            continue
        try:
            (
                idx,
                gene,
                gene_sym,
                chrom,
                strand,
                ex_s,
                ex_e,
                usex_s,
                usex_e,
                dsex_s,
                dsex_e,
                _idx,
                inc_jc_1,
                sk_jc_1,
                inc_jc_2,
                sk_jc_2,
                inc_len,
                sk_len,
                pv,
                fdr,
                inclvl_1,
                inclvl_2,
                delta_incl,
            ) = line.strip("\n").split("\t")
            pv = float(pv)
            if pv > pvalue:
                continue
        except ValueError:
            (
                idx,
                gene,
                gene_sym,
                chrom,
                strand,
                ex_s,
                ex_e,
                usex_s,
                usex_e,
                dsex_s,
                dsex_e,
            ) = line.strip("\n").split("\t")

        ex_s, usex_s, dsex_s = int(ex_s), int(usex_s), int(dsex_s)
        ex_e, usex_e, dsex_e = int(ex_e), int(usex_e), int(dsex_e)
        # Starts are 0-based, ends aren't (or they are exclusive)
        ex_s += 1
        usex_s += 1
        dsex_s += 1

        chrom = chrom[3:]

        intron1 = (int(usex_e), int(ex_s))
        intron2 = (int(ex_e), int(dsex_s))
        intron3 = (intron1[0], intron2[1])

        k = f"{chrom}:{intron1[0]}-{intron1[1]}-{intron2[0]}-{intron2[1]}"
        if novel:
            k = tuple(
                [
                    chrom,
                    get_interval_s(chrom, intron1[0], intron1[1]),
                    get_interval_s(chrom, intron1[1], intron2[0]),
                    get_interval_s(chrom, intron2[0], intron2[1]),
                ]
            )
        EVENTS.add(k)
    return EVENTS


def parse_rmats_a3(fpath, novel, pvalue=0.05):
    EVENTS = set()
    for line in open(fpath):
        if line.startswith("ID"):
            continue
        # Exon position are 0-based
        try:
            (
                idx,
                gene,
                gene_sym,
                chrom,
                strand,
                lex_s,
                lex_e,
                sex_s,
                sex_e,
                ex_s,
                ex_e,
                _idx,
                inc_jc_1,
                sk_jc_1,
                inc_jc_2,
                sk_jc_2,
                inc_len,
                sk_len,
                pv,
                fdr,
                inclvl_1,
                inclvl_2,
                delta_incl,
            ) = line.strip("\n").split("\t")
            pv = float(pv)
            if pv > pvalue:
                continue
        except ValueError:
            (
                idx,
                gene,
                gene_sym,
                chrom,
                strand,
                ex_s,
                ex_e,
                usex_s,
                usex_e,
                dsex_s,
                dsex_e,
            ) = line.strip("\n").split("\t")

        # Starts are 0-based, ends aren't (or they are exclusive)
        ex_s, lex_s, sex_s = int(ex_s), int(lex_s), int(sex_s)
        ex_e, lex_e, sex_e = int(ex_e), int(lex_e), int(sex_e)
        ex_s += 1
        sex_s += 1
        sex_s += 1

        chrom = chrom[3:]

        # strand +
        longer_intron = (int(ex_e), int(sex_s))
        shorter_intron = (int(ex_e), int(lex_s))
        if strand == "-":
            longer_intron = (int(sex_e), int(ex_s))
            shorter_intron = (int(lex_e), int(ex_s))

        k = ""
        if strand == "+":
            # assert longer_intron[0] == shorter_intron[0]
            k = f"{chrom}:{longer_intron[0]}-{shorter_intron[1]+1}-{longer_intron[1]-1}"
        else:
            # assert longer_intron[1] == shorter_intron[1]
            k = f"{chrom}:{longer_intron[0]}-{shorter_intron[0]}-{longer_intron[1]}"
        if novel:
            k = tuple(
                [
                    chrom,
                    get_interval_s(chrom, longer_intron[0], longer_intron[1]),
                    get_interval_s(chrom, shorter_intron[0], shorter_intron[1]),
                ]
            )
        EVENTS.add(k)
    return EVENTS


def parse_rmats_a5(fpath, novel, pvalue=0.05):
    EVENTS = set()
    for line in open(fpath):
        if line.startswith("ID"):
            continue
        # Exon position are 0-based
        try:
            (
                idx,
                gene,
                gene_sym,
                chrom,
                strand,
                lex_s,
                lex_e,
                sex_s,
                sex_e,
                ex_s,
                ex_e,
                _idx,
                inc_jc_1,
                sk_jc_1,
                inc_jc_2,
                sk_jc_2,
                inc_len,
                sk_len,
                pv,
                fdr,
                inclvl_1,
                inclvl_2,
                delta_incl,
            ) = line.strip("\n").split("\t")
            pv = float(pv)
            if pv > pvalue:
                continue
        except ValueError:
            (
                idx,
                gene,
                gene_sym,
                chrom,
                strand,
                ex_s,
                ex_e,
                usex_s,
                usex_e,
                dsex_s,
                dsex_e,
            ) = line.strip("\n").split("\t")

        # Starts are 0-based, ends aren't (or they are exclusive)
        ex_s, lex_s, sex_s = int(ex_s), int(lex_s), int(sex_s)
        ex_e, lex_e, sex_e = int(ex_e), int(lex_e), int(sex_e)
        ex_s += 1
        sex_s += 1
        sex_s += 1

        chrom = chrom[3:]

        pv = float(pv)
        if pv > pvalue:
            continue

        # strand +
        longer_intron = (int(sex_e), int(ex_s))
        shorter_intron = (int(lex_e), int(ex_s))
        if strand == "-":
            longer_intron = (int(ex_e), int(sex_s))
            shorter_intron = (int(ex_e), int(lex_s))

        k = ""
        if strand == "+":
            k = f"{chrom}:{longer_intron[0]}-{shorter_intron[0]}-{longer_intron[1]}"
        else:
            k = f"{chrom}:{longer_intron[0]}-{shorter_intron[1]+1}-{longer_intron[1]-1}"
        if novel:
            k = tuple(
                [
                    chrom,
                    get_interval_s(chrom, longer_intron[0], longer_intron[1]),
                    get_interval_s(chrom, shorter_intron[0], shorter_intron[1]),
                ]
            )
        EVENTS.add(k)
    return EVENTS


def parse_rmats_ri(fpath, novel, pvalue=0.05):
    EVENTS = set()
    for line in open(fpath):
        if line.startswith("ID"):
            continue
        # Exon position are 0-based
        try:
            (
                idx,
                gene,
                gene_sym,
                chrom,
                strand,
                ex_s,  # retained exon
                ex_e,
                fex_s,  # first exon
                fex_e,
                sex_s,  # second exon
                sex_e,
                _idx,
                inc_jc_1,
                sk_jc_1,
                inc_jc_2,
                sk_jc_2,
                inc_len,
                sk_len,
                pv,
                fdr,
                inclvl_1,
                inclvl_2,
                delta_incl,
            ) = line.strip("\n").split("\t")
            pv = float(pv)
            if pv > pvalue:
                continue
        except ValueError:
            (
                idx,
                gene,
                gene_sym,
                chrom,
                strand,
                ex_s,
                ex_e,
                usex_s,
                usex_e,
                dsex_s,
                dsex_e,
            ) = line.strip("\n").split("\t")

        # Starts are 0-based, ends aren't (or they are exclusive)
        ex_s, fex_s, sex_s = int(ex_s), int(fex_s), int(sex_s)
        ex_e, fex_e, sex_e = int(ex_e), int(fex_e), int(sex_e)
        ex_s += 1
        sex_s += 1
        sex_s += 1
        fex_s += 1

        chrom = chrom[3:]

        pv = float(pv)
        if pv > pvalue:
            continue

        assert ex_s == fex_s and ex_e == sex_e
        # strand +
        k = f"{chrom}:{ex_s}-{fex_e}-{sex_s-1}-{ex_e}"
        if novel:
            k = tuple(
                [
                    chrom,
                    get_interval_s(chrom, ex_s, ex_e),
                    get_interval_s(chrom, fex_e, sex_s - 1),
                ]
            )
        EVENTS.add(k)
    return EVENTS


def main_anno(truth_path, rmats_prefix, min_pvalue):
    truth, truth_w = parse_truth(truth_path, False)

    rmats = {x: set() for x in ETYPES}
    rmats["ES"] = parse_rmats_se(rmats_prefix + "/SE.MATS.JC.txt", False, min_pvalue)
    rmats["A3"] = parse_rmats_a3(rmats_prefix + "/A3SS.MATS.JC.txt", False, min_pvalue)
    rmats["A5"] = parse_rmats_a5(rmats_prefix + "/A5SS.MATS.JC.txt", False, min_pvalue)
    rmats["IR"] = parse_rmats_ri(rmats_prefix + "/RI.MATS.JC.txt", False, min_pvalue)

    print("Event", "C", "T", "TP", "FP", "FN", "P", "R", "F1", sep=",")
    for etype in ETYPES:
        # print(rmats[etype])
        # print(truth[etype])
        TP = len(rmats[etype] & truth[etype])
        # print("TP:", rmats[etype] & truth[etype])
        FP = len(rmats[etype] - truth[etype])
        # print("FP:", rmats[etype] - truth[etype])
        FN = len(truth[etype] - rmats[etype])
        # print("FN:", truth[etype] - rmats[etype])
        # TODO: can we relax this too? (see compare_whippet.py)
        P = TP / (TP + FP) if TP + FP != 0 else 0
        R = TP / (TP + FN) if TP + FN != 0 else 0
        F1 = 2 * (P * R) / (P + R) if (P + R) != 0 else 0
        P = round(P, 3)
        R = round(R, 3)
        F1 = round(F1, 3)
        print(
            etype, len(rmats[etype]), len(truth[etype]), TP, FP, FN, P, R, F1, sep=","
        )

        for e in rmats[etype] - truth[etype]:
            print("FP", e, file=sys.stderr)
        for e in truth[etype] - rmats[etype]:
            print("FN", e, file=sys.stderr)


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
            pass
        i += hit
    # print(t1, t2, i)
    return i


def main_novel(truth_path, rmats_prefix, min_pvalue):
    truth, truth_w = parse_truth(truth_path, True)

    rmats = {x: set() for x in ETYPES}
    rmats["ES"] = parse_rmats_se(rmats_prefix + "/SE.MATS.JC.txt", True, min_pvalue)
    rmats["ES"].union(
        parse_rmats_se(
            rmats_prefix + "/fromGTF.novelSpliceSite.SE.txt", True, min_pvalue
        )
    )
    rmats["ES"].union(
        parse_rmats_se(rmats_prefix + "/fromGTF.novelJunction.SE.txt", True, min_pvalue)
    )
    rmats["A3"] = parse_rmats_a3(rmats_prefix + "/A3SS.MATS.JC.txt", True, min_pvalue)
    rmats["A3"].union(
        parse_rmats_se(
            rmats_prefix + "/fromGTF.novelSpliceSite.A3SS.txt", True, min_pvalue
        )
    )
    rmats["A3"].union(
        parse_rmats_se(
            rmats_prefix + "/fromGTF.novelJunction.A3SS.txt", True, min_pvalue
        )
    )
    rmats["A5"] = parse_rmats_a5(rmats_prefix + "/A5SS.MATS.JC.txt", True, min_pvalue)
    rmats["A5"].union(
        parse_rmats_se(
            rmats_prefix + "/fromGTF.novelSpliceSite.A5SS.txt", True, min_pvalue
        )
    )
    rmats["A5"].union(
        parse_rmats_se(
            rmats_prefix + "/fromGTF.novelJunction.A5SS.txt", True, min_pvalue
        )
    )
    rmats["IR"] = parse_rmats_ri(rmats_prefix + "/RI.MATS.JC.txt", True, min_pvalue)
    rmats["IR"].union(
        parse_rmats_se(
            rmats_prefix + "/fromGTF.novelSpliceSite.RI.txt", True, min_pvalue
        )
    )
    rmats["IR"].union(
        parse_rmats_se(rmats_prefix + "/fromGTF.novelJunction.RI.txt", True, min_pvalue)
    )

    print("Event", "C", "T", "TP", "FP", "FN", "P", "R", "F1", sep=",")
    for etype in ETYPES:
        TP = 0
        for e1 in rmats[etype]:
            hit = False
            for e2 in truth[etype]:
                if e1[0] != e2[0]:
                    # different chrom
                    continue
                if (
                    relaxed_intersect(e1, e2) >= len(e1) / 2
                ):  # /2 to check for half junctions, 2 for SE, 1 for others
                    hit = True
                    break
            if hit:
                TP += 1
            else:
                print(etype, e1, file=sys.stderr)
        FP = len(rmats[etype]) - TP
        FN = len(truth[etype]) - TP
        P = TP / (TP + FP) if TP + FP != 0 else 0
        R = TP / (TP + FN) if TP + FN != 0 else 0
        F1 = 2 * (P * R) / (P + R) if (P + R) != 0 else 0
        P = round(P, 3)
        R = round(R, 3)
        F1 = round(F1, 3)
        print(
            etype, len(rmats[etype]), len(truth[etype]), TP, FP, FN, P, R, F1, sep=","
        )

        for e1 in truth[etype]:
            hit = False
            for e2 in rmats[etype]:
                if e1[0] != e2[0]:
                    # different chrom
                    continue
                if (
                    relaxed_intersect(e1, e2) >= len(e2) / 2
                ):  # /2 to check for half junctions, 2 for SE, 1 for others
                    hit = True
                    break
            if not hit:
                print("FN", etype, e1, file=sys.stderr)


if __name__ == "__main__":
    mode = sys.argv[1]
    truth_path = sys.argv[2]
    rmats_path = sys.argv[3]
    min_pvalue = 10
    if mode == "anno":
        main_anno(truth_path, rmats_path, min_pvalue)
    elif mode == "novel":
        main_novel(truth_path, rmats_path, min_pvalue)
