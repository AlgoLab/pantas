#!/usr/bin/env python3
import sys
import os
import itertools
import numpy as np
import math

from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn3
import seaborn as sns
from scipy.stats import pearsonr
import pandas as pd

ETYPES = ["ES", "IR", "A3", "A5"]
EMAP_RMATS = {"SE": "ES", "RI": "IR", "A3SS": "A3", "A5SS": "A5"}
EMAP_WHIPPET = {"CE": "ES", "RI": "IR", "AD": "A5", "AA": "A3"}


def get_interval(region):
    if region == ".":
        return region
    s, e = [int(x) for x in region.split(":")[1].split("-")]
    return s, e


def get_interval_s(c, s, e):
    return f"{c}:{s}-{e}"


def parse_pantas(pantas_path):
    print("Parsing pantas")
    pantas = {x: set() for x in ETYPES}
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
            continue
        k = ""
        if etype == "ES":
            i2 = get_interval(i2)
            i3 = get_interval(i3)
            k = f"{chrom}:{i2[0]}-{i2[1]}-{i3[0]}-{i3[1]}"
        elif etype[0] == "A":
            i1 = get_interval(i1)
            i2 = get_interval(i2)
            if i1[0] == i2[0]:
                k = f"{chrom}:{i1[0]}-{min(i1[1], i2[1])}-{max(i1[1], i2[1])}"
            else:
                k = f"{chrom}:{min(i1[0], i2[0])}-{max(i1[0], i2[0])}-{i1[1]}"
        elif etype == "IR":
            i1 = get_interval(i1)
            i2 = get_interval(i2)
            k = f"{chrom}:{i2[0]}-{i1[0]}-{i1[1]}-{i2[1]}"

        pantas[etype].add((k, dpsi))
    return pantas


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
        EVENTS.add((k, delta_incl))
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
            k = f"{chrom}:{longer_intron[0]}-{shorter_intron[1]}-{longer_intron[1]}"
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
        EVENTS.add((k, delta_incl))
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
            k = f"{chrom}:{longer_intron[0]}-{shorter_intron[1]}-{longer_intron[1]}"
        if novel:
            k = tuple(
                [
                    chrom,
                    get_interval_s(chrom, longer_intron[0], longer_intron[1]),
                    get_interval_s(chrom, shorter_intron[0], shorter_intron[1]),
                ]
            )
        EVENTS.add((k, delta_incl))
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
        EVENTS.add((k, delta_incl))
    return EVENTS


def parse_rmats(rmats_prefix, min_pvalue=10):
    print("Parsing rmats\t prefix: ", rmats_prefix, file=sys.stderr)
    rmats = {x: set() for x in ETYPES}
    rmats["ES"] = parse_rmats_se(rmats_prefix + "/SE.MATS.JC.txt", False, min_pvalue)
    rmats["A3"] = parse_rmats_a3(rmats_prefix + "/A3SS.MATS.JC.txt", False, min_pvalue)
    rmats["A5"] = parse_rmats_a5(rmats_prefix + "/A5SS.MATS.JC.txt", False, min_pvalue)
    rmats["IR"] = parse_rmats_ri(rmats_prefix + "/RI.MATS.JC.txt", False, min_pvalue)
    return rmats


def parse_whippet(whippet_path, min_dpsi=-1, min_prob=-1):
    print("Parsing whippet", file=sys.stderr)
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
        etype = EMAP_WHIPPET[etype]
        whippet[etype].add((region, dpsi))
    return whippet


def main(argv):
    pantas_path = argv[0]
    if pantas_path.endswith("/"):
        pantas_path = pantas_path[:-1]
    rmats_path = argv[1]
    if rmats_path.endswith("/"):
        rmats_path = rmats_path[:-1]
    whippet_path = argv[2]
    if whippet_path.endswith("/"):
        whippet_path = whippet_path[:-1]
    output_dir = argv[3]
    if output_dir.endswith("/"):
        output_dir = output_dir[:-1]
    w = argv[4]
    output = f"{output_dir}"

    p_value = 10
    min_dpsi = -1
    min_prob = -1
    if len(argv) > 5:
        p_value = float(argv[5])
        if len(argv) > 6:
            min_dpsi = float(argv[6])
            min_prob = float(argv[7])

    pantas = parse_pantas(pantas_path)
    rmats = parse_rmats(rmats_path, p_value)
    whippet = parse_whippet(whippet_path, min_dpsi, min_prob)
    # print(rmats)
    # print(whippet)

    wh = set([x[0] for x in whippet["ES"]])
    rm = set([x[0] for x in rmats["ES"]])
    pan = set([x[0] for x in pantas["ES"]])
    plt.title(f"ES")
    venn3([wh, rm, pan], set_labels=("Whippet", "rMATS", "Pantas2"))
    # plt.savefig(f"{output}/SE_venn3_{w}.png")
    plt.clf()

    wh = set([x[0] for x in whippet["A3"]])
    rm = set([x[0] for x in rmats["A3"]])
    pan = set([x[0] for x in pantas["A3"]])
    venn3([wh, rm, pan], set_labels=("Whippet", "rMATS", "Pantas2"))
    plt.title(f"A3")
    # plt.savefig(f"{output}/A3_venn3_{w}.png")
    plt.clf()

    wh = set([x[0] for x in whippet["A5"]])
    rm = set([x[0] for x in rmats["A5"]])
    pan = set([x[0] for x in pantas["A5"]])
    venn3([wh, rm, pan], set_labels=("Whippet", "rMATS", "Pantas2"))
    plt.title(f"A5")
    # plt.savefig(f"{output}/A5_venn3_{w}.png")
    plt.clf()

    wh = set([x[0] for x in whippet["IR"]])
    rm = set([x[0] for x in rmats["IR"]])
    pan = set([x[0] for x in pantas["IR"]])
    venn3([wh, rm, pan], set_labels=("Whippet", "rMATS", "Pantas2"))
    plt.title(f"IR")
    # plt.savefig(f"{output}/RI_venn3_{w}.png")
    plt.clf()

    wh = set(
        itertools.chain(
            set([x[0] for x in whippet["ES"]]),
            set([x[0] for x in whippet["A3"]]),
            set([x[0] for x in whippet["A5"]]),
            set([x[0] for x in whippet["IR"]]),
        )
    )
    rm = set(
        itertools.chain(
            set([x[0] for x in rmats["ES"]]),
            set([x[0] for x in rmats["A3"]]),
            set([x[0] for x in rmats["A5"]]),
            set([x[0] for x in rmats["IR"]]),
        )
    )
    pan = set(
        itertools.chain(
            set([x[0] for x in pantas["ES"]]),
            set([x[0] for x in pantas["A3"]]),
            set([x[0] for x in pantas["A5"]]),
            set([x[0] for x in pantas["IR"]]),
        )
    )
    wh_psi = set(
        itertools.chain(
            set([x[1] for x in whippet["ES"]]),
            set([x[1] for x in whippet["A3"]]),
            set([x[1] for x in whippet["A5"]]),
            set([x[1] for x in whippet["IR"]]),
        )
    )
    rm_psi = set(
        itertools.chain(
            set([x[1] for x in rmats["ES"]]),
            set([x[1] for x in rmats["A3"]]),
            set([x[1] for x in rmats["A5"]]),
            set([x[1] for x in rmats["IR"]]),
        )
    )
    pan_psi = set(
        itertools.chain(
            set([x[1] for x in pantas["ES"]]),
            set([x[1] for x in pantas["A3"]]),
            set([x[1] for x in pantas["A5"]]),
            set([x[1] for x in pantas["IR"]]),
        )
    )
    venn3([wh, rm, pan], set_labels=("Whippet", "rMATS", "Pantas2"))
    plt.title(f"ES/RI/A3/A5")
    # plt.savefig(f"{output}/tot_venn3_{w}.png")
    plt.clf()
    # print(pan)
    # print(rm)
    # print(wh)
    data = [pan_psi, wh_psi]
    ## create pandas dataframe with columns type, event, values
    ## type is one of ES, RI, A3, A5
    ## event is the event name
    ## value is the psi value for pantas, psi value for whippet, psi value for rmats
    data = {}

    for key in ["ES", "A3", "A5", "IR"]:
        for event in pantas[key]:
            data[event] = {
                "type": key,
                "event": event[0],
                "pantas": event[1],
                "whippet": math.nan,
                "rmats": math.nan,
            }
        for event in whippet[key]:
            if event[0] in data:
                data[event[0]]["whippet"] = event[1]
            else:
                data[event[0]] = {
                    "type": key,
                    "event": event[0],
                    "pantas": math.nan,
                    "whippet": event[1],
                    "rmats": math.nan,
                }
        for event in rmats[key]:
            if event[0] in data:
                data[event[0]]["rmats"] = event[1]
            else:
                data[event[0]] = {
                    "type": key,
                    "event": event[0],
                    "pantas": math.nan,
                    "whippet": math.nan,
                    "rmats": event[1],
                }
    df = pd.DataFrame(data.values())
    df.to_csv(f"{output}/res_{w}.csv", index=False)
    # print(df)
    g = sns.jointplot(
        data=df,
        x="pantas",
        y="whippet",
        kind="scatter",
        xlim=(-1.05, 1.05),
        ylim=(-1.05, 1.05),
    )
    plt.tight_layout()
    # plt.savefig(f"{output}/whippet_pantas2_{w}.png")


if __name__ == "__main__":
    main(sys.argv[1:])
