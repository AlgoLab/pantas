#!/usr/bin/env python3

#!/usr/bin/env python3
import sys
import os
import itertools
import numpy as np
import math
import copy
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn3
import seaborn as sns
from scipy.stats import pearsonr
import pandas as pd

ETYPES = ["ES", "IR", "A3", "A5"]
EMAP_RMATS = {"SE": "ES", "RI": "IR", "A3SS": "A3", "A5SS": "A5"}
EMAP_WHIPPET = {"CE": "ES", "RI": "IR", "AD": "A5", "AA": "A3"}


def update_column(row, col_name):
    if not pd.isna(row[col_name]):
        return row["event"]
    else:
        return pd.NA


def pairs(l):
    c = []

    for i in range(len(l)):
        for j in range(i + 1, len(l)):
            c.append((l[i], l[j]))

    return c


def get_interval(region):
    if region == ".":
        return region
    s, e = [int(x) for x in region.split(":")[1].split("-")]
    return s, e


def get_interval_s(c, s, e):
    return f"{c}:{s}-{e}"


def check_whippet(whippet, truth, relax=4):
    hits = {}
    p = []
    for e in ETYPES:
        hits[e] = {}
        if relax > 0:
            for c in whippet[e]:
                c1, c2 = get_interval(c[0])
                hit = False
                ##print("whippet", c, c1, c2)
                for t in truth[e]:
                    i1, i2, i3 = t[2]
                    ##print(i1, i2, i3)
                    if e == "ES":
                        ee = (i1[1], i2[0])
                    elif e == "IR":
                        ee = (i1[0] + 1, i1[1] - 1)
                    else:
                        if i1[0] == i2[0]:
                            ee = (min(i1[1], i2[1]), max(i1[1], i2[1]))
                        else:  # i1[1] = i2[1]
                            ee = (min(i1[0], i2[0]) + 1, max(i1[0], i2[0]))
                    if abs(ee[0] - c1) <= relax and abs(ee[1] - c2) <= relax:
                        hit = True
                        # print("hit", c, ee)
                        break
                if hit:
                    hits[e][t[0]] = c[1]
                    p.append(c[0])
    # print("whippet hits", hits)
    return hits, p


def relaxed_intersect(t1, t2, relax=4):
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


def check_rmats(rmats, truth, relax=4):
    hits = {}
    p = []
    for etype in ETYPES:
        hits[etype] = {}
        if relax > 0:
            for e1 in rmats[etype]:
                hit = False
                c1 = e1[0].split(":")[0]
                r1 = e1[0].split(":")[1].split("-")

                for e2 in truth[etype]:
                    c2 = e2[0].split(":")[0]

                    if c1[0] != c2[0]:
                        continue
                    r2 = e2[0].split(":")[1].split("-")

                    # print(r1, r2)
                    h = 0
                    for t1, t2 in zip(r1, r2):
                        if abs(int(t1) - int(t2)) <= relax:
                            h = h + 1
                        else:
                            break
                    if h == len(r1):
                        hit = True
                        # print("hit", e1, e2)
                        break
                if hit:
                    # print(e1)
                    hits[etype][e2[0]] = e1[1]
                    p.append(e1[0])

    # print("rmats, hits", hits)
    return hits, p


def parse_pantas(pantas_path):
    print(f"Parsing pantas\t file = {pantas_path}", file=sys.stderr)
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
        if etype == "Event":
            continue
        i1_o = (
            int(i1.split(":")[1].split("-")[0]),
            int(i1.split(":")[1].split("-")[1]),
        )
        i2_o = (
            int(i2.split(":")[1].split("-")[0]),
            int(i2.split(":")[1].split("-")[1]),
        )
        if i3 != ".":
            i3_o = (
                int(i3.split(":")[1].split("-")[0]),
                int(i3.split(":")[1].split("-")[1]),
            )
        else:
            i3_o = (math.nan, math.nan)

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

        pantas[etype].add((k, dpsi, (i1_o, i2_o, i3_o)))
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
    print(f"Parsing rmats\tprefix: {rmats_prefix}", file=sys.stderr)
    rmats = {x: set() for x in ETYPES}
    rmats["ES"] = parse_rmats_se(rmats_prefix + "/SE.MATS.JC.txt", False, min_pvalue)
    rmats["A3"] = parse_rmats_a3(rmats_prefix + "/A3SS.MATS.JC.txt", False, min_pvalue)
    rmats["A5"] = parse_rmats_a5(rmats_prefix + "/A5SS.MATS.JC.txt", False, min_pvalue)
    rmats["IR"] = parse_rmats_ri(rmats_prefix + "/RI.MATS.JC.txt", False, min_pvalue)
    return rmats


def parse_whippet(whippet_path, min_dpsi=-1, min_prob=-1):
    print(f"Parsing whippet\t file = {whippet_path}", file=sys.stderr)
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
    pantas_path = snakemake.params.pantas_dir
    if pantas_path.endswith("/"):
        pantas_path = pantas_path[:-1]
    rmats_path = snakemake.params.rmats_dir
    if rmats_path.endswith("/"):
        rmats_path = rmats_path[:-1]
    whippet_path = snakemake.input.w
    if whippet_path.endswith("/"):
        whippet_path = whippet_path[:-1]
    output_dir = snakemake.params.res_dir
    if output_dir.endswith("/"):
        output_dir = output_dir[:-1]
    Ws = snakemake.params.Ws
    output = f"{output_dir}"

    p_value = snakemake.params.p_value
    min_dpsi = snakemake.params.min_dpsi
    min_prob = snakemake.params.min_prob
    relax = snakemake.params.relax
    pantas = {}
    for w in Ws:
        pantas[w] = parse_pantas(f"{pantas_path}/quant.w{w}.csv")
    rmats = parse_rmats(rmats_path, p_value)
    whippet = parse_whippet(whippet_path, min_dpsi, min_prob)

    data = {}
    # pantas_keys = list(pantas.keys())
    p_d = copy.deepcopy(pantas[Ws[0]])
    key_names = [f"pantas_{w}" for w in Ws]
    for key in ["ES", "A3", "A5", "IR"]:
        for event in pantas[Ws[0]][key]:
            data[event[0]] = {
                "type": key,
                "event": event[0],
                "whippet": math.nan,
                "rmats": math.nan,
                f"pantas_{Ws[0]}": float(event[1]),
            }

    for i, w in enumerate(Ws[1:]):
        # print(w)
        for key in ["ES", "A3", "A5", "IR"]:
            for event in pantas[w][key]:
                # print(event[0])
                if event[0] in data.keys():
                    data[event[0]][f"pantas_{w}"] = float(event[1])
                else:
                    tmp_dict = {
                        "type": key,
                        "event": event[0],
                        "whippet": math.nan,
                        "rmats": math.nan,
                        f"pantas_{w}": float(event[1]),
                    }
                    for j in Ws[0 : i + 1]:
                        tmp_dict[f"pantas_{j}"] = math.nan
                    data[event[0]] = tmp_dict
                    p_d[key].add((event[0], event[1], event[2]))
                    # print("no",data[event[0]])
    # print(data)
    # print(p_d)
    # pantas_data = copy.deepcopy(data)
    mask_whippet, p_w = check_whippet(whippet, p_d, relax)
    mask_rmats, p_r = check_rmats(rmats, p_d, relax)
    for event, row in data.items():
        e = row["type"]
        if event in mask_whippet[e].keys():
            data[event]["whippet"] = float(mask_whippet[e][event])
        if event in mask_rmats[e].keys():
            data[event]["rmats"] = float(mask_rmats[e][event])
    for key in ["ES", "A3", "A5", "IR"]:
        for event in whippet[key]:
            if event[0] not in p_w:
                tmp_dict = {
                    "type": key,
                    "event": event[0],
                    "whippet": float(event[1]),
                    "rmats": math.nan,
                }
                for j in Ws:
                    tmp_dict[f"pantas_{j}"] = math.nan
                data[event[0]] = tmp_dict

        for event in rmats[key]:
            if event[0] not in p_w:
                tmp_dict = {
                    "type": key,
                    "event": event[0],
                    "whippet": math.nan,
                    "rmats": float(event[1]),
                }
                for j in Ws:
                    tmp_dict[f"pantas_{j}"] = math.nan
                data[event[0]] = tmp_dict
    df = pd.DataFrame(data.values())

    df.to_csv(f"{output}/res.csv", index=False)
    print(df)
    for w in Ws:
        p = f"pantas_{w}"
        g = sns.jointplot(
            data=df,
            x=p,
            y="rmats",
            hue="type",
            kind="scatter",
            xlim=(-1.05, 1.05),
            ylim=(-1.05, 1.05),
        )
        plt.tight_layout()
        plt.savefig(f"{output}/pantas2_{w}_rmats.png")
        plt.clf()
        g = sns.jointplot(
            data=df,
            x=p,
            y="whippet",
            hue="type",
            kind="scatter",
            xlim=(-1.05, 1.05),
            ylim=(-1.05, 1.05),
        )
        plt.tight_layout()
        plt.savefig(f"{output}/pantas2_{w}_whippet.png")
        plt.clf()

    for (w1, w2) in pairs(Ws):
        g = sns.jointplot(
            data=df,
            x=f"pantas_{w1}",
            y=f"pantas_{w2}",
            hue="type",
            kind="scatter",
            xlim=(-1.05, 1.05),
            ylim=(-1.05, 1.05),
        )
        plt.tight_layout()
        plt.savefig(f"{output}/pantas_{w1}_pantas2_{w2}.png")
        plt.clf()

    g = sns.jointplot(
        data=df,
        x="rmats",
        y="whippet",
        hue="type",
        kind="scatter",
        xlim=(-1.05, 1.05),
        ylim=(-1.05, 1.05),
    )
    plt.tight_layout()
    plt.savefig(f"{output}/rmats_whippet.png")
    plt.clf()

    df_mask = df.copy()
    for col in df_mask.columns:
        if col != "type" and col != "event":
            df_mask[col] = df_mask.apply(update_column, col_name=col, axis=1)
    print(df_mask)
    df_mask.to_csv(f"{output}/res_mask.csv", index=False)
    for e in ETYPES:
        tmp_df = df_mask[df_mask["type"] == e]
        rmats_set = set(tmp_df["rmats"])
        whippet_set = set(tmp_df["whippet"])
        for w in Ws:
            pantas_set = set(tmp_df[f"pantas_{w}"])
            venn3(
                [rmats_set, whippet_set, pantas_set],
                ("rmats", "whippet", f"pantas_{w}"),
            )
            plt.tight_layout()
            plt.savefig(f"{output}/venn_{e}_rmats_whippet_pantas_{w}.png")
            plt.clf()
        for (w1, w2) in pairs(Ws):
            pantas_set1 = set(tmp_df[f"pantas_{w1}"])
            pantas_set2 = set(tmp_df[f"pantas_{w2}"])
            venn2(
                [pantas_set1, pantas_set2],
                (f"pantas_{w1}", f"pantas_{w2}"),
            )
            plt.tight_layout()
            plt.savefig(f"{output}/venn_{e}_pantas_{w1}_pantas_{w2}.png")
            plt.clf()

    tmp_df = df_mask
    rmats_set = set(tmp_df["rmats"])
    whippet_set = set(tmp_df["whippet"])
    for w in Ws:
        pantas_set = set(tmp_df[f"pantas_{w}"])
        venn3(
            [rmats_set, whippet_set, pantas_set],
            ("rmats", "whippet", f"pantas_{w}"),
        )
        plt.tight_layout()
        plt.savefig(f"{output}/venn_full_rmats_whippet_pantas_{w}.png")
        plt.clf()
    for (w1, w2) in pairs(Ws):
        pantas_set1 = set(tmp_df[f"pantas_{w1}"])
        pantas_set2 = set(tmp_df[f"pantas_{w2}"])
        venn2(
            [pantas_set1, pantas_set2],
            (f"pantas_{w1}", f"pantas_{w2}"),
        )
        plt.tight_layout()
        plt.savefig(f"{output}/venn_full_pantas_{w1}_pantas_{w2}.png")
        plt.clf()


if __name__ == "__main__":
    main(sys.argv[1:])

