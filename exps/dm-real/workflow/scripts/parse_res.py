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
import eparser
import math

ETYPES = ["ES", "IR", "A3", "A5"]
EMAP_RMATS = {"SE": "ES", "RI": "IR", "A3SS": "A3", "A5SS": "A5"}
EMAP_WHIPPET = {"CE": "ES", "RI": "IR", "AD": "A5", "AA": "A3"}
## TODO metti come variabili snakemake
FILTER = snakemake.params.min_dpsi
MIN_EVENT_COV = snakemake.params.min_coverage


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
        for c in whippet[e]:
            hit = False
            for t in truth[e]:
                hit = eparser.eq_event(c, t, relax)
                if hit:
                    t_name = f"{t.etype}_{t.chrom}_{t.event_j[0]}_{t.event_j[1]}"
                    hits[e][t_name] = c.dpsi
                    c_name = f"{c.etype}_{c.chrom}_{c.event_j[0]}_{c.event_j[1]}"
                    p.append(c_name)
    return hits, p


def parse_pantas(path):
    event_pantas = {x: [] for x in ETYPES}
    for line in open(path, "r"):
        if line.startswith("etype"):
            continue
        line = line.strip()
        _e = line.split(",")
        e = eparser.EventPantas(*_e)
        if math.isnan(e.psi_c1) or math.isnan(e.psi_c2):
            continue
        if abs(e.dpsi) < FILTER:
            continue
        event_pantas[e.etype].append(e)
    return event_pantas


def parse_rmats(path):
    event_rmats = {x: [] for x in ETYPES}
    for line in open(path, "r"):
        if line.startswith("etype"):
            continue
        line = line.strip()
        _e = line.split(",")
        e = eparser.EventPantas(*_e)
        if math.isnan(e.psi_c1) or math.isnan(e.psi_c2):
            continue
        if abs(e.dpsi) < FILTER:
            continue
        event_rmats[e.etype].append(e)
    return event_rmats


def parse_whippet(path):
    event_whippet = {x: [] for x in ETYPES}
    for line in open(path, "r"):
        if line.startswith("Gene"):
            # header
            continue
        line = line.strip()
        _e = line.split("\t")
        _e[4] = EMAP_WHIPPET.get(_e[4], _e[4])
        if _e[4] not in ETYPES:
            continue
        e = eparser.EventWhippet(*_e, "anno")
        if math.isnan(e.psi_c1) or math.isnan(e.psi_c2):
            continue
        if abs(e.dpsi) < FILTER:
            continue
        event_whippet[e.etype].append(e)
    return event_whippet


def main(argv):
    pantas_path = snakemake.params.pantas_dir
    if pantas_path.endswith("/"):
        pantas_path = pantas_path[:-1]
    rmats_path = snakemake.input.r
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
    rmats = parse_rmats(rmats_path)
    whippet = parse_whippet(whippet_path)
    data = {}
    pantas_keys = list(pantas.keys())
    p_d = copy.deepcopy(pantas[Ws[0]])
    key_names = [f"pantas_{w}" for w in Ws]
    for key in ["ES", "A3", "A5", "IR"]:
        for event in pantas[Ws[0]][key]:
            # print(event)
            e_name = (
                f"{event.etype}_{event.chrom}_{event.event_j[0]}_{event.event_j[1]}"
            )
            data[e_name] = {
                "type": key,
                "event": e_name,
                "whippet": math.nan,
                "rmats": math.nan,
                f"pantas_{Ws[0]}": event.dpsi,
            }
    for i, w in enumerate(Ws[1:]):
        # print(w)
        for key in ["ES", "A3", "A5", "IR"]:
            for event in pantas[w][key]:
                # print(event[0])
                e_name = (
                    f"{event.etype}_{event.chrom}_{event.event_j[0]}_{event.event_j[1]}"
                )
                if e_name in data.keys():
                    data[e_name][f"pantas_{w}"] = event.dpsi
                else:
                    tmp_dict = {
                        "type": key,
                        "event": e_name,
                        "whippet": math.nan,
                        "rmats": math.nan,
                        f"pantas_{w}": event.dpsi,
                    }
                    for j in Ws[0 : i + 1]:
                        tmp_dict[f"pantas_{j}"] = math.nan
                    data[e_name] = tmp_dict
                    p_d[key].append(event)

    for key in ["ES", "A3", "A5", "IR"]:
        for event in rmats[key]:
            e_name = (
                f"{event.etype}_{event.chrom}_{event.event_j[0]}_{event.event_j[1]}"
            )
            if e_name in data.keys():
                data[e_name][f"rmats"] = event.dpsi
            else:
                tmp_dict = {
                    "type": key,
                    "event": e_name,
                    "whippet": math.nan,
                    "rmats": event.dpsi,
                }
                for j in Ws:
                    tmp_dict[f"pantas_{j}"] = math.nan
                data[event[0]] = tmp_dict
                p_d[key].append(event)

    mask_whippet, p_w = check_whippet(whippet, p_d, relax)
    for event, row in data.items():
        e = row["type"]
        if event in mask_whippet[e].keys():
            data[event]["whippet"] = mask_whippet[e][event]
    for key in ["ES", "A3", "A5", "IR"]:
        for event in whippet[key]:
            e_name = (
                f"{event.etype}_{event.chrom}_{event.event_j[0]}_{event.event_j[1]}"
            )
            if e_name not in p_w:
                tmp_dict = {
                    "type": key,
                    "event": e_name,
                    "whippet": event.dpsi,
                    "rmats": math.nan,
                }
                for j in Ws:
                    tmp_dict[f"pantas_{j}"] = math.nan
                data[e_name] = tmp_dict

    df = pd.DataFrame(data.values())
    df.to_csv(f"{output}/res.csv", index=False)
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
