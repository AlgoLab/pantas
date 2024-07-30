#!/usr/bin/env python3

#!/usr/bin/env python3
import sys
import os
import itertools
import numpy as np
import math
import copy
from matplotlib import pyplot as plt
from venn import venn
import seaborn as sns
import matplotlib.gridspec as gridspec
import seaborn as sns; sns.set()
import SeabornFig2Grid as sfg
from scipy.stats import pearsonr
import pandas as pd
import eparser
import math
from matplotlib.patches import Rectangle
import matplotlib.colors 

ETYPES = ["ES", "IR", "A3", "A5"]
EMAP_RMATS = {"SE": "ES", "RI": "IR", "A3SS": "A3", "A5SS": "A5"}
EMAP_WHIPPET = {"CE": "ES", "RI": "IR", "AD": "A5", "AA": "A3"}

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
        del _e[2]
        del _e[5:11]
        e = eparser.EventPantas(*_e)
        if math.isnan(e.psi_c1) or math.isnan(e.psi_c2):
            continue
        if abs(e.dpsi) < FILTER:
            continue
        if e.etype == "IR":
            e.dpsi = 1 - e.dpsi
        event_pantas[e.etype].append(e)
    return event_pantas


def parse_rmats(path):
    filt = {x: {} for x in ETYPES}
    event_rmats = {x: [] for x in ETYPES}
    pos = 0
    for line in open(path, "r"):
        if line.startswith("etype"):
            continue
        line = line.strip()
        _e = line.split(",")
     
        #e = eparser.EventPantas(*_e)
        e = eparser.EventRmats(*_e)
        if math.isnan(e.psi_c1) or math.isnan(e.psi_c2):
            continue
        if abs(e.dpsi) < FILTER:
            continue
        if e not in filt[e.etype].keys():
            event_rmats[e.etype].append(e)
            filt[e.etype][e] = [len(event_rmats[e.etype])-1, float(e.pv)]
            pos = pos + 1
        else:
            if float(e.pv) < filt[e.etype][e][1]: 
                filt[e.etype][e][1] = float(e.pv)
                event_rmats[e.etype][filt[e][0]] = e
        
    return event_rmats

def parse_suppa(path):
    filt = {x: {} for x in ETYPES}
    event_suppa = {x: [] for x in ETYPES}
    pos = 0
    for line in open(path, "r"):
        if line.startswith("etype"):
            continue
        line = line.strip()
        _e = line.split(",")
        e = eparser.EventRmats(*_e)
        e.dpsi = -e.dpsi
        if abs(e.dpsi) < FILTER:
            continue
        if e not in filt[e.etype].keys():
            event_suppa[e.etype].append(e)
            filt[e.etype][e] = [len(event_suppa[e.etype])-1, float(e.pv)]
            pos = pos + 1
        else:
            if float(e.pv) < filt[e.etype][e][1]: 
                filt[e.etype][e][1] = float(e.pv)
                event_suppa[e.etype][filt[e][0]] = e
    return event_suppa

def parse_whippet(path):
    filt = {x: {} for x in ETYPES}
    event_whippet = {x: [] for x in ETYPES}
    pos = 0
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
        #event_whippet[e.etype].append(e)
        if e not in filt[e.etype].keys():
            event_whippet[e.etype].append(e)
            filt[e.etype][e] = [len(event_whippet[e.etype])-1, float(e.pv)]
            pos = pos + 1
        else:
            if float(e.pv) > filt[e.etype][e][1]: # here pv is min prob
                filt[e.etype][e][1] = float(e.pv)
                event_whippet[e.etype][filt[e][0]] = e
        
        
    return event_whippet


def main(argv):
    plt.rc('axes', labelsize=18) #fontsize of the x and y labels
    plt.rc('legend', fontsize=14) #fontsize of the legend
    pantas_path = snakemake.params.pantas_dir
    if pantas_path.endswith("/"):
        pantas_path = pantas_path[:-1]
    rmats_path = snakemake.input.r
    if rmats_path.endswith("/"):
        rmats_path = rmats_path[:-1]
    whippet_path = snakemake.input.w
    if whippet_path.endswith("/"):
        whippet_path = whippet_path[:-1]
    suppa_path = snakemake.input.s
    if suppa_path.endswith("/"):
        suppa_path = suppa_path[:-1]

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
        pantas[w] = parse_pantas(f"{pantas_path}/quant-remap.w{w}.csv")
    rmats = parse_rmats(rmats_path)
    suppa = parse_suppa(suppa_path)
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
                "SUPPA2": math.nan,
                "rMATS": math.nan,
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
                        "SUPPA2": math.nan,
                        "rMATS": math.nan,
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
                data[e_name][f"rMATS"] = event.dpsi
            else:
                tmp_dict = {
                    "type": key,
                    "event": e_name,
                    "whippet": math.nan,
                    "SUPPA2": math.nan,
                    "rMATS": event.dpsi,
                }
                for j in Ws:
                    tmp_dict[f"pantas_{j}"] = math.nan
                data[e_name] = tmp_dict
                p_d[key].append(event)
    for key in ["ES", "A3", "A5", "IR"]:
        for event in suppa[key]:
            e_name = (
                f"{event.etype}_{event.chrom}_{event.event_j[0]}_{event.event_j[1]}"
            )
            if e_name in data.keys():
                data[e_name][f"SUPPA2"] = event.dpsi
            else:
                tmp_dict = {
                    "type": key,
                    "event": e_name,
                    "whippet": math.nan,
                    "SUPPA2": event.dpsi,
                    "rMATS": math.nan,
                }
                for j in Ws:
                    tmp_dict[f"pantas_{j}"] = math.nan
                data[e_name] = tmp_dict
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
                    "SUPPA2": math.nan,
                    "rMATS": math.nan,
                }
                for j in Ws:
                    tmp_dict[f"pantas_{j}"] = math.nan
                data[e_name] = tmp_dict

    df = pd.DataFrame(data.values())
    df.to_csv(f"{output}/res.csv", index=False)
    df_mask = df.copy()
    for col in df_mask.columns:
        if col != "type" and col != "event":
            df_mask[col] = df_mask.apply(update_column, col_name=col, axis=1)
    df_mask.to_csv(f"{output}/res_mask.csv", index=False)
    df=df.dropna(how='any')
    #print(df)
    if len(Ws) == 1:
        fig = plt.figure(figsize=(13,8))
        gs = gridspec.GridSpec(2, 3)
        p = f"pantas_{Ws[0]}"
        sns.set(style="white", color_codes=True)
        g1 = sns.JointGrid(
            data=df,
            x=p,
            y="rMATS",
            #hue="type",
            #color="black",
            #kind="kde",
            xlim=(-1.05, 1.05),
            ylim=(-1.05, 1.05),
            #palette={"ES": "#c08c7d", "A3": "#9549ff", "A5": "#beae03", "IR": "#65d7cd"}
        )
        g1.plot(sns.scatterplot, sns.kdeplot, color="black")
        g1.ax_joint.set_xlabel('pantas')
        #sns.move_legend(g1.ax_joint, "upper left", title='Type')
        
        corr, _ = pearsonr(df[p], df["rMATS"])
        corr = round(corr, 3)
        g1.ax_joint.text(s=f"Pearson correlation: {corr:.3f}", x=-0.85, y=-1, fontsize=16)
        #g1.ax_joint.legend_.remove()
        sns.set(style="white", color_codes=True)
        g2 = sns.JointGrid(
            data=df,
            x=p,
            y="whippet",
            #hue="type",
            #kind="kde",
            #color="black",
            xlim=(-1.05, 1.05),
            ylim=(-1.05, 1.05),
            #palette={"ES": "#c08c7d", "A3": "#9549ff", "A5": "#beae03", "IR": "#65d7cd"}
        )
        g2.plot(sns.scatterplot, sns.kdeplot, color="black")
        g2.ax_joint.set_xlabel('pantas')
        corr, _ = pearsonr(df[p], df["whippet"])
        corr = round(corr, 3)
        g2.ax_joint.text(s=f"Pearson correlation: {corr:.3f}", x=-0.85, y=-1, fontsize=16)
        #g2.ax_joint.legend_.remove()
        sns.set(style="white", color_codes=True)
        g3 = sns.JointGrid(
            data=df,
            x=p,
            y="SUPPA2",
            #hue="type",
            #kind="kde",
            #color="black",
            xlim=(-1.05, 1.05),
            ylim=(-1.05, 1.05),
            #palette={"ES": "#c08c7d", "A3": "#9549ff", "A5": "#beae03", "IR": "#65d7cd"}
        )
        g3.plot(sns.scatterplot, sns.kdeplot, color="black")
        g3.ax_joint.set_xlabel('pantas')
        corr, _ = pearsonr(df[p], df["SUPPA2"])
        corr = round(corr, 3)
        g3.ax_joint.text(s=f"Pearson correlation: {corr:.3f}", x=-0.85, y=-1, fontsize=16)
        #g3.ax_joint.legend_.remove()
        sns.set(style="white", color_codes=True)
        g4 = sns.JointGrid(
            data=df,
            x="rMATS",
            y="whippet",
            #hue="type",
            #kind="kde",
            #color="black",
            xlim=(-1.05, 1.05),
            ylim=(-1.05, 1.05),
            #palette={"ES": "#c08c7d", "A3": "#9549ff", "A5": "#beae03", "IR": "#65d7cd"}
        )
        g4.plot(sns.scatterplot, sns.kdeplot, color="black")
        corr, _ = pearsonr(df["rMATS"], df["whippet"])
        corr = round(corr, 3)
        g4.ax_joint.text(s=f"Pearson correlation: {corr:.3f}", x=-0.85, y=-1, fontsize=16)
        #g4.ax_joint.legend_.remove()
        sns.set(style="white", color_codes=True)
        g5 = sns.JointGrid(
            data=df,
            x="rMATS",
            y="SUPPA2",
            #hue="type",
            #kind="kde",
            #color="black",
            xlim=(-1.05, 1.05),
            ylim=(-1.05, 1.05),
            #palette={"ES": "#c08c7d", "A3": "#9549ff", "A5": "#beae03", "IR": "#65d7cd"}
        )
        corr, _ = pearsonr(df["rMATS"], df["SUPPA2"])
        corr = round(corr, 3)
        g5.ax_joint.text(s=f"Pearson correlation: {corr:.3f}", x=-0.85, y=-1, fontsize=16)
        #g5.ax_joint.legend_.remove()
        g5.plot(sns.scatterplot, sns.kdeplot, color="black")
        sns.set(style="white", color_codes=True)
        g6 = sns.JointGrid(
            data=df,
            x="whippet",
            y="SUPPA2",
            #hue="type",
            #kind="scatter",
            #kind="kde",
            #color="black",
            xlim=(-1.05, 1.05),
            ylim=(-1.05, 1.05),
            #palette={"ES": "#c08c7d", "A3": "#9549ff", "A5": "#beae03", "IR": "#65d7cd"}
        )
        g6.plot(sns.scatterplot, sns.kdeplot, color="black")
        corr, _ = pearsonr(df["whippet"], df["SUPPA2"])
        corr = round(corr, 3)
        g6.ax_joint.text(s=f"Pearson correlation: {corr:.3f}", x=-0.85, y=-1, fontsize=16)
        #g6.ax_joint.legend_.remove()
        
        fig = plt.figure(figsize=(15,10))
        gs = gridspec.GridSpec(2, 3)
        mg0 = sfg.SeabornFig2Grid(g1, fig, gs[0])
        mg1 = sfg.SeabornFig2Grid(g2, fig, gs[1])
        mg2 = sfg.SeabornFig2Grid(g3, fig, gs[2])
        mg3 = sfg.SeabornFig2Grid(g4, fig, gs[3])
        mg4 = sfg.SeabornFig2Grid(g5, fig, gs[4])
        mg5 = sfg.SeabornFig2Grid(g6, fig, gs[5])
        gs.tight_layout(fig)
        plt.savefig(f"{output}/full_corr.png",bbox_inches='tight') 
        
        
    for w in Ws:
        p = f"pantas_{w}"
        g = sns.jointplot(
            data=df,
            x=p,
            y="rMATS",
            hue="type",
            kind="scatter",
            xlim=(-1.05, 1.05),
            ylim=(-1.05, 1.05),
        )
        corr, _ = pearsonr(df[p], df["rMATS"])
        corr = round(corr, 3)
        plt.text(s=f"Pearson correlation: {corr}", x=-0.3, y=-1)
        plt.tight_layout()
        plt.savefig(f"{output}/corr_pantas2_{w}_rmats.png")
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
        corr, _ = pearsonr(df[p], df["whippet"])
        corr = round(corr, 3)
        plt.text(s=f"Pearson correlation: {corr}", x=-0.3, y=-1)
        plt.tight_layout()
        plt.savefig(f"{output}/corr_pantas2_{w}_whippet.png")
        plt.clf()
        
        g = sns.jointplot(
            data=df,
            x=p,
            y="SUPPA2",
            hue="type",
            kind="scatter",
            xlim=(-1.05, 1.05),
            ylim=(-1.05, 1.05),
        )
        corr, _ = pearsonr(df[p], df["SUPPA2"])
        corr = round(corr, 3)
        plt.text(s=f"Pearson correlation: {corr}", x=-0.3, y=-1)
        plt.tight_layout()
        plt.savefig(f"{output}/corr_pantas2_{w}_suppa.png")
        plt.clf()

    if len(Ws) > 1:
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
            corr, _ = pearsonr(df[f"pantas_{w1}"], df[f"pantas_{w2}"])
            corr = round(corr, 3)
            plt.text(s=f"Pearson correlation: {corr}", x=-0.3, y=-1)
            plt.tight_layout()
            plt.savefig(f"{output}/corr_pantas_{w1}_pantas2_{w2}.png")
            plt.clf()

    g = sns.jointplot(
        data=df,
        x="rMATS",
        y="whippet",
        hue="type",
        kind="scatter",
        xlim=(-1.05, 1.05),
        ylim=(-1.05, 1.05),
    )
    corr, _ = pearsonr(df["rMATS"], df["whippet"])
    corr = round(corr, 3)
    plt.text(s=f"Pearson correlation: {corr}", x=-0.3, y=-1)
    plt.tight_layout()
    plt.savefig(f"{output}/corr_rmats_whippet.png")
    plt.clf()
    
    g = sns.jointplot(
        data=df,
        x="rMATS",
        y="SUPPA2",
        hue="type",
        kind="scatter",
        xlim=(-1.05, 1.05),
        ylim=(-1.05, 1.05),
    )
    corr, _ = pearsonr(df["rMATS"], df["SUPPA2"])
    corr = round(corr, 3)
    plt.text(s=f"Pearson correlation: {corr}", x=-0.3, y=-1)
    plt.tight_layout()
    plt.savefig(f"{output}/corr_rmats_suppa.png")
    plt.clf()
    
    g = sns.jointplot(
        data=df,
        x="whippet",
        y="SUPPA2",
        hue="type",
        kind="scatter",
        xlim=(-1.05, 1.05),
        ylim=(-1.05, 1.05),
    )
    corr, _ = pearsonr(df["whippet"], df["SUPPA2"])
    corr = round(corr, 3)
    plt.text(s=f"Pearson correlation: {corr}", x=-0.3, y=-1)
    plt.tight_layout()
    plt.savefig(f"{output}/corr_whippet_suppa.png")
    plt.clf()
    
    print(df)
    for e in ETYPES:
        print(e)
        tmp_df = df[df["type"] == e]
        #print(tmp_df)
        for w in Ws:
            p = f"pantas_{w}"
            g = sns.jointplot(
                data=tmp_df,
                x=p,
                y="rMATS",
                hue="type",
                kind="scatter",
                xlim=(-1.05, 1.05),
                ylim=(-1.05, 1.05),
            )
           #print(tmp_df[p], tmp_df["rMATS"])
            corr, _ = pearsonr(tmp_df[p], tmp_df["rMATS"])
            corr = round(corr, 3)
            plt.text(s=f"Pearson correlation: {corr}", x=-0.3, y=-1)
            plt.tight_layout()
            plt.savefig(f"{output}/corr_{e}_pantas2_{w}_rmats.png")
            plt.clf()
            
            g = sns.jointplot(
                data=tmp_df,
                x=p,
                y="whippet",
                hue="type",
                kind="scatter",
                xlim=(-1.05, 1.05),
                ylim=(-1.05, 1.05),
            )
            corr, _ = pearsonr(tmp_df[p], tmp_df["whippet"])
            corr = round(corr, 3)
            plt.text(s=f"Pearson correlation: {corr}", x=-0.3, y=-1)
            plt.tight_layout()
            plt.savefig(f"{output}/corr_{e}_pantas2_{w}_whippet.png")
            plt.clf()
            
            g = sns.jointplot(
                data=tmp_df,
                x=p,
                y="SUPPA2",
                hue="type",
                kind="scatter",
                xlim=(-1.05, 1.05),
                ylim=(-1.05, 1.05),
            )
            corr, _ = pearsonr(tmp_df[p], tmp_df["SUPPA2"])
            corr = round(corr, 3)
            plt.text(s=f"Pearson correlation: {corr}", x=-0.3, y=-1)
            plt.tight_layout()
            plt.savefig(f"{output}/corr_{e}_pantas2_{w}_suppa.png")
            plt.clf()
        
        if len(Ws) > 1:
            for (w1, w2) in pairs(Ws):
                g = sns.jointplot(
                    data=tmp_df,
                    x=f"pantas_{w1}",
                    y=f"pantas_{w2}",
                    hue="type",
                    kind="scatter",
                    xlim=(-1.05, 1.05),
                    ylim=(-1.05, 1.05),
                )
                corr, _ = pearsonr(tmp_df[f"pantas_{w1}"], tmp_df[f"pantas_{w2}"])
                corr = round(corr, 3)
                plt.text(s=f"Pearson correlation: {corr}", x=-0.3, y=-1)
                plt.tight_layout()
                plt.savefig(f"{output}/corr_{e}_pantas_{w1}_pantas2_{w2}.png")
                plt.clf()

        g = sns.jointplot(
            data=tmp_df,
            x="rMATS",
            y="whippet",
            hue="type",
            kind="scatter",
            xlim=(-1.05, 1.05),
            ylim=(-1.05, 1.05),
        )
        corr, _ = pearsonr(tmp_df["rMATS"], tmp_df["whippet"])
        corr = round(corr, 3)
        plt.text(s=f"Pearson correlation: {corr}", x=-0.3, y=-1)
        plt.tight_layout()
        plt.savefig(f"{output}/corr_{e}_rmats_whippet.png")
        plt.clf()
        
        g = sns.jointplot(
            data=tmp_df,
            x="rMATS",
            y="SUPPA2",
            hue="type",
            kind="scatter",
            xlim=(-1.05, 1.05),
            ylim=(-1.05, 1.05),
        )
        corr, _ = pearsonr(tmp_df["rMATS"], tmp_df["SUPPA2"])
        corr = round(corr, 3)
        plt.text(s=f"Pearson correlation: {corr}", x=-0.3, y=-1)
        plt.tight_layout()
        plt.savefig(f"{output}/corr_{e}_rmats_suppa.png")
        plt.clf()
        g = sns.jointplot(
            data=tmp_df,
            x="whippet",
            y="SUPPA2",
            hue="type",
            kind="scatter",
            xlim=(-1.05, 1.05),
            ylim=(-1.05, 1.05),
        )
        corr, _ = pearsonr(tmp_df["whippet"], tmp_df["SUPPA2"])
        corr = round(corr, 3)
        plt.text(s=f"Pearson correlation: {corr}", x=-0.3, y=-1)
        plt.tight_layout()
        plt.savefig(f"{output}/corr_{e}_whippet_suppa.png")
        plt.clf()

    
    for e in ETYPES:
        tmp_df = df_mask[df_mask["type"] == e]
        rmats_set = set(tmp_df["rMATS"])
        whippet_set = set(tmp_df["whippet"])
        suppa_set = set(tmp_df["SUPPA2"])
        for w in Ws:
            pantas_set = set(tmp_df[f"pantas_{w}"])
            dic_data = {
                "rMATS": rmats_set,
                "whippet": whippet_set,
                "SUPPA2": suppa_set,
                f"pantas_{w}": pantas_set,
            }
            venn(dic_data)
            plt.tight_layout()
            plt.savefig(f"{output}/venn_{e}_rmats_whippet_suppa_pantas_{w}.png")
            plt.clf()
        if len(Ws) > 1:
            pantas_set_full = {}
            for key in key_names:
                pantas_set_full[key] = set(tmp_df[key])

            venn(pantas_set_full)
            plt.tight_layout()
            plt.savefig(f"{output}/venn_{e}_pantas.png")
            plt.clf()
        # for (w1, w2) in pairs(Ws):
        #     pantas_set1 = set(tmp_df[f"pantas_{w1}"])
        #     pantas_set2 = set(tmp_df[f"pantas_{w2}"])
        #     dic_data = {
        #         f"pantas_{w1}": pantas_set1,
        #         f"pantas_{w2}": pantas_set2,
        #     }
        #     venn(dic_data)

        #     plt.tight_layout()
        #     plt.savefig(f"{output}/venn_{e}_pantas_{w1}_pantas_{w2}.png")
        #     plt.clf()

    tmp_df = df_mask
    rmats_set = set(tmp_df["rMATS"].dropna())
    whippet_set = set(tmp_df["whippet"].dropna())
    suppa_set = set(tmp_df["SUPPA2"].dropna())
    if len(Ws) == 1: 
        pantas_set = set(tmp_df[f"pantas_{Ws[0]}"].dropna())
       	legends = []
        for t in [f"pantas_{Ws[0]}", "rMATS", "whippet", "SUPPA2"]:
	        n = tmp_df[t].count()
	        if t == f"pantas_{Ws[0]}":
	            t = "pantas"
	        legends.append(f"{t}: {n}")

        custom_lines = [
	        Rectangle(
	            (0, 0),
	            1,
	            1,
	            facecolor=sns.color_palette()[0],
	            linewidth=1,
	            edgecolor="black",
	        ),
	        Rectangle(
	            (0, 0),
	            1,
	            1,
	            facecolor=sns.color_palette()[1],
	            linewidth=1,
	            edgecolor="black",
	        ),
	        Rectangle(
	            (0, 0),
	            1,
	            1,
	            facecolor=sns.color_palette()[2],
	            linewidth=1,
	            edgecolor="black",
	        ),
	        Rectangle(
	            (0, 0),
	            1,
	            1,
	            facecolor=sns.color_palette()[3],
	            linewidth=1,
	            edgecolor="black",
	        ),
        ]
        print(len(pantas_set))
        dic_data = {
            "pantas": pantas_set,
            "rMATS": rmats_set,
            "whippet": whippet_set,
            "SUPPA2": suppa_set,      
        }
        fig, ax1 = plt.subplots(1, 1, figsize=(5,5))
        venn(
            dic_data,
            fontsize=13,
            legend_loc=None,  # "upper left",
            cmap=sns.color_palette(),
            ax=ax1,
        )
        ax1.legend(
            custom_lines,
            legends,
            title="Tool: #Events",
            loc="lower center",
            bbox_to_anchor=(0.5, -0.1),
            ncol=2,
        )
        plt.tight_layout()
        plt.savefig(f"{output}/venn_full_rmats_whippet_suppa_pantas.png",bbox_inches='tight')
        plt.clf()
        
    for w in Ws:
        pantas_set = set(tmp_df[f"pantas_{w}"])
        dic_data = {
            "rMATS": rmats_set,
            "whippet": whippet_set,
            "SUPPA2": suppa_set,
            f"pantas_{w}": pantas_set,
        }
        venn(dic_data)
        plt.tight_layout()
        plt.savefig(f"{output}/venn_full_rmats_whippet_suppa_pantas_{w}.png")
        plt.clf()

    if len(Ws) > 1:
        pantas_set_full = {}
        for key in key_names:
            pantas_set_full[key] = set(tmp_df[key])
        venn(pantas_set_full)
        plt.tight_layout()
        plt.savefig(f"{output}/venn_full_pantas.png")
        plt.clf()


if __name__ == "__main__":
    main(sys.argv[1:])
