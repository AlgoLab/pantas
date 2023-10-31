import sys
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from matplotlib.patches import Rectangle
from venn import venn

sns.set()


font = {"size": 13}
matplotlib.rc("font", **font)


def get_interval(region):
    return [int(x) for x in region.split(":")[1].split("-")]


def parse_pantas(fpath):
    events = {}
    for line in open(fpath):
        if line.startswith("etype"):
            continue
        etype, novel, chrom, _, strand, _, i1, i2, _, _, psi1, psi2, dpsi = line.strip(
            "\n"
        ).split(",")
        if etype != "ES":
            continue
        dpsi = float(dpsi)
        s1, e1 = get_interval(i1)
        s2, e2 = get_interval(i2)
        k = f"{chrom}:{e1+1}-{s2-1}"
        events[k] = (
            events[k] + [(-float(dpsi), novel)]
            if k in events
            else [(-float(dpsi), novel)]
        )
    return events


def parse_rmats(fpath):
    events = {}
    for line in open(fpath):
        if line.startswith("ID"):
            continue
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
        delta_incl = float(delta_incl)
        pv = float(pv)
        ex_s, ex_e = int(ex_s), int(ex_e)
        k = f"{chrom}:{ex_s+1}-{ex_e}"
        events[k] = (
            events[k] + [(-delta_incl, float(pv))]
            if k in events
            else [(-delta_incl, float(pv))]
        )
    return events


def parse_whippet(fpath):
    events = {}
    for line in open(fpath):
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
            p,
            compl,
            entr,
        ) = line.strip("\t \n").split("\t")
        dpsi = float(dpsi)
        p = float(p)
        if etype != "CE":
            continue
        events[region] = (
            events[region] + [(-dpsi, p)] if region in events else [(-dpsi, p)]
        )
    return events


def parse_suppa(suppa_dpsi):
    events = {}
    for i, line in enumerate(open(suppa_dpsi)):
        if i == 0:
            continue
        idx, dpsi, pvalue = line.strip("\n").split("\t")
        dpsi, pvalue = float(dpsi), float(pvalue)
        gene, rest = idx.split(";")
        etype, chrom, *positions, strand = rest.split(":")
        if not chrom.startswith("chr"):
            chrom = "chr" + chrom
        if etype == "SE":
            ab, cd = positions
            intron1 = tuple(int(x) for x in ab.split("-"))
            intron2 = tuple(int(x) for x in cd.split("-"))
            k = f"{chrom}:{intron1[1]}-{intron2[0]}"
            # if k in events and pvalue > events[k][1]:
            #     continue
            events[k] = (
                events[k] + [(dpsi, pvalue)] if k in events else [(dpsi, pvalue)]
            )
    return events


def parse_truth(fpath):
    truth_pos = {}
    truth_neg = set()
    for line in open(fpath):
        if line.startswith("POS"):
            _, chrom, exs, exe, gidx, dpsi = line.strip("\n").split("\t")
            truth_pos[f"{chrom}:{exs}-{exe}"] = float(dpsi)
        else:
            _, chrom, exs, exe, gidx = line.strip("\n").split("\t")
            truth_neg.add(f"{chrom}:{exs}-{exe}")
    return truth_pos, truth_neg


def main(args):
    # Input parsing
    truth, negatives = parse_truth(args.TRUTH)
    events = {}
    events["pantas"] = parse_pantas(args.PANTAS)
    events["rMATS"] = parse_rmats(args.RMATS)
    events["whippet"] = parse_whippet(args.WHIPPET)
    events["SUPPA2"] = parse_suppa(args.SUPPA)

    # Filtering and dataframe preparation
    print("Truth:", len(truth))
    truth = {
        k: v
        for k, v in truth.items()
        if abs(v) >= args.delta and abs(v) <= 1 - args.delta
    }
    print(f"Filtered truth with delta={args.delta}:", len(truth))
    df = []
    df_neg = []
    for t, Es in events.items():
        TPs = set(Es.keys()) & set(truth.keys())
        for k in TPs:
            best_dpsi = -1
            best_conf = -1
            best_diff = 2
            for dpsi, conf in Es[k]:
                if abs(dpsi) < args.delta or abs(dpsi) > 1 - args.delta:
                    continue
                if t == "pantas":
                    pass
                elif t == "rMATS":
                    if conf > args.pvalue:
                        continue
                elif t == "whippet":
                    if conf < args.prob:
                        continue
                elif t == "SUPPA2":
                    if conf > args.pvalue:
                        continue
                if dpsi - truth[k] > best_diff:
                    continue
                best_dpsi = dpsi
                best_conf = conf
                best_diff = dpsi - truth[k]
            if best_diff == 2:
                continue
            df.append(
                [
                    t,
                    k,
                    best_dpsi,
                    best_conf,
                    truth[k],
                    abs(best_dpsi - truth[k]),
                ]
            )
        FPs = set(Es.keys()) & set(negatives)
        for k in FPs:
            add_flag = False
            for dpsi, conf in Es[k]:
                if abs(dpsi) < args.delta or abs(dpsi) > 1 - args.delta:
                    continue
                if t == "pantas":
                    pass
                elif t == "rMATS":
                    if conf > args.pvalue:
                        continue
                elif t == "whippet":
                    if conf < args.prob:
                        continue
                elif t == "SUPPA2":
                    if conf > args.pvalue:
                        continue
                add_flag = True
                break
            if add_flag:
                df_neg.append(
                    [
                        t,
                        k,
                        dpsi,
                        conf,
                    ]
                )
    df = pd.DataFrame(df, columns=["Tool", "Event", "dPSI", "P", "RTPCR", "X"])
    df_neg = pd.DataFrame(df_neg, columns=["Tool", "Event", "dPSI", "P"])

    # Sets for each tool (for easier intersection)
    pantas = set(df[df["Tool"] == "pantas"]["Event"])
    rmats = set(df[df["Tool"] == "rMATS"]["Event"])
    whippet = set(df[df["Tool"] == "whippet"]["Event"])
    suppa2 = set(df[df["Tool"] == "SUPPA2"]["Event"])
    pantas_neg = set(df_neg[df_neg["Tool"] == "pantas"]["Event"])
    rmats_neg = set(df_neg[df_neg["Tool"] == "rMATS"]["Event"])
    whippet_neg = set(df_neg[df_neg["Tool"] == "whippet"]["Event"])
    suppa2_neg = set(df_neg[df_neg["Tool"] == "SUPPA2"]["Event"])
    pantas_all = set(events["pantas"].keys())
    rmats_all = set(events["rMATS"].keys())
    whippet_all = set(events["whippet"].keys())
    suppa2_all = set(events["SUPPA2"].keys())

    # Negative results
    print(
        "pantas",
        len(pantas_neg & negatives),
        "/",
        len(pantas_all & negatives),
        len(pantas_neg & negatives) / len(pantas_all & negatives),
    )
    print(
        "rMATS",
        len(rmats_neg & negatives),
        "/",
        len(rmats_all & negatives),
        len(rmats_neg & negatives) / len(rmats_all & negatives),
    )
    print(
        "whippet",
        len(whippet_neg & negatives),
        "/",
        len(whippet_all & negatives),
        len(whippet_neg & negatives) / len(whippet_all & negatives),
    )
    print(
        "SUPPA2",
        len(suppa2_neg & negatives),
        "/",
        len(suppa2_all & negatives),
        len(suppa2_neg & negatives) / len(suppa2_all & negatives),
    )

    # Do we want only common events?
    if args.common:
        shared = pantas & rmats & whippet & suppa2
        df = df[df["Event"].isin(shared)]
        pantas &= shared
        rmats &= shared
        whippet &= shared
        suppa2 &= shared

    # Compute correlation and prepare legend
    legends = []
    xticks = []
    for t in events:
        n = len(df[df["Tool"] == t])
        corr, _ = pearsonr(
            df[df["Tool"] == t].sort_values(by="Event")["RTPCR"],
            df[df["Tool"] == t].sort_values(by="Event")["dPSI"],
        )
        print(df[df["Tool"] == t]["X"].describe())
        corr = round(corr, 3)
        # legends.append(f"{t}: {n} events (P={corr:.3f})")
        # legends.append(f"{t}: {n} ({corr:.3f})")
        legends.append(f"{t}: {n}")
        xticks.append(f"{t}\n(r={corr:.3f})")

    # Print events not found by pantas + Some other stuff
    print((rmats | whippet | suppa2) - pantas)
    for k in (suppa2 | whippet | rmats) - pantas:
        if k in whippet:
            print(k, events["whippet"][k])
        elif k in rmats:
            print(k, events["rMATS"][k])
        else:
            print(k, events["SUPPA2"][k])

    print("All:", len((rmats | whippet | suppa2 | pantas) & set(truth.keys())))

    # Start plotting
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 7))

    # First plot: ECDF with legend
    # sns.ecdfplot(data=df, x="X", hue="Tool", ax=ax1)
    # ax1.set_xlabel("|ΔPSI - RTPCR|")
    # ax1.set_ylabel("Cumulative Proportion")

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
    # ax1.set_xlim(-0.01, 1.01)
    # ax1.set_ylim(-0.01, 1.01)

    # Second plot: VENN of TPs
    # ins_ax = ax1.inset_axes([0.2, -0.09, 0.85, 0.85])
    venn_dict = {"pantas": pantas, "rMATS": rmats, "whippet": whippet, "SUPPA2": suppa2}
    venn(
        venn_dict,
        fontsize=13,
        legend_loc=None,  # "upper left",
        cmap=sns.color_palette(),
        ax=ax1,
    )
    # ax1.get_legend().remove()
    ax1.legend(
        custom_lines,
        legends,
        title="Tool: #Events (Pearson)",
        loc="lower center",
        bbox_to_anchor=(0.5, -0.1),
        ncol=2,
    )

    # Third plot: box + strip
    sns.boxplot(
        data=df,
        x="Tool",
        y="X",
        hue="Tool",
        fliersize=7,
        flierprops={"marker": "x"},
        showfliers=False,
        # notch=True,
        # showcaps=False,
        # fill=False,
        ax=ax2,
    )
    sns.stripplot(
        data=df, x="Tool", y="X", hue="Tool", linewidth=1, edgecolor="black", ax=ax2
    )
    ax2.set_xticklabels(xticks)
    ax2.set_ylabel("|ΔPSI - RTPCR|")
    ax2.set_ylim(-0.01, 0.7)

    # Set subplots titles
    ax1.set_title("(a)")  # Number of events reported by each tool")
    ax2.set_title("(b)")  # Difference between predicted and RT-PCR ∆ψ")
    # ax3.set_title("(c)")

    # Plot
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        prog="",
        description="",
    )
    parser.add_argument("TRUTH", help="")
    parser.add_argument("PANTAS", help="")
    parser.add_argument("RMATS", help="")
    parser.add_argument("WHIPPET", help="")
    parser.add_argument("SUPPA", help="")

    parser.add_argument(
        "--strict",
        help="Use strict filtering",
        dest="strict",
        action="store_true",
        required=False,
    )
    parser.add_argument(
        "--common",
        help="Analyze only AS events in the intersection of the 4 tools",
        dest="common",
        action="store_true",
        required=False,
    )
    parser.add_argument(
        "-d",
        "--delta",
        help="DeltaPSI filtering value",
        dest="delta",
        type=float,
        default=-1,
        required=False,
    )
    parser.add_argument(
        "-v",
        "--pvalue",
        help="pvalue filtering value",
        dest="pvalue",
        type=float,
        default=-1,
        required=False,
    )
    parser.add_argument(
        "-p",
        "--prob",
        help="Probability filtering value",
        dest="prob",
        type=float,
        default=-1,
        required=False,
    )
    args = parser.parse_args()
    if args.strict:
        if args.delta == -1:
            args.delta = 0.05
        if args.prob == -1:
            args.prob = 0.9
        if args.pvalue == -1:
            args.pvalue = 0.05
    else:
        if args.delta == -1:
            args.delta = 0
        if args.pvalue == -1:
            args.pvalue = 2
        if args.prob == -1:
            args.prob = -1
    main(args)
