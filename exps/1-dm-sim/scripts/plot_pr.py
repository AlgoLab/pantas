import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sns.set()

plt.rcParams.update({"font.size": 15})
# sns.set_style("whitegrid")
colors = [sns.color_palette("bright")[1]] + [
    sns.color_palette("dark")[i] for i in [0, 2, 6]
]


def main():
    fpath = sys.argv[1]
    data = []
    truth = {}
    for line in open(fpath):
        if line.startswith("p-supp"):
            continue
        w, tool, etype, dpsi, c, tp, fn, fp, p, r, f1, *tot = line.strip("\n").split(
            ","
        )
        c = int(c)
        if c in [0, 2]:
            continue
        tp, fn = int(tp), int(fn)
        if c not in truth:
            truth[c] = {}
        if etype not in truth[c]:
            truth[c][etype] = tp + fn
        else:
            assert truth[c][etype] == tp + fn
        data.append([tool, etype, c, float(p), float(r)])

    etypes = ["ES", "A3", "A5", "IR"]  # sorted(truth[0].keys())
    print("True Support (ω)", ",".join(etypes), sep=",")
    for c in truth:
        print(c, end="")
        for etype in etypes:
            print("," + str(truth[c][etype]), end="")
        print("")

    df = pd.DataFrame(
        data, columns=["Tool", "E.Type", "True Support (ω)", "Precision", "Recall"]
    )

    fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(7, 7))
    axes = axes.flatten()
    for ax, etype in zip(axes, ["ES", "IR", "A3", "A5"]):
        sns.lineplot(
            df[df["E.Type"] == etype],
            x="Precision",
            y="Recall",
            hue="Tool",
            palette=colors,
            legend=None,
            estimator=None,
            linewidth=2,
            alpha=0.2,
            sort=False,
            ax=ax,
        )
        if etype == "IR":
            sns.lineplot(
                df[df["E.Type"] == etype],
                x="Precision",
                y="Recall",
                hue="Tool",
                palette=colors,
                style="True Support (ω)",
                # legend=None,
                markers=True,
                dashes=False,
                markersize=13,
                ax=ax,
            )
        else:
            sns.lineplot(
                df[df["E.Type"] == etype],
                x="Precision",
                y="Recall",
                hue="Tool",
                palette=colors,
                style="True Support (ω)",
                legend=None,
                markers=True,
                dashes=False,
                markersize=13,
                ax=ax,
            )
        ax.set_title(etype)
        ax.set_xlim(0, 1.01)
        ax.set_ylim(0, 1.01)

    # plt.tight_layout()
    plt.subplots_adjust(
        bottom=0.07, right=0.99, top=0.95, left=0.05, wspace=0.03, hspace=0.1
    )
    plt.show()  # savefig(fpath + ".png")


if __name__ == "__main__":
    main()
