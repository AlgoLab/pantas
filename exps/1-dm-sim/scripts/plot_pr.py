import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# plt.rcParams.update({"font.size": 17})
sns.set_style("whitegrid")
colors = [sns.color_palette("bright")[1]] + [
    sns.color_palette("dark")[i] for i in [0, 2, 6]
]


def main():
    fpath = sys.argv[1]
    data = []
    for line in open(fpath):
        if line.startswith("p-supp"):
            continue
        w, tool, etype, dpsi, c, tp, fn, fp, p, r, f1 = line.strip("\n").split(",")
        data.append([tool, etype, int(c), float(p), float(r)])
    df = pd.DataFrame(
        data, columns=["Tool", "E.Type", "TrueCov", "Precision", "Recall"]
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
        if etype == "ES":
            sns.lineplot(
                df[df["E.Type"] == etype],
                x="Precision",
                y="Recall",
                hue="Tool",
                palette=colors,
                style="TrueCov",
                markers=True,
                dashes=False,
                markersize=10,
                ax=ax,
            )
        else:
            sns.lineplot(
                df[df["E.Type"] == etype],
                x="Precision",
                y="Recall",
                hue="Tool",
                palette=colors,
                style="TrueCov",
                legend=None,
                markers=True,
                dashes=False,
                markersize=10,
                ax=ax,
            )
        ax.set_title(etype)
        ax.set_xlim(0, 1.01)
        ax.set_ylim(0, 1.01)

    plt.tight_layout()
    plt.show()  # savefig(fpath + ".png")


if __name__ == "__main__":
    main()
