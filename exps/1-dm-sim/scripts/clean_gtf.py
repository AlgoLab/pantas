import sys
from intervaltree import Interval, IntervalTree


def main():
    gtf = sys.argv[1]

    trees = {}
    for line in open(gtf):
        if line.startswith("#"):
            continue
        chrom, _, t, s, e = line.split("\t")[0:5]
        s = int(s)
        e = int(e) + 1
        if t != "gene":
            continue
        if chrom not in trees:
            trees[chrom] = IntervalTree()
        trees[chrom][s:e] = len(trees[chrom][s:e]) > 0

    p = False
    for line in open(gtf):
        if line.startswith("#"):
            continue
        chrom, _, t, s, e = line.split("\t")[0:5]
        s = int(s)
        e = int(e) + 1
        if t == "gene":
            pflag = len(trees[chrom][s:e]) == 1
        if pflag:
            if not list(trees[chrom][s:e])[0].data:
                print(line, end="")


if __name__ == "__main__":
    main()
