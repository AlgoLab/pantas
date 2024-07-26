import sys
from intervaltree import Interval, IntervalTree

# FIXME: we do not store any information on the "original" position of vertices in the reference path
def main():
    gfa_path = sys.argv[1]

    k = 0

    ref_name = ""
    ref_path = []
    tree = IntervalTree()
    for line in open(gfa_path):
        if line.startswith("P"):
            _, tidx, nodes, _ = line.strip("\n").split("\t")
            if not tidx.startswith("FBtr"):  # FIXME: hardcoded
                ref_name = tidx
                ref_path = [int(x[:-1]) for x in nodes.split(",")]
                continue
            plus = nodes[-1] == "+"
            nodes = [int(x[:-1]) for x in nodes.split(",")]
            min_idx, max_idx = min(nodes), max(nodes)
            # ---
            if not plus:
                nodes.reverse()
            assert all(b >= a for a, b in zip(nodes[:-1], nodes[1:]))
            # ---
            tree[min_idx - k : max_idx + k + 1] = 1

    print(f"We have {len(tree)} unique transcripts", file=sys.stderr)
    tree.merge_overlaps()
    print(f"We have {len(tree)} genic regions", file=sys.stderr)

    for i in tree:
        print(i.begin, i.end - 1, file=sys.stderr)

    for line in open(gfa_path):
        if line.startswith("S"):
            _, idx, _ = line.strip("\n").split("\t")
            idx = int(idx)
            if len(tree[idx]) > 0:
                print(line, end="")
        elif line.startswith("L"):
            _, idx1, _, idx2, _, _ = line.strip("\n").split("\t")
            idx1, idx2 = int(idx1), int(idx2)
            if len(tree[idx1]) > 0 and len(tree[idx2]) > 0:
                print(line, end="")
        elif line.startswith("P"):
            _, tidx, _, _ = line.strip("\n").split("\t")
            if tidx.startswith("FBtr"):  # FIXME: hardcoded
                print(line, end="")

    ref_subpath = []
    i = 0
    for n in ref_path:
        if len(tree[n]) > 0:
            ref_subpath.append(n)
        else:
            if len(ref_subpath) > 0:
                print(
                    "P",
                    f"{ref_name}.{i}",
                    ",".join([str(x) + "+" for x in ref_subpath]),
                    "*",
                    sep="\t",
                )
                i += 1
                ref_subpath = []


if __name__ == "__main__":
    main()
