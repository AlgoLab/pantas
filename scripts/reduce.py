import sys
from intervaltree import Interval, IntervalTree


def main():
    # assuming a single ref path (e.g., 1 chromosome)
    gfa_path = sys.argv[1]
    ogfa_path = sys.argv[2]

    k = 0  # FIXME: if we use a value>0, we must be sure to extend to a node in the ref path

    nodes_l = {}
    ref_path_name = ""
    ref_path = []
    tree = IntervalTree()
    for line in open(gfa_path):
        if line.startswith("S"):
            _, idx, seq = line.strip("\n").split("\t")
            nodes_l[int(idx)] = len(seq)
        elif line.startswith("P"):
            _, tidx, nodes, _ = line.strip("\n").split("\t")
            if not tidx.endswith("_R1"):
                ref_path_name = tidx
                ref_path = [int(x[:-1]) for x in nodes.split(",")]
                continue
            plus = nodes[-1] == "+"
            nodes = [int(x[:-1]) for x in nodes.split(",")]
            if not plus:
                nodes.reverse()
            assert all(b >= a for a, b in zip(nodes[:-1], nodes[1:]))
            tree[nodes[0] - k : nodes[-1] + k + 1] = 1

    # Get reference position for each node on the reference path
    ref_positions = {}
    p = 0
    for n in ref_path:
        ref_positions[n] = p
        p += nodes_l[n]

    print(len(tree))
    tree.merge_overlaps()
    print(len(tree))
    intervals = []
    for i in tree:
        intervals.append((i.begin, i.end - 1))
    intervals.sort()
    starts = [x[0] for x in intervals[1:]]
    ends = [x[1] for x in intervals[:-1]]
    print(len(intervals))

    gfa = open(ogfa_path, "w")
    gfarefpath = open(ogfa_path + ".refpath", "w")
    curr_node_idx = -1
    for line in open(gfa_path):
        if line.startswith("S"):
            _, idx, seq = line.strip("\n").split("\t")
            idx = int(idx)
            curr_node_idx = max(curr_node_idx, idx)
            hit = len(tree[idx]) > 0
            if hit:
                if idx in ref_positions:
                    line = line[:-1] + "\t" + "RP:i:" + str(ref_positions[idx]) + "\n"
                print(line, file=gfa, end="")
        elif line.startswith("L"):
            _, idx1, strand1, idx2, strand2, cigar = line.strip("\n").split("\t")
            idx1, idx2 = int(idx1), int(idx2)
            hit = len(tree[idx1]) > 0 and len(tree[idx2]) > 0
            if hit:
                print(line, file=gfa, end="")
        elif line.startswith("P"):
            _, tidx, nodes, _ = line.strip("\n").split("\t")
            if tidx.endswith("_R1"):
                print(line, file=gfa, end="")

    curr_node_idx += 1  # just to be sure
    print(curr_node_idx)
    dummy_seq = "N" * 32
    ref_path_new = []
    p = 0
    for r in ref_path:
        if len(tree[r]) > 0:
            ref_path_new.append(r)
            if r in ends:
                print("S", curr_node_idx, dummy_seq, file=gfa, sep="\t")
                print("L", r, "+", curr_node_idx, "+", "0M", file=gfa, sep="\t")
                ref_path_new.append(curr_node_idx)
            elif r in starts:
                print("L", curr_node_idx, "+", r, "+", "0M", file=gfa, sep="\t")
                curr_node_idx += 1
        else:
            pass  # print(r, "X")
    print(
        "P",
        ref_path_name,
        ",".join([f"{x}+" for x in ref_path_new]),
        "*",
        file=gfa,
        sep="\t",
    )
    print(
        ref_path_name,
        ",".join(
            [f"{ref_positions[x]}" if x in ref_positions else "." for x in ref_path_new]
        ),
        file=gfarefpath,
        sep="\t",
    )
    gfa.close()
    gfarefpath.close()


if __name__ == "__main__":
    main()
