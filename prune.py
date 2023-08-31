#!/usr/bin/env python3

import sys
import argparse


def main(args):
    # Find min and max nodes used in the gene (assuming IDs to be topological sorted - as vg)
    nodes_to_keep = set()
    nnodes = 0
    for line in open(args.GFA):
        if line.startswith("S"):
            nnodes += 1
        elif line.startswith("P"):
            _, pname, nodes, _ = line.split("\t")
            if not pname.startswith(args.tprefix):
                continue
            nodes = [int(n[:-1]) for n in nodes.split(",")]
            min_idx = min(nodes) - args.w
            max_idx = max(nodes) + args.w
            nodes_to_keep |= set(range(min_idx, max_idx + 1))
    print(f"Pruning {len(nodes_to_keep)}/{nnodes} nodes..", file=sys.stderr)

    # Reiterate over gfa and print only subgraph
    for line in open(args.GFA):
        if line.startswith("H"):
            print(line, end="")
        elif line.startswith("S"):
            idx = int(line.split("\t")[1])
            if idx in nodes_to_keep:
                print(line, end="")
        elif line.startswith("L"):
            idx1 = int(line.split("\t")[1])
            idx2 = int(line.split("\t")[3])
            if idx1 in nodes_to_keep and idx2 in nodes_to_keep:
                print(line, end="")
        elif line.startswith("P"):
            f1, pname, nodes, lf = line.strip("\n").split("\t")
            if pname.startswith(args.tprefix) or pname.startswith(
                "_alt"
            ):  # FIXME: hardcoded
                print(line, end="")
            else:
                paths = []
                curr_path = []
                nodes = [int(n[:-1]) for n in nodes.split(",")]
                for node in nodes:
                    if node in nodes_to_keep:
                        curr_path.append(node)
                    else:
                        if curr_path != []:
                            paths.append(curr_path)
                            curr_path = []
                for i, path in enumerate(paths, 1):
                    print(
                        "P",
                        f"{pname}_{i}",
                        ",".join([f"{x}+" for x in path]),
                        "*",
                        sep="\t",
                    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="Pruning utility",
        description="",
    )
    parser.add_argument("GFA", help="Spliced pangenome in GFA format")
    parser.add_argument(
        "-w",
        help="Extend subgraphs by W nodes (default: 0)",
        dest="w",
        type=int,
        required=False,
        default=0,
    )
    parser.add_argument(
        "-t",
        dest="tprefix",
        help="Transcript prefix (default: ENST)",
        required=False,
        default="ENST",
    )
    args = parser.parse_args()
    main(args)
