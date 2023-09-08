import sys


def main():
    gfa_path = sys.argv[1]
    pruned_gfa_path = sys.argv[2]

    # Double pass to print empty GFA in case of failure (and not incomplete one) - maybe we can avoid this

    edges = set()
    for line in open(pruned_gfa_path):
        if line.startswith("L"):
            _, n1, _, n2, _, _ = line.split("\t")
            n1 = int(n1)
            n2 = int(n2)
            edges.add((n1, n2))

    paths = []
    for line in open(gfa_path):
        if line.startswith("P"):
            nodes = line.split("\t")[2].split(",")
            strand = nodes[-1][-1]
            nodes = [int(x[:-1]) for x in nodes]
            if strand == "-":
                nodes.reverse()
            for n1, n2 in zip(nodes[:-1], nodes[1:]):
                assert (n1, n2) in edges

    for line in open(pruned_gfa_path):
        print(line, end="")
    for line in open(gfa_path):
        if line.startswith("P"):
            print(line, end="")


if __name__ == "__main__":
    main()
