import sys


def main():
    # assuming the first segment of each GFA to be 1
    # + assuming the GFA to be topological sorted (maybe we do not need this) - vg rna outputs a sorted graph anyway
    gfa_fns = sys.argv[1:]

    shift = 0
    max_idx = 0
    print("H", "VN:Z:1.1", sep="\t")
    for gfa_fn in gfa_fns:
        for line in open(gfa_fn):
            if line.startswith("H"):
                continue
            elif line.startswith("S"):
                _, idx, seq, *rest = line.strip("\n").split("\t")
                idx = int(idx) + shift
                max_idx = max(idx, max_idx)
                if len(rest) == 0:
                    print("S", idx, seq, sep="\t")
                else:
                    print("S", idx, seq, "\t".join(rest), sep="\t")
            elif line.startswith("L"):
                _, idx1, strand1, idx2, strand2, *rest = line.strip("\n").split("\t")
                idx1 = int(idx1) + shift
                idx2 = int(idx2) + shift
                # we are sure to have something in rest
                print("L", idx1, strand1, idx2, strand2, "\t".join(rest), sep="\t")
            elif line.startswith("P"):
                _, pname, path, *rest = line.strip("\n").split("\t")
                path = [str(int(x[:-1]) + shift) + "+" for x in path.split(",")]
                # we are sure to have something in rest
                print("P", pname, ",".join(path), "\t".join(rest), sep="\t")
        print(
            f"Dumped {gfa_fn} shifting by {shift}. New shift: {max_idx}",
            file=sys.stderr,
        )
        shift = max_idx


if __name__ == "__main__":
    main()
