# python alignments_augmentation_from_gaf.py  alignment.gaf output.path input.gfa > output.gfa


import sys
import json
import re


def main(argv):
    gaf_file = argv[0]
    output_file = argv[1]
    gfa_file = argv[2]
    # gfa_file_out = argv[3]
    weights = {}
    revs = {}
    print("Building paths and weights", file=sys.stderr)
    with open(gaf_file) as f, open(output_file, "w") as out:
        last_read = ""
        for line in f:
            tokens = line.strip().split()
            if tokens[5] == "*":
                continue
            read_name = tokens[0]
            path = tokens[5]
            cigar_re = re.search(
                r"cs:.*?(?=\s|$)", " ".join(item for item in tokens[12:])
            )
            if cigar_re:
                cigar = cigar_re.group(0).replace("cs:Z:", "")
            else:
                cigar = "*"
            if path[0] == ">":
                nodes = path.split(">")[1:]
            else:
                nodes = path.split("<")[1:]
                print(nodes)
                nodes.reverse()
            print(nodes)
            for n1, n2 in zip(nodes, nodes[1:]):
                if (n1, n2) in weights.keys():
                    weights[(n1, n2)] = weights[(n1, n2)] + 1
                else:
                    weights[(n1, n2)] = 1
            if read_name == last_read:
                out.write(f"{path}\t{cigar}\n")
            else:
                out.write(f"{read_name}\n{path}\t{cigar}\n")
            last_read = read_name

    print("Annotating GFA", file=sys.stderr)
    with open(gfa_file, "r") as f:

        for line in f:
            line = line.strip()
            if not line.startswith("L"):
                print(line)
            else:
                if len(line) == 1:
                    continue
                tokens = line.split()
                w = weights.pop((tokens[1], tokens[3]), 0)
                print(f"{line}\tRC:i:{w}")

    for k, v in weights.items():
        print(f"L\t{k[0]}\t+\t{k[1]}\t+\t*\tRC:i:{v},ID:Z:N")


if __name__ == "__main__":
    main(sys.argv[1:])
