import sys


def main():
    gaf_path = sys.argv[1]
    gfa_path = sys.argv[2]
    ref_path = sys.argv[3]

    print("@HD", "VN:1.5", sep="\t")

    print("@SQ", "SN:chr13", "LN:114364328", sep="\t")
    print("@SQ", "SN:chr6", "LN:170805979", sep="\t")
    print("@SQ", "SN:chrX", "LN:156040895", sep="\t")
    print("@SQ", "SN:chr21", "LN:46709983", sep="\t")
    print("@SQ", "SN:chr16", "LN:90338345", sep="\t")
    print("@SQ", "SN:chr7", "LN:159345973", sep="\t")
    print("@SQ", "SN:chr1", "LN:248956422", sep="\t")
    print("@SQ", "SN:chr9", "LN:138394717", sep="\t")
    print("@SQ", "SN:chr14", "LN:107043718", sep="\t")
    print("@SQ", "SN:chr20", "LN:64444167", sep="\t")
    print("@SQ", "SN:chr17", "LN:83257441", sep="\t")
    print("@SQ", "SN:chr4", "LN:190214555", sep="\t")
    print("@SQ", "SN:chr11", "LN:135086622", sep="\t")
    print("@SQ", "SN:chr15", "LN:101991189", sep="\t")
    print("@SQ", "SN:chr12", "LN:133275309", sep="\t")
    print("@SQ", "SN:chr5", "LN:181538259", sep="\t")
    print("@SQ", "SN:chr10", "LN:133797422", sep="\t")
    print("@SQ", "SN:chr2", "LN:242193529", sep="\t")
    print("@SQ", "SN:chr8", "LN:145138636", sep="\t")
    print("@SQ", "SN:chr3", "LN:198295559", sep="\t")

    ref_positions = {}
    nodes = {}
    nodes_l = {}
    nodes_to_path = {}
    for line in open(ref_path):
        pname, rpos = line.strip("\n").split("\t")
        ref_positions[pname] = [int(x) if x != "." else x for x in rpos.split(",")]
    for line in open(gfa_path):
        if line.startswith("S"):
            _, idx, seq, *attrs = line.strip("\n").split("\t")
            nodes_l[int(idx)] = len(seq)
        elif line.startswith("P"):
            _, pname, Ns, lf = line.strip("\n").split("\t")
            if not pname.endswith("_R1"):
                for i, node in enumerate([int(n[:-1]) for n in Ns.split(",")]):
                    nodes[node] = ref_positions[pname][i]
                    nodes_to_path[node] = pname

    skipped = 0
    skipped_n = 0
    total = 0
    for line in open(gaf_path):
        rname, rl, rs, re, strand, path, pl, ps, pe, _, _, mapq, *attrs = line.strip(
            "\n"
        ).split("\t")
        # if rname != "SRR1513329.971":
        #     continue
        if path == "*":
            continue
        total += 1
        print(line, end="", file=sys.stderr)
        print(rname, path, mapq, file=sys.stderr)
        if "<" in path:
            path = [int(x) for x in path[1:].split("<")]
            path.reverse()
        else:
            path = [int(x) for x in path[1:].split(">")]
        print("...", nodes_l[path[0]], file=sys.stderr)
        if path[0] not in nodes_to_path or nodes[path[0]] == ".":
            skipped += 1
            continue
        rpos = [nodes[x] if x in nodes else -1 for x in path]
        Ls = [nodes_l[x] for x in path]
        last_p = rpos[0]
        cigar = []
        for n, l, p in zip(path, Ls, rpos):
            print(n, l, p, last_p, file=sys.stderr)
            if p == -1 or p == ".":
                p = last_p + l - 1
            print(n, l, p, last_p, file=sys.stderr)
            if p != last_p:
                n = p - last_p
                if n <= 0:
                    skipped_n += 1
                    continue
                cigar.append((n, "N"))
            cigar.append((l, "M"))
            last_p = p + l
        print(cigar, file=sys.stderr)

        compact_cigar = [cigar[0]]
        for l, op in cigar[1:]:
            if op == compact_cigar[-1][1]:
                compact_cigar[-1] = (compact_cigar[-1][0] + l, op)
            else:
                compact_cigar.append((l, op))
        cigar_s = "".join([f"{l}{op}" for l, op in compact_cigar])
        print(cigar_s, file=sys.stderr)
        print(
            rname,
            0 if strand == "+" else 16,
            nodes_to_path[path[0]],
            rpos[0],
            mapq,
            cigar_s,
            "*",
            0,
            0,
            "*",
            "*",
            sep="\t",
        )
        print("", file=sys.stderr)
    print(f"Skipped {skipped} - {skipped_n} over {total} alignments.", file=sys.stderr)


if __name__ == "__main__":
    main()
