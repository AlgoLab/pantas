import gfautils
from Bio import SeqIO
import sys


def debug_align(s1, s2):
    print("Sequence len:", len(s1))
    l1 = "path: "
    l2 = "seq : "
    m = min(len(s1), len(s2))
    for i in range(m):
        if s1[i] == s2[i]:
            l1 += f"\033[1;32m{s1[i]}"
            l2 += f"\033[1;32m{s2[i]}"
        else:
            l1 += f"\033[1;31m{s1[i]}"
            l2 += f"\033[1;31m{s2[i]}"
    if len(s1) > m:
        for i in range(m, len(s1)):
            l1 += f"\033[1;34m{s1[i]}"
            l2 += "\033[1;34m-"
    if len(s2) > m:
        for i in range(m, len(s2)):
            l1 += "\033[1;31m*"
            l2 += f"\033[1;31m{s2[i]}"
    l1 += f"\033[0m"
    l2 += f"\033[0m"
    print(l1)
    print(l2)


def main(args):
    gfa = gfautils.GFA(args.GFA)

    for record in SeqIO.parse(args.GFFW, "fasta"):
        if not f"{record.name}_R1" in gfa.paths:
            print(f"Skipping {record.name}. Missing from GFA", file=sys.stderr)
            continue
        if not args.IS:
            if args.debug:
                debug_align(gfa.pseq(f"{record.name}_R1"), record.seq)
            assert gfa.pseq(f"{record.name}_R1") == record.seq

        desc = record.description.split()
        segs = (
            desc[[desc.index(l) for l in desc if l.startswith("segs:")][0]]
            .split(":")[1]
            .split(",")
        )
        intsegs = []

        for seg in segs:
            s, e = seg.split("-")
            s = int(s)
            e = int(e)
            intsegs.append((s, e))

        # print(intsegs)
        curr_seg = 0
        cum_len = 0
        is_reverse = gfa.paths[f"{record.name}_R1"].is_reverse
        for ix, nid_curr in enumerate(_pp := gfa.pnodes_list(f"{record.name}_R1")):
            cum_len += gfa.nlen(nid_curr)
            gfa.nodes[nid_curr].addexon(f"{record.name}.{curr_seg+1}")

            if cum_len == intsegs[curr_seg][1]:
                curr_seg += 1
                if ix == len(_pp) - 1:
                    continue

                _lkey = (
                    [nid_curr, _pp[ix + 1]]
                    if not is_reverse
                    else [_pp[ix + 1], nid_curr]
                )
                gfa.link(*_lkey).addjunction(f"{record.name}.{curr_seg}.{curr_seg+1}")

        assert curr_seg == len(intsegs)

    gfa.print()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Check that GFA Path are identical to input FAs"
    )
    parser.add_argument("GFA", type=str, help="Path to GFA")
    parser.add_argument("GFFW", type=str, help="Path to gffread out")
    parser.add_argument("--debug", action="store_true", help="Print alignment debug")
    parser.add_argument(
        "--IS",
        action="store_true",
        help="Ignore sequence check. Assume they are correct.",
    )
    args = parser.parse_args()

    main(args)
