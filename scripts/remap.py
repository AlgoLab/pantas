import re


def parse_region(r):
    if r == "?":
        return "?", -1, -1, True
    imprecise = r.endswith("?")
    if imprecise:
        r = r[:-1]
    chrom = r.split(":")[0]
    s, e = r.split(":")[1].split("-")
    return chrom, int(s), int(e), imprecise


def get_region_size(r):
    chrom = r.split(":")[0]
    s, e = r.split(":")[1].split("-")
    return int(e) - int(s)


def get_reference_transcript(hts):
    t = [x for x in hts.split("|") if x.split("_")[-1][0] == "R"]
    if len(t) == 0:
        return "?"
    else:
        return t[0]


def main(args):
    # Iterate over quant csv to get transcripts we need
    transcripts = {}
    for line in open(args.CSV):
        if line.startswith("etype"):
            continue
        junction1_names, junction2_names, junction3_names = line.strip("\n").split(",")[
            5:8
        ]

        junction1_name = get_reference_transcript(junction1_names)
        junction2_name = get_reference_transcript(junction2_names)
        junction3_name = get_reference_transcript(junction3_names)
        for t in [junction1_name, junction2_name, junction3_name]:
            if t == "?":
                continue
            t = "_".join(t.split(".")[0].split("_")[:-1])
            transcripts[t] = []

    # Parse GTF to get exon positions
    for line in open(args.GTF):
        if line.startswith("#"):
            continue
        line = line.strip("\n").split("\t")

        if line[2] != "exon":
            continue
        tidx = (
            re.search('transcript_id "[A-Za-z0-9_]+";', line[-1])
            .group(0)
            .split('"')[-2]
        )
        if tidx in transcripts:
            s, e = int(line[3]), int(line[4])
            transcripts[tidx].append((s, e))
    for exons in transcripts.values():
        exons.sort()

    print(
        "etype",
        "annotation_type",
        "haplotype_type",
        "chrom",
        "gene",
        "strand",
        "junction1_name",
        "junction2_name",
        "junction3_name",
        "junction1_nodes",
        "junction2_nodes",
        "junction3_nodes",
        "junction1_positions",
        "junction2_positions",
        "junction3_positions",
        "W1",
        "W2",
        "psi_c1",
        "psi_c2",
        "dpsi",
        sep=",",
    )
    # Reiterate over quant CSV to add reference positions
    for line in open(args.CSV):
        if line.startswith("etype"):
            continue
        (
            etype,
            annotation_type,
            chrom,
            gene,
            strand,
            junction1_names,
            junction2_names,
            junction3_names,
            junction1_nodes,
            junction2_nodes,
            junction3_nodes,
            W1,
            W2,
            psi_c1,
            psi_c2,
            dpsi,
        ) = line.strip("\n").split(",")

        junction1_name = get_reference_transcript(junction1_names)
        junction2_name = get_reference_transcript(junction2_names)
        junction3_name = get_reference_transcript(junction3_names)

        positions = []
        for t in [junction1_name, junction2_name, junction3_name]:
            if t == "?":
                positions.append("?")
                continue
            tidx = "_".join(t.split(".")[0].split("_")[:-1])
            imprecise = (
                t.split(".")[0].split("_")[-1][0] == "H"
            )  # we are on an haplotype-aware transcript
            exons = [int(n) for n in t.split(".")[1:]]
            assert len(exons) <= 2
            if len(exons) == 0:
                positions.append("?")
            elif len(exons) == 1:
                # Intron retention - We report the full exon
                e = exons[0]
                # start is start of exon
                s = transcripts[tidx][e - 1][0]
                # end is end of exon
                e = transcripts[tidx][e - 1][1]
                assert s <= e  # we can have introns of length 1
                positions.append(f"{chrom}:{s}-{e}" + ("?" if imprecise else ""))
            else:
                # Splicing junction
                e1, e2 = exons
                # start is end of first exon + 1
                s = transcripts[tidx][e1 - 1][1] + 1
                # end is start of second exon - 1
                e = transcripts[tidx][e2 - 1][0] - 1
                assert s <= e  # we can have introns of length 1
                positions.append(f"{chrom}:{s}-{e}" + ("?" if imprecise else ""))

        if annotation_type == "annotated":
            # exon skippings are good. We want inclusion1, inclusion2, skipping junction
            # alternative splice sites, we have to put the shorter junction first
            # intron retention, we have only retained intron. So we put it first
            if etype[0] == "A":
                if positions[0] == "?":
                    positions[0] = positions[1]
                    positions[1] = "?"
                    junction1_name = junction2_name
                    junction2_name = "?"
                    positions[0] = positions[1]
                    positions[1] = "?"
                    junction1_nodes = junction2_nodes
                    junction2_nodes = "."
                    W1, W2 = W2, W1
                    psi_c1, psi_c2 = psi_c2, psi_c1
                    dpsi = -float(dpsi)
                elif positions[1] == "?":
                    pass
                else:
                    p0 = parse_region(positions[0])
                    p1 = parse_region(positions[1])
                    if p0[2] - p0[1] > p1[2] - p1[1]:
                        positions[0], positions[1] = positions[1], positions[0]
                        junction1_name, junction2_name = junction2_name, junction1_name
                        junction1_nodes, junction2_nodes = (
                            junction2_nodes,
                            junction1_nodes,
                        )
                        W1, W2 = W2, W1
                        psi_c1, psi_c2 = psi_c2, psi_c1
                        dpsi = -float(dpsi)
            elif etype == "IR":
                junction1_name = junction2_name
                junction2_name = "?"
                positions[0] = positions[1]
                positions[1] = "?"
                junction1_nodes = junction2_nodes
                junction2_nodes = "."
                W1, W2 = W2, W1
                psi_c1, psi_c2 = psi_c2, psi_c1
                dpsi = 1.0-float(dpsi)

        htype = "reference"
        if annotation_type == "annotated":
            if etype == "ES" and (
                positions[0] == "?" or positions[1] == "?" or positions[2] == "?"
            ):
                htype = "haplotype"
            if etype[0] == "A" and (positions[0] == "?" or positions[1] == "?"):
                htype = "haplotype"
            if etype == "IR" and positions[0] == "?":
                htype = "haplotype"
        else:
            if etype == "ES" and (
                positions[0] == "?" and positions[1] == "?" and positions[2] == "?"
            ):
                htype = "haplotype"
            if etype[0] == "A" and (positions[0] == "?" and positions[1] == "?"):
                htype = "haplotype"
            if etype == "IR" and (positions[0] == "?" and positions[1] == "?"):
                htype = "haplotype"

        if annotation_type == "novel":
            if any(
                [
                    get_region_size(p) < args.min_intron_size
                    for p in positions
                    if p != "?"
                ]
            ):
                continue

        print(
            etype,
            annotation_type,
            htype,
            chrom,
            gene,
            strand,
            junction1_name,
            junction2_name,
            junction3_name,
            junction1_nodes,
            junction2_nodes,
            junction3_nodes,
            positions[0],
            positions[1],
            positions[2],
            W1,
            W2,
            psi_c1,
            psi_c2,
            dpsi,
            sep=",",
        )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        prog="Remap",
        description="",
    )
    parser.add_argument("CSV", help="Quantification in CSV format")
    parser.add_argument("GTF", help="Annotation in GTF format")
    parser.add_argument(
        "-i",
        dest="min_intron_size",
        help="Minimum intron size for a novel junction (default: 100)",
        type=int,
        default=100,
    )
    args = parser.parse_args()
    main(args)
