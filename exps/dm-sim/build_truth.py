import sys


def parse_event_annotation(fpath):
    """
    es: gs/ge is skipped exon ((gs,ge) in read counts as exon)
    ir: gs/ge is retained intron ((gs-1,ge+1) in read counts as junction)
    a5: gs/ge is extension (one is start/end of junction in read counts depending on strand)
    a3: gs/ge is extension (one is start/end of junction in read counts depending on strand)
    """
    events = {}
    for line in open(fpath):
        if line.startswith("event"):
            continue
        etype, tvar, templ, gs, ge, ts, te = line.strip("\n").split("\t")
        gs, ge = int(gs), int(ge)
        assert gs < ge
        events[tvar] = (gs, ge)
    return events


def main():
    # FIXME: hardcoded to 2 conditions, 1 replicate
    # TODO: get diff?
    event_annotation_path = sys.argv[1]
    counts = sys.argv[2]

    events = parse_event_annotation(event_annotation_path)
    templates = {}
    alternates = {}
    strands = {}
    for line in open(counts):
        if line.startswith("seqnames"):
            continue
        (
            chrom,
            gs,
            ge,
            strand,
            feature,
            gene_id,
            transcript_id,
            exon_number,
            tr_start,
            tr_end,
            read_count_1,
            read_count_2,
        ) = line.strip("\n").split(",")
        gs, ge = int(gs), int(ge)
        read_count_1, read_count_2 = int(read_count_1), int(read_count_2)
        if gene_id not in templates:
            templates[gene_id] = []
            strands[gene_id] = strand
            alternates[gene_id] = {}
        if transcript_id.endswith("template"):
            templates[gene_id].append(
                (feature, chrom, gs, ge, read_count_1, read_count_2)
            )
        else:
            if transcript_id not in alternates[gene_id]:
                alternates[gene_id][transcript_id] = []
            alternates[gene_id][transcript_id].append(
                (feature, chrom, gs, ge, read_count_1, read_count_2)
            )

    for gene_id, template in templates.items():
        strand = strands[gene_id]
        for transcript_id, alternate in alternates[gene_id].items():
            if transcript_id.endswith("es"):
                skipped_exon = events[transcript_id]
                chrom = None
                j1, j2, jj = None, None, None  # junction 1 and 2 without skip, jj jumps
                for t, chrom, s, e, rc1, rc2 in template:
                    if t == "junction":
                        if e == skipped_exon[0]:
                            j1 = (s, e, rc1, rc2)
                        elif s == skipped_exon[1]:
                            j2 = (s, e, rc1, rc2)
                for t, chrom, s, e, rc1, rc2 in alternate:
                    if t == "junction":
                        if s == j1[0] and e == j2[1]:
                            jj = (s, e, rc1, rc2)
                assert j1 != None and j2 != None and jj != None
                w1, w2, w3 = j1[2], j2[2], jj[2]
                try:
                    psi1 = ((w1 + w2) / 2) / ((w1 + w2) / 2 + w3)
                except ZeroDivisionError:
                    psi1 = "NaN"
                w1, w2, w3 = j1[3], j2[3], jj[3]
                try:
                    psi2 = ((w1 + w2) / 2) / ((w1 + w2) / 2 + w3)
                except ZeroDivisionError:
                    psi2 = "NaN"
                # Inclusion isoform is canonical (numerator)
                print(
                    "ES",
                    chrom,
                    gene_id,
                    strand,
                    f"{chrom}:{j1[0]}-{j1[1]}",
                    f"{chrom}:{j2[0]}-{j2[1]}",
                    f"{chrom}:{jj[0]}-{jj[1]}",
                    f"{j1[2]}/{j2[2]}/{jj[2]}",
                    f"{j1[3]}/{j2[3]}/{jj[3]}",
                    psi1,
                    psi2,
                    sep=",",
                )
            elif (strand == "+" and transcript_id.endswith("a5")) or (
                strand == "-" and transcript_id.endswith("a3")
            ):
                extension = events[transcript_id]
                chrom = None
                sj, lj = (
                    None,
                    None,
                )  # shorter (from template) and longer (from alternate) junctions
                for t, chrom, s, e, rc1, rc2 in template:
                    if t == "junction":
                        if s == extension[1]:
                            sj = (s, e, rc1, rc2)
                for t, chrom, s, e, rc1, rc2 in alternate:
                    if t == "junction":
                        if (
                            s == extension[0] - 1
                        ):  # +1 since we want the same as in the annotation
                            lj = (s, e, rc1, rc2)
                assert sj != None and lj != None
                assert sj[1] == lj[1]
                try:
                    psi1 = sj[2] / (sj[2] + lj[2])
                except ZeroDivisionError:
                    psi1 = "NaN"
                try:
                    psi2 = sj[3] / (sj[3] + lj[3])
                except ZeroDivisionError:
                    psi2 = "NaN"
                # Isoform with longer exon is canonical
                print(
                    "A3" if transcript_id.endswith("a3") else "A5",
                    chrom,
                    gene_id,
                    strand,
                    f"{chrom}:{sj[0]}-{sj[1]}",
                    f"{chrom}:{lj[0]}-{lj[1]}",
                    ".",
                    f"{sj[2]}/{lj[2]}",
                    f"{sj[3]}/{lj[3]}",
                    psi1,
                    psi2,
                    sep=",",
                )
            elif (strand == "-" and transcript_id.endswith("a5")) or (
                strand == "+" and transcript_id.endswith("a3")
            ):
                extension = events[transcript_id]
                chrom = None
                sj, lj = (
                    None,
                    None,
                )  # shorter (from template) and longer (from alternate) junctions
                for t, chrom, s, e, rc1, rc2 in template:
                    if t == "junction":
                        if e == extension[0]:
                            sj = (s, e, rc1, rc2)
                for t, chrom, s, e, rc1, rc2 in alternate:
                    if t == "junction":
                        if (
                            e == extension[1] + 1
                        ):  # +1 since we want the same as in the annotation
                            lj = (s, e, rc1, rc2)
                assert sj != None and lj != None
                assert sj[0] == lj[0]
                try:
                    psi1 = sj[2] / (sj[2] + lj[2])
                except ZeroDivisionError:
                    psi1 = "NaN"
                try:
                    psi2 = sj[3] / (sj[3] + lj[3])
                except ZeroDivisionError:
                    psi2 = "NaN"
                # Isoform with longer exon is canonical
                print(
                    "A3" if transcript_id.endswith("a3") else "A5",
                    chrom,
                    gene_id,
                    strand,
                    f"{chrom}:{sj[0]}-{sj[1]}",
                    f"{chrom}:{lj[0]}-{lj[1]}",
                    ".",
                    f"{sj[2]}/{lj[2]}",
                    f"{sj[3]}/{lj[3]}",
                    psi1,
                    psi2,
                    sep=",",
                )
            elif transcript_id.endswith("ir"):
                retained_intron = events[transcript_id]
                chrom = None
                sj, exon = (
                    None,
                    (0, 0, 0, 0),
                )  # splice junction and retained intron
                for t, chrom, s, e, rc1, rc2 in template:
                    if t == "junction":
                        if s == retained_intron[0] - 1 and e == retained_intron[1] + 1:
                            sj = (s, e, rc1, rc2)
                for t, chrom, s, e, rc1, rc2 in alternate:
                    if (
                        t == "exon"
                        and s < retained_intron[0]
                        and retained_intron[1] < e
                    ):
                        exon = (s, e, exon[2], exon[3])
                    if (
                        t == "-exon"
                        and s == retained_intron[0]
                        and e == retained_intron[1]
                    ):
                        # Keep longer exon, use coverage of retained portion
                        exon = (exon[0], exon[1], rc1, rc2)

                assert sj != None and exon != (0, 0, 0, 0)
                assert exon[0] < sj[0] and sj[0] < sj[1] and sj[1] < exon[1]
                try:
                    psi1 = sj[2] / (sj[2] + exon[2])
                except ZeroDivisionError:
                    psi1 = "NaN"
                try:
                    psi2 = sj[3] / (sj[3] + exon[3])
                except ZeroDivisionError:
                    psi2 = "NaN"
                # Isoform with splice junction is canonical
                print(
                    "IR",
                    chrom,
                    gene_id,
                    strand,
                    f"{chrom}:{sj[0]}-{sj[1]}",
                    f"{chrom}:{exon[0]}-{exon[1]}",
                    ".",
                    f"{sj[2]}/{exon[2]}",
                    f"{sj[3]}/{exon[3]}",
                    psi1,
                    psi2,
                    sep=",",
                )


if __name__ == "__main__":
    main()
