import sys
import re
import gffutils
import pysam


def split(i):
    return (
        i.split(":")[0],
        int(i.split(":")[1].split("-")[0]),
        int(i.split(":")[1].split("-")[1]),
    )


def main():
    events_path = sys.argv[1]
    gtf_path = sys.argv[2]
    bam_paths = sys.argv[3:]
    bams = [pysam.AlignmentFile(bam, "rb") for bam in bam_paths]

    gtf = None
    try:
        gtf = gffutils.FeatureDB(gtf_path + ".db")
    except ValueError:
        gtf = gffutils.create_db(
            gtf_path,
            gtf_path + ".db",
            disable_infer_genes=True,
            disable_infer_transcripts=True,
        )

    events = set()
    for line in open(events_path):
        events.add(line.strip("\n"))

    found_skipped_events = {}
    for gene in gtf.features_of_type("gene"):
        chrom = gene.seqid
        skipped_exon = ""
        skipping_introns = set()
        for transcript in gtf.children(
            gene, featuretype="transcript", order_by="start"
        ):
            exons = list(
                (e.start, e.end)
                for e in gtf.children(transcript, featuretype="exon", order_by="start")
            )
            for e1, e2, e3 in zip(exons[:-2], exons[1:-1], exons[2:]):
                k1 = f"{chrom}:{e1[0]}-{e1[1]}"
                k2 = f"{chrom}:{e2[0]}-{e2[1]}"
                k3 = f"{chrom}:{e3[0]}-{e3[1]}"
                # print(k2, events)
                if k2 in events:
                    if k2 not in found_skipped_events:
                        found_skipped_events[k2] = []
                    found_skipped_events[k2].append((k1, k2, k3))
    events = set()
    for k, Ks in found_skipped_events.items():
        for k1, k2, k3 in Ks:
            assert k == k2
            k1, k2, k3 = split(k1), split(k2), split(k3)
            events.add((k1[0], k1[2], k2[1], k2[2], k3[1]))

    p = re.compile(r"[0-9]+N")
    for event in events:
        chrom, c1, c2, c3, c4 = event
        for i, bam in enumerate(bams):
            introns = {c2 - c1 - 1: 0, c4 - c3 - 1: 0, c4 - c1 - 1: 0}
            for aln in bam.fetch(chrom, c1, c4):
                if "N" not in aln.cigarstring:
                    continue
                for m in re.findall(p, aln.cigarstring):
                    m = int(m[:-1])
                    if m in introns:
                        introns[m] += 1
            print(
                f"{event[0]}:{event[2]}-{event[3]}",
                bam_paths[i],
                " ".join([str(v) for v in introns.values()]),
            )


if __name__ == "__main__":
    main()
