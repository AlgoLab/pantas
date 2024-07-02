import sys
import gffutils

from compare import parse_truth


def main():
    truth_path = sys.argv[1]
    gtf_path = sys.argv[2]

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

    truth, _ = parse_truth(truth_path)
    print("Truth:", len(truth))
    truth = {k: v for k, v in truth.items() if abs(v) >= 0.05 and abs(v) <= 1 - 0.05}
    print(f"Filtered truth with delta=0.05:", len(truth))
    true_skipped_exons = set(truth.keys())

    print(len(true_skipped_exons))

    found_skipped_exons = set()
    for gene in gtf.features_of_type("gene"):
        chrom = gene.seqid
        skipped_exon = ""
        skipping_introns = set()
        for transcript in gtf.children(
            gene, featuretype="transcript", order_by="start"
        ):
            exons = list(gtf.children(transcript, featuretype="exon", order_by="start"))
            introns = [(ex1.end, ex2.start) for ex1, ex2 in zip(exons[:-1], exons[1:])]
            for i1, i2 in zip(introns[:-1], introns[1:]):
                if f"{chrom}:{i1[1]}-{i2[0]}" in true_skipped_exons:
                    # print(f"{chrom}:{i1[1]}-{i2[0]}")
                    skipping_introns.add((i1[0], i2[1]))
                    skipped_exon = f"{chrom}:{i1[1]}-{i2[0]}"

        if len(skipping_introns) == 0:
            continue
        novel = True
        for transcript in gtf.children(
            gene, featuretype="transcript", order_by="start"
        ):
            exons = list(gtf.children(transcript, featuretype="exon", order_by="start"))
            introns = set(
                [(ex1.end, ex2.start) for ex1, ex2 in zip(exons[:-1], exons[1:])]
            )
            if len(skipping_introns & introns) != 0:
                novel = False
                break
        found_skipped_exons.add(skipped_exon)
        print(skipped_exon, novel)

    # Manually checked on IGV
    assert len(true_skipped_exons - found_skipped_exons) == 1
    print(list(true_skipped_exons - found_skipped_exons)[0], "False")


if __name__ == "__main__":
    main()
