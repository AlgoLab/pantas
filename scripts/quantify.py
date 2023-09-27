import sys

ETYPES = ["ES", "IR", "A3", "A5"]


def get_interval(region):
    if region == ".":
        return region
    s, e = [int(x) for x in region.split(":")[1].split("-")]
    return s, e


class Event:
    def __init__(
        self, etype, novel, chrom, gene, strand, i1, i2, i3, n1, n2, n3, w1, w2, w3
    ):
        self.chrom = chrom
        self.etype = etype
        self.novel = novel
        self.strand = strand
        self.gene = gene
        self.intron1_r = i1
        self.intron2_r = i2
        self.intron3_r = i3
        self.intron1 = get_interval(self.intron1_r)
        self.intron2 = get_interval(self.intron2_r)
        self.intron3 = get_interval(self.intron3_r)
        self.nodes1 = n1
        self.nodes2 = n2
        self.nodes3 = n3
        self.w1 = w1
        self.w2 = w2
        self.w3 = w3
        self.sort()
        self.psi = ["NaN"]
        self.compute_psi()
        self.dpsi = "NaN"

    # Sort info
    def sort(self):
        if self.etype == "ES":
            pass  # 1 is skipping, 2-3 are inclusion
        elif self.etype in ["A3", "A5"]:
            # force 1 to be shorter intron
            if self.intron1[1] - self.intron1[0] > self.intron2[1] - self.intron2[0]:
                self.intron1, self.intron2 = self.intron2, self.intron1
                self.intron1_r, self.intron2_r = self.intron2_r, self.intron1_r
                self.w1, self.w2 = self.w2, self.w1
                self.nodes1, self.nodes2 = self.nodes2, self.nodes1
        elif self.etype == "IR":
            pass  # 1 is retaining junction, 2 is exon

    def compute_psi(self):
        psi = "NaN"
        if self.etype == "ES":
            if self.w1 + self.w2 + self.w3 != 0:
                psi = ((self.w2 + self.w3) / 2) / (self.w1 + (self.w2 + self.w3) / 2)
        elif self.etype in ["A3", "A5"]:
            if self.w1 + self.w2 != 0:
                psi = self.w1 / (self.w1 + self.w2)
        elif self.etype == "IR":
            # w1 is splicing, w2 is full exon (canonical)
            if self.w1 + self.w2 != 0:
                psi = self.w2 / (self.w1 + self.w2)
        self.psi = [psi]

    def __repr__(self):
        return f"{self.etype} {self.intron1_r} {self.intron2_r} {self.intron3_r}"

    def __eq__(self, o):
        return self.__hash__() == o.__hash__()

    def __hash__(self):
        return hash(self.__repr__())


# Given splicedfa from gffread, build a map between pantas' internal IDs and intervals on the reference genome
def get_referencebased_coordinates(splicedfa_path):
    exons = {}
    introns = {}
    for line in open(splicedfa_path):
        if not line.startswith(">"):
            continue
        line = line[1:]
        transcript, loc, exs, segs = line.split(" ")
        # print(transcript, loc, exs, segs)

        chrom = loc.split(":")[1].split("|")[0]
        strand = loc[-1]

        exs = [
            (int(s), int(e))
            for s, e in list(map(lambda x: x.split("-"), exs.split(":")[1].split(",")))
        ]
        if strand == "-":
            exs.reverse()
        exons[f"{transcript}.1"] = f"{chrom}:{exs[0][0]}-{exs[0][1]}"
        for i, (ex1, ex2) in enumerate(zip(exs[:-1], exs[1:]), 1):
            s, e = ex1[1], ex2[0]
            if strand == "-":
                s, e = ex2[1], ex1[0]
            intron = f"{chrom}:{s}-{e}"
            introns[f"{transcript}.{i}.{i+1}"] = intron
            exons[f"{transcript}.{i+1}"] = f"{chrom}:{ex2[0]}-{ex2[1]}"
    return exons, introns


def parse_pantas(fpath, exons, introns):
    pantas = {x: set() for x in ETYPES}
    for line in open(fpath):
        (
            etype,
            novel,
            chrom,
            gene,
            strand,
            i1,
            n1,
            w1,
            i2,
            n2,
            w2,
            i3,
            n3,
            w3,
        ) = line.strip("\n").split(",")
        w1, w2, w3 = int(w1), int(w2), int(w3) if w3 != "." else -1
        if etype == "IR":
            # one is exon, other is intron
            # we force i1 to be intron
            if i1 not in introns:
                i1, i2, i3 = introns[i2], exons[i1], "."
            else:
                i1, i2, i3 = introns[i1], exons[i2], "."
        else:
            # SE/A3/A5
            i1, i2, i3 = introns[i1], introns[i2], introns[i3] if i3 in introns else "."
        pantas[etype].add(
            Event(etype, novel, chrom, gene, strand, i1, i2, i3, n1, n2, n3, w1, w2, w3)
        )
    return pantas


def main():
    # FIXME: hardcoded to 1 replicate and 2 conditions, class Event should be ready btw
    splicedfa_path = sys.argv[1]
    pantas2_1_path = sys.argv[2]
    pantas2_2_path = sys.argv[3]

    exons, introns = get_referencebased_coordinates(splicedfa_path)

    events_1 = parse_pantas(pantas2_1_path, exons, introns)
    events_2 = parse_pantas(pantas2_2_path, exons, introns)

    # TODO: merge replicates from same condition

    # Merge two conditions and compute dPSI when possible
    for et, E1 in events_1.items():
        E2 = events_2[et]
        for e1 in E1:
            hit = False
            for e2 in E2:
                if e1 == e2:
                    psi1 = e1.psi
                    psi2 = e2.psi
                    if all([x == "NaN" for x in psi1]) or all(
                        [x == "NaN" for x in psi2]
                    ):
                        dpsi = "NaN"
                    else:
                        dpsi = sum(psi1) / len(psi1) - sum(psi2) / len(psi2)
                    hit = True
                    print(
                        e1.etype,
                        e1.chrom,
                        e1.gene,
                        e1.strand,
                        e1.intron1_r,
                        # e1.nodes1,
                        # e1.w1,
                        e1.intron2_r,
                        # e1.nodes2,
                        # e1.w2,
                        e1.intron3_r,
                        # e1.nodes3,
                        # e1.w3,
                        "/".join([str(x) for x in psi1]),
                        "/".join([str(x) for x in psi2]),
                        dpsi,
                    )
                    break
            if not hit:
                # we have the event in condition 1, no dPSI
                print(
                    e1.etype,
                    e1.chrom,
                    e1.gene,
                    e1.strand,
                    e1.intron1_r,
                    # e1.nodes1,
                    # e1.w1,
                    e1.intron2_r,
                    # e1.nodes2,
                    # e1.w2,
                    e1.intron3_r,
                    # e1.nodes3,
                    # e1.w3,
                    "/".join([str(x) for x in e1.psi]),
                    "NaN",
                    "NaN",
                )
        for e2 in E2:
            if e2 not in E1:
                # we have the event in condition 2, no dPSI
                print(
                    e2.etype,
                    e2.chrom,
                    e2.gene,
                    e2.strand,
                    e2.intron1_r,
                    # e2.nodes1,
                    # e2.w1,
                    e2.intron2_r,
                    # e2.nodes2,
                    # e2.w2,
                    e2.intron3_r,
                    # e2.nodes3,
                    # e2.w3,
                    "NaN",
                    "/".join([str(x) for x in e2.psi]),
                    "NaN",
                )


if __name__ == "__main__":
    main()
