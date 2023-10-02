import sys

ETYPES = ["ES", "CE", "IR", "A3", "A5"]


def get_interval(region):
    if region in [".", "?"]:
        return region
    s, e = [int(x) for x in region.split(":")[1].split("-")]
    return s, e


class Event:
    def __init__(
        self,
        etype,
        novel,
        chrom,
        gene,
        strand,
        i1,
        i2,
        i3,
        n1,
        n2,
        n3,
        w1,
        w2,
        w3,
        rep,
        nreps,
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
        self.w1 = [0] * nreps
        self.w1[rep] = w1
        self.w2 = [0] * nreps
        self.w2[rep] = w2
        self.w3 = [0] * nreps
        self.w3[rep] = w3
        self.sort()
        self.psi = ["NaN"] * nreps
        self.compute_psi(rep)

    # Sort info
    def sort(self):
        if self.etype == "ES":
            pass  # 1 is skipping, 2-3 are inclusion
        elif self.etype in ["A3", "A5"]:
            # force 1 to be shorter intron
            if self.intron1_r == "?" or self.intron2_r == "?":
                pass
            else:
                if (
                    self.intron1[1] - self.intron1[0]
                    > self.intron2[1] - self.intron2[0]
                ):
                    self.intron1, self.intron2 = self.intron2, self.intron1
                    self.intron1_r, self.intron2_r = self.intron2_r, self.intron1_r
                    self.w1, self.w2 = self.w2, self.w1
                    self.nodes1, self.nodes2 = self.nodes2, self.nodes1
        elif self.etype == "IR":
            pass  # 1 is retaining junction, 2 is exon

    def compute_psi(self, rep):
        psi = "NaN"
        w1, w2, w3 = self.w1[rep], self.w2[rep], self.w3[rep]
        if self.etype == "ES":
            if w1 + w2 + w3 != 0:
                psi = ((w2 + w3) / 2) / (w1 + (w2 + w3) / 2)
        elif self.etype in ["A3", "A5"]:
            if w1 + w2 != 0:
                psi = w1 / (w1 + w2)
        elif self.etype == "IR":
            if w1 + w2 != 0:
                psi = w2 / (w1 + w2)
        self.psi[rep] = psi

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


def parse_pantas(fpath, rep, nreps, exons, introns):
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
        w1, w2, w3 = int(w1), int(w2), int(w3) if w3 != "." else "."
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
        if etype not in pantas:
            print(f"Skipping {etype}..", file=sys.stderr)
        else:
            pantas[etype].add(
                Event(
                    etype,
                    novel,
                    chrom,
                    gene,
                    strand,
                    i1,
                    i2,
                    i3,
                    n1,
                    n2,
                    n3,
                    w1,
                    w2,
                    w3,
                    rep,
                    nreps,
                )
            )
    return pantas


def main():
    splicedfa_path = sys.argv[1]
    pantas_paths = sys.argv[2:]
    c1_paths = pantas_paths[: int(len(pantas_paths) / 2)]
    c2_paths = pantas_paths[int(len(pantas_paths) / 2) :]

    exons, introns = get_referencebased_coordinates(splicedfa_path)
    introns["."] = "."
    introns["?"] = "?"

    events_1 = {x: set() for x in ETYPES}
    for i, fpath in enumerate(c1_paths):
        events = parse_pantas(fpath, i, len(c1_paths), exons, introns)
        for etype in events:
            for new_e in events[etype]:
                if new_e not in events_1[etype]:
                    events_1[etype].add(new_e)
                else:
                    for old_e in events_1[etype]:
                        if old_e == new_e:
                            old_e.psi[i] = new_e.psi[i]
                            old_e.w1[i] = new_e.w1[i]
                            old_e.w2[i] = new_e.w2[i]
                            old_e.w3[i] = new_e.w3[i]
                            break
    events_2 = {x: set() for x in ETYPES}
    for i, fpath in enumerate(c2_paths):
        events = parse_pantas(fpath, i, len(c2_paths), exons, introns)
        for etype in events:
            for new_e in events[etype]:
                if new_e not in events_2[etype]:
                    events_2[etype].add(new_e)
                else:
                    for old_e in events_2[etype]:
                        if old_e == new_e:
                            old_e.psi[i] = new_e.psi[i]
                            old_e.w1[i] = new_e.w1[i]
                            old_e.w2[i] = new_e.w2[i]
                            old_e.w3[i] = new_e.w3[i]
                            break

    print(
        "Event",
        "Novel",
        "Chrom",
        "Gene",
        "Strand",
        "Region1",
        "Region2",
        "Region3",
        # e1.nodes1,
        # e1.nodes2,
        # e1.nodes3,
        "W1",
        "W2",
        "psi1",
        "psi2",
        "dPSI",
        sep=",",
    )
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
                        e1.novel,
                        e1.chrom,
                        e1.gene,
                        e1.strand,
                        e1.intron1_r,
                        e1.intron2_r,
                        e1.intron3_r,
                        # e1.nodes1,
                        # e1.nodes2,
                        # e1.nodes3,
                        "-".join([str(w) for w in e1.w1])
                        + "/"
                        + "-".join([str(w) for w in e1.w2])
                        + "/"
                        + "-".join([str(w) for w in e1.w3]),
                        "-".join([str(w) for w in e2.w1])
                        + "/"
                        + "-".join([str(w) for w in e2.w2])
                        + "/"
                        + "-".join([str(w) for w in e2.w3]),
                        "/".join([str(x) for x in psi1]),
                        "/".join([str(x) for x in psi2]),
                        dpsi,
                        sep=",",
                    )
                    break
            if not hit:
                # we have the event in condition 1, no dPSI
                print(
                    e1.etype,
                    e1.novel,
                    e1.chrom,
                    e1.gene,
                    e1.strand,
                    e1.intron1_r,
                    e1.intron2_r,
                    e1.intron3_r,
                    # e1.nodes1,
                    # e1.nodes2,
                    # e1.nodes3,
                    "-".join([str(w) for w in e1.w1])
                    + "/"
                    + "-".join([str(w) for w in e1.w2])
                    + "/"
                    + "-".join([str(w) for w in e1.w3]),
                    "-".join([str(w) for w in e2.w1])
                    + "/"
                    + "-".join([str(w) for w in e2.w2])
                    + "/"
                    + "-".join([str(w) for w in e2.w3]),
                    "/".join([str(x) for x in e1.psi]),
                    "/".join(["NaN" for x in e1.psi]),
                    "NaN",
                    sep=",",
                )
        for e2 in E2:
            if e2 not in E1:
                # we have the event in condition 2, no dPSI
                print(
                    e2.etype,
                    e2.novel,
                    e2.chrom,
                    e2.gene,
                    e2.strand,
                    e2.intron1_r,
                    e2.intron2_r,
                    e2.intron3_r,
                    # e2.nodes1,
                    # e2.nodes2,
                    # e2.nodes3,
                    "-".join([str(w) for w in e1.w1])
                    + "/"
                    + "-".join([str(w) for w in e1.w2])
                    + "/"
                    + "-".join([str(w) for w in e1.w3]),
                    "-".join([str(w) for w in e2.w1])
                    + "/"
                    + "-".join([str(w) for w in e2.w2])
                    + "/"
                    + "-".join([str(w) for w in e2.w3]),
                    "/".join(["NaN" for x in e2.psi]),
                    "/".join([str(x) for x in e2.psi]),
                    "NaN",
                    sep=",",
                )


if __name__ == "__main__":
    main()
