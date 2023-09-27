import sys

ETYPES = ["ES", "IR", "A3", "A5"]


def get_chrom(region):
    return region.split(":")[0]


def get_interval(region):
    if region == ".":
        return region
    s, e = [int(x) for x in region.split(":")[1].split("-")]
    return s, e


class Event:
    def from_line(self, line):
        (
            self.etype,
            self.chrom,
            self.gene,
            self.strand,
            self.intron1,
            self.intron2,
            self.intron3,
            self.wrep1,
            self.wrep2,
            self.psi1,
            self.psi2,
        ) = line.strip("\n").split(",")
        self.intron1 = get_interval(self.intron1)
        self.intron2 = get_interval(self.intron2)
        self.intron3 = get_interval(self.intron3)
        self.sort()

    def from_args(self, etype, intron1, intron2, intron3):
        self.etype = etype
        self.chrom = get_chrom(intron1)
        self.gene = ""
        self.strand = ""
        self.intron1 = get_interval(intron1)
        self.intron2 = get_interval(intron2)
        self.intron3 = get_interval(intron3)
        self.wrep1 = "."
        self.wrep2 = "."
        self.psi1 = "."
        self.psi2 = "."
        self.sort()

    def sort(self):
        if self.intron3 != ".":
            self.intron1, self.intron2, self.intron3 = sorted(
                [self.intron1, self.intron2, self.intron3], key=lambda x: x[1] - x[0]
            )
        else:
            self.intron1, self.intron2 = sorted(
                [self.intron1, self.intron2], key=lambda x: x[1] - x[0]
            )

    def __repr__(self):
        return f"{self.etype} {self.chrom} {self.intron1} {self.intron2} {self.intron3}"

    def __eq__(self, o):
        return self.__hash__() == o.__hash__()

    def __hash__(self):
        return hash(self.__repr__())


def compare(e1, e2):
    return e1 == e2


def main():
    truth_path = sys.argv[1]
    splicedfa_path = sys.argv[2]
    pantas2_1_path = sys.argv[3]
    pantas2_2_path = sys.argv[4]

    truth = {}
    for line in open(truth_path):
        e = Event()
        e.from_line(line)
        if e.etype not in truth:
            truth[e.etype] = set()
        truth[e.etype].add(e)

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
            key = f"{transcript}.{i}.{i+1}"
            s, e = ex1[1], ex2[0]
            if strand == "-":
                s, e = ex2[1], ex1[0]
            intron = f"{chrom}:{s}-{e}"
            introns[key] = intron
            exons[f"{transcript}.{i+1}"] = f"{chrom}:{ex2[0]}-{ex2[1]}"

    pantas_1 = {}
    for line in open(pantas2_1_path):
        etype, novel, i1, w1, i2, w2, i3, w3 = line.strip("\n").split(",")
        etype = etype[:-1] if etype.startswith("A") else etype
        if etype == "ES":
            i1, i2, i3 = introns[i1], introns[i2], introns[i3]
        elif etype == "IR":
            # one is exon, other is intron
            # we force i1 to be intron
            if i1 not in introns:
                i1, i2, i3 = introns[i2], exons[i1], "."
            else:
                i1, i2, i3 = introns[i1], exons[i2], "."
        else:
            # A3/A5
            i1, i2, i3 = introns[i1], introns[i2], "."

        if etype not in pantas_1:
            pantas_1[etype] = set()
        e = Event()
        e.from_args(etype, i1, i2, i3)
        pantas_1[e.etype].add(e)

    for e in ETYPES:
        TP = len(pantas_1[e] & truth[e])
        FP = len(pantas_1[e] - truth[e])
        FN = len(truth[e] - pantas_1[e])
        print(e, TP, FP, FN, TP / (TP + FP), TP / (TP + FN))


if __name__ == "__main__":
    main()
