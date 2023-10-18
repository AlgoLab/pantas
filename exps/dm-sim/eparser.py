import sys
import re
from abc import abstractmethod, ABC

ETYPES = ["ES", "IR", "A3", "A5", "CE"]


def parse_region(string: str) -> list[int]:
    if string == ".":
        return string
    reg = re.match(r"(?P<chr>[\w\d]+):(?P<start>\d+)-(?P<end>\d+)", string)
    if not reg:
        print(
            f"Unable to read region {string}. Ignoring it",
            file=sys.stderr,
        )
        sys.exit(1)
    return [int(reg.group("start")), int(reg.group("end"))]


def fix_region(reg: list[int]) -> list[int]:
    return [reg[0] + 1, reg[1] - 1]


def build_region(regions):
    if regions == ".":
        return "."
    elif type(regions[0]) == int:
        return f"{regions[0]}-{regions[1]}"
    else:
        return ",".join([f"{r[0]}-{r[1]}" for r in regions])


class Event(ABC):
    def __init__(
        self,
        event_type: str,
        annotation_type: str,
        chrom: str,
        gene: str,
        strand: str,
        junction1_refpos: str,
        junction2_refpos: str,
        junction3_refpos: str,
        W1: str,
        W2: str,
        psi_c1: str,
        psi_c2: str,
        dpsi: str,
    ):
        self.chrom = chrom
        self.etype = event_type
        self.annotation_type = annotation_type
        self.strand = strand
        self.gene = gene
        self.junction1_refpos = junction1_refpos
        self.junction2_refpos = junction2_refpos
        self.junction3_refpos = junction3_refpos
        self.psi_c1 = float(psi_c1)
        self.psi_c2 = float(psi_c2)
        self.dpsi = float(dpsi)

        self.event_j = ""
        self.canonic_j = ""

        self.build_conditions()

    @abstractmethod
    def build_conditions():
        pass

    def __repr__(self) -> str:
        return f"({self.etype}) E: {self.event_j} [{self.psi_c1}] --- T {self.canonic_j} [{self.psi_c2}]"

    def to_csv(self) -> str:
        return ",".join(
            map(
                str,
                [
                    self.etype,
                    self.annotation_type,
                    self.chrom,
                    self.gene,
                    self.strand,
                    f"{self.chrom}:{build_region(self.event_j)}",
                    f"{self.chrom}:{build_region(self.canonic_j)}",
                    self.psi_c1,
                    self.psi_c2,
                    self.dpsi,
                ],
            )
        )


class EventPantas(Event):
    def __init__(
        self,
        event_type: str,
        annotation_type: str,
        chrom: str,
        gene: str,
        strand: str,
        junction1_refpos: str,
        junction2_refpos: str,
        junction3_refpos: str,
        W1: str,
        W2: str,
        psi_c1: str,
        psi_c2: str,
        dpsi: str,
    ):
        super().__init__(
            event_type,
            annotation_type,
            chrom,
            gene,
            strand,
            junction1_refpos,
            junction2_refpos,
            junction3_refpos,
            W1,
            W2,
            psi_c1,
            psi_c2,
            dpsi,
        )

    def build_conditions(self):
        match self.etype:
            case "ES":
                self.event_j = parse_region(self.junction1_refpos)
                self.canonic_j = [
                    parse_region(self.junction2_refpos),
                    parse_region(self.junction3_refpos),
                ]

            case "A5":
                self.event_j = parse_region(self.junction1_refpos)
                self.canonic_j = parse_region(self.junction2_refpos)

            case "A3":
                self.event_j = parse_region(self.junction1_refpos)
                self.canonic_j = parse_region(self.junction2_refpos)

            case "IR":
                self.event_j = parse_region(self.junction1_refpos)
                self.canonic_j = parse_region(self.junction2_refpos)

            case "CE":
                self.event_j = [
                    parse_region(self.junction2_refpos),
                    parse_region(self.junction3_refpos),
                ]
                self.canonic_j = parse_region(self.junction1_refpos)


class EventTruth(Event):
    def __init__(
        self,
        event_type: str,
        annotation_type: str,
        chrom: str,
        gene: str,
        strand: str,
        junction1_refpos: str,
        junction2_refpos: str,
        junction3_refpos: str,
        W1: str,
        W2: str,
        psi_c1: str,
        psi_c2: str,
        dpsi: str,
    ):
        super().__init__(
            event_type,
            annotation_type,
            chrom,
            gene,
            strand,
            junction1_refpos,
            junction2_refpos,
            junction3_refpos,
            W1,
            W2,
            psi_c1,
            psi_c2,
            dpsi,
        )
        self.rc_c1 = list(map(int, W1.split("/")))
        self.rc_c2 = list(map(int, W2.split("/")))
        if event_type == "ES":
            self.event_cov_c1 = self.rc_c1[2]
            self.event_cov_c2 = self.rc_c2[2]
        else:
            self.event_cov_c1 = self.rc_c1[1]
            self.event_cov_c2 = self.rc_c2[1]

        self.min_event_cov = min(self.event_cov_c1, self.event_cov_c2)

    def build_conditions(self):
        match self.etype:
            case "ES":
                self.event_j = fix_region(parse_region(self.junction3_refpos))
                self.canonic_j = [
                    fix_region(parse_region(self.junction1_refpos)),
                    fix_region(parse_region(self.junction2_refpos)),
                ]

            case "A5":
                self.event_j = fix_region(parse_region(self.junction2_refpos))
                self.canonic_j = fix_region(parse_region(self.junction1_refpos))

            case "A3":
                self.event_j = fix_region(parse_region(self.junction2_refpos))
                self.canonic_j = fix_region(parse_region(self.junction1_refpos))

            case "IR":
                self.event_j = fix_region(parse_region(self.junction1_refpos))
                # self.canonic_j = fix_region(parse_region(self.junction2_refpos))
                self.canonic_j = "."

            case "CE":
                # TODO: fix
                self.event_j = [
                    parse_region(self.junction2_refpos),
                    parse_region(self.junction3_refpos),
                ]
                self.canonic_j = parse_region(self.junction1_refpos)


def eq_event(e1: Event, e2: Event, relax: int = 0, reverse: bool = False):
    if e1.etype != e2.etype:
        return False
    if e1.gene != e2.gene:
        return False
    e1_canonic_j = e1.canonic_j
    e1_event_j = e1.event_j
    e2_canonic_j = e2.canonic_j
    e2_event_j = e2.event_j
    if reverse:
        e1_canonic_j = e1.event_j
        e1_event_j = e1.canonic_j
        e2_canonic_j = e2.event_j
        e2_event_j = e2.canonic_j
    if e1.etype == "CE":
        dt0 = abs(e1_canonic_j[0] - e2_canonic_j[0]) <= relax
        dt1 = abs(e1_canonic_j[1] - e2_canonic_j[1]) <= relax
        de00 = abs(e1_event_j[0][0] - e2_event_j[0][0]) <= relax
        de01 = abs(e1_event_j[0][1] - e2_event_j[0][1]) <= relax
        de10 = abs(e1_event_j[1][0] - e2_event_j[1][0]) <= relax
        de11 = abs(e1_event_j[1][1] - e2_event_j[1][1]) <= relax
        return dt0 & dt1 & de00 & de01 & de10 & de11
    elif e1.etype == "ES":
        print("ev", e1_event_j, e2_event_j)
        print("can", e1_canonic_j, e2_canonic_j)
        de0 = abs(e1_event_j[0] - e2_event_j[0]) <= relax
        de1 = abs(e1_event_j[1] - e2_event_j[1]) <= relax
        dt00 = abs(e1_canonic_j[0][0] - e2_canonic_j[0][0]) <= relax
        dt01 = abs(e1_canonic_j[0][1] - e2_canonic_j[0][1]) <= relax
        dt10 = abs(e1_canonic_j[1][0] - e2_canonic_j[1][0]) <= relax
        dt11 = abs(e1_canonic_j[1][1] - e2_canonic_j[1][1]) <= relax
        return de0 & de1 & dt00 & dt01 & dt10 & dt11
    elif e1.etype == "IR":
        # print(e1_event_j, e2_event_j)
        de0 = abs(e1_event_j[0] - e2_event_j[0]) <= relax
        de1 = abs(e1_event_j[1] - e2_event_j[1]) <= relax
        return de0 & de1
    else:
        dt0 = abs(e1_canonic_j[0] - e2_canonic_j[0]) <= relax
        dt1 = abs(e1_canonic_j[1] - e2_canonic_j[1]) <= relax
        de0 = abs(e1_event_j[0] - e2_event_j[0]) <= relax
        de1 = abs(e1_event_j[1] - e2_event_j[1]) <= relax
        return dt0 & dt1 & de0 & de1


# def get_interval(region):
#     if region == ".":
#         return region
#     s, e = [int(x) if x != "?" else -1 for x in region.split(":")[1].split("-")]
#     return s, e

# def parse_truth(truth_path, novel):
#     truth = {x: set() for x in ETYPES}
#     truth_w = {x: {} for x in ETYPES}
#     truth_psi = {x: dict() for x in ETYPES}

#     for line in open(truth_path):
#         etype, chrom, gene, strand, i1, i2, i3, W1, W2, psi1, psi2 = line.strip(
#             "\n"
#         ).split(",")
#         if psi1 == "NaN" or psi2 == "NaN":
#             continue
#         k = None
#         if novel:
#             if W1[-2:] == "/0" or W2[-2:] == "/0":
#                 continue
#             if i3 == ".":
#                 k = (chrom, i1, i2)
#             else:
#                 k = (chrom, i1, i2, i3)
#         else:
#             i1 = get_interval(i1)
#             i2 = get_interval(i2)
#             if etype == "ES":
#                 k = f"{chrom}:{i1[0]}-{i1[1]}-{i2[0]}-{i2[1]}"
#             elif etype[0] == "A":
#                 if i1[0] == i2[0]:
#                     k = f"{chrom}:{i1[0]}-{min(i1[1], i2[1])}-{max(i1[1], i2[1])}"
#                 else:
#                     k = f"{chrom}:{min(i1[0], i2[0])}-{max(i1[0], i2[0])}-{i1[1]}"
#             elif etype == "IR":
#                 k = f"{chrom}:{i1[0]}-{i1[1]}"
#         truth[etype].add(k)
#         truth_w[etype][k] = (W1, W2)
#         assert not k in truth_psi[etype]
#         truth_psi[etype][k] = (float(psi1), float(psi2))

#     return truth, truth_w, truth_psi
