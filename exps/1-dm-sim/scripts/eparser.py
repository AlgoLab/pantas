import sys
import re
from abc import abstractmethod, ABC

ETYPES = ["ES", "IR", "A3", "A5", "CE"]


def parse_region(string: str) -> list[int]:
    if string == "." or string == "?":
        return "."
    if string.endswith("?"):
        string = string[:-1]
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
    elif regions[0] == ".":
        return f"{regions[1][0]}-{regions[1][1]}"
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
        htype: str,
        chrom: str,
        gene: str,
        strand: str,
        junction1_name: str,
        junction2_name: str,
        junction3_name: str,
        junction1_nodes: str,
        junction2_nodes: str,
        junction3_nodes: str,
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
                self.event_j = parse_region(self.junction3_refpos)
                self.canonic_j = [
                    parse_region(self.junction1_refpos),
                    parse_region(self.junction2_refpos),
                ]

            case "A5":
                self.event_j = parse_region(self.junction2_refpos)
                self.canonic_j = parse_region(self.junction1_refpos)

            case "A3":
                self.event_j = parse_region(self.junction2_refpos)
                self.canonic_j = parse_region(self.junction1_refpos)

            case "IR":
                self.event_j = parse_region(self.junction2_refpos)
                self.canonic_j = parse_region(self.junction1_refpos)
                if self.event_j == ".":
                    self.event_j = self.canonic_j
                    self.canonic_j = "."

            case "CE":
                self.event_j = [
                    parse_region(self.junction2_refpos),
                    parse_region(self.junction3_refpos),
                ]
                self.canonic_j = parse_region(self.junction1_refpos)

class EventRmats(Event):
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
                if self.event_j == ".":
                    self.event_j = self.canonic_j
                    self.canonic_j = "."

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

        self.min_event_cov = [0]
        if event_type == "ES":
            self.min_event_cov = [self.rc_c1[2], self.rc_c2[2]]
        elif event_type == "IR":
            # this works for annotated events
            # CHECKME for novels
            self.min_event_cov = [self.rc_c1[0], self.rc_c2[0]]
        else:
            self.min_event_cov = [self.rc_c1[1], self.rc_c2[1]]
        self.min_event_cov = min(self.min_event_cov)

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
                self.canonic_j = fix_region(parse_region(self.junction2_refpos))
            case "CE":
                # TODO: fix
                self.event_j = [
                    parse_region(self.junction2_refpos),
                    parse_region(self.junction3_refpos),
                ]
                self.canonic_j = parse_region(self.junction1_refpos)


class EventWhippet(Event):
    def __init__(
        self,
        gene: str,
        t1: str,
        junction1_refpos: str,
        strand: str,
        event_type: str,
        psi_c1: str,
        psi_c2: str,
        dpsi: str,
        t2: str,
        t3: str,
        t4: str,
        annotation_type: str,
    ):
        super().__init__(
            event_type,
            annotation_type,
            junction1_refpos.split(":")[0],
            gene,
            strand,
            junction1_refpos,
            ".",
            ".",
            "W1",
            "W2",
            psi_c1,
            psi_c2,
            dpsi,
        )

    def build_conditions(self):
        match self.etype:
            case "ES":
                self.event_j = parse_region(self.junction1_refpos)
                # self.canonic_j = [math.nan, math.nan]
                self.canonic_j = "."
            case "A5":
                self.event_j = parse_region(self.junction1_refpos)
                self.canonic_j = "."

            case "A3":
                self.event_j = parse_region(self.junction1_refpos)
                self.canonic_j = "."

            case "IR":
                reg = parse_region(self.junction1_refpos)
                self.event_j = [reg[0] - 1, reg[1] + 1]
                self.canonic_j = "."

            case "CE":
                self.event_j = parse_region(self.junction1_refpos)
                self.canonic_j = "."


def eq_event_anno(e1: Event, e2: Event, relax: int = 0):
    if e1.etype != e2.etype:
        return False
    if e1.gene != e2.gene:
        return False
    if type(e1).__name__ == "EventWhippet" and type(e2).__name__ != "EventWhippet":
        e1_canonic_j = e2.canonic_j
        e1_event_j = e1.event_j
        e2_canonic_j = e2.canonic_j
        e2_event_j = e2.event_j
    elif type(e1).__name__ != "EventWhippet" and type(e2).__name__ == "EventWhippet":
        e1_canonic_j = e1.canonic_j
        e1_event_j = e1.event_j
        e2_canonic_j = e1.canonic_j
        e2_event_j = e2.event_j
    else:
        e1_canonic_j = e1.canonic_j
        e1_event_j = e1.event_j
        e2_canonic_j = e2.canonic_j
        e2_event_j = e2.event_j

    if e1.etype == "CE":
        dt0 = abs(e1_canonic_j[0] - e2_canonic_j[0]) <= relax
        dt1 = abs(e1_canonic_j[1] - e2_canonic_j[1]) <= relax
        de00 = abs(e1_event_j[0][0] - e2_event_j[0][0]) <= relax
        de01 = abs(e1_event_j[0][1] - e2_event_j[0][1]) <= relax
        de10 = abs(e1_event_j[1][0] - e2_event_j[1][0]) <= relax
        de11 = abs(e1_event_j[1][1] - e2_event_j[1][1]) <= relax

        return dt0 & dt1 & de00 & de01 & de10 & de11
    elif e1.etype == "ES":
        if type(e1).__name__ == "EventWhippet" and type(e2).__name__ != "EventWhippet":
            dt00 = abs(e1_event_j[0] - 1 - e2_canonic_j[0][1]) <= relax
            dt01 = abs(e1_event_j[1] + 1 - e2_canonic_j[1][0]) <= relax
            return dt00 & dt01
        elif (
            type(e1).__name__ != "EventWhippet" and type(e2).__name__ == "EventWhippet"
        ):
            dt00 = abs(e2_event_j[0] - 1 - e1_canonic_j[0][1]) <= relax
            dt01 = abs(e2_event_j[1] + 1 - e1_canonic_j[1][0]) <= relax
            return dt00 & dt01
        else:
            de0 = abs(e1_event_j[0] - e2_event_j[0]) <= relax
            de1 = abs(e1_event_j[1] - e2_event_j[1]) <= relax
            dt00 = abs(e1_canonic_j[0][0] - e2_canonic_j[0][0]) <= relax
            dt01 = abs(e1_canonic_j[0][1] - e2_canonic_j[0][1]) <= relax
            dt10 = abs(e1_canonic_j[1][0] - e2_canonic_j[1][0]) <= relax
            dt11 = abs(e1_canonic_j[1][1] - e2_canonic_j[1][1]) <= relax

            return de0 & de1 & dt00 & dt01 & dt10 & dt11
    elif e1.etype == "IR":
        if type(e1).__name__ == "EventWhippet" and type(e2).__name__ != "EventWhippet":
            dt00 = abs(e1_event_j[0] + 1 - e2_event_j[0]) <= relax
            dt11 = abs(e1_event_j[1] - 1 - e2_event_j[1]) <= relax
            return dt00 & dt11
        elif (
            type(e1).__name__ != "EventWhippet" and type(e2).__name__ == "EventWhippet"
        ):
            dt00 = abs(e2_event_j[0] - e1_event_j[0] + 1) <= relax
            dt11 = abs(e2_event_j[1] - e1_event_j[1] - 1) <= relax
            return dt00 & dt11
        else:
            de0 = abs(e1_event_j[0] - e2_event_j[0]) <= relax
            de1 = abs(e1_event_j[1] - e2_event_j[1]) <= relax
            return de0 & de1
    else:
        if type(e1).__name__ == "EventWhippet" and type(e2).__name__ != "EventWhippet":
            dt0 = abs(e1_event_j[0] - e2_event_j[0]) <= relax
            dt1 = abs(e1_event_j[1] + 1 - e2.canonic_j[0]) <= relax
            return dt0 & dt1

        elif (
            type(e1).__name__ != "EventWhippet" and type(e2).__name__ == "EventWhippet"
        ):
            dt0 = abs(e2_event_j[0] - e1_event_j[0]) <= relax
            dt1 = abs(e2_event_j[1] + 1 - e1.canonic_j[0]) <= relax
            return dt0 & dt1

        else:
            dt0 = abs(e1_canonic_j[0] - e2_canonic_j[0]) <= relax
            dt1 = abs(e1_canonic_j[1] - e2_canonic_j[1]) <= relax
            de0 = abs(e1_event_j[0] - e2_event_j[0]) <= relax
            de1 = abs(e1_event_j[1] - e2_event_j[1]) <= relax
            return dt0 & dt1 & de0 & de1

# First event is always truth!
def eq_event_novel(e1: Event, e2: Event, print_flag: bool = False, relax: int = 0):
    # this works only for pantas and rmats. No whippet support. SUPPA does not report novel
    if e1.etype != e2.etype:
        return False
    if e1.gene != e2.gene:
        return False
    # if print_flag:
    #     print(e1.etype)
    #     print(e1.canonic_j, e1.event_j)
    #     print(e2.canonic_j, e2.event_j)
    #     print("")

    if e1.etype == "CE":
        assert False, "We have a novel cassete exon!"
    elif e1.etype == "ES":
        return e1.canonic_j == e2.canonic_j
    elif e1.etype == "IR":
        if e1.canonic_j == "." or e1.event_j == ".":
            e1_j = e1.canonic_j if e1.event_j == "." else e1.event_j
            return e1_j == e2.event_j
        elif e2.canonic_j == "." or e2.event_j == ".":
            e2_j = e2.canonic_j if e2.event_j == "." else e2.event_j
            return e2_j == e1.event_j
        else:
            assert False, "Compare novel IR, why are we here?"
    else:
        e2_j = e2.canonic_j if e2.event_j == "." else e2.event_j
        return e1.canonic_j == e2_j or e1.event_j == e2_j

def eq_event(e1: Event, e2: Event, novel: bool, print_flag=False):
    if novel:
        return eq_event_novel(e1, e2, print_flag=print_flag, relax=0)
    else:
        return eq_event_anno(e1, e2, relax=0)