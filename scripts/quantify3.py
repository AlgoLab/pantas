import sys
import re
from functools import partial
from statistics import mean

ETYPES = ["ES", "CE", "IR", "A3", "A5"]


def parse_region(string: str) -> tuple:
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


def calc_psi(event_cov: int, canonic_cov: int, mode: str = "canonic"):
    den = canonic_cov + event_cov
    if mode == "canonic":
        num = canonic_cov
    elif mode == "event":
        num = event_cov
    else:
        raise TypeError

    return float(num) / den if den != 0 else -1


class Event:
    def __init__(
        self,
        event_type: str,
        annotation_type: str,
        chrom: str,
        gene: str,
        strand: str,
        junction1_name: str,
        junction1_nodes: str,
        junction1_refpos: str,
        junction1_coverage: str,
        junction2_name: str,
        junction2_nodes: str,
        junction2_refpos: str,
        junction2_coverage: str,
        junction3_name: str,
        junction3_nodes: str,
        junction3_refpos: str,
        junction3_coverage: str,
        min_junction_len: int = 3,
    ):
        self.chrom = chrom
        self.etype = event_type
        self.annotation_type = annotation_type
        self.strand = strand
        self.gene = gene
        self.junction1_name = junction1_name
        self.junction1_nodes = junction1_nodes
        self.junction1_refpos = junction1_refpos
        self.junction1_coverage = (
            int(junction1_coverage) if junction1_coverage != "." else -1
        )
        self.junction2_name = junction2_name
        self.junction2_nodes = junction2_nodes
        self.junction2_refpos = junction2_refpos
        self.junction2_coverage = (
            int(junction2_coverage) if junction2_coverage != "." else -1
        )
        self.junction3_name = junction3_name
        self.junction3_nodes = junction3_nodes
        self.junction3_refpos = junction3_refpos
        self.junction3_coverage = (
            int(junction3_coverage) if junction3_coverage != "." else -1
        )

        self.valid = True
        self.replicates = []

        match self.etype:
            case "ES":
                self.event_cov = self.junction1_coverage
                self.event_j = parse_region(self.junction1_refpos)
                self.canonic_cov = (
                    self.junction2_coverage + self.junction3_coverage
                ) // 2
                self.canonic_j = [
                    parse_region(self.junction2_refpos),
                    parse_region(self.junction3_refpos),
                ]
                self.canonic_nodes = [self.junction2_nodes, self.junction3_nodes]

                self.csv_j1 = junction1_refpos
                self.csv_j2 = junction2_refpos
                self.csv_j3 = junction3_refpos

                if self.event_j[1] - self.event_j[0] < min_junction_len:
                    self.valid = False
            case "A5":
                if self.strand == "+":
                    self.event_cov = self.junction1_coverage
                    self.event_j = parse_region(self.junction1_refpos)
                    self.canonic_cov = self.junction2_coverage
                    self.canonic_j = parse_region(self.junction2_refpos)
                    self.canonic_nodes = self.junction2_nodes

                    self.csv_j1 = junction1_refpos
                    self.csv_j2 = junction2_refpos
                    self.csv_j3 = junction3_refpos
                else:
                    self.event_cov = self.junction2_coverage
                    self.event_j = parse_region(self.junction2_refpos)
                    self.canonic_cov = self.junction1_coverage
                    self.canonic_j = parse_region(self.junction1_refpos)
                    self.canonic_nodes = self.junction1_nodes

                    self.csv_j1 = junction2_refpos
                    self.csv_j2 = junction1_refpos
                    self.csv_j3 = junction3_refpos

                if self.event_j[1] - self.event_j[0] < min_junction_len:
                    self.valid = False
            case "A3":
                if self.strand == "+":
                    self.event_cov = self.junction2_coverage
                    self.event_j = parse_region(self.junction2_refpos)
                    self.canonic_cov = self.junction1_coverage
                    self.canonic_j = parse_region(self.junction1_refpos)
                    self.canonic_nodes = self.junction1_nodes

                    self.csv_j1 = junction2_refpos
                    self.csv_j2 = junction1_refpos
                    self.csv_j3 = junction3_refpos
                else:
                    self.event_cov = self.junction1_coverage
                    self.event_j = parse_region(self.junction1_refpos)
                    self.canonic_cov = self.junction2_coverage
                    self.canonic_j = parse_region(self.junction2_refpos)
                    self.canonic_nodes = self.junction2_nodes

                    self.csv_j1 = junction1_refpos
                    self.csv_j2 = junction2_refpos
                    self.csv_j3 = junction3_refpos

                if self.event_j[1] - self.event_j[0] < min_junction_len:
                    self.valid = False

            case "IR":
                if annotation_type == "novel" and self.junction2_name == "?":
                    self.event_cov = self.junction2_coverage
                    self.event_j = parse_region(self.junction2_refpos)
                    self.canonic_cov = self.junction1_coverage
                    self.canonic_j = parse_region(self.junction1_refpos)
                    self.canonic_nodes = self.junction1_nodes
                    if self.canonic_j[1] - self.canonic_j[0] < min_junction_len:
                        self.valid = False

                    self.csv_j1 = junction2_refpos
                    self.csv_j2 = junction1_refpos
                    self.csv_j3 = junction3_refpos
                else:
                    self.event_cov = self.junction1_coverage
                    self.event_j = parse_region(self.junction1_refpos)
                    self.canonic_cov = self.junction2_coverage
                    self.canonic_j = parse_region(self.junction2_refpos)
                    self.canonic_nodes = self.junction2_nodes

                    if self.event_j[1] - self.event_j[0] < min_junction_len:
                        self.valid = False

                    self.csv_j1 = junction1_refpos
                    self.csv_j2 = junction2_refpos
                    self.csv_j3 = junction3_refpos
            case "CE":
                self.event_cov = (
                    self.junction2_coverage + self.junction3_coverage
                ) // 2
                self.event_j = [
                    parse_region(self.junction2_refpos),
                    parse_region(self.junction3_refpos),
                ]
                self.canonic_cov = self.junction1_coverage
                self.canonic_j = parse_region(self.junction1_refpos)
                self.canonic_nodes = self.junction1_nodes

                if self.event_j[0][1] - self.event_j[0][0] < min_junction_len:
                    self.valid = False
                if self.event_j[1][1] - self.event_j[1][0] < min_junction_len:
                    self.valid = False

                self.csv_j1 = junction2_refpos
                self.csv_j2 = junction3_refpos
                self.csv_j3 = junction1_refpos

        self.add_replicate(self.event_cov, self.canonic_cov)

    def __repr__(self) -> str:
        return f"({self.etype}) E: {self.event_j} [{self.event_cov}] --- T {self.canonic_j} [{self.canonic_cov}]"

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
                    self.csv_j1,
                    # self.junction1_coverage,
                    self.csv_j2,
                    # self.junction2_coverage,
                    self.csv_j3,
                    # self.junction3_coverage,
                ],
            )
        )

    def add_replicate(self, event_cov: int, canonical_cov: int):
        self.replicates.append([event_cov, canonical_cov])

    def psi(self, mode: str = "canonic"):
        _psis = [calc_psi(*x, mode) for x in self.replicates]
        if all([x == -1 for x in _psis]):
            return -1
        return mean([x for x in _psis if x != -1])

    def get_event_cov(self):
        return int(mean([x[0] for x in self.replicates]))

    def get_canonic_cov(self):
        return int(mean([x[1] for x in self.replicates]))


def eq_event(e1: Event, e2: Event, relax: int = 0):
    if e1.annotation_type == "annotated" and e2.annotation_type == "annotated":
        relax = 0
    if e1.etype != e2.etype:
        return False
    if e1.gene != e2.gene:
        return False
    if e1.etype == "CE":
        dt0 = abs(e1.canonic_j[0] - e2.canonic_j[0]) <= relax
        dt1 = abs(e1.canonic_j[1] - e2.canonic_j[1]) <= relax
        de00 = abs(e1.event_j[0][0] - e2.event_j[0][0]) <= relax
        de01 = abs(e1.event_j[0][1] - e2.event_j[0][1]) <= relax
        de10 = abs(e1.event_j[1][0] - e2.event_j[1][0]) <= relax
        de11 = abs(e1.event_j[1][1] - e2.event_j[1][1]) <= relax
        return dt0 & dt1 & de00 & de01 & de10 & de11
    elif e1.etype == "ES":
        de0 = abs(e1.event_j[0] - e2.event_j[0]) <= relax
        de1 = abs(e1.event_j[1] - e2.event_j[1]) <= relax
        dt00 = abs(e1.canonic_j[0][0] - e2.canonic_j[0][0]) <= relax
        dt01 = abs(e1.canonic_j[0][1] - e2.canonic_j[0][1]) <= relax
        dt10 = abs(e1.canonic_j[1][0] - e2.canonic_j[1][0]) <= relax
        dt11 = abs(e1.canonic_j[1][1] - e2.canonic_j[1][1]) <= relax
        return de0 & de1 & dt00 & dt01 & dt10 & dt11
    elif e1.etype == "IR":
        if e1.event_j == e2.event_j == ".":
            dt0 = abs(e1.canonic_j[0] - e2.canonic_j[0]) <= relax
            dt1 = abs(e1.canonic_j[1] - e2.canonic_j[1]) <= relax
            return dt0 & dt1
        elif e1.canonic_j == e2.canonic_j == ".":
            de0 = abs(e1.event_j[0] - e2.event_j[0]) <= relax
            de1 = abs(e1.event_j[1] - e2.event_j[1]) <= relax
            return de0 & de1
        else:
            return False
    else:
        dt0 = abs(e1.canonic_j[0] - e2.canonic_j[0]) <= relax
        dt1 = abs(e1.canonic_j[1] - e2.canonic_j[1]) <= relax
        de0 = abs(e1.event_j[0] - e2.event_j[0]) <= relax
        de1 = abs(e1.event_j[1] - e2.event_j[1]) <= relax
        return dt0 & dt1 & de0 & de1


def parse_pantas(fpath: str, min_junction_len: int = 3) -> dict[str, Event]:
    pantas = {x: [] for x in ETYPES}
    for line in open(fpath):
        line = line.strip()
        # print(line, file=sys.stderr)
        if "-?" in line or "?-" in line:
            continue
        _e = line.split(",")
        e = Event(*_e, min_junction_len=min_junction_len)
        if e.valid:
            pantas[e.etype].append(e)
        # else:
        #     print("NOTVALID", e.to_csv())
    return pantas


def load_gfa(fpath: str, gfa: dict, cond: int) -> None:
    for line in open(fpath, "r"):
        if line.startswith("L") and not "RC:i:0" in line:
            line = line.strip()
            _, start, _, end, _, _, *ann = line.split()
            rc = [int(x.split(":")[-1]) for x in ann if "RC" in x]
            assert len(rc) == 1
            rc = rc[0]
            key = f"{start}>{end}"
            if not key in gfa:
                gfa[key] = [0, 0]
                gfa[key][cond - 1] = rc
            else:
                gfa[key][cond - 1] += rc


def main(args):
    RELAX = args.relax

    GFA = dict()
    if args.gfac1:
        for _gfa in args.gfac1:
            load_gfa(_gfa, GFA, 1)
    if args.gfac2:
        for _gfa in args.gfac2:
            load_gfa(_gfa, GFA, 2)

    events_1 = {x: [] for x in ETYPES}
    for i, fpath in enumerate(args.c1):
        _ei = parse_pantas(fpath, min_junction_len=args.min_junction_len)
        for etype in ETYPES:
            # print(i, etype, len(events_1[etype]))
            for _ee in _ei[etype]:
                eqs = [
                    x
                    for x in events_1[etype]
                    if eq_event(_ee, x, relax=RELAX * max(1, (len(args.c1) - 1)))
                ]
                if len(eqs) > 0:
                    # assert len(eqs) == 1
                    eqs[0].add_replicate(_ee.event_cov, _ee.canonic_cov)
                else:
                    events_1[etype].append(_ee)
            # print(i, etype, len(events_1[etype]))

    events_2 = {x: [] for x in ETYPES}
    for i, fpath in enumerate(args.c2):
        _ei = parse_pantas(fpath, min_junction_len=args.min_junction_len)
        for etype in ETYPES:
            # print(i, etype, len(events_2[etype]))
            for _ee in _ei[etype]:
                eqs = [
                    x
                    for x in events_2[etype]
                    if eq_event(_ee, x, relax=RELAX * max(1, (len(args.c1) - 1)))
                ]
                if len(eqs) > 0:
                    # assert len(eqs) == 1
                    eqs[0].add_replicate(_ee.event_cov, _ee.canonic_cov)
                else:
                    events_2[etype].append(_ee)
            # print(i, etype, len(events_2[etype]))

    print(
        "etype",
        "annotation_type",
        "chrom",
        "gene",
        "strand",
        "junction1_refpos",
        # "junction1_coverage",
        "junction2_refpos",
        # "junction2_coverage",
        "junction3_refpos",
        "W1",
        "W2",
        "psi_c1",
        "psi_c2",
        "dpsi",
        sep=",",
    )

    # Get delta psi for two conditions
    for etype in ETYPES:
        for e1 in events_1[etype]:
            # eqs = partial(eq_event, e1, relax=0)
            # print(e1)
            eqs = [
                x
                for x in events_2[etype]
                if eq_event(e1, x, relax=RELAX * max(1, (len(args.c1) - 1)))
            ]
            if len(eqs) > 0:
                # assert len(eqs) == 1
                psi1 = e1.psi()
                psi2 = eqs[0].psi()
                dpsi = max(0, psi1) - max(0, psi2)
                if psi1 == -1 and psi2 == -1:
                    dpsi = -1
                if abs(dpsi) < args.mindpsi:
                    continue
                print(
                    e1.to_csv(),
                    f"{e1.get_canonic_cov()}/{e1.get_event_cov()}",
                    f"{eqs[0].get_canonic_cov()}/{eqs[0].get_event_cov()}",
                    psi1,
                    psi2,
                    dpsi,
                    sep=",",
                )
            else:
                if not e1.psi() == -1:
                    if e1.get_event_cov() < args.minrc:
                        continue
                    psi2 = "NaN"
                    dpsi = "NaN"
                    w2 = "."
                    if args.gfac2:
                        s = sum(ws := [GFA.get(k, [0, 0])[1] for k in e1.canonic_nodes])
                        if s > 0:
                            psi2 = 1
                            dpsi = 1 - e1.psi()
                        w2 = f"{mean(ws)}/0"
                    print(
                        e1.to_csv(),
                        f"{e1.get_canonic_cov()}/{e1.get_event_cov()}",
                        w2,
                        e1.psi(),
                        psi2,
                        dpsi,
                        sep=",",
                    )
        for e2 in events_2[etype]:
            eqs = [
                x
                for x in events_1[etype]
                if eq_event(e2, x, relax=RELAX * max(1, (len(args.c1) - 1)))
            ]
            if len(eqs) == 0 and not e2.psi() == -1:
                if e2.get_event_cov() < args.minrc:
                    continue
                psi1 = "NaN"
                dpsi = "NaN"
                w1 = "."
                if args.gfac1:
                    s = sum(ws := [GFA.get(k, [0, 0])[0] for k in e2.canonic_nodes])
                    if s > 0:
                        psi1 = 1
                        dpsi = 1 - e2.psi()
                        w1 = f"{mean(ws)}/0"
                print(
                    e2.to_csv(),
                    w1,
                    f"{e2.get_canonic_cov()}/{e2.get_event_cov()}",
                    psi1,
                    e2.psi(),
                    dpsi,
                    sep=",",
                )

    # for etype in ETYPES:
    #     for es1 in events_1_matched[etype]:
    #         es1eqs = []
    #         for e1 in es1:
    #             eqs = [x for x in events_2_matched[etype] if any([eq_event(e1, y, relax=RELAX) for y in x])]
    #             print(eqs)
    #             es1eqs += eqs
    #         if len(es1eqs) > 0:
    #             print(es1eqs)
    #             print(
    #                 es1[0].to_csv(),
    #                 mps1:=mean([x.psi() for x in es1]),
    #                 mps2:=mean([x.psi() for x in es1eqs]),
    #                 mps1-mps2,
    #                 sep=","
    #             )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        prog="Quantify",
        description="",
    )
    parser.add_argument(
        "-c1",
        help="Replicates of condition 1",
        dest="c1",
        type=str,
        required=True,
        nargs="+",
    )
    parser.add_argument(
        "-c2",
        help="Replicates of condition 2",
        dest="c2",
        type=str,
        required=True,
        nargs="+",
    )
    parser.add_argument(
        "--relax",
        dest="relax",
        help="Relaxation of reference positions matching (Default: 0)",
        type=int,
        default=0,
    )
    parser.add_argument(
        "--minj",
        dest="min_junction_len",
        help="Minimum number of junction length to be valid (Default: 3)",
        type=int,
        default=3,
    )
    parser.add_argument(
        "--mindpsi",
        dest="mindpsi",
        help="Minimum value of delta-psi (absolute value) to be valid (Default: 0.0)",
        type=float,
        default=0.0,
    )
    parser.add_argument(
        "--minrc",
        dest="minrc",
        help="Minimum value of read count of the event to be valid (Default: 0)",
        type=int,
        default=0,
    )
    parser.add_argument(
        "--gfac1",
        help="Annotated spliced pangenome of condition 1",
        dest="gfac1",
        type=str,
        required=False,
        nargs="*",
    )
    parser.add_argument(
        "--gfac2",
        help="Annotated spliced pangenome of condition 2",
        dest="gfac2",
        type=str,
        required=False,
        nargs="*",
    )
    args = parser.parse_args()

    if not len(args.c1) == len(args.c2):
        print("Provide the same number of replicates for each condition")
    else:
        main(args)
