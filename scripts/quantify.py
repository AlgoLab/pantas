import sys
import re
from functools import partial
from statistics import mean

ETYPES = ["ES", "CE", "IR", "A3", "A5"]


def parse_nodes(string: str) -> list:
    return [
        int(x) for x in string.split(">") if x != ".." and x != "?"
    ]  # we have ".." in long IR


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
        junction2_name: str,
        junction3_name: str,
        junction1_nodes: str,
        junction1_coverage: str,
        junction2_nodes: str,
        junction2_coverage: str,
        junction3_nodes: str,
        junction3_coverage: str,
        replicate: int,
    ):
        self.chrom = chrom
        self.etype = event_type
        self.annotation_type = annotation_type
        self.strand = strand
        self.gene = gene

        junction1_coverage = (
            float(junction1_coverage) if junction1_coverage != "." else -1
        )
        junction2_coverage = (
            float(junction2_coverage) if junction2_coverage != "." else -1
        )
        junction3_coverage = (
            float(junction3_coverage) if junction3_coverage != "." else -1
        )

        self.replicates = []

        # canonic_jname and event_jname are sets since we can have the same "junction"
        # on multiple haplotype-aware transcripts. We use the update function to
        # update the sets. When printing (to_csv()), we prioritize transcripts on reference "_R1"

        if self.etype == "ES":
            self.event_cov = junction1_coverage
            self.event_j = parse_nodes(junction1_nodes)
            self.canonic_cov = (junction2_coverage + junction3_coverage) // 2
            self.canonic_j = [
                parse_nodes(junction2_nodes),
                parse_nodes(junction3_nodes),
            ]
            self.event_nodes = [junction1_nodes]
            self.canonic_nodes = sorted(
                [junction2_nodes, junction3_nodes]
            )  # CHECKME: this should work since ids are topological sorted
            self.event_jname = [junction1_name]
            self.canonic_jname = [junction2_name, junction3_name]
        elif self.etype == "A5":
            if self.strand == "+":
                self.event_cov = junction1_coverage
                self.event_j = parse_nodes(junction1_nodes)
                self.canonic_cov = junction2_coverage
                self.canonic_j = parse_nodes(junction2_nodes)
                self.event_nodes = [junction1_nodes, "."]
                self.canonic_nodes = [junction2_nodes]
                self.event_jname = [junction1_name, "."]
                self.canonic_jname = [junction2_name]
            else:
                self.event_cov = junction2_coverage
                self.event_j = parse_nodes(junction2_nodes)
                self.canonic_cov = junction1_coverage
                self.canonic_j = parse_nodes(junction1_nodes)
                self.event_nodes = [junction2_nodes, "."]
                self.canonic_nodes = [junction1_nodes]
                self.event_jname = [junction2_name, "."]
                self.canonic_jname = [junction1_name]
        elif self.etype == "A3":
            if self.strand == "+":
                self.event_cov = junction2_coverage
                self.event_j = parse_nodes(junction2_nodes)
                self.canonic_cov = junction1_coverage
                self.canonic_j = parse_nodes(junction1_nodes)
                self.event_nodes = [junction2_nodes, "."]
                self.canonic_nodes = [junction1_nodes]
                self.event_jname = [junction2_name, "."]
                self.canonic_jname = [junction1_name]
            else:
                self.event_cov = junction1_coverage
                self.event_j = parse_nodes(junction1_nodes)
                self.canonic_cov = junction2_coverage
                self.canonic_j = parse_nodes(junction2_nodes)
                self.event_nodes = [junction1_nodes, "."]
                self.canonic_nodes = [junction2_nodes]
                self.event_jname = [junction1_name, "."]
                self.canonic_jname = [junction2_name]
        elif self.etype == "IR":
            if annotation_type == "novel" and junction2_name == "?":
                self.event_cov = junction2_coverage
                self.event_j = parse_nodes(junction2_nodes)
                self.canonic_cov = junction1_coverage
                self.canonic_j = parse_nodes(junction1_nodes)
                self.event_nodes = [junction2_nodes, "."]
                self.canonic_nodes = [junction1_nodes]
                self.event_jname = [junction2_name, "."]
                self.canonic_jname = [junction1_name]
            else:
                self.event_cov = junction1_coverage
                self.event_j = parse_nodes(junction1_nodes)
                self.canonic_cov = junction2_coverage
                self.canonic_j = parse_nodes(junction2_nodes)
                self.event_nodes = [junction1_nodes, "."]
                self.canonic_nodes = [junction2_nodes]
                self.event_jname = [junction1_name, "."]
                self.canonic_jname = [junction2_name]
        elif self.etype == "CE":
            self.event_cov = (junction2_coverage + junction3_coverage) // 2
            self.event_j = [
                parse_nodes(junction2_nodes),
                parse_nodes(junction3_nodes),
            ]
            self.canonic_cov = junction1_coverage
            self.canonic_j = parse_nodes(junction1_nodes)
            self.event_nodes = sorted(
                [junction2_nodes, junction3_nodes]
            )  # CHECKME: this should work since ids are topological sorted
            self.canonic_nodes = [junction1_nodes]
            self.event_jname = [junction2_name, junction3_name]
            self.canonic_jname = [junction1_name]
        self.add_replicate(replicate, self.event_cov, self.canonic_cov)

    def __repr__(self) -> str:
        return f"({self.etype}) E: {self.event_j} [{self.event_cov}] --- T {self.canonic_j} [{self.canonic_cov}]"

    def to_csv(self) -> str:
        # we need to select the junction names. We iterate to search for reference transcript. Otherwise, first haplotype-aware transcript
        # canonic_jname = []
        # for JNs in self.canonic_jname:
        #     jname = ""
        #     for jn in JNs:
        #         if jn.split(".")[0].endswith("R1"):
        #             jname = jn
        #             break
        #     if jname == "":
        #         jname = next(iter(JNs))
        #     canonic_jname.append(jname)

        # event_jname = []
        # for JNs in self.event_jname:
        #     jname = ""
        #     for jn in JNs:
        #         if jn.split(".")[0].endswith("R1"):
        #             jname = jn
        #             break
        #     if jname == "":
        #         jname = next(iter(JNs))
        #     event_jname.append(jname)

        return ",".join(
            map(
                str,
                [
                    self.etype,
                    self.annotation_type,
                    self.chrom,
                    self.gene,
                    self.strand,
                    ",".join(self.canonic_jname),
                    ",".join(self.event_jname),
                    ",".join(self.canonic_nodes),
                    ",".join(self.event_nodes),
                ],
            )
        )

    # def update(self, canonic_jname: list, event_jname: list):
    #     for i, _ in enumerate(self.canonic_jname):
    #         self.canonic_jname[i] |= canonic_jname[i]
    #     for i, _ in enumerate(self.event_jname):
    #         self.event_jname[i] |= event_jname[i]

    def add_replicate(self, replicate: int, event_cov: int, canonical_cov: int):
        if replicate >= len(self.replicates):
            # we need to add new replicate(s) - more if we have missing replicates
            while replicate + 1 != len(self.replicates):
                self.replicates.append([0, 0])
        # then just add current replicate
        self.replicates[replicate][0] = event_cov
        self.replicates[replicate][1] = canonical_cov

    def psi(self, mode: str = "canonic"):
        _psis = [calc_psi(*x, mode) for x in self.replicates]
        if all([x == -1 for x in _psis]):
            return -1
        return mean([x for x in _psis if x != -1])

    def get_event_cov(self):
        return int(mean([x[0] for x in self.replicates]))

    def get_canonic_cov(self):
        return int(mean([x[1] for x in self.replicates]))


def eq_event(e1: Event, e2: Event):
    if e1.etype != e2.etype or e1.chrom != e2.chrom or e1.gene != e2.gene:
        return False
    if e1.etype[0] == "A":
        return e1.canonic_j == e2.canonic_j and e1.event_j == e2.event_j
    elif e1.etype == "CE":
        # TODO
        print("CE not done yet", file=sys.stderr)
        sys.exit(1)
    elif e1.etype == "ES":
        return e1.canonic_j == e2.canonic_j and e1.event_j == e2.event_j
    elif e1.etype == "IR":
        # FIXME: canonic_j can follow a different path due to variations. We may do some sort of jaccard
        return e1.event_j == e2.event_j
    # FIXME: we may just do this for all event types
    # return e1.canonic_j == e2.canonic_j and e1.event_j == e2.event_j


def parse_pantas(fpath: str, i: int) -> dict[str, Event]:
    pantas = {x: [] for x in ETYPES}
    for line in open(fpath):
        line = line.strip()
        _e = line.split(",")
        if _e[0] == "event_type":
            continue
        e = Event(*_e, i)
        pantas[e.etype].append(e)
    return pantas


def main(args):
    events_1 = {x: [] for x in ETYPES}
    for i, fpath in enumerate(args.c1):
        _ei = parse_pantas(fpath, i)
        for etype in ETYPES:
            for _ee in _ei[etype]:
                if _ee.canonic_cov < args.w or _ee.event_cov < args.w:
                    continue
                eqs = [x for x in events_1[etype] if eq_event(_ee, x)]
                if len(eqs) > 0:
                    assert len(eqs) == 1
                    # we already found this event in this replicate or a previous one
                    # so we add the replicate (if new). This check is done in the add_replicate method thanks to i
                    eqs[0].add_replicate(i, _ee.event_cov, _ee.canonic_cov)
                    # we then update the junction names of the event
                    # eqs[0].update(_ee.canonic_jname, _ee.event_jname)
                else:
                    events_1[etype].append(_ee)

    events_2 = {x: [] for x in ETYPES}
    for i, fpath in enumerate(args.c2):
        _ei = parse_pantas(fpath, i)
        for etype in ETYPES:
            for _ee in _ei[etype]:
                if _ee.canonic_cov < args.w or _ee.event_cov < args.w:
                    continue
                eqs = [x for x in events_2[etype] if eq_event(_ee, x)]
                if len(eqs) > 0:
                    assert len(eqs) == 1
                    eqs[0].add_replicate(i, _ee.event_cov, _ee.canonic_cov)
                    # eqs[0].update(_ee.canonic_jname, _ee.event_jname)
                else:
                    events_2[etype].append(_ee)

    print(
        "etype",
        "annotation_type",
        "chrom",
        "gene",
        "strand",
        "junction1_name", # canonic
        "junction2_name", # canonic or event
        "junction3_name", # event if 2 canonic
        "junction1_nodes",
        "junction2_nodes",
        "junction3_nodes",
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
            eqs = [x for x in events_2[etype] if eq_event(e1, x)]
            if len(eqs) > 0:
                assert len(eqs) == 1
                psi1 = e1.psi()
                psi2 = eqs[0].psi()
                dpsi = max(0, psi1) - max(0, psi2)
                if psi1 == -1 and psi2 == -1:
                    dpsi = -1
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
                if args.both:
                    continue
                if not e1.psi() == -1:
                    psi2 = "NaN"
                    dpsi = "NaN"
                    w2 = "."
                    print(
                        e1.to_csv(),
                        f"{e1.get_canonic_cov()}/{e1.get_event_cov()}",
                        w2,
                        e1.psi(),
                        psi2,
                        dpsi,
                        sep=",",
                    )
        if not args.both:
            for e2 in events_2[etype]:
                eqs = [x for x in events_1[etype] if eq_event(e2, x)]
                if len(eqs) == 0 and not e2.psi() == -1:
                    psi1 = "NaN"
                    dpsi = "NaN"
                    w1 = "."
                    print(
                        e2.to_csv(),
                        w1,
                        f"{e2.get_canonic_cov()}/{e2.get_event_cov()}",
                        psi1,
                        e2.psi(),
                        dpsi,
                        sep=",",
                    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        prog="Quantify",
        description="",
    )
    parser.add_argument(
        "--c1",
        help="Replicates of condition 1",
        dest="c1",
        type=str,
        required=True,
        nargs="+",
    )
    parser.add_argument(
        "--c2",
        help="Replicates of condition 2",
        dest="c2",
        type=str,
        required=True,
        nargs="+",
    )

    parser.add_argument(
        "--both",
        dest="both",
        help="Report only events present in both conditions",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "-w",
        dest="w",
        help="Minimum value of read count of the event to be valid (Default: 1)",
        type=int,
        default=1,
    )
    args = parser.parse_args()

    if not len(args.c1) == len(args.c2):
        print("Provide the same number of replicates for each condition")
    else:
        main(args)
