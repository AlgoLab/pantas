import sys
import numpy as np
import eparser
from math import isnan

ETYPES = ["ES", "IR", "A3", "A5", "CE"]
EMAP_WHIPPET = {"CE": "ES", "RI": "IR", "AD": "A5", "AA": "A3"}


def precision_recall_f1(tp: int, fn: int, fp: int) -> list[int]:
    prec = round(float(tp) / (tp + fp) if tp + fp != 0 else 0, 3)
    rec = round(float(tp) / (tp + fn) if tp + fn != 0 else 0, 3)
    f1 = round(2 * float(tp) / (2 * tp + fp + fn) if tp + fp + fn != 0 else 0, 3)
    return [prec, rec, f1]


def is_good(e, DPSI_FILTER, MIN_EVENT_COV, novel=False):
    if abs(e.dpsi) < DPSI_FILTER:
        return False
    if not novel:
        if any([c < MIN_EVENT_COV for c in e.rc_c1 + e.rc_c2]):
            return False
    else:
        return e.min_event_cov >= MIN_EVENT_COV
    return True


def main(args):
    sep = "\t" if args.tabs else ","
    event_truth = {x: [] for x in ETYPES}
    for line in open(args.t, "r"):
        line = line.strip()
        (etype, chrom, gene, strand, j1, j2, j3, w1, w2, psi1, psi2) = line.split(",")
        if etype not in args.events:
            continue
        psi1 = float(psi1)
        psi2 = float(psi2)
        if isnan(psi1) or isnan(psi2):
            continue
        dpsi = max(0, psi1) - max(0, psi2)
        if psi1 == -1 and psi2 == -1:
            dpsi = -1
        e = eparser.EventTruth(
            etype, "truth", chrom, gene, strand, j1, j2, j3, w1, w2, psi1, psi2, dpsi
        )
        if args.mode2 and not is_good(e, args.min_dpsi, args.min_cov, args.novel):
            continue
        event_truth[e.etype].append(e)

    # sys.exit()

    event_pantas = {x: [] for x in ETYPES}
    for line in open(args.p, "r"):
        if line.startswith("etype"):
            continue
        line = line.strip()
        _e = line.split(",")
        if _e[2] == "haplotype":
            continue
        e = eparser.EventPantas(*_e)
        if isnan(e.psi_c1) or isnan(e.psi_c2):
            continue
        if abs(e.dpsi) < args.min_dpsi:
            continue
        event_pantas[e.etype].append(e)

    event_rmats = {x: [] for x in ETYPES}
    if args.r:
        for line in open(args.r, "r"):
            if line.startswith("etype"):
                continue
            line = line.strip()
            _e = line.split(",")
            e = eparser.EventRmats(*_e)
            if isnan(e.psi_c1) or isnan(e.psi_c2):
                continue
            if args.mode2 and abs(e.dpsi) < args.min_dpsi:
                continue
            event_rmats[e.etype].append(e)

    event_whippet = {x: [] for x in ETYPES}
    if args.w:
        for line in open(args.w, "r"):
            if line.startswith("Gene"):
                continue
            line = line.strip()
            _e = line.split("\t")
            _e[4] = EMAP_WHIPPET.get(_e[4], _e[4])
            if _e[4] not in ETYPES:
                continue
            e = eparser.EventWhippet(*_e, "anno")
            if isnan(e.psi_c1) or isnan(e.psi_c2):
                continue
            if args.mode2 and abs(e.dpsi) < args.min_dpsi:
                continue
            event_whippet[e.etype].append(e)

    event_suppa = {x: [] for x in ETYPES}
    if args.s:
        for line in open(args.s, "r"):
            if line.startswith("Gene"):
                # header
                continue
            line = line.strip()
            _e = line.split(",")
            e = eparser.EventRmats(*_e)
            if isnan(e.dpsi):
                continue
            if args.mode2 and abs(e.dpsi) < args.min_dpsi:
                continue
            event_suppa[e.etype].append(e)

    TP_PANTAS = {x: 0 for x in ETYPES}
    FN_PANTAS = {x: 0 for x in ETYPES}
    FP_PANTAS = {x: 0 for x in ETYPES}

    TP_RMATS = {x: 0 for x in ETYPES}
    FN_RMATS = {x: 0 for x in ETYPES}
    FP_RMATS = {x: 0 for x in ETYPES}

    TP_WHIPPET = {x: 0 for x in ETYPES}
    FN_WHIPPET = {x: 0 for x in ETYPES}
    FP_WHIPPET = {x: 0 for x in ETYPES}

    TP_SUPPA = {x: 0 for x in ETYPES}
    FN_SUPPA = {x: 0 for x in ETYPES}
    FP_SUPPA = {x: 0 for x in ETYPES}

    ### Loop to compute TPs and FNs
    for etype in ETYPES:
        if etype not in args.events:
            continue
        for e1 in event_truth[etype]:
            if not args.mode2 and not is_good(
                e1, args.min_dpsi, args.min_cov, args.novel
            ):
                continue
            str_event = e1.to_csv()
            eqsp = [
                x
                for x in event_pantas[etype]
                if eparser.eq_event(e1, x, args.novel, print_flag=args.print)
            ]
            if len(eqsp) > 0:
                # True positives
                TP_PANTAS[etype] += 1
                str_event += ",TP"
                # if args.print:
                #     print("TP", e1.to_csv())
            else:
                # False negatives
                FN_PANTAS[etype] += 1
                str_event += ",FN"
                if args.print:
                    print("FN", e1.to_csv(), file=sys.stderr)

            eqsr = [
                x for x in event_rmats[etype] if eparser.eq_event(e1, x, args.novel)
            ]
            if len(eqsr) > 0:
                # True positives
                assert len(eqsr) == 1
                TP_RMATS[etype] += 1
                str_event += ",TP"
            else:
                # False negatives
                FN_RMATS[etype] += 1
                str_event += ",FN"
                if args.print:
                    print("FN RMATS", e1.to_csv())

            eqsw = [
                x for x in event_whippet[etype] if eparser.eq_event(e1, x, args.novel)
            ]
            if len(eqsw) > 0:
                # True positives
                assert len(eqsw) == 1
                TP_WHIPPET[etype] += 1
                str_event += ",TP"
            else:
                # False negatives
                FN_WHIPPET[etype] += 1
                str_event += ",FN"
                if args.print:
                    print("FN WHIPPET", e1.to_csv())

            eqss = [
                x for x in event_suppa[etype] if eparser.eq_event(e1, x, args.novel)
            ]
            if len(eqss) > 0:
                # True positives
                assert len(eqss) == 1
                TP_SUPPA[etype] += 1
                str_event += ",TP"
            else:
                # False negatives
                FN_SUPPA[etype] += 1
                str_event += ",FN"
                if args.print:
                    print("FN SUPPA2", e1.to_csv())

            # print(str_event)

    ### Loop to compute FPs
    for etype in ETYPES:
        if etype not in args.events:
            continue
        for e2 in event_pantas[etype]:
            if not args.mode2 and abs(e2.dpsi) < args.min_dpsi:
                continue
            # print(len(event_truth[etype]))
            eqs = [
                x
                for x in event_truth[etype]
                if eparser.eq_event(x, e2, args.novel, print_flag=args.print)
            ]
            if len(eqs) == 0:
                # False positives
                FP_PANTAS[etype] += 1
                if args.print:
                    print("FP-PANTAS", e2.to_csv(), file=sys.stderr)

        for e2 in event_rmats[etype]:
            if not args.mode2 and abs(e2.dpsi) < args.min_dpsi:
                continue
            eqs = [x for x in event_truth[etype] if eparser.eq_event(x, e2, args.novel)]
            if len(eqs) == 0:
                # False positives
                FP_RMATS[etype] += 1
                if args.print:
                    print("FP-RMATS", e2.to_csv())

        for e2 in event_whippet[etype]:
            if not args.mode2 and abs(e2.dpsi) < args.min_dpsi:
                continue
            eqs = [x for x in event_truth[etype] if eparser.eq_event(x, e2, args.novel)]
            if len(eqs) == 0:
                # False positives
                FP_WHIPPET[etype] += 1
                if args.print:
                    print("FP-WHIPPET", e2.to_csv(), file=sys.stderr)

        for e2 in event_suppa[etype]:
            if not args.mode2 and abs(e2.dpsi) < args.min_dpsi:
                continue
            eqs = [x for x in event_truth[etype] if eparser.eq_event(x, e2, args.novel)]
            if len(eqs) == 0:
                # False positives
                FP_SUPPA[etype] += 1
                if args.print:
                    print("FP-SUPPA2", e2.to_csv(), file=sys.stderr)

    # print("PANTAS")
    print(
        "p-supp",
        "tool",
        "etype",
        "mindpsi",
        "mincov",
        "TP",
        "FN",
        "FP",
        "Prec",
        "Rec",
        "F1",
        sep=sep,
    )
    for etype in ETYPES:
        if etype not in args.events:
            continue
        print(
            args.supp,
            "pantas",
            etype,
            args.min_dpsi,
            args.min_cov,
            TP_PANTAS[etype],
            FN_PANTAS[etype],
            FP_PANTAS[etype],
            *precision_recall_f1(TP_PANTAS[etype], FN_PANTAS[etype], FP_PANTAS[etype]),
            sep=sep
        )

    # print("RMATS")
    # print("etype", "TP", "FN", "FP", "Prec", "Rec", "F1", sep=",")
    if args.r:
        for etype in ETYPES:
            if etype not in args.events:
                continue
            print(
                0,
                "rMATS",
                etype,
                args.min_dpsi,
                args.min_cov,
                TP_RMATS[etype],
                FN_RMATS[etype],
                FP_RMATS[etype],
                *precision_recall_f1(TP_RMATS[etype], FN_RMATS[etype], FP_RMATS[etype]),
                sep=sep
            )

    if args.w:
        # print("WHIPPET")
        # print("etype", "TP", "FN", "FP", "Prec", "Rec", "F1", sep=",")
        for etype in ETYPES:
            if etype not in args.events:
                continue
            print(
                0,
                "Whippet",
                etype,
                args.min_dpsi,
                args.min_cov,
                TP_WHIPPET[etype],
                FN_WHIPPET[etype],
                FP_WHIPPET[etype],
                *precision_recall_f1(
                    TP_WHIPPET[etype], FN_WHIPPET[etype], FP_WHIPPET[etype]
                ),
                sep=sep
            )
    if args.s:
        # print("SUPPA2")
        # print("etype", "TP", "FN", "FP", "Prec", "Rec", "F1", sep=",")
        for etype in ETYPES:
            if etype not in args.events:
                continue
            print(
                0,
                "SUPPA2",
                etype,
                args.min_dpsi,
                args.min_cov,
                TP_SUPPA[etype],
                FN_SUPPA[etype],
                FP_SUPPA[etype],
                *precision_recall_f1(TP_SUPPA[etype], FN_SUPPA[etype], FP_SUPPA[etype]),
                sep=sep
            )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        prog="Compare",
        description="",
    )
    parser.add_argument(
        "--tabs",
        help="",
        dest="tabs",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--supp",
        help="Support parameters for pantas",
        dest="supp",
        type=int,
        default=0,
    )
    parser.add_argument(
        "-c",
        help="Minimum coverage for junctions",
        dest="min_cov",
        type=int,
        default=5,
    )
    parser.add_argument(
        "-d",
        help="Minimum dPSI for events",
        dest="min_dpsi",
        type=float,
        default=0.05,
    )
    parser.add_argument(
        "-t",
        help="Truth CSV",
        dest="t",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-p",
        help="Pantas CSV",
        dest="p",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-r",
        help="RMATS CSV",
        dest="r",
        type=str,
        required=False,
    )
    parser.add_argument(
        "-w",
        help="Whippet psi file",
        dest="w",
        type=str,
        required=False,
    )
    parser.add_argument(
        "-s",
        help="suppa2 CSV",
        dest="s",
        type=str,
        required=False,
    )
    parser.add_argument(
        "--novel",
        dest="novel",
        help="Perform comparison of novel events (Default: False)",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--print",
        dest="print",
        help="Print explicitly TP/FP/FN (Default: False)",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--events",
        help="Events to call (default: [ES, A3, A5, IR])",
        dest="events",
        nargs="+",
        required=False,
        default=["ES", "A3", "A5", "IR"],
    )
    parser.add_argument(
        "-2",
        dest="mode2",
        help="(Default: False)",
        action="store_true",
        default=False,
    )
    args = parser.parse_args()
    main(args)
