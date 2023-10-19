import sys
import numpy as np
import eparser
from math import isnan

ETYPES = ["ES", "IR", "A3", "A5", "CE"]
EMAP_WHIPPET = {"CE": "ES", "RI": "IR", "AD": "A5", "AA": "A3"}


def main(args):
    FILTER = 0.05
    MIN_EVENT_COV = 10
    event_truth = {x: [] for x in ETYPES}
    for line in open(args.t, "r"):
        line = line.strip()
        (etype, chrom, gene, strand, j1, j2, j3, w1, w2, psi1, psi2) = line.split(",")
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
        # if abs(e.dpsi) < FILTER:
        #     continue
        event_truth[e.etype].append(e)
        # print(e.to_csv(), e.rc_c1, e.rc_c2, e.event_cov_c1, e.event_cov_c2)

    # sys.exit()

    event_pantas = {x: [] for x in ETYPES}
    for line in open(args.p, "r"):
        if line.startswith("etype"):
            # header
            continue
        line = line.strip()
        _e = line.split(",")
        e = eparser.EventPantas(*_e)
        if isnan(e.psi_c1) or isnan(e.psi_c2):
            continue
        # if abs(e.dpsi) < FILTER:
        #     continue
        event_pantas[e.etype].append(e)

    event_rmats = {x: [] for x in ETYPES}
    for line in open(args.r, "r"):
        if line.startswith("etype"):
            # header
            continue
        line = line.strip()
        _e = line.split(",")
        e = eparser.EventPantas(*_e)
        if isnan(e.psi_c1) or isnan(e.psi_c2):
            continue
        # if abs(e.dpsi) < FILTER:
        #     continue
        event_rmats[e.etype].append(e)

    event_whippet = {x: [] for x in ETYPES}
    for line in open(args.w, "r"):
        if line.startswith("Gene"):
            # header
            continue
        line = line.strip()
        _e = line.split("\t")
        _e[4] = EMAP_WHIPPET.get(_e[4], _e[4])
        if _e[4] not in ETYPES:
            continue
        e = eparser.EventWhippet(*_e, "anno")
        if isnan(e.psi_c1) or isnan(e.psi_c2):
            continue
        if abs(e.dpsi) < FILTER:
            continue
        event_whippet[e.etype].append(e)

    TP_PANTAS = {x: 0 for x in ETYPES}
    FN_PANTAS = {x: 0 for x in ETYPES}
    FP_PANTAS = {x: 0 for x in ETYPES}

    TP_RMATS = {x: 0 for x in ETYPES}
    FN_RMATS = {x: 0 for x in ETYPES}
    FP_RMATS = {x: 0 for x in ETYPES}

    TP_WHIPPET = {x: 0 for x in ETYPES}
    FN_WHIPPET = {x: 0 for x in ETYPES}
    FP_WHIPPET = {x: 0 for x in ETYPES}

    for etype in ETYPES:
        for e1 in event_truth[etype]:
            # print(e1)
            if e1.min_event_cov < MIN_EVENT_COV:
                continue
            if abs(e1.dpsi) < FILTER:
                continue
            str_event = e1.to_csv()

            eqsp = [
                x
                for x in event_pantas[etype]
                if eparser.eq_event(e1, x, relax=args.relax)
            ]
            if len(eqsp) > 0:
                # True positives
                assert len(eqsp) == 1
                TP_PANTAS[etype] += 1
                str_event += ",TP"
            else:
                # False negatives
                FN_PANTAS[etype] += 1
                str_event += ",FN"
                # print("FN", e1.to_csv())

            eqsr = [
                x
                for x in event_rmats[etype]
                if eparser.eq_event(e1, x, relax=args.relax)
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
                # print("FN", e1.to_csv())

            eqsw = [
                x
                for x in event_whippet[etype]
                if eparser.eq_event(e1, x, relax=args.relax)
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

            # print(str_event)

    for etype in ETYPES:
        for e2 in event_pantas[etype]:
            if abs(e2.dpsi) < FILTER:
                continue
            # print(len(event_truth[etype]))
            eqs = [
                x
                for x in event_truth[etype]
                if eparser.eq_event(e2, x, relax=args.relax)
            ]
            if len(eqs) == 0:
                # False positives
                FP_PANTAS[etype] += 1
                # print("FP-PANTAS", e2.to_csv())
        for e2 in event_rmats[etype]:
            if abs(e2.dpsi) < FILTER:
                continue
            eqs = [
                x
                for x in event_truth[etype]
                if eparser.eq_event(e2, x, relax=args.relax)
            ]
            if len(eqs) == 0:
                # False positives
                FP_RMATS[etype] += 1
                # print("FP-RMATS", e2.to_csv())

        for e2 in event_whippet[etype]:
            if abs(e2.dpsi) < FILTER:
                continue
            eqs = [
                x
                for x in event_truth[etype]
                if eparser.eq_event(e2, x, relax=args.relax)
            ]
            if len(eqs) == 0:
                # False positives
                FP_WHIPPET[etype] += 1
                # print("FP-WHIPPET", e2.to_csv())

    print("etype", "TP", "FN", "FP", sep=",")
    for etype in ETYPES:
        print(etype, TP_PANTAS[etype], FN_PANTAS[etype], FP_PANTAS[etype], sep=",")

    for etype in ETYPES:
        print(etype, TP_RMATS[etype], FN_RMATS[etype], FP_RMATS[etype], sep=",")

    for etype in ETYPES:
        print(etype, TP_WHIPPET[etype], FN_WHIPPET[etype], FP_WHIPPET[etype], sep=",")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        prog="Compare",
        description="",
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
        required=True,
    )
    parser.add_argument(
        "-w",
        help="Whippet psi file",
        dest="w",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--relax",
        dest="relax",
        help="Relaxation of reference positions matching (Default: 0)",
        type=int,
        default=0,
    )
    args = parser.parse_args()
    main(args)
