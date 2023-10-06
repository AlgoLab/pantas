import itertools
import re
import sys
from functools import partial
import logging
from math import floor

eprint = logging.debug


def collapse_linkcounts(lc: list):
    count = sum([x[1] for x in lc])
    # weighted sum of positions
    r = [(x[0] * x[1]) / count for x in lc]
    pos = int(floor(sum(r)))
    # eprint(f"{lc=} {count=} {pos=}")
    return [pos, count]


def build_attrs(fields: str, d: int):
    attrs = dict()
    for f in fields:
        name, _, value = f.split(":")
        if name == "LN":
            attrs[name] = int(value)
        elif name == "RC" or name == "NC":
            attrs[name] = int(value)
        elif name == "IL" or name == "OL":
            _v = [x for x in value.split(",")]
            _v = [list(map(int, x.split("."))) for x in _v]
            eprint(f"{_v=}")
            if len(_v) >= 2:
                minv = min(_v, key=lambda x: x[0])
                k1 = [minv.copy()]
                maxv = max(_v, key=lambda x: x[0])
                k2 = [maxv.copy()]
                eprint(f"{minv=} {maxv=}")

                if abs(minv[0] - maxv[0]) < d:
                    # collapse all
                    _v = [collapse_linkcounts(_v)]
                else:
                    for x in _v:
                        if x == minv or x == maxv:
                            continue
                        dmin = abs(x[0] - minv[0])
                        dmax = abs(x[0] - maxv[0])

                        if dmin < dmax:
                            k1.append(x)
                        else:
                            k2.append(x)

                        # eprint(f"{x=}, {dmin=}, {dmax=}")
                    k1 = collapse_linkcounts(k1)
                    k2 = collapse_linkcounts(k2)
                    _v = [k1, k2]
                # eprint(f"{k1=} {k2=}")

            attrs[name] = _v
            attrs[f"MAX{name}"] = max(_v, key=lambda x: x[1])[0]
        else:
            attrs[name] = value.split(",")
    return attrs


def get_refpos(
    segments: dict,
    start: str,
    end: str,
):
    if "RP" in segments[start] and "RP" in segments[end]:
        add_start = segments[start].get("LN")
        add_end = 0
        # eprint(f"{add_start=}")
        # eprint(f"{add_end=}")
        return (
            f"{segments[start]['RP'] + add_start + 1}"
            + "-"
            + f"{segments[end]['RP'] + add_end}"
        )
    else:
        return "?-?"


def get_refpos_node(segments: dict, nid: str, key: str, jn_w: int = -1):
    if "RP" in segments[nid]:
        if key == "LN":
            add = segments[nid]["LN"] + 1
            # return f"{segments[nid]['RP'] + segments[nid]['LN'] + 1}"
        elif key == "OL":
            assert jn_w > 0
            add = (
                min(
                    [(x[0], abs(jn_w - x[1])) for x in segments[nid]["OL"]]
                    if "OL" in segments[nid]
                    else [(segments[nid]["LN"], 0)],
                    key=lambda y: y[1],
                )[0]
                + 1
            )
        elif key == "MAXOL":
            add = segments[nid].get("MAXOL", segments[nid]["LN"]) + 1
        elif key == "IL":
            assert jn_w > 0
            add = min(
                [(x[0], abs(jn_w - x[1])) for x in segments[nid]["IL"]]
                if "IL" in segments[nid]
                else [(0, 1)],
                key=lambda y: y[1],
            )[0]
        elif key == "MAXIL":
            add = segments[nid].get("MAXIL", 0)
        elif key == "RP":
            add = 0
        eprint(f"[get_refpos_node]: {nid}= {segments[nid]}")
        eprint(f"[get_refpos_node]: {add=}")
        return f"{segments[nid]['RP'] + add}"
    else:
        return "?"


def get_outgoing_nodes(
    links: dict, nid: str, segments: dict = None, rc: int = -1
) -> list:
    ret = [k[1] for k in links.keys() if k[0] == nid]
    if segments:
        ret = [x for x in ret if segments[x]["NC"] > rc]
    return ret


def get_incoming_nodes(
    links: dict, nid: str, segments: dict = None, rc: int = -1
) -> list:
    ret = [k[0] for k in links.keys() if k[1] == nid]
    if segments:
        ret = [x for x in ret if segments[x]["NC"] > rc]
    return ret


def get_outgoing_links(links: dict, nid: str) -> list:
    return [k for k in links.keys() if k[0] == nid]


def get_incoming_links(links: dict, nid: str) -> list:
    return [k for k in links.keys() if k[1] == nid]


def get_set_exons(nodes: dict, nid: str) -> set:
    return set(nodes[nid]["EX"]) if "EX" in nodes[nid] else set()


def get_transcript_from_exons(exons) -> map:
    return map(lambda x: ".".join(x.split(".")[:-1]), exons)


def get_path_transcript(path: dict, pid: str, start=None, end=None):
    _key = pid
    if not pid.endswith("_R1"):
        _key = f"{pid}_R1"
    p = path[_key]["path"]
    start = p.index(start) if start else 0
    end = p.index(end) if end else len(p)
    _s = min(start, end)
    _e = max(start, end) + 1
    return p[_s:_e]


def check_junction(ix_j: tuple, segments: dict, links: dict, window: int, rc: int):
    next_n0 = get_outgoing_nodes(links, ix_j[0], segments=segments, rc=rc)
    prev_n1 = get_incoming_nodes(links, ix_j[1], segments=segments, rc=rc)

    eprint(f"pre {next_n0=}")
    eprint(f"pre {prev_n1=}")

    _intron_next = set(next_n0) - set(ix_j)
    _intron_prev = set(prev_n1) - set(ix_j)
    i = 0
    eprint(f"{i=} {_intron_next=}")
    eprint(f"{i=} {_intron_prev=}")

    _subpath_n = []
    _subpath_p = []
    _subpath_count = 0
    if len(_intron_next) > 0 and len(_intron_prev) > 0:
        _max_n = max([(x, segments[x]["NC"]) for x in _intron_next], key=lambda x: x[1])
        _max_p = max([(x, segments[x]["NC"]) for x in _intron_prev], key=lambda x: x[1])

        _subpath_n.append(_max_n[0])
        _subpath_p.append(_max_p[0])
        _subpath_count += _max_n[1] + _max_p[1]

        eprint(f"{i=} {_subpath_n=}")
        eprint(f"{i=} {_subpath_p=}")
    else:
        return False

    while i < window:
        i += 1
        _intron_next = [
            get_outgoing_nodes(links, x, segments=segments, rc=rc) for x in _intron_next
        ]
        _intron_prev = [
            get_incoming_nodes(links, x, segments=segments, rc=rc) for x in _intron_prev
        ]
        # flatten lists
        _intron_next = [x for y in _intron_next for x in y]
        _intron_prev = [x for y in _intron_prev for x in y]
        _intron_next = set(_intron_next) - set(ix_j)
        _intron_prev = set(_intron_prev) - set(ix_j)

        eprint(f"{i=} {_intron_next=}")
        eprint(f"{i=} {_intron_prev=}")

        if len(_intron_next & _intron_prev) > 0:
            _max_n = max(
                [(x, segments[x]["NC"]) for x in _intron_next & _intron_prev],
                key=lambda x: x[1],
            )
            _subpath_n.append(_max_n[0])
            _subpath_count += _max_n[1]
            i = window
            break

        if len(_intron_next) > 0 and len(_intron_prev) > 0:
            _max_n = max(
                [(x, segments[x]["NC"]) for x in _intron_next], key=lambda x: x[1]
            )
            _max_p = max(
                [(x, segments[x]["NC"]) for x in _intron_prev], key=lambda x: x[1]
            )

            _subpath_n.append(_max_n[0])
            _subpath_p.append(_max_p[0])
            _subpath_count += _max_n[1] + _max_p[1]

            eprint(f"{i=} {_subpath_n=}")
            eprint(f"{i=} {_subpath_p=}")

            if len(set(_subpath_n) & set(_subpath_p)) > 1:
                i = window
                break

        else:
            break

    if i == window:
        return True
    return False


def main(args):
    gfaS = dict()
    gfaL = dict()
    gfaP = dict()
    junctions = list()
    noveljunctions = list()
    refpath = []
    for line in open(args.GFA, "r"):
        line = line.strip()
        if line.startswith("S"):
            _, nid, seq, *fields = line.split()
            gfaS[nid] = build_attrs(fields, args.d)
            # TODO: uncomment if needed
            # gfaS[nid]['seq'] = seq
        elif line.startswith("P"):
            _, pid, p, _ = line.split()
            if "+," in p[:-1]:
                gfaP[pid] = {"path": p[:-1].split("+,")}
                gfaP[pid]["reverse"] = False
            else:
                gfaP[pid] = {"path": p[:-1].split("-,")}
                gfaP[pid]["reverse"] = True
            if not "_R1" in pid:
                refpath = gfaP[pid]["path"]
        elif line.startswith("L"):
            (
                _,
                nid_from,
                _,
                nid_to,
                _,
                overlap,
                *fields,
            ) = line.split()
            gfaL[(nid_from, nid_to)] = build_attrs(fields, args.d)
            # TODO: uncomment if needed
            # gfaL[(nid_from, nid_to)]['overlap'] = overlap
            if "JN" in gfaL[(nid_from, nid_to)]:
                junctions.append((nid_from, nid_to))
            if "ID" in gfaL[(nid_from, nid_to)]:
                noveljunctions.append((nid_from, nid_to))

    curr = 0
    for n in refpath:
        gfaS[n]["RP"] = curr
        curr += gfaS[n]["LN"]

    transcript2gene = dict()
    genestrand = dict()
    genechr = dict()
    for line in open(args.GTF):
        if line.startswith("#"):
            continue
        line = line.strip("\n").split("\t")

        if line[2] in [
            "mRNA",
            "transcript",
            "miRNA",
            "ncRNA",
            "pre_miRNA",
            "snoRNA",
            "pseudogene",  # FIXME: do we want all these?
        ]:
            # FIXME: we may need to add more chars to the regex
            gidx = (
                re.search('gene_id "[A-Za-z0-9_]+";', line[-1]).group(0).split('"')[-2]
            )
            tidx = (
                re.search('transcript_id "[A-Za-z0-9_]+";', line[-1])
                .group(0)
                .split('"')[-2]
            )
            transcript2gene[tidx] = gidx
            genestrand[gidx] = line[6]
            genechr[gidx] = line[0]

    if args.header:
        print(
            "event_type",
            "annotated/novel",
            "chr",
            "gene",
            "strand",
            "junction1_name",
            "junction1_nodes",
            "junction1_refpos",
            "junction1_coverage",
            "junction2_name",
            "junction2_nodes",
            "junction2_refpos",
            "junction2_coverage",
            "junction3_name",
            "junction3_nodes",
            "junction3_refpos",
            "junction3_coverage",
            sep=",",
        )

    # Check all the junctions
    def check_nonnovel():
        for ix_j in junctions:
            junc = gfaL[ix_j]
            if junc["RC"] > args.rc:
                _trjunc = set(map(lambda x: ".".join(x.split(".")[:-2]), junc["JN"]))
                eprint(f"[Checking junction {ix_j}]: {junc}, {_trjunc}")

                exons_n0 = get_set_exons(gfaS, ix_j[0])
                exons_n1 = get_set_exons(gfaS, ix_j[1])
                transcripts_n0 = set(
                    map(lambda x: ".".join(x.split(".")[:-1]), exons_n0)
                )
                transcripts_n1 = set(
                    map(lambda x: ".".join(x.split(".")[:-1]), exons_n1)
                )

                if (
                    len(cap := ((set(transcripts_n0) & set(transcripts_n1)) - _trjunc))
                    > 0
                ):
                    # cap contains all the trascripts that:
                    # 1. visit n0 and n1
                    # 2. are not part of the junction (n0, n1)

                    eprint(f"{cap=}")

                    # Checking for non-novel ES
                    for _tr in cap:
                        _fex0 = list(filter(lambda x: x.startswith(_tr), exons_n0))
                        _fex1 = list(filter(lambda x: x.startswith(_tr), exons_n1))
                        assert len(_fex0) == len(_fex1) == 1

                        _tex0 = int(_fex0[0].split(".")[-1])
                        _tex1 = int(_fex1[0].split(".")[-1])

                        if abs(_tex0 - _tex1) > 1:
                            # Checking that the trascripts belong to the same gene
                            for _j, _te in itertools.product(junc["JN"], [_tr]):
                                _tj = ".".join(_j.split(".")[:-2])
                                if transcript2gene[_tj] == transcript2gene[_te]:
                                    # Find the junctions of this exon
                                    _exno_min = min(_tex0, _tex1)
                                    _exno_max = max(_tex0, _tex1)
                                    _es_j1_name = f"{_tr}.{_exno_min}.{_exno_min+1}"
                                    _es_j2_name = f"{_tr}.{_exno_max-1}.{_exno_max}"

                                    if genestrand[transcript2gene[_tj]] == "-":
                                        _es_j1_name = f"{_tr}.{_exno_max-1}.{_exno_max}"
                                        _es_j2_name = f"{_tr}.{_exno_min}.{_exno_min+1}"

                                    _n_j1 = [x for x in junctions if x[0] == ix_j[0]]
                                    _n_j2 = [x for x in junctions if x[1] == ix_j[1]]
                                    _es_j1 = [
                                        x for x in _n_j1 if _es_j1_name in gfaL[x]["JN"]
                                    ]
                                    _es_j2 = [
                                        x for x in _n_j2 if _es_j2_name in gfaL[x]["JN"]
                                    ]

                                    if len(_es_j1) == 1 and len(_es_j2) == 1:
                                        print(
                                            "ES",
                                            "annotated",
                                            genechr[transcript2gene[_tj]],
                                            transcript2gene[_tj],
                                            genestrand[transcript2gene[_tj]],
                                            _j,
                                            ">".join(ix_j),
                                            f"{genechr[transcript2gene[_tr]]}:{get_refpos(gfaS, *ix_j)}",
                                            junc["RC"],
                                            _es_j1_name,
                                            ">".join(_es_j1[0]),
                                            f"{genechr[transcript2gene[_tr]]}:{get_refpos(gfaS, *_es_j1[0])}",
                                            gfaL[_es_j1[0]]["RC"],
                                            _es_j2_name,
                                            ">".join(_es_j2[0]),
                                            f"{genechr[transcript2gene[_tr]]}:{get_refpos(gfaS, *_es_j2[0])}",
                                            gfaL[_es_j2[0]]["RC"],
                                            sep=",",
                                        )

                    # Checking for non-novel A5
                    # this is A5 on + / A3 on -
                    for n in get_outgoing_nodes(gfaL, ix_j[0]):
                        if (
                            len(
                                cap_a5_ex := (
                                    (get_set_exons(gfaS, n) & exons_n0) - exons_n1
                                )
                            )
                            > 0
                        ):
                            cap_a5 = set(
                                map(lambda x: ".".join(x.split(".")[:-1]), cap_a5_ex)
                            )
                            cap_a5 = cap_a5 & cap
                            if len(cap_a5) > 0:
                                # Checking that the trascripts belong to the same gene
                                for _j, _te in itertools.product(junc["JN"], cap_a5):
                                    _tj = ".".join(_j.split(".")[:-2])
                                    if transcript2gene[_tj] == transcript2gene[_te]:
                                        _a_j = [x for x in junctions if x[1] == ix_j[1]]
                                        _a_j = list(
                                            filter(
                                                lambda x: any(
                                                    [
                                                        y.startswith(_te)
                                                        for y in gfaL[x]["JN"]
                                                    ]
                                                ),
                                                _a_j,
                                            )
                                        )

                                        if len(_a_j) == 1:
                                            _a_j = _a_j[0]
                                            _a_j_name = [
                                                x
                                                for x in gfaL[_a_j]["JN"]
                                                if x.startswith(_te)
                                            ]
                                            assert len(_a_j_name) == 1
                                            print(
                                                "A5"
                                                if genestrand[transcript2gene[_tj]]
                                                == "+"
                                                else "A3",
                                                "annotated",
                                                genechr[transcript2gene[_tj]],
                                                transcript2gene[_tj],
                                                genestrand[transcript2gene[_tj]],
                                                _j,
                                                ">".join(ix_j),
                                                f"{genechr[transcript2gene[_tr]]}:{get_refpos(gfaS, *ix_j)}",
                                                junc["RC"],
                                                _a_j_name[0],
                                                ">".join(_a_j),
                                                f"{genechr[transcript2gene[_tr]]}:{get_refpos(gfaS, *_a_j)}",
                                                gfaL[_a_j]["RC"],
                                                ".",
                                                ".",
                                                ".",
                                                ".",
                                                sep=",",
                                            )

                    # Checking for non-novel A3
                    # this is A3 on + / A5 on -
                    for n in get_incoming_nodes(gfaL, ix_j[1]):
                        if (
                            len(
                                cap_a3_ex := (
                                    (get_set_exons(gfaS, n) & exons_n1) - exons_n0
                                )
                            )
                            > 0
                        ):
                            cap_a3 = set(
                                map(lambda x: ".".join(x.split(".")[:-1]), cap_a3_ex)
                            )
                            cap_a3 = cap & cap_a3
                            if len(cap_a3) > 0:
                                # Checking that the trascripts belong to the same gene
                                for _j, _te in itertools.product(junc["JN"], cap_a3):
                                    _tj = ".".join(_j.split(".")[:-2])
                                    if transcript2gene[_tj] == transcript2gene[_te]:
                                        _a_j = [x for x in junctions if x[0] == ix_j[0]]
                                        _a_j = list(
                                            filter(
                                                lambda x: any(
                                                    [
                                                        y.startswith(_te)
                                                        for y in gfaL[x]["JN"]
                                                    ]
                                                ),
                                                _a_j,
                                            )
                                        )
                                        if len(_a_j) == 1:
                                            _a_j = _a_j[0]
                                            _a_j_name = [
                                                x
                                                for x in gfaL[_a_j]["JN"]
                                                if x.startswith(_te)
                                            ]
                                            assert len(_a_j_name) == 1

                                            # NOTE: In this case to keep ordering consistent with the reference
                                            # the order of trascript is inverted.
                                            # The first one is the one the event is considered to
                                            # and the second is the one containing the event
                                            print(
                                                "A3"
                                                if genestrand[transcript2gene[_tj]]
                                                == "+"
                                                else "A5",
                                                "annotated",
                                                genechr[transcript2gene[_tj]],
                                                transcript2gene[_tj],
                                                genestrand[transcript2gene[_tj]],
                                                _a_j_name[0],
                                                ">".join(_a_j),
                                                f"{genechr[transcript2gene[_tr]]}:{get_refpos(gfaS, *_a_j)}",
                                                gfaL[_a_j]["RC"],
                                                _j,
                                                ">".join(ix_j),
                                                f"{genechr[transcript2gene[_tr]]}:{get_refpos(gfaS, *ix_j)}",
                                                junc["RC"],
                                                ".",
                                                ".",
                                                ".",
                                                ".",
                                                sep=",",
                                            )

                    # Checking for non-novel IR
                    next_n0 = get_outgoing_nodes(gfaL, ix_j[0])
                    ex_next_n0 = [get_set_exons(gfaS, x) for x in next_n0]
                    prev_n1 = get_incoming_nodes(gfaL, ix_j[1])
                    ex_prev_n1 = [get_set_exons(gfaS, x) for x in prev_n1]
                    cap_ir = exons_n0.intersection(exons_n1, *ex_next_n0, *ex_prev_n1)

                    if len(cap_ir) > 0:
                        # Checking that the trascripts belong to the same gene
                        for _j, _e in itertools.product(junc["JN"], cap_ir):
                            _tj = ".".join(_j.split(".")[:-2])
                            _te = ".".join(_e.split(".")[:-1])
                            if transcript2gene[_tj] == transcript2gene[_te]:
                                # get retained nodes to get their read count
                                _subpath = get_path_transcript(
                                    gfaP, _te, start=ix_j[0], end=ix_j[1]
                                )
                                assert len(_subpath) > 2
                                # excluding the junction nodes
                                _subpath = _subpath[1:-1]
                                _count_sum = 0
                                for _in in _subpath:
                                    _count_sum += gfaS[_in].get("NC", 0)

                                print(
                                    "IR",
                                    "annotated",
                                    genechr[transcript2gene[_tj]],
                                    transcript2gene[_tj],
                                    genestrand[transcript2gene[_tj]],
                                    _j,
                                    ">".join(ix_j),
                                    f"{genechr[transcript2gene[_tr]]}:{get_refpos(gfaS, *ix_j)}",
                                    junc["RC"],
                                    _e,
                                    ">".join(_subpath),
                                    ".",  # CHECKME: this should not be needed
                                    _count_sum // len(_subpath),
                                    ".",
                                    ".",
                                    ".",
                                    ".",
                                    sep=",",
                                )
            eprint("-" * 15)

    if not args.annotated:
        check_nonnovel()

    def check_novel():
        def from_single_novel_junctions():
            # Check all novel junctions
            for ix_j in noveljunctions:
                junc = gfaL[ix_j]
                if junc["RC"] > args.rc:
                    _trjunc = set()
                    eprint(f"[Checking junction {ix_j}]: {junc}, {_trjunc}")

                    exons_n0 = get_set_exons(gfaS, ix_j[0])
                    eprint(f"{exons_n0=}")
                    exons_n1 = get_set_exons(gfaS, ix_j[1])
                    eprint(f"{exons_n1=}")
                    transcripts_n0 = set(
                        map(lambda x: ".".join(x.split(".")[:-1]), exons_n0)
                    )
                    eprint(f"{transcripts_n0=}")
                    transcripts_n1 = set(
                        map(lambda x: ".".join(x.split(".")[:-1]), exons_n1)
                    )
                    eprint(f"{transcripts_n1=}")

                    if (
                        len(
                            cap := (
                                (set(transcripts_n0) & set(transcripts_n1)) - _trjunc
                            )
                        )
                        > 0
                    ):
                        # cap contains all the trascripts that:
                        # 1. visit n0 and n1

                        eprint(f"{cap=}")

                        # Checking for novel ES
                        for _tr in cap:
                            _fex0 = list(filter(lambda x: x.startswith(_tr), exons_n0))
                            _fex1 = list(filter(lambda x: x.startswith(_tr), exons_n1))
                            assert len(_fex0) == len(_fex1) == 1

                            _tex0 = int(_fex0[0].split(".")[-1])
                            _tex1 = int(_fex1[0].split(".")[-1])

                            if abs(_tex0 - _tex1) > 1:
                                _exno_min = min(_tex0, _tex1)
                                _exno_max = max(_tex0, _tex1)
                                _es_j1_name = f"{_tr}.{_exno_min}.{_exno_min+1}"
                                _es_j2_name = f"{_tr}.{_exno_max-1}.{_exno_max}"

                                if genestrand[transcript2gene[_tr]] == "-":
                                    _es_j1_name = f"{_tr}.{_exno_max-1}.{_exno_max}"
                                    _es_j2_name = f"{_tr}.{_exno_min}.{_exno_min+1}"

                                _n_j1 = [x for x in junctions if x[0] == ix_j[0]]
                                _n_j2 = [x for x in junctions if x[1] == ix_j[1]]

                                _es_j1 = [
                                    x for x in _n_j1 if _es_j1_name in gfaL[x]["JN"]
                                ]
                                _es_j2 = [
                                    x for x in _n_j2 if _es_j2_name in gfaL[x]["JN"]
                                ]

                                print(
                                    "ES",
                                    "novel",
                                    genechr[transcript2gene[_tr]],
                                    transcript2gene[_tr],
                                    genestrand[transcript2gene[_tr]],
                                    "?",  # _j,
                                    ">".join(ix_j),
                                    f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, ix_j[0], 'LN')}-{get_refpos_node(gfaS, ix_j[1], 'RP')}",
                                    junc["RC"],
                                    _es_j1_name,
                                    ">".join(_es_j1[0]),
                                    f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, _es_j1[0][0], 'LN')}-{get_refpos_node(gfaS, _es_j1[0][1], 'RP')}",
                                    gfaL[_es_j1[0]]["RC"],
                                    _es_j2_name,
                                    ">".join(_es_j2[0]),
                                    f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, _es_j2[0][0], 'LN')}-{get_refpos_node(gfaS, _es_j2[0][1], 'RP')}",
                                    gfaL[_es_j2[0]]["RC"],
                                    sep=",",
                                )

                        # Checking for novel A5+ before / A3- after
                        for _tr in cap:
                            _fex0 = list(filter(lambda x: x.startswith(_tr), exons_n0))
                            _fex1 = list(filter(lambda x: x.startswith(_tr), exons_n1))
                            assert len(_fex0) == len(_fex1) == 1

                            eprint(f"{_fex0=}")
                            eprint(f"{_fex1=}")

                            _tex0 = int(_fex0[0].split(".")[-1])
                            _tex1 = int(_fex1[0].split(".")[-1])
                            if abs(_tex0 - _tex1) == 1:
                                ex_next_n0 = set().union(
                                    *[
                                        get_set_exons(gfaS, x)
                                        for x in get_outgoing_nodes(gfaL, ix_j[0])
                                    ]
                                )
                                eprint(f"check A5b {ex_next_n0=}")
                                if _fex0[0] in ex_next_n0:
                                    _a_j = [x for x in junctions if x[1] == ix_j[1]]
                                    _a_j = list(
                                        filter(
                                            lambda x: any(
                                                [
                                                    y.startswith(_tr)
                                                    for y in gfaL[x]["JN"]
                                                ]
                                            ),
                                            _a_j,
                                        )
                                    )
                                    if len(_a_j) == 1:
                                        _a_j = _a_j[0]
                                        _a_j_name = [
                                            x
                                            for x in gfaL[_a_j]["JN"]
                                            if x.startswith(_tr)
                                        ]
                                        assert len(_a_j_name) == 1

                                        # A5b+ / A3a-
                                        eprint("A5b: A5b+ / A3a-")
                                        print(
                                            "A5"
                                            if genestrand[transcript2gene[_tr]] == "+"
                                            else "A3",
                                            "novel",
                                            genechr[transcript2gene[_tr]],
                                            transcript2gene[_tr],
                                            genestrand[transcript2gene[_tr]],
                                            "?",  # _j,
                                            ">".join(ix_j),
                                            f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, ix_j[0], 'OL', junc['RC'])}-{get_refpos_node(gfaS, ix_j[1], 'RP')}",
                                            junc["RC"],
                                            _a_j_name[0],
                                            ">".join(_a_j),
                                            f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, _a_j[0], 'LN')}-{get_refpos_node(gfaS, _a_j[1], 'RP')}",
                                            gfaL[_a_j]["RC"],
                                            ".",
                                            ".",
                                            ".",
                                            ".",
                                            sep=",",
                                        )

                                ex_prev_n1 = set().union(
                                    *[
                                        get_set_exons(gfaS, x)
                                        for x in get_incoming_nodes(gfaL, ix_j[1])
                                    ]
                                )
                                eprint(f"check A3a {ex_prev_n1=}")
                                if _fex1[0] in ex_prev_n1:
                                    _a_j = [x for x in junctions if x[0] == ix_j[0]]
                                    _a_j = list(
                                        filter(
                                            lambda x: any(
                                                [
                                                    y.startswith(_tr)
                                                    for y in gfaL[x]["JN"]
                                                ]
                                            ),
                                            _a_j,
                                        )
                                    )
                                    if len(_a_j) == 1:
                                        _a_j = _a_j[0]
                                        _a_j_name = [
                                            x
                                            for x in gfaL[_a_j]["JN"]
                                            if x.startswith(_tr)
                                        ]
                                        assert len(_a_j_name) == 1
                                        eprint(f"{_fex1=}")

                                        # NOTE: In this case to keep ordering consistent with the reference
                                        # the order of trascript is inverted.
                                        # The first one is the one the event is considered to
                                        # and the second is the one containing the event

                                        # A3a+ / A5b-
                                        eprint("A3a: A3a+ / A5b-")
                                        print(
                                            "A3"
                                            if genestrand[transcript2gene[_tr]] == "+"
                                            else "A5",
                                            "novel",
                                            genechr[transcript2gene[_tr]],
                                            transcript2gene[_tr],
                                            genestrand[transcript2gene[_tr]],
                                            _a_j_name[0],
                                            ">".join(_a_j),
                                            f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, _a_j[0], 'LN')}-{get_refpos_node(gfaS, _a_j[1], 'RP')}",
                                            gfaL[_a_j]["RC"],
                                            "?",  # _j,
                                            ">".join(ix_j),
                                            f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, ix_j[0], 'LN')}-{get_refpos_node(gfaS, ix_j[1], 'IL', junc['RC'])}",
                                            junc["RC"],
                                            ".",
                                            ".",
                                            ".",
                                            ".",
                                            sep=",",
                                        )

                        # Checking for novel IR reverse
                        next_n0 = get_outgoing_nodes(
                            gfaL, ix_j[0], segments=gfaS, rc=args.rc
                        )
                        ex_next_n0 = [get_set_exons(gfaS, x) for x in next_n0]
                        prev_n1 = get_incoming_nodes(
                            gfaL, ix_j[1], segments=gfaS, rc=args.rc
                        )
                        ex_prev_n1 = [get_set_exons(gfaS, x) for x in prev_n1]
                        cap_ir = set().intersection(
                            exons_n0, exons_n1, *ex_next_n0, *ex_prev_n1
                        )
                        eprint(f"EX {cap_ir=}")

                        if len(cap_ir) > 0:
                            for ex_ir in cap_ir:
                                _tr = ".".join(ex_ir.split(".")[:-1])
                                _subpath = get_path_transcript(
                                    gfaP, _tr, start=ix_j[0], end=ix_j[1]
                                )
                                assert len(_subpath) > 2
                                # excluding the junction nodes
                                _subpath = _subpath[1:-1]
                                _count_sum = 0
                                for _in in _subpath:
                                    _count_sum += gfaS[_in].get("NC", 0)
                                eprint("IRr")

                                # if genestrand[transcript2gene[_tr]] == "+":
                                _refpos = f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, ix_j[0], 'LN')}-{get_refpos_node(gfaS, ix_j[1], 'RP')}"
                                # else:
                                #     _refpos = f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, ix_j[1], 'LN')}-{get_refpos_node(gfaS, ix_j[0], 'RP')}"

                                print(
                                    "IR",
                                    "novel",
                                    genechr[transcript2gene[_tr]],
                                    transcript2gene[_tr],
                                    genestrand[transcript2gene[_tr]],
                                    "?",  # _j,
                                    ">".join(ix_j),
                                    _refpos,
                                    junc["RC"],
                                    ex_ir,
                                    ">".join(_subpath),
                                    ".",  # CHECKME: this should not be needed
                                    _count_sum // len(_subpath),
                                    ".",
                                    ".",
                                    ".",
                                    ".",
                                    sep=",",
                                )

                    # Check A3 - before
                    # if len(exons_n1) == 0:
                    if (
                        len(set(get_transcript_from_exons(exons_n1)) & transcripts_n0) == 0
                    ):
                        # n1 is an intron, check if there is a junction
                        # from n0 to somewhere else

                        nX_j = [x for x in junctions if x[0] == ix_j[0]]
                        nX = [x[1] for x in nX_j if gfaS[x[1]].get("NC", 0) > args.rc]
                        if len(nX) > 0:
                            eprint(f"{nX=}")
                            exons_nX = [get_set_exons(gfaS, x) for x in nX]
                            eprint(f"{exons_nX=}")
                            transcripts_nX = set(
                                *[list(get_transcript_from_exons(x)) for x in exons_nX]
                            )
                            eprint(f"{transcripts_nX=}")

                            for _tr in transcripts_nX & transcripts_n0:
                                _fex0 = list(
                                    filter(lambda x: x.startswith(_tr), exons_n0)
                                )
                                _fexX = list(
                                    filter(lambda x: x.startswith(_tr), *exons_nX)
                                )
                                assert len(_fex0) == len(_fex1) == 1

                                _fnX = [
                                    x
                                    for x in nX
                                    if any(
                                        map(
                                            lambda y: y.startswith(_tr),
                                            get_set_exons(gfaS, x),
                                        )
                                    )
                                ]
                                eprint(f"{_fnX=}")

                                _tex0 = int(_fex0[0].split(".")[-1])
                                _texX = int(_fexX[0].split(".")[-1])

                                if abs(_tex0 - _texX) == 1:
                                    _a_j = list(
                                        filter(
                                            lambda x: any(
                                                [
                                                    y.startswith(_tr)
                                                    for y in gfaL[x]["JN"]
                                                ]
                                            ),
                                            nX_j,
                                        )
                                    )
                                    if len(_a_j) == 1 and any(
                                        check_junction(
                                            [ix_j[1], _x], gfaS, gfaL, args.irw, args.rc
                                        )
                                        for _x in _fnX
                                    ):
                                        eprint(f"{_a_j=} {gfaL[_a_j[0]]=}")
                                        _a_j = _a_j[0]
                                        _a_j_name = [
                                            x
                                            for x in gfaL[_a_j]["JN"]
                                            if x.startswith(_tr)
                                        ]
                                        assert len(_a_j_name) == 1
                                        eprint("A3b: A3b+ / A5a-")
                                        print(
                                            "A3"
                                            if genestrand[transcript2gene[_tr]] == "+"
                                            else "A5",
                                            "novel",
                                            genechr[transcript2gene[_tr]],
                                            transcript2gene[_tr],
                                            genestrand[transcript2gene[_tr]],
                                            _a_j_name[0],
                                            ">".join(_a_j),
                                            f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, _a_j[0], 'LN')}-{get_refpos_node(gfaS, _a_j[1], 'RP')}",
                                            gfaL[_a_j]["RC"],
                                            "?",  # _j,
                                            ">".join(ix_j),
                                            f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, ix_j[0], 'LN')}-{get_refpos_node(gfaS, ix_j[1], 'IL', junc['RC'])}",
                                            junc["RC"],
                                            ".",
                                            ".",
                                            ".",
                                            ".",
                                            sep=",",
                                        )

                    # Chek A5 - after
                    # if len(exons_n0) == 0:
                    if (
                        len(set(get_transcript_from_exons(exons_n0)) & transcripts_n1) == 0
                    ):
                        # n0 is an intron, check if there is a junction
                        # to n1 from somewhere else

                        nX_j = [x for x in junctions if x[1] == ix_j[1]]
                        nX = [x[0] for x in nX_j if gfaS[x[0]].get("NC", 0) > args.rc]
                        if len(nX) > 0:
                            eprint(f"{nX=}")
                            exons_nX = [get_set_exons(gfaS, x) for x in nX]
                            eprint(f"{exons_nX=}")
                            transcripts_nX = set(
                                *[list(get_transcript_from_exons(x)) for x in exons_nX]
                            )
                            eprint(f"{transcripts_nX=}")

                            for _tr in transcripts_nX & transcripts_n1:
                                eprint(f"{_tr=}")
                                _fex0 = list(
                                    filter(lambda x: x.startswith(_tr), exons_n1)
                                )
                                _fexX = list(
                                    filter(lambda x: x.startswith(_tr), *exons_nX)
                                )
                                eprint(f"{_fex0=}")
                                eprint(f"{_fexX=}")
                                assert len(_fex0) == len(_fexX) == 1

                                _fnX = [
                                    x
                                    for x in nX
                                    if any(
                                        map(
                                            lambda y: y.startswith(_tr),
                                            get_set_exons(gfaS, x),
                                        )
                                    )
                                ]
                                eprint(f"{_fnX=}")

                                _tex0 = int(_fex0[0].split(".")[-1])
                                _texX = int(_fexX[0].split(".")[-1])

                                if abs(_tex0 - _texX) == 1:
                                    _a_j = list(
                                        filter(
                                            lambda x: any(
                                                [
                                                    y.startswith(_tr)
                                                    for y in gfaL[x]["JN"]
                                                ]
                                            ),
                                            nX_j,
                                        )
                                    )
                                    if len(_a_j) == 1 and any(
                                        check_junction(
                                            [_x, ix_j[0]], gfaS, gfaL, args.irw, args.rc
                                        )
                                        for _x in _fnX
                                    ):
                                        eprint(f"{_a_j=}  {gfaL[_a_j[0]]=}")
                                        _a_j = _a_j[0]
                                        _a_j_name = [
                                            x
                                            for x in gfaL[_a_j]["JN"]
                                            if x.startswith(_tr)
                                        ]
                                        assert len(_a_j_name) == 1
                                        eprint("A5a: A5a+ / A3b-")
                                        print(
                                            "A5"
                                            if genestrand[transcript2gene[_tr]] == "+"
                                            else "A3",
                                            "novel",
                                            genechr[transcript2gene[_tr]],
                                            transcript2gene[_tr],
                                            genestrand[transcript2gene[_tr]],
                                            "?",  # _j,
                                            ">".join(ix_j),
                                            f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, ix_j[0], 'OL', junc['RC'])}-{get_refpos_node(gfaS, ix_j[1], 'RP')}",
                                            junc["RC"],
                                            _a_j_name[0],
                                            ">".join(_a_j),
                                            f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, _a_j[0], 'LN')}-{get_refpos_node(gfaS, _a_j[1], 'RP')}",
                                            gfaL[_a_j]["RC"],
                                            ".",
                                            ".",
                                            ".",
                                            ".",
                                            sep=",",
                                        )

                    eprint("-" * 15)

        from_single_novel_junctions()

        # Check potential CE
        # Known junctions (n0 > n1) that have novel junctions (n0 > nX) and (nY > n1)
        def from_novel_inside_nonnovel():
            for ix_j in junctions:
                junc = gfaL[ix_j]
                if junc["RC"] > args.rc:
                    _trjunc = set(
                        map(lambda x: ".".join(x.split(".")[:-2]), junc["JN"])
                    )
                    # nX = [x[1] for x in noveljunctions if x[0] == ix_j[0]]
                    # nY = [x[0] for x in noveljunctions if x[1] == ix_j[1]]
                    nX = [
                        x
                        for x in noveljunctions
                        if x[0] == ix_j[0] and gfaS[x[0]].get("NC", 0) > args.rc
                    ]
                    nY = [
                        x
                        for x in noveljunctions
                        if x[1] == ix_j[1] and gfaS[x[1]].get("NC", 0) > args.rc
                    ]

                    if len(nX) > 0 and len(nY) > 0:
                        eprint(f"[Checking junction {ix_j}]: {junc}, {_trjunc}")
                        eprint(
                            f"n0>nX= "
                            + f"{[(x, gfaL[x]) for x in noveljunctions if x[0] == ix_j[0]]}"
                        )
                        eprint(
                            f"nY>n1= "
                            + f"{[(x, gfaL[x]) for x in noveljunctions if x[1] == ix_j[1]]}"
                        )
                        eprint(f"{nX=}")
                        eprint(f"{nY=}")
                        _enx = get_set_exons(gfaS, ix_j[0])
                        eprint(f"exons_nx: {_enx}")
                        _eny = get_set_exons(gfaS, ix_j[1])
                        eprint(f"exons_ny: {_eny}")

                        for _nx, _ny in itertools.product(nX, nY):
                            eprint(f"pair: {_nx} - {_ny}")

                            _tnx = set(get_transcript_from_exons(_enx))
                            eprint(f"TR_nx: {_tnx}")
                            _tny = set(get_transcript_from_exons(_eny))
                            eprint(f"TR_ny: {_tny}")

                            for _tr in _tnx & _tny:
                                _fex0 = list(filter(lambda x: x.startswith(_tr), _enx))
                                _fex1 = list(filter(lambda x: x.startswith(_tr), _eny))
                                assert len(_fex0) == len(_fex1) == 1

                                _tex0 = int(_fex0[0].split(".")[-1])
                                _tex1 = int(_fex1[0].split(".")[-1])

                                if abs(_tex0 - _tex1) == 1:
                                    # TODO: get sequence if necessary

                                    # _ce_trs = set(
                                    #     get_transcript_from_exons(
                                    #         get_set_exons(gfaS, _nx[1])
                                    #     )
                                    # ) | set(
                                    #     get_transcript_from_exons(
                                    #         get_set_exons(gfaS, _ny[0])
                                    #     )
                                    # )
                                    # for x in _ce_trs:
                                    #     try:
                                    #         _subpath = get_path_transcript(
                                    #             gfaP, x, start=_nx[1], end=_ny[0]
                                    #         )
                                    #         seq_ce = "".join([gfaS[x]["seq"] for x in _subpath])
                                    #         break
                                    #     except:
                                    #         continue
                                    # else:
                                    #     # this might be because:
                                    #     # 1. the nodes are intronic or
                                    #     # 2. the assumption of nX and nY for the CE is not valid:
                                    #     #   e.g. there is more than one jump, it might jump onto different genes/transcripts
                                    #     #       maybe we can check and differentiate the two cases and treat them differently
                                    #     seq_ce = "?"

                                    print(
                                        "CE",
                                        "novel",
                                        genechr[transcript2gene[_tr]],
                                        transcript2gene[_tr],
                                        genestrand[transcript2gene[_tr]],
                                        # Intron on annotation
                                        f"{_tr}.{min(_tex0, _tex1)}.{max(_tex0, _tex1)}",
                                        ">".join(ix_j),
                                        f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, ix_j[0], 'LN')}-{get_refpos_node(gfaS, ix_j[1], 'RP')}",
                                        junc["RC"],
                                        # Cassette junction 1
                                        "?",
                                        ">".join(_nx),
                                        # CHECKME: maybe this is not alway true, but yes
                                        f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, _nx[0], 'LN')}-{get_refpos_node(gfaS, _nx[1], 'IL', junc['RC'])}",
                                        gfaL[_nx]["RC"],
                                        # Cassette junction 2
                                        "?",
                                        ">".join(_ny),
                                        # CHECKME: maybe this is not alway true, but yes
                                        f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, _ny[0], 'OL', junc['RC'])}-{get_refpos_node(gfaS, _ny[1], 'RP')}",
                                        gfaL[_ny]["RC"],
                                        sep=",",
                                    )

                        eprint("-" * 15)

        from_novel_inside_nonnovel()

        # checking for IR
        for ix_j in junctions:
            junc = gfaL[ix_j]
            _trjunc = set(map(lambda x: ".".join(x.split(".")[:-2]), junc["JN"]))
            eprint(f"[IIR Checking junction {ix_j}]: {junc}, {_trjunc}")

            next_n0 = get_outgoing_nodes(gfaL, ix_j[0], segments=gfaS, rc=args.rc)
            prev_n1 = get_incoming_nodes(gfaL, ix_j[1], segments=gfaS, rc=args.rc)

            eprint(f"pre {next_n0=}")
            eprint(f"pre {prev_n1=}")

            # filter only intronic
            # _intron_next = set(
            #     filter(lambda x: len(
            #         set(get_transcript_from_exons(get_set_exons(gfaS, x))) - _trjunc
            #         ) == 0, next_n0)
            # )
            # _intron_prev = set(
            #     filter(lambda x: len(
            #         set(get_transcript_from_exons(get_set_exons(gfaS, x))) - _trjunc
            #         ) == 0, prev_n1)
            # )

            _intron_next = set(next_n0) - set(ix_j)
            _intron_prev = set(prev_n1) - set(ix_j)
            i = 0
            eprint(f"{i=} {_intron_next=}")
            eprint(f"{i=} {_intron_prev=}")

            _subpath_n = []
            _subpath_p = []
            _subpath_count = 0
            if len(_intron_next) > 0 and len(_intron_prev) > 0:
                _max_n = max(
                    [(x, gfaS[x]["NC"]) for x in _intron_next], key=lambda x: x[1]
                )
                _max_p = max(
                    [(x, gfaS[x]["NC"]) for x in _intron_prev], key=lambda x: x[1]
                )

                _subpath_n.append(_max_n[0])
                _subpath_p.append(_max_p[0])
                _subpath_count += _max_n[1] + _max_p[1]

                eprint(f"{i=} {_subpath_n=}")
                eprint(f"{i=} {_subpath_p=}")
            else:
                continue

            _subpath_total = False

            while i < args.irw:
                i += 1
                _intron_next = [
                    get_outgoing_nodes(gfaL, x, segments=gfaS, rc=args.rc)
                    for x in _intron_next
                ]
                _intron_prev = [
                    get_incoming_nodes(gfaL, x, segments=gfaS, rc=args.rc)
                    for x in _intron_prev
                ]
                # flatten lists
                _intron_next = [x for y in _intron_next for x in y]
                _intron_prev = [x for y in _intron_prev for x in y]
                _intron_next = set(_intron_next) - set(ix_j)
                _intron_prev = set(_intron_prev) - set(ix_j)

                eprint(f"{i=} {_intron_next=}")
                eprint(f"{i=} {_intron_prev=}")

                # _intron_next = set(
                #     filter(
                #         lambda x: len(
                #             set(get_transcript_from_exons(get_set_exons(gfaS, x))) - _trjunc
                #             ) == 0,
                #         _intron_next,
                #     )
                # )
                # _intron_prev = set(
                #     filter(
                #         lambda x: len(
                #             set(get_transcript_from_exons(get_set_exons(gfaS, x))) - _trjunc
                #                       ) == 0,
                #         _intron_prev,
                #     )
                # )
                eprint(f"{i=} {_intron_next=}")
                eprint(f"{i=} {_intron_prev=}")

                if len(_intron_next & _intron_prev) > 0:
                    _subpath_total = True
                    _max_n = max(
                        [(x, gfaS[x]["NC"]) for x in _intron_next & _intron_prev],
                        key=lambda x: x[1],
                    )
                    _subpath_n.append(_max_n[0])
                    _subpath_count += _max_n[1]
                    i = args.irw
                    break

                if len(_intron_next) > 0 and len(_intron_prev) > 0:
                    _max_n = max(
                        [(x, gfaS[x]["NC"]) for x in _intron_next], key=lambda x: x[1]
                    )
                    _max_p = max(
                        [(x, gfaS[x]["NC"]) for x in _intron_prev], key=lambda x: x[1]
                    )

                    _subpath_n.append(_max_n[0])
                    _subpath_p.append(_max_p[0])
                    _subpath_count += _max_n[1] + _max_p[1]

                    eprint(f"{i=} {_subpath_n=}")
                    eprint(f"{i=} {_subpath_p=}")

                    if len(set(_subpath_n) & set(_subpath_p)) > 1:
                        _subpath_total = True
                        i = args.irw
                        break

                else:
                    break

            if i == args.irw:
                _subpath = _subpath_n + _subpath_p[::-1]
                _subpath_name = ">".join(_subpath_n + ["?"] + _subpath_p[::-1])
                eprint(f"{i=} {_subpath_n=}")
                eprint(f"{i=} {_subpath_p=}")
                eprint(f"{i=} {_subpath=}")
                # Must be done before collapsing because collapsed nodes
                # are counted twice
                _subpath_avg = _subpath_count // len(_subpath)
                if _subpath_total:
                    _subpath = list(dict.fromkeys(_subpath))
                    _subpath_name = ">".join(_subpath)
                eprint(f"{i=} {_subpath=}")

                for _j in junc["JN"]:
                    _tr = ".".join(_j.split(".")[:-2])

                    # if genestrand[transcript2gene[_tr]] == "+":
                    _refpos = f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, ix_j[0], 'LN')}-{get_refpos_node(gfaS, ix_j[1], 'RP')}"
                    # else:
                    #     _refpos = f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, ix_j[1], 'LN')}-{get_refpos_node(gfaS, ix_j[0], 'RP')}"

                    eprint("IR")
                    print(
                        "IR",
                        "novel",
                        genechr[transcript2gene[_tr]],
                        transcript2gene[_tr],
                        genestrand[transcript2gene[_tr]],
                        _j,
                        ">".join(ix_j),
                        _refpos,
                        junc["RC"],
                        "?",  # ex_ir,
                        _subpath_name,
                        ".",
                        _subpath_avg,
                        ".",
                        ".",
                        ".",
                        ".",
                        sep=",",
                    )

    if args.novel:
        check_novel()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        prog="Caller",
        description="",
    )
    parser.add_argument("GFA", help="Spliced pangenome in GFA format")
    parser.add_argument("GTF", help="Annotation in GTF format")
    parser.add_argument(
        "--rc",
        help="Minimum read count (default: -1)",
        dest="rc",
        type=int,
        required=False,
        default=-1,
    )
    parser.add_argument(
        "--novel",
        dest="novel",
        help="Call novel events (default: False)",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--no-annotated",
        dest="annotated",
        help="Do not call known annotated events (default: False)",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--w",
        dest="irw",
        help="Intronic search window for novel events, larger values reduce FP but increase time (default: 5)",
        type=int,
        default=5,
    )
    parser.add_argument(
        "--d",
        dest="d",
        help="Maximum distance for OL/IL 2-clustering position collapsing (default: 3)",
        type=int,
        default=3,
    )
    parser.add_argument(
        "--debug",
        dest="debug",
        help="Debug (default: False)",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--header",
        dest="header",
        help="Print CSV header (default: False)",
        action="store_true",
        default=False,
    )
    args = parser.parse_args()
    if args.debug:
        # logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
        eprint = print
    main(args)
