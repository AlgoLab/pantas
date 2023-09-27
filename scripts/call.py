import itertools
import re
import sys
from functools import partial

eprint = partial(print, file=sys.stderr)
# eprint = print


def build_attrs(fields: str):
    attrs = dict()
    for f in fields:
        name, _, value = f.split(":")
        if name == "LN":
            attrs[name] = int(value)
        elif name == "RC":
            attrs[name] = int(value)
        else:
            attrs[name] = value.split(",")
    return attrs


def get_outgoing_nodes(links: dict, nid: str) -> list:
    return [k[1] for k in links.keys() if k[0] == nid]


def get_incoming_nodes(links: dict, nid: str) -> list:
    return [k[0] for k in links.keys() if k[1] == nid]


def get_outgoing_links(links: dict, nid: str) -> list:
    return [k for k in links.keys() if k[0] == nid]


def get_incoming_links(links: dict, nid: str) -> list:
    return [k for k in links.keys() if k[1] == nid]


def get_set_exons(nodes: dict, nid: str) -> set:
    return set(nodes[nid]["EX"]) if "EX" in nodes[nid] else set()


def get_transcript_from_exons(exons) -> map:
    return map(lambda x: ".".join(x.split(".")[:-1]), exons)


# def traverse_se(links: dict, start: str, end: str):
#     out_curr = get_outgoing_nodes(links, start)
#     print(f"[traverse s:{start}-{end}] curr: {start} out: {out_curr}")
#     if end == start:
#         return
#     for n in out_curr:
#         traverse_se(links, n, end)


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


def main(args):
    gfaS = dict()
    gfaL = dict()
    gfaP = dict()
    junctions = list()
    noveljunctions = list()
    for line in open(args.GFA, "r"):
        line = line.strip()
        if line.startswith("S"):
            _, nid, seq, *fields = line.split()
            gfaS[nid] = build_attrs(fields)
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
            gfaL[(nid_from, nid_to)] = build_attrs(fields)
            # TODO: uncomment if needed
            # gfaL[(nid_from, nid_to)]['overlap'] = overlap
            if "JN" in gfaL[(nid_from, nid_to)]:
                junctions.append((nid_from, nid_to))
            if "ID" in gfaL[(nid_from, nid_to)]:
                noveljunctions.append((nid_from, nid_to))

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

    # Check all the junctions
    def check_nonnovel():
        for ix_j in junctions:
            junc = gfaL[ix_j]
            if junc["RC"] > args.rc:
                _trjunc = set(map(lambda x: ".".join(x.split(".")[:-2]), junc["JN"]))
                eprint(f"[Checking junction {ix_j}]: {junc}, {_trjunc}")

                # exons_n0 = set(gfaS[ix_j[0]]["EX"])
                # exons_n1 = set(gfaS[ix_j[1]]["EX"])
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

                    eprint("cap:", cap)

                    # Checking for non-novel ES https://hackmd.io/DoQzt8ceThOwyIdUvQZN3w#Exon-skipping
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
                                        if args.format == "junctions":
                                            print(
                                                "ES",
                                                "annotated",
                                                genechr[transcript2gene[_tj]],
                                                transcript2gene[_tj],
                                                genestrand[transcript2gene[_tj]],
                                                _j,
                                                junc["RC"],
                                                _es_j1_name,
                                                gfaL[_es_j1[0]]["RC"],
                                                _es_j2_name,
                                                gfaL[_es_j2[0]]["RC"],
                                                sep=",",
                                            )
                                        elif args.format == "nodes":
                                            print(
                                                "ES",
                                                "annotated",
                                                genechr[transcript2gene[_tj]],
                                                transcript2gene[_tj],
                                                genestrand[transcript2gene[_tj]],
                                                ">".join(ix_j),
                                                junc["RC"],
                                                ">".join(_es_j1[0]),
                                                gfaL[_es_j1[0]]["RC"],
                                                ">".join(_es_j2[0]),
                                                gfaL[_es_j2[0]]["RC"],
                                                sep=",",
                                            )
                                    # TODO: CHECKME: continue?

                    # Checking for non-novel A5 https://hackmd.io/DoQzt8ceThOwyIdUvQZN3w#Alternative-5%E2%80%99
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
                                            if args.format == "junctions":
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
                                                    junc["RC"],
                                                    _a_j_name[0],
                                                    gfaL[_a_j]["RC"],
                                                    ".",
                                                    ".",
                                                    sep=",",
                                                )
                                            elif args.format == "nodes":
                                                print(
                                                    "A5"
                                                    if genestrand[transcript2gene[_tj]]
                                                    == "+"
                                                    else "A3",
                                                    "annotated",
                                                    genechr[transcript2gene[_tj]],
                                                    transcript2gene[_tj],
                                                    genestrand[transcript2gene[_tj]],
                                                    ">".join(ix_j),
                                                    junc["RC"],
                                                    ">".join(_a_j),
                                                    gfaL[_a_j]["RC"],
                                                    ".",
                                                    ".",
                                                    sep=",",
                                                )
                                            # TODO: CHECKME: continue?

                    # Checking for non-novel A3 https://hackmd.io/DoQzt8ceThOwyIdUvQZN3w#Alternative-3%E2%80%99
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
                                        _a_j = [
                                            gfaL[x]
                                            for x in junctions
                                            if x[0] == ix_j[0]
                                        ]
                                        _a_j = list(
                                            filter(
                                                lambda x: any(
                                                    [y.startswith(_te) for y in x["JN"]]
                                                ),
                                                _a_j,
                                            )
                                        )
                                        if len(_a_j) == 1:
                                            _a_j = _a_j[0]
                                            _a_j_name = [
                                                x
                                                for x in _a_j["JN"]
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
                                                _a_j["RC"],
                                                _j,
                                                junc["RC"],
                                                ".",
                                                ".",
                                                sep=",",
                                            )
                                            # TODO: CHECKME: continue?

                    # Checking for non-novel IR https://hackmd.io/DoQzt8ceThOwyIdUvQZN3w#Intron-retention
                    next_n0 = get_outgoing_nodes(gfaL, ix_j[0])
                    ex_next_n0 = [get_set_exons(gfaS, x) for x in next_n0]
                    prev_n1 = get_incoming_nodes(gfaL, ix_j[1])
                    ex_prev_n1 = [get_set_exons(gfaS, x) for x in prev_n1]
                    cap_ir = exons_n0.intersection(exons_n1, *ex_next_n0, *ex_prev_n1)

                    if len(cap_ir) > 0:
                        # cap_ir = set(map(lambda x: ".".join(x.split(".")[:-1]), cap_ir))
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
                                    # TODO: change this to actual RC once we have it
                                    _count_sum += sum(
                                        [
                                            int(x.split(".")[1])
                                            for x in gfaS[_in].get("IL", ["0.0"])
                                        ]
                                    )

                                print(
                                    "IR",
                                    "annotated",
                                    genechr[transcript2gene[_tj]],
                                    transcript2gene[_tj],
                                    genestrand[transcript2gene[_tj]],
                                    _j,
                                    junc["RC"],
                                    _e,
                                    _count_sum // len(_subpath),
                                    ".",
                                    ".",
                                    sep=",",
                                )
                                # TODO: CHECKME: continue?
            eprint("-" * 15)

    if not args.annotated:
        check_nonnovel()

    def check_novel():
        def from_single_novel_junctions():
            # Check all novel junctions
            for ix_j in noveljunctions:
                junc = gfaL[ix_j]
                if junc["RC"] > args.rc:
                    # _trjunc = set(map(lambda x: ".".join(x.split(".")[:-2]), junc["JN"]))
                    _trjunc = set()
                    eprint(f"[Checking junction {ix_j}]: {junc}, {_trjunc}")

                    exons_n0 = get_set_exons(gfaS, ix_j[0])
                    eprint("exons_n0:", exons_n0)
                    exons_n1 = get_set_exons(gfaS, ix_j[1])
                    eprint("exons_n1:", exons_n1)
                    transcripts_n0 = set(
                        map(lambda x: ".".join(x.split(".")[:-1]), exons_n0)
                    )
                    eprint("transcripts_n0:", transcripts_n0)
                    transcripts_n1 = set(
                        map(lambda x: ".".join(x.split(".")[:-1]), exons_n1)
                    )
                    eprint("transcripts_n1:", transcripts_n1)

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
                        # 2. are not part of the junction (n0, n1)

                        eprint("cap:", cap)

                        # Checking for novel ES
                        for _tr in cap:
                            _fex0 = list(filter(lambda x: x.startswith(_tr), exons_n0))
                            _fex1 = list(filter(lambda x: x.startswith(_tr), exons_n1))
                            assert len(_fex0) == len(_fex1) == 1

                            _tex0 = int(_fex0[0].split(".")[-1])
                            _tex1 = int(_fex1[0].split(".")[-1])

                            if abs(_tex0 - _tex1) > 1:
                                print("****************** NOVEL ES!!!", _tr)
                                # Checking that the trascripts belong to the same gene
                                # for _tj, _te in itertools.product(_trjunc, [_tr]):
                                #     if transcript2gene[_tj] == transcript2gene[_te]:
                                #         print("****************** NOVEL ES!!!", _tj, _te)
                                #         # TODO: CHECKME: continue?

                        # Checking for novel A5+ before / A3- after
                        # Checking for novel A3+ after
                        for _tr in cap:
                            _fex0 = list(filter(lambda x: x.startswith(_tr), exons_n0))
                            _fex1 = list(filter(lambda x: x.startswith(_tr), exons_n1))
                            assert len(_fex0) == len(_fex1) == 1

                            _tex0 = int(_fex0[0].split(".")[-1])
                            _tex1 = int(_fex1[0].split(".")[-1])
                            if abs(_tex0 - _tex1) == 1:
                                # TODO: maybe this check is in not necessary
                                ex_next_n0 = set().union(
                                    *[
                                        get_set_exons(gfaS, x)
                                        for x in get_outgoing_nodes(gfaL, ix_j[0])
                                    ]
                                )
                                if _fex0[0] in ex_next_n0:
                                    print(
                                        "****************** NOVEL A5b!!!",
                                        _tr,
                                        min(_tex0, _tex1),
                                        ">",
                                        max(_tex0, _tex1),
                                    )
                                    print(
                                        "A5b+"
                                        if genestrand[transcript2gene[_tr]] == "+"
                                        else "A3a-"
                                    )

                                ex_prev_n1 = set().union(
                                    *[
                                        get_set_exons(gfaS, x)
                                        for x in get_incoming_nodes(gfaL, ix_j[1])
                                    ]
                                )
                                if _fex1[0] in ex_prev_n1:
                                    print(
                                        "****************** NOVEL A3a!!!",
                                        _tr,
                                        min(_tex0, _tex1),
                                        ">",
                                        max(_tex0, _tex1),
                                    )
                                    # print("A3a+" if genestrand[transcript2gene[_tr]] == '+' else "A5b-")

                        # Checking for novel IR reverse
                        next_n0 = get_outgoing_nodes(gfaL, ix_j[0])
                        ex_next_n0 = [get_set_exons(gfaS, x) for x in next_n0]
                        prev_n1 = get_incoming_nodes(gfaL, ix_j[1])
                        ex_prev_n1 = [get_set_exons(gfaS, x) for x in prev_n1]
                        cap_ir = exons_n0.intersection(
                            exons_n1, *ex_next_n0, *ex_prev_n1
                        )
                        eprint("cap_ir EX:", cap_ir)

                        if len(cap_ir) > 0:
                            cap_ir = set(
                                map(lambda x: ".".join(x.split(".")[:-1]), cap_ir)
                            )
                            print("****************** NOVEL IRr!!!", cap_ir)
                            # Checking that the trascripts belong to the same gene
                            # for _tj, _te in itertools.product(_trjunc, cap_ir):
                            #     if transcript2gene[_tj] == transcript2gene[_te]:
                            #         print("****************** NOVEL IRr!!!", _tj, _te)
                            #         # TODO: CHECKME: continue?

                    # Check A3 - before
                    if len(exons_n1) == 0:
                        # n1 is an intron, check if there is a junction
                        # from n0 to somewhere else

                        nX = [x[1] for x in junctions if x[0] == ix_j[0]]
                        if len(nX) > 0:
                            # eprint("========================== CHECK A3b")
                            eprint("nX:", nX)
                            exons_nX = [get_set_exons(gfaS, x) for x in nX]
                            eprint("exons_nx:", exons_nX)
                            transcripts_nX = set(
                                *[list(get_transcript_from_exons(x)) for x in exons_nX]
                            )
                            eprint("TR_nx:", transcripts_nX)

                            for _tr in transcripts_nX & transcripts_n0:
                                # eprint("_tr", _tr)
                                _fex0 = list(
                                    filter(lambda x: x.startswith(_tr), exons_n0)
                                )
                                _fexX = list(
                                    filter(lambda x: x.startswith(_tr), *exons_nX)
                                )
                                assert len(_fex0) == len(_fex1) == 1

                                _tex0 = int(_fex0[0].split(".")[-1])
                                _texX = int(_fexX[0].split(".")[-1])

                                if abs(_tex0 - _texX) == 1:
                                    print(
                                        "****************** NOVEL A3b!!!",
                                        _tr,
                                        min(_tex0, _texX),
                                        ">",
                                        max(_tex0, _texX),
                                    )

                    # Chek A5 - after
                    if len(exons_n0) == 0:
                        # n0 is an intron, check if there is a junction
                        # to n1 from somewhere else

                        nX = [x[0] for x in junctions if x[1] == ix_j[1]]
                        if len(nX) > 0:
                            eprint("nX:", nX)
                            exons_nX = [get_set_exons(gfaS, x) for x in nX]
                            eprint("exons_nx:", exons_nX)
                            transcripts_nX = set(
                                *[list(get_transcript_from_exons(x)) for x in exons_nX]
                            )
                            eprint("TR_nx:", transcripts_nX)

                            for _tr in transcripts_nX & transcripts_n0:
                                eprint("_tr", _tr)
                                _fex0 = list(
                                    filter(lambda x: x.startswith(_tr), exons_n0)
                                )
                                _fexX = list(
                                    filter(lambda x: x.startswith(_tr), *exons_nX)
                                )
                                eprint("fex0", _fex0)
                                eprint("fex1", _fex1)
                                assert len(_fex0) == len(_fex1) == 1

                                # TODO: this has not been tested. Is not present in the example I am using
                                _tex0 = int(_fex0[0].split(".")[-1])
                                _texX = int(_fexX[0].split(".")[-1])

                                if abs(_tex0 - _texX) == 1:
                                    print(
                                        "****************** NOVEL A5a!!!",
                                        _tr,
                                        min(_tex0, _texX),
                                        ">",
                                        max(_tex0, _texX),
                                    )

                    eprint("-" * 15)

        from_single_novel_junctions()

        # TODO: check IR

        # Check potential CE
        # Known junctions (n0 > n1) that have novel junctions (n0 > nX) and (nY > n1)
        def from_novel_inside_nonnovel():
            for ix_j in junctions:
                junc = gfaL[ix_j]
                if junc["RC"] > args.rc:
                    _trjunc = set(
                        map(lambda x: ".".join(x.split(".")[:-2]), junc["JN"])
                    )
                    nX = [x[1] for x in noveljunctions if x[0] == ix_j[0]]
                    nY = [x[0] for x in noveljunctions if x[1] == ix_j[1]]

                    if len(nX) > 0 and len(nY) > 0:
                        eprint(f"[Checking junction {ix_j}]: {junc}, {_trjunc}")
                        eprint(
                            "n0>nX",
                            [(x, gfaL[x]) for x in noveljunctions if x[0] == ix_j[0]],
                        )
                        eprint(
                            "nY>n1",
                            [(x, gfaL[x]) for x in noveljunctions if x[1] == ix_j[1]],
                        )
                        eprint("nX:", nX)
                        eprint("nY:", nY)
                        for _nx, _ny in itertools.product(nX, nY):
                            eprint("pair:", _nx, _ny)
                            _enx = get_set_exons(gfaS, ix_j[0])
                            eprint("exons_nx:", _enx)
                            _eny = get_set_exons(gfaS, ix_j[1])
                            eprint("exons_ny:", _eny)

                            _tnx = set(get_transcript_from_exons(_enx))
                            eprint("TR_nx:", _tnx)
                            _tny = set(get_transcript_from_exons(_eny))
                            eprint("TR_ny:", _tny)

                            for _tr in _tnx & _tny:
                                _fex0 = list(filter(lambda x: x.startswith(_tr), _enx))
                                _fex1 = list(filter(lambda x: x.startswith(_tr), _eny))
                                assert len(_fex0) == len(_fex1) == 1

                                _tex0 = int(_fex0[0].split(".")[-1])
                                _tex1 = int(_fex1[0].split(".")[-1])

                                if abs(_tex0 - _tex1) == 1:
                                    print("****************** NOVEL CE!!!", _tr)
                                    # # Checking that the trascripts belong to the same gene
                                    # for _tj, _te in itertools.product(_trjunc, [_tr]):
                                    #     if transcript2gene[_tj] == transcript2gene[_te]:
                                    #         print("****************** NOVEL CE!!!", _tj, _te)
                                    #         # TODO: CHECKME: continue?

                        eprint("-" * 15)

        from_novel_inside_nonnovel()

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
        "--format",
        help="Minimum read count (default: junctions)",
        dest="format",
        choices=["junctions", "nodes"],
        required=False,
        default="junctions",
    )
    args = parser.parse_args()
    main(args)
