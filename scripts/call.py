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


def get_outgoing_nodes(links: dict, nid: str) ->list:
    return [k[1] for k in links.keys() if k[0] == nid]


def get_incoming_nodes(links: dict, nid: str) ->list:
    return [k[0] for k in links.keys() if k[1] == nid]


def get_outgoing_links(links: dict, nid: str) ->list:
    return [k for k in links.keys() if k[0] == nid]


def get_incoming_links(links: dict, nid: str) ->list:
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


def main():
    gfaS = dict()
    gfaL = dict()
    gfaP = dict()
    junctions = list()
    noveljunctions = list()
    gfa = sys.argv[1]
    for line in open(gfa, "r"):
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

    gtf_path = sys.argv[2]
    transcript2gene = {}
    for line in open(gtf_path):
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

    # Check all the junctions
    def check_nonnovel():
        for ix_j in junctions:
            junc = gfaL[ix_j]
            if junc["RC"] > -10:
                _trjunc = set(map(lambda x: ".".join(x.split(".")[:-2]), junc["JN"]))
                eprint(f"[Checking junction {ix_j}]: {junc}, {_trjunc}")

                # exons_n0 = set(gfaS[ix_j[0]]["EX"])
                # exons_n1 = set(gfaS[ix_j[1]]["EX"])
                exons_n0 = get_set_exons(gfaS, ix_j[0])
                exons_n1 = get_set_exons(gfaS, ix_j[1])
                transcripts_n0 = set(map(lambda x: ".".join(x.split(".")[:-1]), exons_n0))
                transcripts_n1 = set(map(lambda x: ".".join(x.split(".")[:-1]), exons_n1))

                if len(cap := ((set(transcripts_n0) & set(transcripts_n1)) - _trjunc)) > 0:
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
                            for _tj, _te in itertools.product(_trjunc, [_tr]):
                                if transcript2gene[_tj] == transcript2gene[_te]:
                                    print("****************** KNOWN ES!!!", _tj, _te)
                                    # TODO: CHECKME: continue?
                                

                    # Checking for non-novel A5 https://hackmd.io/DoQzt8ceThOwyIdUvQZN3w#Alternative-5%E2%80%99
                    for n in get_outgoing_nodes(gfaL, ix_j[0]):
                        if (
                            len(cap_a5 := ((get_set_exons(gfaS, n) & exons_n0) - exons_n1))
                            > 0
                        ):
                            cap_a5 = set(map(lambda x: ".".join(x.split(".")[:-1]), cap_a5))
                            cap_a5 = cap_a5 & cap
                            if len(cap_a5) > 0:
                                # Checking that the trascripts belong to the same gene
                                for _tj, _te in itertools.product(_trjunc, cap_a5):
                                    if transcript2gene[_tj] == transcript2gene[_te]:
                                        print("****************** KNOWN A5!!!", _tj, _te)
                                        # TODO: CHECKME: continue?

                    # Checking for non-novel A3 https://hackmd.io/DoQzt8ceThOwyIdUvQZN3w#Alternative-3%E2%80%99
                    for n in get_incoming_nodes(gfaL, ix_j[1]):
                        if (
                            len(cap_a3 := ((get_set_exons(gfaS, n) & exons_n1) - exons_n0))
                            > 0
                        ):
                            cap_a3 = set(map(lambda x: ".".join(x.split(".")[:-1]), cap_a3))
                            cap_a3 = cap & cap_a3
                            if len(cap_a3) > 0:
                                # Checking that the trascripts belong to the same gene
                                for _tj, _te in itertools.product(_trjunc, cap_a3):
                                    if transcript2gene[_tj] == transcript2gene[_te]:
                                        print("****************** KNOWN A3!!!", _tj, _te)
                                        # TODO: CHECKME: continue?

                    # Checking for non-novel IR https://hackmd.io/DoQzt8ceThOwyIdUvQZN3w#Intron-retention
                    next_n0 = get_outgoing_nodes(gfaL, ix_j[0])
                    ex_next_n0 = [get_set_exons(gfaS, x) for x in next_n0]
                    prev_n1 = get_incoming_nodes(gfaL, ix_j[1])
                    ex_prev_n1 = [get_set_exons(gfaS, x) for x in prev_n1]
                    cap_ir = exons_n0.intersection(exons_n1, *ex_next_n0, *ex_prev_n1)

                    if len(cap_ir) > 0:
                        cap_ir = set(map(lambda x: ".".join(x.split(".")[:-1]), cap_ir))
                        # Checking that the trascripts belong to the same gene
                        for _tj, _te in itertools.product(_trjunc, cap_ir):
                            if transcript2gene[_tj] == transcript2gene[_te]:
                                print("****************** KNOWN IR!!!", _tj, _te)
                                # TODO: CHECKME: continue?
            eprint("-" * 15)
    check_nonnovel()

    def check_novel():
        def from_single_novel_junctions():
            # Check all novel junctions
            for ix_j in noveljunctions:
                junc = gfaL[ix_j]
                if junc["RC"] > -10:
                    # _trjunc = set(map(lambda x: ".".join(x.split(".")[:-2]), junc["JN"]))
                    _trjunc = set()
                    eprint(f"[Checking junction {ix_j}]: {junc}, {_trjunc}")

                    exons_n0 = get_set_exons(gfaS, ix_j[0]) ; eprint("exons_n0:", exons_n0)
                    exons_n1 = get_set_exons(gfaS, ix_j[1]) ; eprint("exons_n1:", exons_n1)
                    transcripts_n0 = set(map(lambda x: ".".join(x.split(".")[:-1]), exons_n0)) ; eprint("transcripts_n0:", transcripts_n0)
                    transcripts_n1 = set(map(lambda x: ".".join(x.split(".")[:-1]), exons_n1)) ; eprint("transcripts_n1:", transcripts_n1)

                    if len(cap := ((set(transcripts_n0) & set(transcripts_n1)) - _trjunc)) > 0:
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

                        # Checking for novel A3 after / A5 before
                        for n in get_incoming_nodes(gfaL, ix_j[1]):
                            if (
                                len(cap_a3 := ((get_set_exons(gfaS, n) & exons_n1) - exons_n0))
                                > 0
                            ):
                                eprint("cap_a3 EX:", cap_a3)
                                cap_a3 = set(map(lambda x: ".".join(x.split(".")[:-1]), cap_a3))
                                cap_a3 = cap & cap_a3
                                eprint("cap_a3 & cap:", cap_a3, len(cap_a3))
                                if len(cap_a3) > 0:
                                    print("****************** NOVEL A3a / A5b!!!", cap_a3)
                                    # Checking that the trascripts belong to the same gene
                                    # for _tj, _te in itertools.product(_trjunc, cap_a3):
                                    #     print(_tj, _te)
                                    #     if transcript2gene[_tj] == transcript2gene[_te]:
                                    #         print("****************** NOVEL A3!!!", _tj, _te)
                                    #         # TODO: CHECKME: continue?

                        # Checking for novel IR reverse
                        next_n0 = get_outgoing_nodes(gfaL, ix_j[0])
                        ex_next_n0 = [get_set_exons(gfaS, x) for x in next_n0]
                        prev_n1 = get_incoming_nodes(gfaL, ix_j[1])
                        ex_prev_n1 = [get_set_exons(gfaS, x) for x in prev_n1]
                        cap_ir = exons_n0.intersection(exons_n1, *ex_next_n0, *ex_prev_n1)
                        eprint("cap_ir EX:", cap_ir)

                        if len(cap_ir) > 0:
                            cap_ir = set(map(lambda x: ".".join(x.split(".")[:-1]), cap_ir))
                            print("****************** NOVEL IRr!!!", cap_ir)
                            # Checking that the trascripts belong to the same gene
                            # for _tj, _te in itertools.product(_trjunc, cap_ir):
                            #     if transcript2gene[_tj] == transcript2gene[_te]:
                            #         print("****************** NOVEL IRr!!!", _tj, _te)
                            #         # TODO: CHECKME: continue?
                    eprint("-" * 15)

        from_single_novel_junctions()

        # Check potential CE and IR 
        # Known junctions (n0 > n1) that have novel junctions (n0 > nX) and (nY > n1)
        def from_novel_inside_nonnovel():
            for ix_j in junctions:
                junc = gfaL[ix_j]
                if junc["RC"] > -10:
                    _trjunc = set(map(lambda x: ".".join(x.split(".")[:-2]), junc["JN"]))
                    nX = [x[1] for x in noveljunctions if x[0] == ix_j[0]]
                    nY = [x[0] for x in noveljunctions if x[1] == ix_j[1]]

                    if len(nX) > 0 and len(nY) > 0:
                        eprint(f"[Checking junction {ix_j}]: {junc}, {_trjunc}")
                        eprint("n0>nX", [(x, gfaL[x]) for x in noveljunctions if x[0] == ix_j[0]])
                        eprint("nY>n1", [(x, gfaL[x]) for x in noveljunctions if x[1] == ix_j[1]])
                        eprint("nX:", nX)
                        eprint("nY:", nY)
                        for _nx, _ny in itertools.product(nX, nY):
                            eprint("pair:", _nx, _ny)
                            _enx = get_set_exons(gfaS, ix_j[0]) ; eprint("exons_nx:", _enx)
                            _eny = get_set_exons(gfaS, ix_j[1]) ; eprint("exons_ny:", _eny)

                            _tnx = set(get_transcript_from_exons(_enx)) ; eprint("TR_nx:", _tnx)
                            _tny = set(get_transcript_from_exons(_eny)) ; eprint("TR_ny:", _tny)

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

    check_novel()

            # print(f"[Checking junction {ix_j}]: {junc}")
            # for ol in outgoing_l0:
            #     print(outgoing_l0, gfaL[ol])
            #     transcripts_n1 = (
            #         # set(gfaS[ol[1]]["EX"]) if "EX" in gfaS[ol[1]] else set()
            #         set(map(lambda x: ".".join(x.split(".")[:-1]), gfaS[ol[1]]["EX"]))
            #         if "EX" in gfaS[ol[1]]
            #         else set()
            #     )
            #     if len(cap := (set(transcripts_n0) & set(transcripts_n1))) > 0:
            #         exons_n1 = set(gfaS[ol[1]]["EX"]) if "EX" in gfaS[ol[1]] else set()

            #         # Checking for non-novel ES https://hackmd.io/DoQzt8ceThOwyIdUvQZN3w#Exon-skipping
            #         for _tr in cap:
            #             # print(f"[check transcript {_tr}]", exons_n0, exons_n1)
            #             _fex0 = list(filter(lambda x: x.startswith(_tr), exons_n0))
            #             _fex1 = list(filter(lambda x: x.startswith(_tr), exons_n1))
            #             assert(len(_fex0) == len(_fex1) == 1)
            #             # print(f"[check transcript {_tr}] filtered", _fex0, _fex1)
            #             _tex0 = int(_fex0[0].split('.')[-1])
            #             _tex1 = int(_fex1[0].split('.')[-1])

            #             if abs(_tex0 - _tex1) > 1 and not _tr in _trjunc:
            #                 print("****************** KNOWN EXON SKIPPING!!!", _tr)

            #             # Checking for non-novel A5 https://hackmd.io/DoQzt8ceThOwyIdUvQZN3w#Alternative-5%E2%80%99
            #             # elif not _tr in _trjunc:
            #             elif _tex0 - _tex1 == 0:
            #                 print("****************** KNOWN A5!!!", _tr)
            #             # non_junctions = list(filter(lambda x: not "JN" in x, ol))
            #             # print("nonj", non_junctions)

            #         # traverse_se(gfaL, ix_j[0], ix_j[1])
            #     else:
            #         # print("* No known exon skipping")
            #         pass



if __name__ == "__main__":
    main()
