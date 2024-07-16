import sys
import itertools
import re
from math import floor, ceil
import datetime


def dopass(*_):
    pass


eprint = dopass


def fix_tr_(t: str) -> str:
    # we could use some sort of regex here. If last is [RH][0-9]+, then remove it
    if t.split("_")[-1][0] in ["R", "H"]:
        return "_".join(t.split("_")[:-1])
    else:
        return t


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
            # eprint(f"{_v=}")
            if len(_v) >= 2:
                minv = min(_v, key=lambda x: x[0])
                k1 = [minv.copy()]
                maxv = max(_v, key=lambda x: x[0])
                k2 = [maxv.copy()]
                # eprint(f"{minv=} {maxv=}")

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
            # assert jn_w > 0
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
            # assert jn_w > 0
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


def get_outgoing_nodes_old(
    links: dict, nid: str, segments: dict = None, rc: int = -1
) -> list:
    ret = [k[1] for k in links.keys() if k[0] == nid]
    if segments:
        ret = [x for x in ret if segments[x]["NC"] > rc]
    return ret


def get_outgoing_nodes(segments: dict, nid: str, rc: int = -1) -> list:
    ret = segments[nid]["O"]
    if rc > 0:
        ret = [x for x in ret if segments[x]["NC"] > rc]
    return ret


def get_incoming_nodes_old(
    links: dict, nid: str, segments: dict = None, rc: int = -1
) -> list:
    ret = [k[0] for k in links.keys() if k[1] == nid]
    if segments:
        ret = [x for x in ret if segments[x]["NC"] > rc]
    return ret


def get_incoming_nodes(segments: dict, nid: str, rc: int = -1) -> list:
    ret = segments[nid]["I"]
    if rc > 0:
        ret = [x for x in ret if segments[x]["NC"] > rc]
    return ret


# def get_outgoing_links(links: dict, nid: str) -> list:
#     return [k for k in links.keys() if k[0] == nid]


# def get_incoming_links(links: dict, nid: str) -> list:
#     return [k for k in links.keys() if k[1] == nid]


def get_set_exons(nodes: dict, nid: str) -> set:
    return set(nodes[nid]["EX"]) if "EX" in nodes[nid] else set()


def get_transcript_from_exons(exons) -> map:
    return map(lambda x: ".".join(x.split(".")[:-1]), exons)


def get_haplotranscripts(transcripts) -> dict:
    HTs = dict()
    for ht in transcripts:
        t, h = "_".join(ht.split("_")[:-1]), ht.split("_")[-1]
        HTs[t] = HTs[t] | set([h]) if t in HTs else set([h])
    return HTs


def haplotranscripts_to_str(haplotranscripts: dict) -> str:
    ht = []
    for k, Vs in haplotranscripts.items():
        for v in Vs:
            ht.append(f"{k}_{v}")
    return "|".join(ht)


def get_haplotranscripts_from_junction(transcripts: set) -> dict:
    HTs = dict()
    for transcript in transcripts:
        ht = ".".join(transcript.split(".")[:-2])
        t, h = "_".join(ht.split("_")[:-1]), ht.split("_")[-1]
        HTs[t] = HTs[t] | set([h]) if t in HTs else set([h])
    return HTs


def get_haplotranscripts_from_exon(exon: str) -> dict:
    HTs = dict()
    ht = ".".join(exon.split(".")[:-1])
    t, h = "_".join(ht.split("_")[:-1]), ht.split("_")[-1]
    HTs[t] = set([h])
    return HTs

def get_haplotranscripts_from_exons(exons: set) -> dict:
    HTs = dict()
    for exon in exons:
        ht = ".".join(exon.split(".")[:-1])
        t, h = "_".join(ht.split("_")[:-1]), ht.split("_")[-1]
        HTs[t] = HTs[t] | set([h]) if t in HTs else set([h])
    return HTs


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
    return check_junction_retall(ix_j, segments, links, window, rc)[0]


def check_junction_retall(
    ix_j: tuple, segments: dict, links: dict, window: int, rc: int
):
    next_n0 = get_outgoing_nodes(segments, ix_j[0], rc=rc)
    prev_n1 = get_incoming_nodes(segments, ix_j[1], rc=rc)
    is_complete = False

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
        return False, 0, None, False

    while i < window:
        i += 1
        _intron_next = [get_outgoing_nodes(segments, x, rc=rc) for x in _intron_next]
        _intron_prev = [get_incoming_nodes(segments, x, rc=rc) for x in _intron_prev]
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
            is_complete = True
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
                is_complete = True
                break

        else:
            break

    if i == window:
        return True, _subpath_count, (_subpath_n, _subpath_p), is_complete
    return False, 0, None, False


def main(args):
    gfaS = dict()
    gfaL = dict()
    gfaP = dict()

    junctions = set()
    noveljunctions = set()

    print(
        datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "Parsing GFA..",
        file=sys.stderr,
    )
    for line in open(args.GFA, "r"):
        line = line.strip()
        if line.startswith("S"):
            _, nid, seq, *fields = line.split()
            gfaS[nid] = build_attrs(fields, args.d)
            gfaS[nid]["LN"] = len(seq)  # Done to avoid LN in GFA
            # TODO: uncomment if needed
            # gfaS[nid]['seq'] = seq
            gfaS[nid]["I"] = []
            gfaS[nid]["O"] = []
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
                _,  # overlap
                *fields,
            ) = line.split()
            gfaL[(nid_from, nid_to)] = build_attrs(fields, args.d)
            gfaS[nid_from]["O"].append(nid_to)
            gfaS[nid_to]["I"].append(nid_from)
            # NOTE: uncomment if needed
            # gfaL[(nid_from, nid_to)]['overlap'] = overlap
            if "JN" in gfaL[(nid_from, nid_to)]:
                junctions.add((nid_from, nid_to))
            if "ID" in gfaL[(nid_from, nid_to)]:
                noveljunctions.add((nid_from, nid_to))

    print(
        datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "Parsing GTF..",
        file=sys.stderr,
    )
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
            "chrom",
            "gene",
            "strand",
            "transcripts1",
            "transcripts2",
            "transcripts3",
            "nodes1",
            "coverage1",
            "nodes2",
            "coverage2",
            "nodes3",
            "coverage3",
            sep=",",
        )

    # TODO: we may need adjacency lists

    # Check all the junctions
    def check_nonnovel():
        # used = set() # CHECKME: do we want this? can we have a junction to be part of multiple events? probably yes
        for _j in junctions:
            # if ix_j in used:
            #     continue
            if gfaL[_j]["RC"] < args.rca:
                continue
            _ht = get_haplotranscripts_from_junction(gfaL[_j]["JN"])
            _genes = set(transcript2gene[t] for t in _ht)
            if len(_genes) > 1:
                # FIXME: this could be a quite strong assumption
                continue
            _gene = next(iter(_genes))
            _exons1 = get_set_exons(gfaS, _j[0])
            _exons2 = get_set_exons(gfaS, _j[1])

            # we want exons on same gene only
            # FIXME: this can be done way better
            _exons1 = set(e for e in _exons1 if len(_genes & set(transcript2gene[t] for t in get_haplotranscripts_from_exon(e))) > 0)
            _exons2 = set(e for e in _exons2 if len(_genes & set(transcript2gene[t] for t in get_haplotranscripts_from_exon(e))) > 0)

            assert len(_exons1) > 0 and len(_exons2) > 0

            # Find the outgoing junctions of the head of the junction we are checking
            Js1 = set(x for x in junctions if x[0] == _j[0]) - set([_j])
            # Find the incoming junctions of the tail of the junction we are checking
            Js2 = set(x for x in junctions if x[1] == _j[1]) - set([_j])
            
            # filter by weigth # CHECKME: do we want this?
            # Js1 = set(filter(lambda x: gfaL[x]["RC"] >= args.rca, Js1))
            # Js2 = set(filter(lambda x: gfaL[x]["RC"] >= args.rca, Js2))

            # filter by same gene
            Js1 = set(
                filter(
                    lambda x: len(
                        _genes
                        & set(
                            transcript2gene[t]
                            for t in get_haplotranscripts_from_junction(gfaL[x]["JN"])
                        )
                    )
                    > 0,
                    Js1,
                )
            )
            Js2 = set(
                filter(
                    lambda x: len(
                        _genes
                        & set(
                            transcript2gene[t]
                            for t in get_haplotranscripts_from_junction(gfaL[x]["JN"])
                        )
                    )
                    > 0,
                    Js2,
                )
            )

            if "ES" in args.events:
                for j1, j2 in itertools.product(Js1, Js2):
                    haplotranscripts1 = get_haplotranscripts_from_junction(
                        gfaL[j1]["JN"]
                    )
                    haplotranscripts2 = get_haplotranscripts_from_junction(
                        gfaL[j2]["JN"]
                    )

                    haplotranscripts_inclusion = (set(haplotranscripts1) & set(
                        haplotranscripts2
                    )) - set(_ht)
                    if len(haplotranscripts_inclusion) == 0:
                        # no transcript
                        continue

                    # exons1_1 and exons2_2 are already computed for the junction
                    # exons1_1 = get_set_exons(gfaS, j1[0])
                    exons1_2 = get_set_exons(gfaS, j1[1])
                    exons2_1 = get_set_exons(gfaS, j2[0])
                    # exons2_2 = get_set_exons(gfaS, j2[1])

                    # this does not work for multiple exons skipping. commenting
                    # if len(exons1_2) & len(exons2_1) == 0:
                    #     continue
                    if len((exons1_2 | exons2_1) - (_exons1 | _exons2)) == 0:
                        # no change in exon
                        continue

                    # TODO: we could report an event only if on "same" haplotype
                    print(
                        "ES",
                        "annotated",
                        genechr[_gene],
                        _gene,
                        genestrand[_gene],
                        "|".join(gfaL[_j]["JN"]),
                        "|".join([x for x in gfaL[j1]["JN"] if "_".join(x.split("_")[:-1]) in haplotranscripts_inclusion]),
                        "|".join([x for x in gfaL[j2]["JN"] if "_".join(x.split("_")[:-1]) in haplotranscripts_inclusion]),
                        ">".join(_j),
                        gfaL[_j]["RC"],
                        ">".join(j1),
                        gfaL[j1]["RC"],
                        ">".join(j2),
                        gfaL[j2]["RC"],
                        sep=",",
                    )

            if "SS" in args.events:
                # this is A5 on + / A3 on -
                if len(Js2) != 0:
                    for n in get_outgoing_nodes(gfaS, _j[0]):
                        if n == _j[1]:
                            # we don't want the junction we are on
                            continue
                        exons = get_set_exons(gfaS, n) & _exons1
                        if len(exons) == 0:
                            # the exon does not continue
                            continue
                        for j2 in Js2:
                            if j2[0] == _j[0]:
                                # we don't want the junction we are on
                                continue
                            j2_exons = get_set_exons(gfaS, j2[0]) & exons
                            if len(j2_exons) == 0:
                                # exon continues but does not continue with the exon we want wrt to current j2
                                continue
                            ht = get_haplotranscripts_from_exons(j2_exons)
                            # CHECKME: I think this is always False
                            if len(set(ht) - set(_ht)) == 0:
                                # no change in transcript
                                continue
                            assert len(set(gfaL[_j]["JN"]) & set(gfaL[j2]["JN"])) == 0
                            print(
                                "A5" if genestrand[_gene] == "+" else "A3",
                                "annotated",
                                genechr[_gene],
                                _gene,
                                genestrand[_gene],
                                "|".join(gfaL[_j]["JN"]),
                                "|".join(gfaL[j2]["JN"]), # TODO: do we have to filter which transcripts we want? I don't think so. _j and j2 are always annotated with disjoint sets of transcripts 
                                ".",
                                ">".join(_j),
                                gfaL[_j]["RC"],
                                ">".join(j2),
                                gfaL[j2]["RC"],
                                ".",
                                ".",
                                sep=",",
                            )
                # this is A3 on + / A5 on -
                if len(Js1) != 0:
                    for n in get_incoming_nodes(gfaS, _j[1]):
                        if n == _j[0]:
                            # we don't want the junction we are on
                            continue
                        exons = get_set_exons(gfaS, n) & _exons2
                        if len(exons) == 0:
                            # the exon does not continue
                            continue
                        for j1 in Js1:
                            if j1[1] == _j[1]:
                                # we don't want the junction we are on
                                continue
                            j1_exons = get_set_exons(gfaS, j1[1]) & exons
                            if len(j1_exons) == 0:
                                # exon continues but does not continue with the exon we want wrt to current j1
                                continue
                            ht = get_haplotranscripts_from_exons(j1_exons)
                            # CHECKME: I think this is always False
                            if len(set(ht) - set(_ht)) == 0:
                                # no change in transcript
                                continue
                            print(
                                "A3" if genestrand[_gene] == "+" else "A5",
                                "annotated",
                                genechr[_gene],
                                _gene,
                                genestrand[_gene],
                                "|".join(gfaL[_j]["JN"]),
                                "|".join(gfaL[j1]["JN"]), # TODO: filter which transcripts we want
                                ".",
                                ">".join(_j),
                                gfaL[_j]["RC"],
                                ">".join(j1),
                                gfaL[j1]["RC"],
                                ".",
                                ".",
                                sep=",",
                            )

            if "IR" in args.events:
                nnext = get_outgoing_nodes(gfaS, _j[0])
                nprev = get_incoming_nodes(gfaS, _j[1])

                exons = _exons1 & _exons2
                retained_transcripts = {}
                subpath = []
                for n1, n2 in itertools.product(nnext, nprev):
                    if n1 == _j[1] or n2 == _j[0]:
                        # we don't want "same junction" _j
                        continue
                    if n2 < n1:
                        # assuming topological sorting, if we have an exon, we can't have this
                        continue
                    i_exons = get_set_exons(gfaS, n1) & get_set_exons(gfaS, n2) & exons
                    if len(i_exons) == 0:
                        continue

                    e = next(
                        iter(i_exons)
                    )  # CHECKME: we choose one exon to follow. This should work
                    # FIXME: improve this if needed
                    n = n1
                    subpath = [n]
                    while n != n2:
                        nn = -1
                        for nn in get_outgoing_nodes(gfaS, n):
                            # nn must be on same exon AND be smaller than the vertex we want to reach (this should hold if we assume topological sorting)
                            if nn <= n2 and e in get_set_exons(gfaS, nn):
                                break
                        # Here I am assuming that if i_exons is not empty, then there must be a path between n1 and n2 (since there is an exon for sure)
                        if nn == -1:
                            print(
                                "Error while reconstrucing IR path",
                                n1,
                                n2,
                                n,
                                get_outgoing_nodes(gfaS, n),
                                subpath,
                                flush=True,
                            )
                        assert nn != -1
                        subpath.append(nn)
                        n = nn
                    retained_transcripts = get_haplotranscripts_from_exons(i_exons)
                    # CHECKME: we need just one exon, since all exons we can find should produce the same path
                    break

                if len(retained_transcripts) == 0:
                    continue
                print(
                    "IR",
                    "annotated",
                    genechr[_gene],
                    _gene,
                    genestrand[_gene],
                    "|".join(gfaL[_j]["JN"]),
                    ".", # TODO: recover exon if we want
                    ".",
                    ">".join(_j),
                    gfaL[_j]["RC"],
                    ">".join(subpath),
                    ceil(sum([gfaS[x]["NC"] for x in subpath]) / len(subpath)),
                    ".",
                    ".",
                    sep=",",
                )

    if not args.annotated:
        check_nonnovel()

    # def check_novel():
    #     def from_single_novel_junctions():
    #         # Check all novel junctions
    #         for ix_j in noveljunctions:
    #             junc = gfaL[ix_j]
    #             if junc["RC"] >= args.rc:
    #                 _trjunc = set()
    #                 eprint(f"[Checking junction {ix_j}]: {junc}, {_trjunc}")

    #                 exons_n0 = get_set_exons(gfaS, ix_j[0])
    #                 eprint(f"{exons_n0=}")
    #                 exons_n1 = get_set_exons(gfaS, ix_j[1])
    #                 eprint(f"{exons_n1=}")
    #                 transcripts_n0 = set(
    #                     map(lambda x: ".".join(x.split(".")[:-1]), exons_n0)
    #                 )
    #                 eprint(f"{transcripts_n0=}")
    #                 transcripts_n1 = set(
    #                     map(lambda x: ".".join(x.split(".")[:-1]), exons_n1)
    #                 )
    #                 eprint(f"{transcripts_n1=}")

    #                 if (
    #                     len(
    #                         cap := (
    #                             (set(transcripts_n0) & set(transcripts_n1)) - _trjunc
    #                         )
    #                     )
    #                     > 0
    #                 ):
    #                     # cap contains all the trascripts that:
    #                     # 1. visit n0 and n1

    #                     eprint(f"{cap=}")

    #                     # Checking for novel ES
    #                     if "ES" in args.events:
    #                         for _tr in cap:
    #                             _tr = fix_tr_(_tr)
    #                             _fex0 = list(
    #                                 filter(lambda x: x.startswith(_tr), exons_n0)
    #                             )
    #                             _fex1 = list(
    #                                 filter(lambda x: x.startswith(_tr), exons_n1)
    #                             )
    #                             # assert len(_fex0) == len(_fex1) == 1
    #                             assert len(_fex0) > 0 and len(_fex1) > 0

    #                             for _fex0_i, _fex1_j in itertools.product(_fex0, _fex1):
    #                                 _tex0 = int(_fex0_i.split(".")[-1])
    #                                 _tex1 = int(_fex1_j.split(".")[-1])

    #                                 if abs(_tex0 - _tex1) > 1:
    #                                     _n_j1 = [
    #                                         x for x in junctions if x[0] == ix_j[0]
    #                                     ]
    #                                     _n_j2 = [
    #                                         x for x in junctions if x[1] == ix_j[1]
    #                                     ]
    #                                     eprint(f"{_n_j1=}")
    #                                     eprint(f"{_n_j2=}")

    #                                     _es_j1 = list(
    #                                         filter(
    #                                             lambda x: any(
    #                                                 [
    #                                                     y.startswith(_tr)
    #                                                     for y in gfaL[x]["JN"]
    #                                                 ]
    #                                             ),
    #                                             _n_j1,
    #                                         )
    #                                     )
    #                                     _es_j2 = list(
    #                                         filter(
    #                                             lambda x: any(
    #                                                 [
    #                                                     y.startswith(_tr)
    #                                                     for y in gfaL[x]["JN"]
    #                                                 ]
    #                                             ),
    #                                             _n_j2,
    #                                         )
    #                                     )
    #                                     eprint(f"{_es_j1=}")
    #                                     eprint(f"{_es_j2=}")

    #                                     _es_j1_name = [x for x in gfaL[_es_j1[0]]["JN"] if x.startswith(_tr)][0]
    #                                     _es_j2_name = [x for x in gfaL[_es_j2[0]]["JN"] if x.startswith(_tr)][0]
    #                                     if genestrand[transcript2gene[_tr]] == "-":
    #                                         # CHECKME: do we need to swap the names? Which one do we want first?
    #                                         pass
    #                                     eprint(f"{_es_j1_name=}")
    #                                     eprint(f"{_es_j2_name=}")

    #                                     if len(_es_j1) > 0 and len(_es_j2) > 0:
    #                                         print(
    #                                             "ES",
    #                                             "novel",
    #                                             genechr[transcript2gene[_tr]],
    #                                             transcript2gene[_tr],
    #                                             genestrand[transcript2gene[_tr]],
    #                                             "?",  # _j,
    #                                             ">".join(ix_j),
    #                                             f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, ix_j[0], 'LN')}-{get_refpos_node(gfaS, ix_j[1], 'RP')}",
    #                                             junc["RC"],
    #                                             _es_j1_name,
    #                                             ">".join(_es_j1[0]),
    #                                             f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, _es_j1[0][0], 'LN')}-{get_refpos_node(gfaS, _es_j1[0][1], 'RP')}",
    #                                             gfaL[_es_j1[0]]["RC"],
    #                                             _es_j2_name,
    #                                             ">".join(_es_j2[0]),
    #                                             f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, _es_j2[0][0], 'LN')}-{get_refpos_node(gfaS, _es_j2[0][1], 'RP')}",
    #                                             gfaL[_es_j2[0]]["RC"],
    #                                             sep=",",
    #                                         )

    #                     # Checking for novel A5+ before / A3- after
    #                     if "SS" in args.events:
    #                         for _tr in cap:
    #                             _tr = fix_tr_(_tr)
    #                             _fex0 = list(
    #                                 filter(lambda x: x.startswith(_tr), exons_n0)
    #                             )
    #                             _fex1 = list(
    #                                 filter(lambda x: x.startswith(_tr), exons_n1)
    #                             )
    #                             # assert len(_fex0) == len(_fex1) == 1
    #                             assert len(_fex0) > 0 and len(_fex1) > 0
    #                             eprint(f"{_fex0=}")
    #                             eprint(f"{_fex1=}")

    #                             for _fex0_i, _fex1_j in itertools.product(_fex0, _fex1):
    #                                 _tex0 = int(_fex0_i.split(".")[-1])
    #                                 _tex1 = int(_fex1_j.split(".")[-1])

    #                                 # _tex0 = int(_fex0[0].split(".")[-1])
    #                                 # _tex1 = int(_fex1[0].split(".")[-1])
    #                                 if abs(_tex0 - _tex1) == 1:
    #                                     ex_next_n0 = set().union(
    #                                         *[
    #                                             get_set_exons(gfaS, x)
    #                                             for x in get_outgoing_nodes(
    #                                                 gfaS, ix_j[0]
    #                                             )
    #                                         ]
    #                                     )
    #                                     eprint(f"check A5b {ex_next_n0=}")
    #                                     if _fex0[0] in ex_next_n0:
    #                                         _a_j = [
    #                                             x for x in junctions if x[1] == ix_j[1]
    #                                         ]
    #                                         _a_j = list(
    #                                             filter(
    #                                                 lambda x: any(
    #                                                     [
    #                                                         y.startswith(_tr)
    #                                                         for y in gfaL[x]["JN"]
    #                                                     ]
    #                                                 ),
    #                                                 _a_j,
    #                                             )
    #                                         )
    #                                         if len(_a_j) == 1:
    #                                             _a_j = _a_j[0]
    #                                             _a_j_name = [
    #                                                 x
    #                                                 for x in gfaL[_a_j]["JN"]
    #                                                 if x.startswith(_tr)
    #                                             ]
    #                                             assert len(_a_j_name) == 1

    #                                             # A5b+ / A3a-
    #                                             eprint("A5b: A5b+ / A3a-")
    #                                             print(
    #                                                 "A5"
    #                                                 if genestrand[transcript2gene[_tr]]
    #                                                 == "+"
    #                                                 else "A3",
    #                                                 "novel",
    #                                                 genechr[transcript2gene[_tr]],
    #                                                 transcript2gene[_tr],
    #                                                 genestrand[transcript2gene[_tr]],
    #                                                 "?",  # _j,
    #                                                 ">".join(ix_j),
    #                                                 f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, ix_j[0], 'OL', junc['RC'])}-{get_refpos_node(gfaS, ix_j[1], 'RP')}",
    #                                                 junc["RC"],
    #                                                 _a_j_name[0],
    #                                                 ">".join(_a_j),
    #                                                 f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, _a_j[0], 'LN')}-{get_refpos_node(gfaS, _a_j[1], 'RP')}",
    #                                                 gfaL[_a_j]["RC"],
    #                                                 ".",
    #                                                 ".",
    #                                                 ".",
    #                                                 ".",
    #                                                 sep=",",
    #                                             )

    #                                     ex_prev_n1 = set().union(
    #                                         *[
    #                                             get_set_exons(gfaS, x)
    #                                             for x in get_incoming_nodes(
    #                                                 gfaS, ix_j[1]
    #                                             )
    #                                         ]
    #                                     )
    #                                     eprint(f"check A3a {ex_prev_n1=}")
    #                                     if _fex1[0] in ex_prev_n1:
    #                                         _a_j = [
    #                                             x for x in junctions if x[0] == ix_j[0]
    #                                         ]
    #                                         _a_j = list(
    #                                             filter(
    #                                                 lambda x: any(
    #                                                     [
    #                                                         y.startswith(_tr)
    #                                                         for y in gfaL[x]["JN"]
    #                                                     ]
    #                                                 ),
    #                                                 _a_j,
    #                                             )
    #                                         )
    #                                         if len(_a_j) == 1:
    #                                             _a_j = _a_j[0]
    #                                             _a_j_name = [
    #                                                 x
    #                                                 for x in gfaL[_a_j]["JN"]
    #                                                 if x.startswith(_tr)
    #                                             ]
    #                                             assert len(_a_j_name) == 1
    #                                             eprint(f"{_fex1=}")

    #                                             # NOTE: In this case to keep ordering consistent with the reference
    #                                             # the order of trascript is inverted.
    #                                             # The first one is the one the event is considered to
    #                                             # and the second is the one containing the event

    #                                             # A3a+ / A5b-
    #                                             eprint("A3a: A3a+ / A5b-")
    #                                             print(
    #                                                 "A3"
    #                                                 if genestrand[transcript2gene[_tr]]
    #                                                 == "+"
    #                                                 else "A5",
    #                                                 "novel",
    #                                                 genechr[transcript2gene[_tr]],
    #                                                 transcript2gene[_tr],
    #                                                 genestrand[transcript2gene[_tr]],
    #                                                 _a_j_name[0],
    #                                                 ">".join(_a_j),
    #                                                 f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, _a_j[0], 'LN')}-{get_refpos_node(gfaS, _a_j[1], 'RP')}",
    #                                                 gfaL[_a_j]["RC"],
    #                                                 "?",  # _j,
    #                                                 ">".join(ix_j),
    #                                                 f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, ix_j[0], 'LN')}-{get_refpos_node(gfaS, ix_j[1], 'IL', junc['RC'])}",
    #                                                 junc["RC"],
    #                                                 ".",
    #                                                 ".",
    #                                                 ".",
    #                                                 ".",
    #                                                 sep=",",
    #                                             )

    #                     # Checking for novel IR reverse
    #                     if "IR" in args.events:
    #                         next_n0 = get_outgoing_nodes(gfaS, ix_j[0])
    #                         ex_next_n0 = [get_set_exons(gfaS, x) for x in next_n0]
    #                         prev_n1 = get_incoming_nodes(gfaS, ix_j[1])
    #                         ex_prev_n1 = [get_set_exons(gfaS, x) for x in prev_n1]

    #                         ex_next_n0 = set().union(*ex_next_n0)
    #                         ex_prev_n1 = set().union(*ex_prev_n1)

    #                         cap_ir = set.intersection(
    #                             exons_n0, exons_n1, ex_next_n0, ex_prev_n1
    #                         )
    #                         eprint(f"EX {cap_ir=}")

    #                         if len(cap_ir) > 0:
    #                             for ex_ir in cap_ir:
    #                                 _tr = ".".join(ex_ir.split(".")[:-1])
    #                                 # FIXME: this is an heuristic in place of having the paths in the graph
    #                                 (
    #                                     _found,
    #                                     _count_sum,
    #                                     _subpaths,
    #                                     _complete,
    #                                 ) = check_junction_retall(
    #                                     [ix_j[0], ix_j[1]],
    #                                     gfaS,
    #                                     gfaL,
    #                                     args.irw,
    #                                     args.rc,
    #                                 )
    #                                 if not _found:
    #                                     continue
    #                                 # FIXME: the subpath is not correct
    #                                 if _complete:
    #                                     _subpath = _subpaths[0] + _subpaths[1]
    #                                     _splen = len(set(_subpath))
    #                                 else:
    #                                     _subpath = _subpaths[0] + [".."] + _subpaths[1]
    #                                     _splen = len(set(_subpath)) - 1

    #                                 # _subpath = get_path_transcript(
    #                                 #     gfaP, _tr, start=ix_j[0], end=ix_j[1]
    #                                 # )
    #                                 # assert len(_subpath) > 2
    #                                 # # excluding the junction nodes
    #                                 # _subpath = _subpath[1:-1]
    #                                 # _count_sum = 0
    #                                 # for _in in _subpath:
    #                                 #     _count_sum += gfaS[_in].get("NC", 0)
    #                                 eprint("IRr")

    #                                 # if genestrand[transcript2gene[_tr]] == "+":
    #                                 _refpos = f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, ix_j[0], 'LN')}-{get_refpos_node(gfaS, ix_j[1], 'RP')}"
    #                                 # else:
    #                                 #     _refpos = f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, ix_j[1], 'LN')}-{get_refpos_node(gfaS, ix_j[0], 'RP')}"

    #                                 print(
    #                                     "IR",
    #                                     "novel",
    #                                     genechr[transcript2gene[_tr]],
    #                                     transcript2gene[_tr],
    #                                     genestrand[transcript2gene[_tr]],
    #                                     "?",  # _j,
    #                                     ">".join(ix_j),
    #                                     _refpos,
    #                                     junc["RC"],
    #                                     ex_ir,
    #                                     ">".join(_subpath),
    #                                     ".",  # CHECKME: this should not be needed
    #                                     _count_sum // _splen,
    #                                     ".",
    #                                     ".",
    #                                     ".",
    #                                     ".",
    #                                     sep=",",
    #                                 )

    #                 # Check A3 - before
    #                 if "SS" in args.events:
    #                     # if len(exons_n1) == 0:
    #                     if (
    #                         len(
    #                             set(get_transcript_from_exons(exons_n1))
    #                             & transcripts_n0
    #                         )
    #                         == 0
    #                     ):
    #                         # n1 is an intron, check if there is a junction
    #                         # from n0 to somewhere else

    #                         nX_j = [x for x in junctions if x[0] == ix_j[0]]
    #                         nX = [x[1] for x in nX_j]
    #                         if len(nX) > 0:
    #                             eprint(f"{nX=}")
    #                             exons_nX = [get_set_exons(gfaS, x) for x in nX]
    #                             eprint(f"{exons_nX=}")

    #                             tr_list_nX = [
    #                                 list(get_transcript_from_exons(x)) for x in exons_nX
    #                             ]

    #                             # INFO: flatting a list of [[], []]
    #                             # in case we have one more level it breaks. But we shouldn't
    #                             if len(exons_nX) > 1:
    #                                 _tmp_ = []
    #                                 for _x_ in exons_nX:
    #                                     _tmp_ += _x_
    #                                 exons_nX = [_tmp_]

    #                             if len(tr_list_nX) > 1:
    #                                 _tmp_ = []
    #                                 for _x_ in tr_list_nX:
    #                                     _tmp_ += _x_
    #                                 tr_list_nX = [_tmp_]

    #                             transcripts_nX = set(*tr_list_nX)

    #                             eprint(f"{transcripts_nX=}")

    #                             for _tr in transcripts_nX & transcripts_n0:
    #                                 _tr = fix_tr_(_tr)
    #                                 _fex0 = list(
    #                                     filter(lambda x: x.startswith(_tr), exons_n0)
    #                                 )
    #                                 _fexX = list(
    #                                     filter(lambda x: x.startswith(_tr), *exons_nX)
    #                                 )
    #                                 eprint(f"{_fex0=}")
    #                                 eprint(f"{_fexX=}")
    #                                 assert len(_fex0) > 0 and len(_fexX) > 0
    #                                 # assert len(_fex0) == len(_fexX) == 1

    #                                 _fnX = [
    #                                     x
    #                                     for x in nX
    #                                     if any(
    #                                         map(
    #                                             lambda y: y.startswith(_tr),
    #                                             get_set_exons(gfaS, x),
    #                                         )
    #                                     )
    #                                 ]
    #                                 eprint(f"{_fnX=}")

    #                                 for _fex0_i, _fexX_j in itertools.product(
    #                                     _fex0, _fexX
    #                                 ):
    #                                     _tex0 = int(_fex0_i.split(".")[-1])
    #                                     _texX = int(_fexX_j.split(".")[-1])

    #                                     if abs(_tex0 - _texX) == 1:
    #                                         _a_j = list(
    #                                             filter(
    #                                                 lambda x: any(
    #                                                     [
    #                                                         y.startswith(_tr)
    #                                                         for y in gfaL[x]["JN"]
    #                                                     ]
    #                                                 ),
    #                                                 nX_j,
    #                                             )
    #                                         )
    #                                         if len(_a_j) == 1 and any(
    #                                             check_junction(
    #                                                 [ix_j[1], _x],
    #                                                 gfaS,
    #                                                 gfaL,
    #                                                 args.irw,
    #                                                 args.rc,
    #                                             )
    #                                             for _x in _fnX
    #                                         ):
    #                                             eprint(f"{_a_j=} {gfaL[_a_j[0]]=}")
    #                                             _a_j = _a_j[0]
    #                                             _a_j_name = [
    #                                                 x
    #                                                 for x in gfaL[_a_j]["JN"]
    #                                                 if x.startswith(_tr)
    #                                             ]
    #                                             assert len(_a_j_name) >= 1

    #                                             eprint("A3b: A3b+ / A5a-")
    #                                             for _ajn in _a_j_name:
    #                                                 print(
    #                                                     "A3"
    #                                                     if genestrand[
    #                                                         transcript2gene[_tr]
    #                                                     ]
    #                                                     == "+"
    #                                                     else "A5",
    #                                                     "novel",
    #                                                     genechr[transcript2gene[_tr]],
    #                                                     transcript2gene[_tr],
    #                                                     genestrand[
    #                                                         transcript2gene[_tr]
    #                                                     ],
    #                                                     _ajn,
    #                                                     ">".join(_a_j),
    #                                                     f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, _a_j[0], 'LN')}-{get_refpos_node(gfaS, _a_j[1], 'RP')}",
    #                                                     gfaL[_a_j]["RC"],
    #                                                     "?",  # _j,
    #                                                     ">".join(ix_j),
    #                                                     f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, ix_j[0], 'LN')}-{get_refpos_node(gfaS, ix_j[1], 'IL', junc['RC'])}",
    #                                                     junc["RC"],
    #                                                     ".",
    #                                                     ".",
    #                                                     ".",
    #                                                     ".",
    #                                                     sep=",",
    #                                                 )

    #                     # Chek A5 - after
    #                     # if len(exons_n0) == 0:
    #                     if (
    #                         len(
    #                             set(get_transcript_from_exons(exons_n0))
    #                             & transcripts_n1
    #                         )
    #                         == 0
    #                     ):
    #                         # n0 is an intron, check if there is a junction
    #                         # to n1 from somewhere else

    #                         nX_j = [x for x in junctions if x[1] == ix_j[1]]
    #                         nX = [x[0] for x in nX_j]
    #                         if len(nX) > 0:
    #                             eprint(f"{nX=}")
    #                             exons_nX = [get_set_exons(gfaS, x) for x in nX]
    #                             eprint(f"{exons_nX=}")
    #                             tr_list_nX = [
    #                                 list(get_transcript_from_exons(x)) for x in exons_nX
    #                             ]

    #                             # INFO: flatting a list of [[], []]
    #                             # in case we have one more level it breaks. But we shouldn't
    #                             if len(exons_nX) > 1:
    #                                 _tmp_ = []
    #                                 for _x_ in exons_nX:
    #                                     _tmp_ += _x_
    #                                 exons_nX = [_tmp_]

    #                             if len(tr_list_nX) > 1:
    #                                 _tmp_ = []
    #                                 for _x_ in tr_list_nX:
    #                                     _tmp_ += _x_
    #                                 tr_list_nX = [_tmp_]

    #                             transcripts_nX = set(*tr_list_nX)
    #                             eprint(f"{transcripts_nX=}")

    #                             for _tr in transcripts_nX & transcripts_n1:
    #                                 eprint(f"{_tr=}")
    #                                 _tr = fix_tr_(_tr)

    #                                 _fex0 = list(
    #                                     filter(lambda x: x.startswith(_tr), exons_n1)
    #                                 )
    #                                 _fexX = list(
    #                                     filter(lambda x: x.startswith(_tr), *exons_nX)
    #                                 )
    #                                 eprint(f"{_fex0=}")
    #                                 eprint(f"{_fexX=}")
    #                                 assert len(_fex0) > 0 and len(_fexX) > 0

    #                                 _fnX = [
    #                                     x
    #                                     for x in nX
    #                                     if any(
    #                                         map(
    #                                             lambda y: y.startswith(_tr),
    #                                             get_set_exons(gfaS, x),
    #                                         )
    #                                     )
    #                                 ]
    #                                 eprint(f"{_fnX=}")

    #                                 for _fex0_i, _fexX_j in itertools.product(
    #                                     _fex0, _fexX
    #                                 ):
    #                                     _tex0 = int(_fex0_i.split(".")[-1])
    #                                     _texX = int(_fexX_j.split(".")[-1])
    #                                     if abs(_tex0 - _texX) == 1:
    #                                         _a_j = list(
    #                                             filter(
    #                                                 lambda x: any(
    #                                                     [
    #                                                         y.startswith(_tr)
    #                                                         for y in gfaL[x]["JN"]
    #                                                     ]
    #                                                 ),
    #                                                 nX_j,
    #                                             )
    #                                         )
    #                                         if len(_a_j) == 1 and any(
    #                                             check_junction(
    #                                                 [_x, ix_j[0]],
    #                                                 gfaS,
    #                                                 gfaL,
    #                                                 args.irw,
    #                                                 args.rc,
    #                                             )
    #                                             for _x in _fnX
    #                                         ):
    #                                             eprint(f"{_a_j=}  {gfaL[_a_j[0]]=}")
    #                                             _a_j = _a_j[0]
    #                                             _a_j_name = [
    #                                                 x
    #                                                 for x in gfaL[_a_j]["JN"]
    #                                                 if x.startswith(_tr)
    #                                             ]
    #                                             assert len(_a_j_name) == 1
    #                                             eprint("A5a: A5a+ / A3b-")
    #                                             print(
    #                                                 "A5"
    #                                                 if genestrand[transcript2gene[_tr]]
    #                                                 == "+"
    #                                                 else "A3",
    #                                                 "novel",
    #                                                 genechr[transcript2gene[_tr]],
    #                                                 transcript2gene[_tr],
    #                                                 genestrand[transcript2gene[_tr]],
    #                                                 "?",  # _j,
    #                                                 ">".join(ix_j),
    #                                                 f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, ix_j[0], 'OL', junc['RC'])}-{get_refpos_node(gfaS, ix_j[1], 'RP')}",
    #                                                 junc["RC"],
    #                                                 _a_j_name[0],
    #                                                 ">".join(_a_j),
    #                                                 f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, _a_j[0], 'LN')}-{get_refpos_node(gfaS, _a_j[1], 'RP')}",
    #                                                 gfaL[_a_j]["RC"],
    #                                                 ".",
    #                                                 ".",
    #                                                 ".",
    #                                                 ".",
    #                                                 sep=",",
    #                                             )

    #                 eprint("-" * 15)

    #     from_single_novel_junctions()

    #     for ix_j in junctions:
    #         junc = gfaL[ix_j]

    #         # Check potential CE
    #         if "SE" in args.events:
    #             # Known junctions (n0 > n1) that have novel junctions (n0 > nX) and (nY > n1)

    #             _trjunc = set(map(lambda x: ".".join(x.split(".")[:-2]), junc["JN"]))
    #             # nX = [x[1] for x in noveljunctions if x[0] == ix_j[0]]
    #             # nY = [x[0] for x in noveljunctions if x[1] == ix_j[1]]
    #             nX = [
    #                 x
    #                 for x in noveljunctions
    #                 if x[0] == ix_j[0] and gfaS[x[0]].get("NC", 0) >= args.rc
    #             ]
    #             nY = [
    #                 x
    #                 for x in noveljunctions
    #                 if x[1] == ix_j[1] and gfaS[x[1]].get("NC", 0) >= args.rc
    #             ]

    #             if len(nX) > 0 and len(nY) > 0:
    #                 eprint(f"[Checking junction {ix_j}]: {junc}, {_trjunc}")
    #                 eprint(
    #                     f"n0>nX= "
    #                     + f"{[(x, gfaL[x]) for x in noveljunctions if x[0] == ix_j[0]]}"
    #                 )
    #                 eprint(
    #                     f"nY>n1= "
    #                     + f"{[(x, gfaL[x]) for x in noveljunctions if x[1] == ix_j[1]]}"
    #                 )
    #                 eprint(f"{nX=}")
    #                 eprint(f"{nY=}")
    #                 _enx = get_set_exons(gfaS, ix_j[0])
    #                 eprint(f"exons_nx: {_enx}")
    #                 _eny = get_set_exons(gfaS, ix_j[1])
    #                 eprint(f"exons_ny: {_eny}")

    #                 for _nx, _ny in itertools.product(nX, nY):
    #                     eprint(f"pair: {_nx} - {_ny}")

    #                     _tnx = set(get_transcript_from_exons(_enx))
    #                     eprint(f"TR_nx: {_tnx}")
    #                     _tny = set(get_transcript_from_exons(_eny))
    #                     eprint(f"TR_ny: {_tny}")

    #                     for _tr in _tnx & _tny:
    #                         _fex0 = list(filter(lambda x: x.startswith(_tr), _enx))
    #                         _fex1 = list(filter(lambda x: x.startswith(_tr), _eny))
    #                         assert len(_fex0) == len(_fex1) == 1

    #                         _tex0 = int(_fex0[0].split(".")[-1])
    #                         _tex1 = int(_fex1[0].split(".")[-1])

    #                         if abs(_tex0 - _tex1) == 1:
    #                             print(
    #                                 "CE",
    #                                 "novel",
    #                                 genechr[transcript2gene[_tr]],
    #                                 transcript2gene[_tr],
    #                                 genestrand[transcript2gene[_tr]],
    #                                 # Intron on annotation
    #                                 f"{_tr}.{min(_tex0, _tex1)}.{max(_tex0, _tex1)}",
    #                                 ">".join(ix_j),
    #                                 f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, ix_j[0], 'LN')}-{get_refpos_node(gfaS, ix_j[1], 'RP')}",
    #                                 junc["RC"],
    #                                 # Cassette junction 1
    #                                 "?",
    #                                 ">".join(_nx),
    #                                 # CHECKME: maybe this is not alway true, but yes
    #                                 f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, _nx[0], 'LN')}-{get_refpos_node(gfaS, _nx[1], 'IL', junc['RC'])}",
    #                                 gfaL[_nx]["RC"],
    #                                 # Cassette junction 2
    #                                 "?",
    #                                 ">".join(_ny),
    #                                 # CHECKME: maybe this is not alway true, but yes
    #                                 f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, _ny[0], 'OL', junc['RC'])}-{get_refpos_node(gfaS, _ny[1], 'RP')}",
    #                                 gfaL[_ny]["RC"],
    #                                 sep=",",
    #                             )

    #             eprint("-" * 15)

    #         # checking for IR
    #         if "IR" in args.events:
    #             _trjunc = set(map(lambda x: ".".join(x.split(".")[:-2]), junc["JN"]))
    #             eprint(f"[IIR Checking junction {ix_j}]: {junc}, {_trjunc}")

    #             next_n0 = get_outgoing_nodes(gfaS, ix_j[0], rc=args.rc)
    #             prev_n1 = get_incoming_nodes(gfaS, ix_j[1], rc=args.rc)

    #             eprint(f"pre {next_n0=}")
    #             eprint(f"pre {prev_n1=}")

    #             _intron_next = set(next_n0) - set(ix_j)
    #             _intron_prev = set(prev_n1) - set(ix_j)
    #             i = 0
    #             eprint(f"{i=} {_intron_next=}")
    #             eprint(f"{i=} {_intron_prev=}")

    #             _subpath_n = []
    #             _subpath_p = []
    #             _subpath_count = 0
    #             if len(_intron_next) > 0 and len(_intron_prev) > 0:
    #                 _max_n = max(
    #                     [(x, gfaS[x]["NC"]) for x in _intron_next], key=lambda x: x[1]
    #                 )
    #                 _max_p = max(
    #                     [(x, gfaS[x]["NC"]) for x in _intron_prev], key=lambda x: x[1]
    #                 )

    #                 _subpath_n.append(_max_n[0])
    #                 _subpath_p.append(_max_p[0])
    #                 _subpath_count += _max_n[1] + _max_p[1]

    #                 eprint(f"{i=} {_subpath_n=}")
    #                 eprint(f"{i=} {_subpath_p=}")
    #             else:
    #                 continue

    #             _subpath_total = False

    #             while i < args.irw:
    #                 i += 1
    #                 _intron_next = [
    #                     get_outgoing_nodes(gfaS, x, rc=args.rc) for x in _intron_next
    #                 ]
    #                 _intron_prev = [
    #                     get_incoming_nodes(gfaS, x, rc=args.rc) for x in _intron_prev
    #                 ]
    #                 # flatten lists
    #                 _intron_next = [x for y in _intron_next for x in y]
    #                 _intron_prev = [x for y in _intron_prev for x in y]
    #                 _intron_next = set(_intron_next) - set(ix_j)
    #                 _intron_prev = set(_intron_prev) - set(ix_j)

    #                 eprint(f"{i=} {_intron_next=}")
    #                 eprint(f"{i=} {_intron_prev=}")

    #                 if len(_intron_next & _intron_prev) > 0:
    #                     _subpath_total = True
    #                     _max_n = max(
    #                         [(x, gfaS[x]["NC"]) for x in _intron_next & _intron_prev],
    #                         key=lambda x: x[1],
    #                     )
    #                     _subpath_n.append(_max_n[0])
    #                     _subpath_count += _max_n[1]
    #                     i = args.irw
    #                     break

    #                 if len(_intron_next) > 0 and len(_intron_prev) > 0:
    #                     _max_n = max(
    #                         [(x, gfaS[x]["NC"]) for x in _intron_next],
    #                         key=lambda x: x[1],
    #                     )
    #                     _max_p = max(
    #                         [(x, gfaS[x]["NC"]) for x in _intron_prev],
    #                         key=lambda x: x[1],
    #                     )

    #                     _subpath_n.append(_max_n[0])
    #                     _subpath_p.append(_max_p[0])
    #                     _subpath_count += _max_n[1] + _max_p[1]

    #                     eprint(f"{i=} {_subpath_n=}")
    #                     eprint(f"{i=} {_subpath_p=}")

    #                     if len(set(_subpath_n) & set(_subpath_p)) > 1:
    #                         _subpath_total = True
    #                         i = args.irw
    #                         break

    #                 else:
    #                     break

    #             if i == args.irw:
    #                 _subpath = _subpath_n + _subpath_p[::-1]
    #                 _subpath_name = ">".join(_subpath_n + ["?"] + _subpath_p[::-1])
    #                 eprint(f"{i=} {_subpath_n=}")
    #                 eprint(f"{i=} {_subpath_p=}")
    #                 eprint(f"{i=} {_subpath=}")
    #                 # Must be done before collapsing because collapsed nodes
    #                 # are counted twice
    #                 _subpath_avg = _subpath_count // len(_subpath)
    #                 if _subpath_total:
    #                     _subpath = list(dict.fromkeys(_subpath))
    #                     _subpath_name = ">".join(_subpath)
    #                 eprint(f"{i=} {_subpath=}")

    #                 for _j in junc["JN"]:
    #                     _tr = ".".join(_j.split(".")[:-2])
    #                     _tr = fix_tr_(_tr)

    #                     _refpos = f"{genechr[transcript2gene[_tr]]}:{get_refpos_node(gfaS, ix_j[0], 'LN')}-{get_refpos_node(gfaS, ix_j[1], 'RP')}"

    #                     eprint("IR")
    #                     print(
    #                         "IR",
    #                         "novel",
    #                         genechr[transcript2gene[_tr]],
    #                         transcript2gene[_tr],
    #                         genestrand[transcript2gene[_tr]],
    #                         _j,
    #                         ">".join(ix_j),
    #                         _refpos,
    #                         junc["RC"],
    #                         "?",  # ex_ir,
    #                         _subpath_name,
    #                         ".",
    #                         _subpath_avg,
    #                         ".",
    #                         ".",
    #                         ".",
    #                         ".",
    #                         sep=",",
    #                     )

    # if args.novel:
    #     check_novel()


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        prog="Caller",
        description="",
    )
    parser.add_argument("GFA", help="Spliced pangenome in GFA format")
    parser.add_argument("GTF", help="Annotation in GTF format")
    parser.add_argument(
        "--rp",
        help='Reduceed spliced pangenome reference paths (default: "")',
        dest="RP",
        type=str,
        required=False,
        default="",
    )
    parser.add_argument(
        "--rc",
        help="Minimum read count (default: 3)",
        dest="rc",
        type=int,
        required=False,
        default=3,
    )
    parser.add_argument(
        "--rca",
        help="Minimum read count for annotated events (default: -1)",
        dest="rca",
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
        "--events",
        help="Events to call (default: [ES, SS, IR])",
        dest="events",
        nargs="+",
        required=False,
        default=["ES", "SS", "IR"],
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
        eprint = print
    main(args)
