import sys
import itertools
import re
from math import floor, ceil
import datetime


def dopass(*_, file=None):
    pass


eprint = dopass


def collapse_linkcounts(lc: list):
    count = sum([x[1] for x in lc])
    # weighted sum of positions
    r = [(x[0] * x[1]) / count for x in lc]
    pos = int(floor(sum(r)))
    # eprint(f"{lc=} {count=} {pos=}")
    return [pos, count]


# FIXME: d was in CLI but in the current version we do not need it
def build_attrs(fields: str, d: int = 3):
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


def get_outgoing_nodes(segments: dict, nid: str, rc: int = -1) -> list:
    ret = segments[nid]["O"]
    if rc > 0:
        ret = [x for x in ret if segments[x]["NC"] > rc]
    return ret


def get_incoming_nodes(segments: dict, nid: str, rc: int = -1) -> list:
    ret = segments[nid]["I"]
    if rc > 0:
        ret = [x for x in ret if segments[x]["NC"] > rc]
    return ret


def get_set_exons(nodes: dict, nid: str) -> set:
    return set(nodes[nid]["EX"]) if "EX" in nodes[nid] else set()


def get_haplotranscripts(transcripts) -> dict:
    HTs = dict()
    for ht in transcripts:
        t, h = "_".join(ht.split("_")[:-1]), ht.split("_")[-1]
        HTs[t] = HTs[t] | set([h]) if t in HTs else set([h])
    return HTs


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


def haplotranscripts_to_str(haplotranscripts: dict) -> str:
    ht = []
    for k, Vs in haplotranscripts.items():
        for v in Vs:
            ht.append(f"{k}_{v}")
    return "|".join(ht)


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
            gfaS[nid] = build_attrs(fields)
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
            gfaL[(nid_from, nid_to)] = build_attrs(fields)
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
        for _j in junctions:
            if args.junction != None and "f{_j[0]}-{_j[1]}" != args.junction:
                continue
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

            eprint(f"Checking annotated junction {_j[0]} -> {_j[1]}", file=sys.stderr)

            # we want exons on same gene only
            # FIXME: this can be done way better
            _exons1 = set(
                e
                for e in _exons1
                if len(
                    _genes
                    & set(transcript2gene[t] for t in get_haplotranscripts_from_exon(e))
                )
                > 0
            )
            _exons2 = set(
                e
                for e in _exons2
                if len(
                    _genes
                    & set(transcript2gene[t] for t in get_haplotranscripts_from_exon(e))
                )
                > 0
            )

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
                eprint("Checking annotated ES", file=sys.stderr)
                for j1, j2 in itertools.product(Js1, Js2):
                    haplotranscripts1 = get_haplotranscripts_from_junction(
                        gfaL[j1]["JN"]
                    )
                    haplotranscripts2 = get_haplotranscripts_from_junction(
                        gfaL[j2]["JN"]
                    )

                    haplotranscripts_inclusion = (
                        set(haplotranscripts1) & set(haplotranscripts2)
                    ) - set(_ht)
                    if len(haplotranscripts_inclusion) > 0:
                        # we have some transcript that includes exons in between the junction

                        # exons1_1 and exons2_2 are already computed for the junction
                        # exons1_1 = get_set_exons(gfaS, j1[0])
                        exons1_2 = get_set_exons(gfaS, j1[1])
                        exons2_1 = get_set_exons(gfaS, j2[0])
                        # exons2_2 = get_set_exons(gfaS, j2[1])

                        # this does not work for multiple exons skipping. commenting
                        # if len(exons1_2) & len(exons2_1) == 0:
                        #     continue
                        if len((exons1_2 | exons2_1) - (_exons1 | _exons2)) > 0:
                            # we have new exons in between the junction
                            # TODO: we could report an event only if on "same" haplotype
                            # but maybe we are already doing so
                            print(
                                "ES",
                                "annotated",
                                genechr[_gene],
                                _gene,
                                genestrand[_gene],
                                "|".join(gfaL[_j]["JN"]),
                                "|".join(
                                    [
                                        x
                                        for x in gfaL[j1]["JN"]
                                        if "_".join(x.split("_")[:-1])
                                        in haplotranscripts_inclusion
                                    ]
                                ),
                                "|".join(
                                    [
                                        x
                                        for x in gfaL[j2]["JN"]
                                        if "_".join(x.split("_")[:-1])
                                        in haplotranscripts_inclusion
                                    ]
                                ),
                                ">".join(_j),
                                gfaL[_j]["RC"],
                                ">".join(j1),
                                gfaL[j1]["RC"],
                                ">".join(j2),
                                gfaL[j2]["RC"],
                                sep=",",
                            )

            if "SS" in args.events:
                eprint("Checking annotated SS", file=sys.stderr)
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
                                "|".join(
                                    gfaL[j2]["JN"]
                                ),  # TODO: do we have to filter which transcripts we want? I don't think so. _j and j2 are always annotated with disjoint sets of transcripts
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
                                "|".join(
                                    gfaL[j1]["JN"]
                                ),  # TODO: filter which transcripts we want
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
                eprint("Checking annotated IR", file=sys.stderr)
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
                        assert nn != -1, "Error while reconstrucing IR path"
                        subpath.append(nn)
                        n = nn
                    retained_transcripts = get_haplotranscripts_from_exons(i_exons)
                    # CHECKME: we need just one exon, since all exons we can find should produce the same path
                    break

                if len(retained_transcripts) > 0:
                    print(
                        "IR",
                        "annotated",
                        genechr[_gene],
                        _gene,
                        genestrand[_gene],
                        "|".join(gfaL[_j]["JN"]),
                        ".",  # TODO: recover exon if we want
                        ".",
                        ">".join(_j),
                        gfaL[_j]["RC"],
                        ">".join(subpath),
                        ceil(
                            sum([gfaS[x]["NC"] if x in gfaS else 0 for x in subpath])
                            / len(subpath)
                        ),
                        ".",
                        ".",
                        sep=",",
                    )

    if not args.annotated:
        check_nonnovel()

    def check_novel():
        # Check all novel junctions
        for _j in noveljunctions:
            if args.junction != None and f"{_j[0]}-{_j[1]}" != args.junction:
                continue
            if gfaL[_j]["RC"] < args.rca:
                continue
            _exons0 = get_set_exons(gfaS, _j[0])
            _exons1 = get_set_exons(gfaS, _j[1])
            if len(_exons0) == 0 and len(_exons1) == 0:
                # at least one end of the junction must be an annotated exon
                continue
            eprint(f"Checking novel junction {_j[0]} -> {_j[1]}", file=sys.stderr)
            _ht0 = get_haplotranscripts_from_exons(_exons0)
            _ht1 = get_haplotranscripts_from_exons(_exons1)

            cap = set(_ht0) & set(_ht1)  # TODO: change name

            _next0 = get_outgoing_nodes(gfaS, _j[0])
            _prev1 = get_incoming_nodes(gfaS, _j[1])

            if "ES" in args.events:
                eprint(f"Checking novel ES", file=sys.stderr)
                if len(_exons0) != 0 and len(_exons1) != 0 and len(_exons0 & _exons1) != len(_exons0):
                    # both ends must be annotated exons that are different
                    nodes1 = [n for n in _next0 if (_j[0], n) in junctions]
                    nodes2 = [p for p in _prev1 if (p, _j[1]) in junctions]
                    if len(nodes1) != 0 and len(nodes2) != 0:
                        # we have junctions to check
                        for n, p in itertools.product(nodes1, nodes2):
                            j1 = (_j[0], n)
                            j2 = (p, _j[1])
                            nht = get_haplotranscripts_from_junction(gfaL[j1]["JN"])
                            pht = get_haplotranscripts_from_junction(gfaL[j2]["JN"])
                            haplotranscripts_inclusion = set(nht) & set(pht)
                            if len(haplotranscripts_inclusion) == 0:
                                # no single transcript covers both nodes
                                continue
                            _genes = set(
                                transcript2gene[t] for t in haplotranscripts_inclusion
                            )
                            if len(_genes) > 1:
                                # CHECKME: do we need this here? We already checked this for the junction
                                # FIXME: this could be a quite strong assumption
                                print(
                                    "Skipping ES due to multiple genes",
                                    file=sys.sterr,
                                )
                            else:
                                _gene = next(iter(_genes))

                                print(
                                    "ES",
                                    "novel",
                                    genechr[_gene],
                                    _gene,
                                    genestrand[_gene],
                                    "?",
                                    "|".join(
                                        [
                                            x
                                            for x in gfaL[j1]["JN"]
                                            if "_".join(x.split("_")[:-1])
                                            in haplotranscripts_inclusion
                                        ]
                                    ),
                                    "|".join(
                                        [
                                            x
                                            for x in gfaL[j2]["JN"]
                                            if "_".join(x.split("_")[:-1])
                                            in haplotranscripts_inclusion
                                        ]
                                    ),
                                    ">".join(_j),
                                    gfaL[_j]["RC"],
                                    ">".join(j1),
                                    gfaL[j1]["RC"],
                                    ">".join(j2),
                                    gfaL[j2]["RC"],
                                    sep=",",
                                )

            if "SS" in args.events:
                if len(cap) != 0:
                    eprint(f"Checking novel SS (1)", file=sys.stderr)
                    # we may have an "internal" SS since both exons share a transcript
                    # -
                    # let's check second exonic node for A3+ or A5-
                    spliced_exons = set()
                    if _j[0] in gfaS[_j[1]]["I"]:
                        # the novel junction breaks the node
                        # FIXME: how can we know where? do we need to "change" the reported junction somehow?
                        spliced_exons = _exons1
                    else:
                        # iterate over parents of second exonic vertex and see if is on some of the transcript we are interested in
                        for p in _prev1:
                            # TODO: check if this is correct. But I'm not sure this is common with novel events
                            spliced_exons = get_set_exons(gfaS, p) & _exons1
                    if len(spliced_exons) > 0:
                        # we first need to find the annotated junctions
                        annotated_js = []
                        for n in _next0:
                            exons_n = get_set_exons(gfaS, n) & spliced_exons
                            # check only those exonic vertices on one of the spliced exons we are interested in
                            if len(exons_n) == 0:
                                continue
                            ht_n = get_haplotranscripts_from_exons(exons_n)
                            if (_j[0], n) in junctions and len(
                                set(ht_n) & cap
                            ) != 0:  # CHECKME: the and may be useless, we already checked the exons. Do we need to check the transcripts?
                                annotated_js.append([(_j[0], n), set(ht_n) & cap])
                        for j1, transcripts in annotated_js:
                            _genes = set(transcript2gene[t] for t in transcripts)
                            if len(_genes) > 1:
                                # CHECKME: do we need this here? We already checked this for the junction
                                # FIXME: this could be a quite strong assumption
                                print(
                                    "Skipping SS due to multiple genes",
                                    file=sys.sterr,
                                )
                            else:
                                _gene = next(iter(_genes))
                                print(
                                    "A3" if genestrand[_gene] == "+" else "A5",
                                    "novel",
                                    genechr[_gene],
                                    _gene,
                                    genestrand[_gene],
                                    "?",
                                    "|".join(gfaL[j1]["JN"]),
                                    ".",
                                    ">".join(_j),
                                    gfaL[_j]["RC"],
                                    ">".join(j1),
                                    gfaL[j1]["RC"],
                                    ".",
                                    ".",
                                    sep=",",
                                )
                    # -
                    # let's check first exonic node for A5+ or A3-
                    # these are the exons what can be spliced by the novel junction
                    spliced_exons = set()
                    if _j[1] in gfaS[_j[0]]["O"]:
                        # the novel junction breaks the node
                        # FIXME: how can we know where? do we need to "change" the reported junction somehow?
                        spliced_exons = _exons0
                    else:
                        # iterate over parents of second exonic vertex and see if is on some of the transcript we are interested in
                        for n in _next0:
                            # TODO: check if this is correct. But I'm not sure this is common with novel events
                            spliced_exons = get_set_exons(gfaS, p) & _exons0
                    if len(spliced_exons) != 0:
                        # we first need to find the annotated junctions
                        annotated_js = []
                        for p in _prev1:
                            exons_p = get_set_exons(gfaS, p) & spliced_exons
                            # check only those exonic vertices on one of the spliced exons we are interested in
                            if len(exons_p) == 0:
                                continue
                            ht_p = get_haplotranscripts_from_exons(exons_p)
                            if (p, _j[1]) in junctions and len(
                                set(ht_p) & cap
                            ) != 0:  # CHECKME: the and may be useless, we already checked the exons. Do we need to check the transcripts?
                                annotated_js.append([(p, _j[1]), set(ht_p) & cap])
                        for j1, transcripts in annotated_js:
                            _genes = set(transcript2gene[t] for t in transcripts)
                            if len(_genes) > 1:
                                # CHECKME: do we need this here? We already checked this for the junction
                                # FIXME: this could be a quite strong assumption
                                print(
                                    "Skipping SS due to multiple genes",
                                    file=sys.sterr,
                                )
                            else:
                                _gene = next(iter(_genes))
                                print(
                                    ("A5" if genestrand[_gene] == "+" else "A3"),
                                    "novel",
                                    genechr[_gene],
                                    _gene,
                                    genestrand[_gene],
                                    "?",
                                    "|".join(gfaL[j1]["JN"]),
                                    ".",
                                    ">".join(_j),
                                    gfaL[_j]["RC"],
                                    ">".join(j1),
                                    gfaL[j1]["RC"],
                                    ".",
                                    ".",
                                    sep=",",
                                )
                # ---
                eprint(f"Checking novel SS (2)", file=sys.stderr)
                # in any case, we may have an intronic SS
                if len(_exons0) > 0 and len(_exons1) == 0:
                    # second vertex is not on exon. So intronic A3+ or A5-
                    exonic_next = set(n for n in _next0 if (_j[0], n) in junctions)
                    # we have an event if we can reach one of the exonic nodes from _j[1]
                    visit = set([_j[1]])
                    pvisitl = 1
                    _i = 0
                    # TODO: we don't need to store the subpath, but we may want it
                    while (
                        len(visit & exonic_next) == 0 and _i < args.isw
                    ):  # FIXME: hardcoded
                        n = visit.pop()
                        pvisitl -= 1
                        visit |= set(get_outgoing_nodes(gfaS, n))
                        if pvisitl == 0:
                            _i += 1
                            pvisitl = len(visit)
                    # TODO: here we are reporting only one event per novel junction. We could do a visit for **each** exonic_next
                    if _i < args.isw:
                        # we report the event since we found a path
                        j1 = (_j[0], next(iter(visit & exonic_next)))
                        _genes = set(
                            transcript2gene[t]
                            for t in get_haplotranscripts_from_junction(gfaL[j1]["JN"])
                        )
                        if len(_genes) > 1:
                            # CHECKME: do we need this here? We already checked this for the junction
                            # FIXME: this could be a quite strong assumption
                            print(
                                "Skipping SS due to multiple genes",
                                file=sys.sterr,
                            )
                        else:
                            _gene = next(iter(_genes))
                            print(
                                "A3" if genestrand[_gene] == "+" else "A5",
                                "novel",
                                genechr[_gene],
                                _gene,
                                genestrand[_gene],
                                "|".join(gfaL[j1]["JN"]),
                                "?",
                                ".",
                                ">".join(j1),
                                gfaL[j1]["RC"],
                                ">".join(_j),
                                gfaL[_j]["RC"],
                                ".",
                                ".",
                                sep=",",
                            )
                elif len(_exons0) == 0 and len(_exons1) > 0:
                    # first vertex is not on exon. So A5+ or A3-
                    exonic_prev = set(p for p in _prev1 if (p, _j[1]) in junctions)
                    # we have an event if we can reach one of the exonic nodes from _j[1]
                    visit = set([_j[0]])
                    pvisitl = 1
                    _i = 0
                    # TODO: we don't need to store the subpath, but we may want it
                    while len(visit & exonic_prev) == 0 and _i < args.isw:
                        n = visit.pop()
                        pvisitl -= 1
                        visit |= set(get_incoming_nodes(gfaS, n))
                        if pvisitl == 0:
                            _i += 1
                            pvisitl = len(visit)
                    # TODO: here we are reporting only one event per novel junction. We could do a visit for **each** exonic_next
                    if _i < args.isw:
                        # we report the event since we found a path
                        j1 = (next(iter(visit & exonic_prev)), _j[1])
                        _genes = set(
                            transcript2gene[t]
                            for t in get_haplotranscripts_from_junction(gfaL[j1]["JN"])
                        )
                        if len(_genes) > 1:
                            # CHECKME: do we need this here? We already checked this for the junction
                            # FIXME: this could be a quite strong assumption
                            print(
                                "Skipping SS due to multiple genes",
                                file=sys.sterr,
                            )
                        else:
                            _gene = next(iter(_genes))
                            print(
                                "A3" if genestrand[_gene] == "+" else "A5",
                                "novel",
                                genechr[_gene],
                                _gene,
                                genestrand[_gene],
                                "|".join(gfaL[j1]["JN"]),
                                "?",
                                ".",
                                ">".join(j1),
                                gfaL[j1]["RC"],
                                ">".join(_j),
                                gfaL[_j]["RC"],
                                ".",
                                ".",
                                sep=",",
                            )

            if "IR" in args.events:
                if len(cap) != 0:
                    eprint(f"Checking novel IR (1)", file=sys.stderr)
                    exons = _exons0 & _exons1
                    if len(exons) > 0:
                        # novel retained intron inside some exon \in exons
                        subpath = [_j[0]]
                        while subpath[-1] != _j[1]:
                            nn = -1
                            for nn in get_outgoing_nodes(gfaS, subpath[-1]):
                                # nn must be on same exon AND be smaller than the vertex we want to reach (this should hold if we assume topological sorting)
                                if nn <= _j[1] and len(exons & get_set_exons(gfaS, nn)):
                                    break
                            # we **must** have a path from the two nodes since they share exon
                            assert nn != -1, "Error while reconstrucing novel IR path"
                            subpath.append(nn)
                        if sum([gfaS[x]["LN"] for x in subpath]) >= args.minintronsize:
                            retained_transcripts = get_haplotranscripts_from_exons(
                                exons
                            )
                            _genes = set(
                                transcript2gene[t] for t in retained_transcripts
                            )
                            if len(_genes) > 1:
                                # CHECKME: do we need this here? We already checked this for the junction
                                # FIXME: this could be a quite strong assumption
                                print(
                                    "Skipping novel IR due to multiple genes",
                                    file=sys.sterr,
                                )
                            else:
                                _gene = next(iter(_genes))
                                print(
                                    "IR",
                                    "novel",
                                    genechr[_gene],
                                    _gene,
                                    genestrand[_gene],
                                    "?",
                                    "|".join(exons),
                                    ".",
                                    ">".join(_j),
                                    gfaL[_j]["RC"],
                                    ">".join(subpath),
                                    ceil(
                                        sum(
                                            [
                                                gfaS[x]["NC"] if x in gfaS else 0
                                                for x in subpath
                                            ]
                                        )
                                        / len(subpath)
                                    ),
                                    ".",
                                    ".",
                                    sep=",",
                                )

        if "IR" in args.events or "ES" in args.events:
            for _j in junctions:
                if args.junction != None and "f{_j[0]}-{_j[1]}" != args.junction:
                    continue
                if gfaL[_j]["RC"] < args.rca:
                    continue
                _ht = get_haplotranscripts_from_junction(gfaL[_j]["JN"])
                _genes = set(transcript2gene[t] for t in _ht)
                if len(_genes) > 1:
                    # FIXME: this could be a quite strong assumption
                    continue
                _gene = next(iter(_genes))

                _exons0 = get_set_exons(gfaS, _j[0])
                _exons1 = get_set_exons(gfaS, _j[1])

                eprint(
                    f"Checking annotated junction {_j[0]} -> {_j[1]}", file=sys.stderr
                )

                # we want exons on same gene only
                # FIXME: this can be done way better
                _exons0 = set(
                    e
                    for e in _exons0
                    if len(
                        _genes
                        & set(
                            transcript2gene[t]
                            for t in get_haplotranscripts_from_exon(e)
                        )
                    )
                    > 0
                )
                _exons1 = set(
                    e
                    for e in _exons1
                    if len(
                        _genes
                        & set(
                            transcript2gene[t]
                            for t in get_haplotranscripts_from_exon(e)
                        )
                    )
                    > 0
                )

                assert len(_exons0) > 0 and len(_exons1) > 0

                # Find the outgoing novel junctions of the head of the junction we are checking
                Js1 = set(x for x in noveljunctions if x[0] == _j[0]) - set([_j])
                # Find the incoming novel junctions of the tail of the junction we are checking
                Js2 = set(x for x in noveljunctions if x[1] == _j[1]) - set([_j])

                # filter by weigth # CHECKME: do we want this?
                # Js1 = set(filter(lambda x: gfaL[x]["RC"] >= args.rca, Js1))
                # Js2 = set(filter(lambda x: gfaL[x]["RC"] >= args.rca, Js2))

                # Cassete exons
                if "ES" in args.events:
                    # print(_j, Js1, Js2,  file=sys.stderr)
                    if len(Js1) > 0 and len(Js2) > 0:
                        eprint(f"Checking novel CE", file=sys.stderr)
                        novel_exons = set()
                        for j1, j2 in itertools.product(Js1, Js2):
                            # find the novel exon
                            if j1[1] <= j2[0]:
                                # this holds if we assume topological sorting
                                # CHECKME: do we need additional conditions?
                                novel_exons.add((j1[1], j2[0]))
                        for es, ee in novel_exons:
                            j1 = (_j[0], es)
                            j2 = (ee, _j[1])
                            print(
                                "CE",
                                "novel",
                                genechr[_gene],
                                _gene,
                                genestrand[_gene],
                                "|".join(gfaL[_j]["JN"]),
                                "?",
                                "?",
                                ">".join(_j),
                                gfaL[_j]["RC"],
                                ">".join(j1),
                                gfaL[j1]["RC"],
                                ">".join(j2),
                                gfaL[j2]["RC"],
                                sep=",",
                            )
                if "IR" in args.events:
                    # Assuming that we may have a variation after/before the exons, we check few edges based on topological sorting
                    if any(
                        [
                            gfaL[(_j[0], str(x))]["RC"] >= args.rca
                            for x in range(int(_j[0]) + 1, int(_j[0]) + 1 + 3)
                            if (_j[0], str(x)) in gfaL
                            and (_j[0], str(x)) not in junctions
                        ]
                    ) and any(
                        [
                            gfaL[(str(x), _j[1])]["RC"] >= args.rca
                            for x in range(int(_j[1]) - 3, int(_j[1]))
                            if (str(x), _j[1]) in gfaL
                            and (str(x), _j[1]) not in junctions
                        ]
                    ):
                        eprint(f"Checking novel IR (2)", file=sys.stderr)
                        # get exonic nodes that are at the real end or start
                        exons0_end = [
                            e
                            for e in _exons0
                            if all(
                                [
                                    e not in get_set_exons(gfaS, x)
                                    for x in get_outgoing_nodes(gfaS, _j[0])
                                ]
                            )
                        ]
                        exons1_start = [
                            e
                            for e in _exons1
                            if all(
                                [
                                    e not in get_set_exons(gfaS, x)
                                    for x in get_incoming_nodes(gfaS, _j[1])
                                ]
                            )
                        ]
                        exon_pairs = [
                            (e0, e1)
                            for (e0, e1) in itertools.product(exons0_end, exons1_start)
                            if len(
                                set(get_haplotranscripts_from_exon(e0))
                                & set(get_haplotranscripts_from_exon(e1))
                            )
                            > 0
                        ]
                        if len(exon_pairs) > 0:
                            # we just get the reference path
                            # FIXME: visit and get best path
                            subpath = [_j[0]]
                            while subpath[-1] != _j[1]:
                                onodes = get_outgoing_nodes(gfaS, subpath[-1])
                                if _j[1] in onodes:
                                    subpath.append(_j[1])
                                else:
                                    subpath.append(min(onodes))
                            if (
                                sum([gfaS[x]["LN"] for x in subpath[1:-1]])
                                >= args.minintronsize
                            ):
                                retained_transcripts = get_haplotranscripts_from_exons(
                                    [ep[0] for ep in exon_pairs]
                                )
                                _genes = set(
                                    transcript2gene[t] for t in retained_transcripts
                                )
                                if len(_genes) > 1:
                                    # CHECKME: do we need this here? We already checked this for the junction
                                    # FIXME: this could be a quite strong assumption
                                    print(
                                        "Skipping novel IR due to multiple genes",
                                        file=sys.sterr,
                                    )
                                else:
                                    _gene = next(iter(_genes))

                                    # we report just one junction on one transcript
                                    jann = "|".join(
                                        [
                                            e1 + "." + e2.split(".")[-1]
                                            for e1, e2 in exon_pairs
                                        ]
                                    )
                                    print(
                                        "IR",
                                        "novel",
                                        genechr[_gene],
                                        _gene,
                                        genestrand[_gene],
                                        jann,
                                        "?",
                                        ".",
                                        ">".join(subpath),
                                        ceil(
                                            sum(
                                                [
                                                    gfaS[x]["NC"] if x in gfaS else 0
                                                    for x in subpath
                                                ]
                                            )
                                            / len(subpath)
                                        ),
                                        ">".join(_j),
                                        gfaL[_j]["RC"],
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
        "--isw",
        dest="isw",
        help="Intronic search window for novel events, larger values reduce FP but increase time (default: 5)",
        type=int,
        default=5,
    )
    parser.add_argument(
        "--minintronsize",
        dest="minintronsize",
        help="Minimum intron size (default: 100)",
        type=int,
        default=100,
    )
    parser.add_argument(
        "--debug",
        dest="debug",
        help="Debug (default: False)",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "--junction",
        dest="junction",
        help="Junction to check, in the form '1-2' (default: None)",
        type=str,
        default=None,
    )
    args = parser.parse_args()
    if args.debug:
        eprint = print
    main(args)
