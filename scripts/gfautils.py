from typing import Any


def fa_complement(fa: str):
    _fa = ""
    for x in fa.lower():
        if x == "a":
            _fa += "T"
        elif x == "c":
            _fa += "G"
        elif x == "g":
            _fa += "C"
        elif x == "t":
            _fa += "A"
        else:
            _fa += "N"
    return _fa


class _gfanode:
    def __init__(self, nid: str, seq: str, fields: list[str]) -> None:
        self.nid = nid
        self.seq = seq
        self.len = len(self.seq)
        self.exons = []
        self.fields = fields

    # def __getattr__(self, __name: str) -> Any:
    #     return self.attrs[__name]

    # def __setattr__(self, __name: str, __value: Any) -> None:
    #     self.attrs[__name] == __value

    def addexon(self, exon: str) -> None:
        self.exons.append(exon)


class _gfapath:
    def __init__(
        self,
        pid: str,
        nodes: list[str],
        overlap: str,
        fields: list[str],
        is_reverse: bool = False,
    ) -> None:
        self.pid = pid
        self.nodes = nodes
        self.overlap = overlap
        self.fields = fields

        self.is_reverse = is_reverse

    def seq(self, nodes: list[_gfanode]):
        _fa = ""
        if not self.is_reverse:
            for n in self.nodes:
                _fa += nodes[n].seq
        else:
            for n in self.nodes:
                _fa += fa_complement(nodes[n].seq[::-1])
        return _fa

    def strnodes(self) -> str:
        _joiner = "+," if not self.is_reverse else "-,"
        return _joiner.join(self.nodes) + _joiner[0]


class _gfalink:
    def __init__(
        self,
        nid_from: str,
        orient_from: str,
        nid_to: str,
        orient_to: str,
        overlap: str,
        fields: list[str],
    ) -> None:
        self.nid_from = nid_from
        self.orient_from = orient_from
        self.nid_to = nid_to
        self.orient_to = orient_to
        self.overlap = overlap
        self.fields = fields

        self.junctions = []

    def addjunction(self, transcript: str) -> None:
        self.junctions.append(transcript)


class GFA:
    def __init__(self, filepath: str = None) -> None:
        self.nodes = {}
        self.paths = {}
        self.links = {}
        self.header = ""

        if filepath:
            for line in open(filepath):
                line = line.strip()
                if line.startswith("S"):
                    _, nid, seq, *fields = line.split()
                    self.add_node(nid, seq, fields)
                elif line.startswith("P"):
                    _, pid, p, overlap, *fields = line.split()
                    assert not ("+," in p[:-1] and "-," in p[:-1])
                    if "+," in p[:-1]:
                        self.add_path(pid, p[:-1].split("+,"), overlap, fields)
                    else:
                        self.add_path(
                            pid, p[:-1].split("-,"), overlap, fields, is_reverse=True
                        )
                elif line.startswith("L"):
                    (
                        _,
                        nid_from,
                        orient_from,
                        nid_to,
                        orient_to,
                        overlap,
                        *fields,
                    ) = line.split()
                    self.add_link(
                        nid_from, orient_from, nid_to, orient_to, overlap, fields
                    )
                elif line.startswith("H"):
                    self.header = line

    def add_node(self, nid: str, seq: str, fields: list) -> None:
        self.nodes[nid] = _gfanode(nid, seq, fields)

    def add_path(
        self,
        pid: str,
        nodes: list[str],
        overlap: str,
        fields: list[str],
        is_reverse: bool = False,
    ) -> None:
        self.paths[pid] = _gfapath(pid, nodes, overlap, fields, is_reverse)

    def add_link(
        self,
        nid_from: str,
        orient_from: str,
        nid_to: str,
        orient_to: str,
        overlap: str,
        fields: list[str],
    ) -> None:
        self.links[(nid_from, nid_to)] = _gfalink(
            nid_from, orient_from, nid_to, orient_to, overlap, fields
        )

    def link(self, nid_from: str, nid_to: str) -> _gfalink:
        return self.links[(nid_from, nid_to)]

    def pseq(self, pid: str) -> str:
        return self.paths[pid].seq(self.nodes)

    def pnodes(self, pid: str) -> str:
        for n in self.paths[pid].nodes:
            yield n

    def nlen(self, nid: str) -> int:
        return self.nodes[nid].len

    def print(self) -> None:
        if len(self.header) > 0:
            print(self.header)
        for nid in self.nodes:
            _node = self.nodes[nid]
            print(
                "S",
                nid,
                _node.seq,
                *_node.fields,
                f"LN:i:{_node.len}",
                sep="\t",
                end="",
            )
            if len(_node.exons) > 0:
                print(f'\tEX:Z:{",".join(_node.exons)}', end="")
            print()

        for lid in self.links:
            _link = self.links[lid]
            print(
                "L",
                _link.nid_from,
                _link.orient_from,
                _link.nid_to,
                _link.orient_to,
                _link.overlap,
                *_link.fields,
                sep="\t",
                end="",
            )
            if len(_link.junctions) > 0:
                print(f'\tJN:Z:{",".join(_link.junctions)}', end="")
            print()

        for pid in self.paths:
            _path = self.paths[pid]
            print(
                "P",
                _path.pid,
                _path.strnodes(),
                _path.overlap,
                sep="\t",
                *_path.fields,
                end="",
            )
            print()
